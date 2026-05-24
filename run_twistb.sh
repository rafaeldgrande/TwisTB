#!/usr/bin/env bash
# run_twistb.sh  —  Full TwisTB workflow for a twisted bilayer TMD
#
# Usage (run from your calculation directory):
#   bash /path/to/TwisTB/run_twistb.sh [OPTIONS]
#
# Expects in the current directory:
#   input.dat          — TwisTB input file (see below)
#   <lammps_file>      — LAMMPS molecular-style data file named in input.dat
#
# input.dat keys (colon-delimited):
#   task            1=band structure, 3=DOS, 6=wavefunction
#   geom_file       LAMMPS data file (raw) or *_TB.dat (already converted)
#   nlayer          2
#   tmdc            44 = WSe2, 42 = MoS2, etc.
#   convention      2
#   knum            k-points per segment (band) or per axis (DOS/Chern)
#   interlayer_int  true/false
#   spin_orbit      true/false
#   num_workers     number of parallel workers (requires PCT)
#   flipped         per-layer flip flag (e.g. 11 for LAMMPS AA-stacked twisted)
#   mini_bz         true/false  (use mini BZ for twisted bilayer)
#
# Environment overrides:
#   MATLAB_BIN      full path to matlab binary
#   TWISTB_ROOT     override repo root detection
#   PYTHON          override python3 interpreter
set -euo pipefail

# ── Locate repo root and working directory ───────────────────────────────────
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
TWISTB_ROOT="${TWISTB_ROOT:-$(cd "$(dirname "$SCRIPT_PATH")" && pwd)}"
WORKDIR="$(pwd)"
TWISTB_SRC="$TWISTB_ROOT/TB/src"
SCRIPTS_DIR="$TWISTB_ROOT/scripts"
NO_PCT_DIR="$TWISTB_ROOT/no_pct"
INPUT_DAT="$WORKDIR/input.dat"

# ── Find MATLAB ───────────────────────────────────────────────────────────────
if [[ -z "${MATLAB_BIN:-}" ]]; then
    if command -v matlab &>/dev/null; then
        MATLAB_BIN="$(command -v matlab)"
    else
        # Common macOS location
        MATLAB_GUESS=$(ls -d /Applications/MATLAB_*.app/bin/matlab 2>/dev/null | sort -r | head -1)
        MATLAB_BIN="${MATLAB_GUESS:-matlab}"
    fi
fi

# ── Find Python ───────────────────────────────────────────────────────────────
if [[ -z "${PYTHON:-}" ]]; then
    CONDA_PY="$HOME/opt/miniconda3/envs/py311env/bin/python3"
    if [[ -x "$CONDA_PY" ]]; then
        PYTHON="$CONDA_PY"
    else
        PYTHON="$(command -v python3)"
    fi
fi

echo "==> TwisTB run script"
echo "    repo root  : $TWISTB_ROOT"
echo "    workdir    : $WORKDIR"
echo "    matlab     : $MATLAB_BIN"
echo "    python     : $PYTHON"

# ── Parse input.dat ───────────────────────────────────────────────────────────
if [[ ! -f "$INPUT_DAT" ]]; then
    echo "ERROR: input.dat not found in $WORKDIR" >&2; exit 1
fi

get_val() {
    local key="$1" default="${2:-}"
    awk -v key="$key" '
        {
            sub(/#.*/, "")
            gsub(/\t/, " ")
            gsub(/^ +| +$/, "")
            if (length($0) == 0) next
            if (index($0, ":")) {
                n = split($0, a, /:/)
                k = a[1]; gsub(/^ +| +$/, "", k)
                v = ""; for (i=2; i<=n; i++) v = v (i>2 ? ":" : "") a[i]
                gsub(/^ +| +$/, "", v)
            } else {
                k = $1; v = $2
            }
            if (k == key) { print v; exit }
        }
    ' "$INPUT_DAT" || echo "$default"
}

GEOM_RAW=$(get_val "geom_file" "")
[[ -z "$GEOM_RAW" ]] && GEOM_RAW=$(get_val "geomfname" "")
if [[ -z "$GEOM_RAW" ]]; then
    echo "ERROR: input.dat must contain 'geom_file' or 'geomfname'" >&2; exit 1
fi

TASK=$(get_val        "task"           "1")
NLAYER=$(get_val      "nlayer"         "2")
TMDC=$(get_val        "tmdc"           "44")
CONVENTION=$(get_val  "convention"     "2")
KNUM=$(get_val        "knum"           "30")
INTERLAYER=$(get_val  "interlayer_int" "true")
SOC=$(get_val         "spin_orbit"     "false")
NUM_WORKERS=$(get_val "num_workers"    "1")
FLIPPED=$(get_val     "flipped"        "11")
MINI_BZ=$(get_val     "mini_bz"        "true")

echo ""
echo "==> Parsed input.dat:"
echo "    geom_file     : $GEOM_RAW"
echo "    task/tmdc     : $TASK / $TMDC"
echo "    nlayer/knum   : $NLAYER / $KNUM"
echo "    flipped       : $FLIPPED"
echo "    spin_orbit    : $SOC"

# ── Convert LAMMPS → TwisTB if needed ────────────────────────────────────────
LAMMPS_FILE="$WORKDIR/$GEOM_RAW"
if [[ ! -f "$LAMMPS_FILE" ]]; then
    echo "ERROR: geometry file not found: $LAMMPS_FILE" >&2; exit 1
fi

BASENAME="${GEOM_RAW%.*}"

if [[ "$GEOM_RAW" == *_TB.dat ]]; then
    TB_FILE="$GEOM_RAW"
    echo ""
    echo "==> $GEOM_RAW is already in TwisTB format — skipping conversion"
else
    TB_FILE="${BASENAME}_TB.dat"
    echo ""
    echo "==> Converting $GEOM_RAW → $TB_FILE"
    "$PYTHON" "$SCRIPTS_DIR/convert_lammps_to_tb.py" "$LAMMPS_FILE" "$WORKDIR/$TB_FILE"
fi

# ── Rewrite input.dat with colon format and all mandatory fields ─────────────
echo ""
echo "==> Writing input.dat"
cat > "$INPUT_DAT" <<EOF
task: $TASK
geom_file: $TB_FILE
nlayer: $NLAYER
tmdc: $TMDC
convention: $CONVENTION
knum: $KNUM
interlayer_int: $INTERLAYER
spin_orbit: $SOC
num_workers: $NUM_WORKERS
flipped: $FLIPPED
mini_bz: $MINI_BZ
EOF

# ── Copy no-PCT serial overrides ─────────────────────────────────────────────
# These must be in WORKDIR so they shadow the PCT versions in TB/src.
# Remove them from WORKDIR if you have a PCT license and want parallel execution.
echo ""
echo "==> Installing serial (no-PCT) overrides from $NO_PCT_DIR"
for f in opt_build_and_diag_H_3rdnn.m parpool.m parcluster.m gcp.m; do
    if [[ ! -f "$WORKDIR/$f" ]]; then
        cp "$NO_PCT_DIR/$f" "$WORKDIR/$f"
        echo "    copied $f"
    else
        echo "    $f already present — not overwriting"
    fi
done

# ── Run MATLAB ────────────────────────────────────────────────────────────────
echo ""
echo "==> Running TwisTB (task=$TASK, knum=$KNUM)"
"$MATLAB_BIN" -nodisplay -nodesktop -nosplash \
    -r "cd('$WORKDIR'); addpath('$WORKDIR'); addpath('$TWISTB_SRC'); main; exit" \
    2>&1 | tee "$WORKDIR/matlab_run.log"

echo ""
echo "==> Done. Output files:"
find "$WORKDIR" -maxdepth 1 \( -name "*_BS.dat" -o -name "*_DOS.dat" -o -name "*.png" \) \
    | sort || echo "    (none found — check matlab_run.log)"
