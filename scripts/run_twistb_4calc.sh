#!/usr/bin/env bash
# run_twistb_4calc.sh — run 4 TwisTB calculations for one (N, M) system:
#   1. val_bands  : band structure near VBM  (lambda_val,  num_eigs_bs)
#   2. cond_bands : band structure near CBM  (lambda_cond, num_eigs_bs)
#   3. val_wfc    : wavefunctions near VBM   (lambda_val,  num_eigs_wfc, eigvecs)
#   4. cond_wfc   : wavefunctions near CBM   (lambda_cond, num_eigs_wfc, eigvecs)
#
# Usage:
#   bash run_twistb_4calc.sh -n N -m M [-l] [-f] [-b]
#     -n N  larger supercell index  (e.g. 4 for 4_3)
#     -m M  smaller supercell index (e.g. 3 for 4_3)
#     -l    local mode (Mac, no PCT, uses no_pct/ overrides)
#     -f    full diagonalization (eig instead of eigs; ignores lambda/num_eigs)
#     -b    bands only (skip val_wfc and cond_wfc)

set -euo pipefail

# ── argument parsing ──────────────────────────────────────────────────────────
N=; M=; LOCAL=0; FULL_DIAG=0; BANDS_ONLY=0
while getopts "n:m:lfb" opt; do
    case $opt in
        n) N=$OPTARG      ;;
        m) M=$OPTARG      ;;
        l) LOCAL=1        ;;
        f) FULL_DIAG=1    ;;
        b) BANDS_ONLY=1   ;;
        *) echo "Usage: $0 -n N -m M [-l] [-f] [-b]" >&2; exit 1 ;;
    esac
done
[[ -z "$N" || -z "$M" ]] && { echo "Error: -n N -m M required" >&2; exit 1; }

echo "============================================================"
echo " TwisTB 4-calculation run: ${N}_${M}"
[[ "$FULL_DIAG" -eq 1 ]]  && echo "  mode: full diagonalization (eig)"
[[ "$BANDS_ONLY" -eq 1 ]] && echo "  scope: bands only (val_bands + cond_bands)"
echo "============================================================"

# ── physics parameters ────────────────────────────────────────────────────────
# VBM ≈ −4.40 eV, CBM ≈ −2.84 eV (gap ≈ 1.57 eV; user estimate ~1.40 eV)
# lambda_val:  −lambda just above VBM → captures near-VBM valence states
# lambda_cond: −lambda just above CBM → captures near-CBM conduction states
# Both must lie strictly inside the gap (2.84 < lambda < 4.40).
LAMBDA_VAL=4.35         # inside gap, 0.05 eV above VBM
LAMBDA_COND=2.95        # inside gap, 0.11 eV above CBM  (≈ lambda_val − 1.40)
NUM_EIGS_BS=10          # neigs=20 per k-pt — dense enough for a band plot
NUM_EIGS_WFC=5          # neigs=10 per k-pt — enough for first 2 val/cond states
KNUM=30                 # 90 total k-points along Γ-M-K-Γ

# ── paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TWISTB_ROOT="$(dirname "$SCRIPT_DIR")"         # TwisTB root = parent of scripts/
TWISTB_SRC="${TWISTB_ROOT}/TB/src"
NO_PCT_DIR="${TWISTB_ROOT}/no_pct"

if [[ "$LOCAL" -eq 1 ]]; then
    MATLAB=/Applications/MATLAB_R2026a.app/bin/matlab
    PYTHON=${HOME}/opt/miniconda3/envs/py311env/bin/python3
    LAMMPS_BASE=/Users/rdelgrande/work/TEMP_DATA/WSe2/classical_ff_calcs/twisted_2LWSe2
    WDIR_BASE=/Users/rdelgrande/work/TEMP_DATA/WSe2/tight_binding
    NUM_WORKERS=1
else
    MATLAB=matlab
    PYTHON=python3
    LAMMPS_BASE=/home/rrodriguesdelgrand/scratch/WSe2/classical_FF_relaxations
    WDIR_BASE=/home/rrodriguesdelgrand/scratch/WSe2/tight_binding
    NUM_WORKERS=${SLURM_CPUS_PER_TASK:-4}
fi

LAMMPS_FILE="${LAMMPS_BASE}/${N}_${M}/twisted_2LWSe2_relaxed.dat"
WDIR="${WDIR_BASE}/${N}_${M}"

mkdir -p "$WDIR"

# ── geometry conversion ───────────────────────────────────────────────────────
TB_FILE="${WDIR}/twisted_2LWSe2_relaxed_TB.dat"
if [[ ! -f "$TB_FILE" ]]; then
    echo "--> Converting LAMMPS → TwisTB geometry"
    "$PYTHON" "${SCRIPT_DIR}/convert_lammps_to_tb.py" \
        "$LAMMPS_FILE" "$TB_FILE" --n "$N" --m "$M"
else
    echo "--> Using existing $TB_FILE"
fi

# ── helper: run one TwisTB calculation ───────────────────────────────────────
# run_calc NAME LAMBDA NUM_EIGS EIGVECS SAVE_WS READ_NEIGHBORS SAVE_NEIGHBORS
run_calc() {
    local name=$1 lambda=$2 num_eigs=$3
    local eigvecs=$4 save_ws=$5
    local read_nb=$6 save_nb=$7

    local calcdir="${WDIR}/${name}"
    mkdir -p "$calcdir"

    # geometry symlink
    ln -sf "$TB_FILE" "${calcdir}/twisted_2LWSe2_relaxed_TB.dat"

    if [[ "$FULL_DIAG" -eq 1 ]]; then
        cat > "${calcdir}/input.dat" << EOF
task: 1
geom_file: twisted_2LWSe2_relaxed_TB.dat
nlayer: 2
tmdc: 44
convention: 2
knum: ${KNUM}
interlayer_int: true
spin_orbit: true
num_workers: ${NUM_WORKERS}
flipped: 11
mini_bz: true
reduced_workspace: false
eigvecs: ${eigvecs}
save_workspace: ${save_ws}
save_neighbors: ${save_nb}
read_neighbors: ${read_nb}
EOF
    else
        cat > "${calcdir}/input.dat" << EOF
task: 1
geom_file: twisted_2LWSe2_relaxed_TB.dat
nlayer: 2
tmdc: 44
convention: 2
knum: ${KNUM}
interlayer_int: true
spin_orbit: true
num_workers: ${NUM_WORKERS}
flipped: 11
mini_bz: true
reduced_workspace: true
num_eigs: ${num_eigs}
lambda: ${lambda}
eigvecs: ${eigvecs}
save_workspace: ${save_ws}
save_neighbors: ${save_nb}
read_neighbors: ${read_nb}
EOF
    fi

    echo ""
    echo "--> ${name}  (lambda=${lambda}, num_eigs=${num_eigs}, eigvecs=${eigvecs})"

    if [[ "$LOCAL" -eq 1 ]]; then
        "$MATLAB" -nodisplay -nodesktop -nosplash \
            -r "cd('${calcdir}'); addpath('${TWISTB_SRC}'); addpath('${NO_PCT_DIR}'); main; exit" \
            2>&1 | tee "${calcdir}/matlab_run.log"
    else
        "$MATLAB" -nodisplay -nodesktop -nosplash \
            -r "cd('${calcdir}'); addpath('${TWISTB_SRC}'); main; exit" \
            2>&1 | tee "${calcdir}/matlab_run.log"
    fi

    local rc=${PIPESTATUS[0]}
    [[ $rc -ne 0 ]] && echo "WARNING: MATLAB exited with status $rc for ${name}" >&2
    echo "--> ${name} done"
}

# ── calculation 1 — valence band structure ────────────────────────────────────
# save_neighbors=true so the neighbor list is cached for calcs 2-4
run_calc "val_bands" "$LAMBDA_VAL" "$NUM_EIGS_BS" "false" "false" "false" "true"

# ── share neighbor list with remaining calcs ─────────────────────────────────
SAVE_DIR_1=$(find "${WDIR}/val_bands" -maxdepth 1 -name "*.save" -type d 2>/dev/null | head -1)
if [[ -n "$SAVE_DIR_1" && -f "${SAVE_DIR_1}/neighbors.mat" ]]; then
    OUTNAME=$(basename "$SAVE_DIR_1" .save)
    if [[ "$BANDS_ONLY" -eq 1 ]]; then
        share_list="cond_bands"
    else
        share_list="cond_bands val_wfc cond_wfc"
    fi
    echo "--> Sharing neighbor list (${OUTNAME}.save/neighbors.mat) with: ${share_list}"
    for calc in $share_list; do
        mkdir -p "${WDIR}/${calc}/${OUTNAME}.save"
        cp "${SAVE_DIR_1}/neighbors.mat" "${WDIR}/${calc}/${OUTNAME}.save/"
    done
    READ_NB=true
else
    echo "WARNING: neighbors.mat not found after calc 1 — remaining calcs will recompute" >&2
    READ_NB=false
fi

# ── calculation 2 — conduction band structure ─────────────────────────────────
run_calc "cond_bands" "$LAMBDA_COND" "$NUM_EIGS_BS" "false" "false" "$READ_NB" "false"

if [[ "$BANDS_ONLY" -eq 0 ]]; then
    # ── calculation 3 — valence wavefunctions ─────────────────────────────────
    run_calc "val_wfc" "$LAMBDA_VAL" "$NUM_EIGS_WFC" "true" "true" "$READ_NB" "false"

    # ── calculation 4 — conduction wavefunctions ──────────────────────────────
    run_calc "cond_wfc" "$LAMBDA_COND" "$NUM_EIGS_WFC" "true" "true" "$READ_NB" "false"
fi

echo ""
echo "============================================================"
if [[ "$BANDS_ONLY" -eq 1 ]]; then
    echo " Band structure calculations done for ${N}_${M} (val_bands + cond_bands)"
else
    echo " All 4 calculations done for ${N}_${M}"
fi
echo "============================================================"
