```
                ______                         ______  ____
               /\__  _\            __         /\__  _\/\  _`\
               \/_/\ \/ __  __  __/\_\    ____\/_/\ \/\ \ \L\ \
                  \ \ \/\ \/\ \/\ \/\ \  /',__\  \ \ \ \ \  _ <'
                   \ \ \ \ \_/ \_/ \ \ \/\__, `\  \ \ \ \ \ \L\ \
                    \ \_\ \___x___/'\ \_\/\____/   \ \_\ \ \____/
                     \/_/\/__//__/   \/_/\/___/     \/_/  \/___/

```

# TwisTB — Fork

> **This is a fork** of the original [TwisTB](https://github.com/Valerio-Vitale/TwisTB) repository
> by Valerio Vitale and Kemal Atalar (University of Trieste / Imperial College London).
> The original code is described in:
>
> Vitale, V., Atalar, K., *et al.* "Flat band properties of twisted transition metal
> dichalcogenide homo- and heterobilayers of MoS₂, MoSe₂, WS₂ and WSe₂"
> *2D Materials* **8**, 045010 (2021).
>
> **Please cite this paper** in any publication arising from use of this code.

## What this fork adds

| Addition | Location | Purpose |
|---|---|---|
| `run_twistb.sh` | repo root | End-to-end workflow driver: parse `input.dat`, convert geometry, run MATLAB |
| `scripts/convert_lammps_to_tb.py` | `scripts/` | Convert LAMMPS molecular data → TwisTB v2.0 geometry format |
| `scripts/convert_asg_to_tb.py` | `scripts/` | Convert ASG 9-line format → TwisTB v2.0 |
| `scripts/plot_bands.py` | `scripts/` | Plot band structure from `*_BS.dat` output |
| `scripts/make_2h_tb.py` | `scripts/` | Generate 2H bilayer geometry from QE positions |
| `TB/src/plot_DOS.m` | `TB/src/` | DOS via Gaussian smearing (was missing from original) |
| `TB/src/plot_wfc_new.m` | `TB/src/` | Fixed sublayer detection for LAMMPS-relaxed structures |
| `no_pct/` | `no_pct/` | Serial overrides for machines without MATLAB PCT license |

---

## Quick Start — LAMMPS-relaxed twisted bilayer TMD

### Prerequisites

- MATLAB (any recent release)
- Python 3 with NumPy (`pip install numpy` or use conda)
- A LAMMPS molecular-style data file of your twisted bilayer structure

### Step 1 — Clone and set up

```bash
git clone https://github.com/rafaeldgrande/TwisTB.git
```

No compilation required.

### Step 2 — Create a calculation directory

```bash
mkdir my_calc && cd my_calc
cp /path/to/my_structure.lammps .
```

### Step 3 — Write `input.dat`

```
task:           1          # 1=band structure, 3=DOS, 6=wavefunction
geom_file:      my_structure.lammps
nlayer:         2
tmdc:           44         # 44=WSe2, 42=MoS2, 43=WS2, 41=MoSe2
convention:     2
knum:           30         # k-points per segment (band) or per axis (DOS)
interlayer_int: true
spin_orbit:     false
num_workers:    1
flipped:        11         # see "flipped convention" section below
mini_bz:        true
```

### Step 4 — Run

```bash
bash /path/to/TwisTB/run_twistb.sh
```

The script will:
1. Parse `input.dat`
2. Convert `my_structure.lammps` → `my_structure_TB.dat` (TwisTB v2.0 format)
3. Rewrite `input.dat` pointing to the converted file
4. Copy serial PCT overrides from `no_pct/` into the working directory (if not already present)
5. Launch MATLAB and run `main.m`
6. Log all output to `matlab_run.log`

### Step 5 — Plot the band structure

```bash
python /path/to/TwisTB/scripts/plot_bands.py my_structure_BS.dat
```

---

## `input.dat` reference

| Key | Default | Description |
|---|---|---|
| `task` | `1` | 1 = band structure, 3 = DOS, 6 = wavefunction |
| `geom_file` | — | Geometry file: raw LAMMPS `.lammps` or converted `*_TB.dat` |
| `nlayer` | `2` | Number of layers |
| `tmdc` | `44` | Material code: 44=WSe₂, 42=MoS₂, 43=WS₂, 41=MoSe₂ |
| `convention` | `2` | SK convention (2 = Fang PRB 2015) |
| `knum` | `30` | k-points per BZ segment (band) or per axis (DOS mesh = knum×knum) |
| `interlayer_int` | `true` | Include interlayer hopping |
| `spin_orbit` | `false` | Include spin-orbit coupling |
| `num_workers` | `1` | Parallel workers (requires PCT license; see `no_pct/README.md`) |
| `flipped` | `11` | Per-layer sign convention flag — **read the section below** |
| `mini_bz` | `true` | Use moiré mini BZ (recommended for twisted bilayers) |

### The `flipped` convention

`flipped` is a per-layer integer string where each digit indicates whether
the tungsten (or Mo) atom sits at the fractional origin (0, 0) of the unit cell:

| Digit | Meaning |
|---|---|
| `1` | W/Mo is at fractional position (0, 0) — SK Hamiltonian sign is *flipped* |
| `0` | W/Mo is at fractional position (1/3, 1/3) — standard sign |

**LAMMPS-generated twisted bilayers** (`generate_twisted_2L_*.py`): both layers
use W at (0, 0) as reference → `flipped: 11`

**2H QE-derived bilayers** (`make_2h_tb.py`): layer 1 has W at (0, 0),
layer 2 has W at (1/3, 1/3) → `flipped: 10`

> **Common mistake**: leaving `flipped: 0` (or `flipped: 00`) for a LAMMPS
> twisted bilayer gives a band gap of ~0.12 eV instead of the correct ~1.9 eV
> for WSe₂. Always verify the gap against the known monolayer value (~1.7–2.0 eV)
> before trusting the twisted-bilayer spectrum.

---

## Directory structure

```
TwisTB/
├── run_twistb.sh          # workflow driver (run from your calc dir)
├── no_pct/                # serial overrides for machines without PCT
│   ├── opt_build_and_diag_H_3rdnn.m
│   ├── parpool.m
│   ├── parcluster.m
│   ├── gcp.m
│   └── README.md
├── scripts/               # Python utilities
│   ├── convert_lammps_to_tb.py
│   ├── convert_asg_to_tb.py
│   ├── plot_bands.py
│   └── make_2h_tb.py
├── ASG/                   # Atomic Structures Generator (original)
└── TB/
    └── src/               # MATLAB tight-binding source (original + fixes)
        ├── main.m
        ├── plot_DOS.m     # added in this fork
        ├── plot_wfc_new.m # fixed in this fork
        └── ...
```

---

## Cluster / SLURM usage (serial, no PCT)

```bash
#!/bin/bash
#SBATCH --job-name=twistb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00

module load matlab

cd $SLURM_SUBMIT_DIR
bash /path/to/TwisTB/run_twistb.sh
```

The `no_pct/` serial overrides are copied automatically by `run_twistb.sh`.
No PCT license is required.

### Parallel execution with PCT

If your cluster has a MATLAB PCT license, delete the `no_pct/` override files
from your working directory and set `num_workers` in `input.dat`:

```
num_workers: 16
```

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
```

---

## Environment variables

| Variable | Default | Purpose |
|---|---|---|
| `MATLAB_BIN` | auto-detected | Full path to the `matlab` binary |
| `TWISTB_ROOT` | directory of `run_twistb.sh` | Override repo root |
| `PYTHON` | `python3` (or conda env if found) | Python interpreter for geometry conversion |


---

# Original TwisTB documentation

## Authors (original code)
* Valerio Vitale — MATLAB codes for generation of initial structures
* Valerio Vitale and Kemal Atalar — MATLAB codes for tight-binding calculations

## ASG — Atomic Structures Generator

This module generates atomic structures of twisted homo- and hetero-bilayers of 2D
materials with a hexagonal unit cell (graphene, hBN, TMDs).
Input parameters are specified in the `main.m` file inside the relevant ASG subdirectory.

Outputs:
1. `<rootname>.xyz` / `<rootname>.xsf` for visualisation
2. `lammps_positions.<rootname>.dat` for use as LAMMPS input

For twisted homo-bilayers the twist angle is specified by a pair of integers (n, m);
for hetero-bilayers a target angle and lattice parameter ratio are used.
See the header of `main.m` in the relevant ASG subdirectory for full parameter descriptions.

## TB — Tight-Binding

See the header of `TB/src/read_input.m` for a full description of `input.dat` keywords.

## Funding

This project has received funding from the European Union's Horizon 2020 research
and innovation programme under the Marie Skłodowska-Curie grant agreement no. 101067977.
