# no_pct — Serial override for machines without MATLAB PCT

This directory contains a serial replacement for `opt_build_and_diag_H_3rdnn.m`
and stub implementations of the Parallel Computing Toolbox functions
(`parpool`, `parcluster`, `gcp`).

## How it works

MATLAB resolves function names by searching paths in order.  `run_twistb.sh`
calls MATLAB with:

```
addpath(WORKDIR)           ← searched first  → picks up serial overrides
addpath(TwisTB/TB/src)     ← searched second → rest of TwisTB
```

The serial `opt_build_and_diag_H_3rdnn.m` shadows the parallel version in
`TB/src`, so no PCT license is needed.

## Files

| File | Purpose |
|------|---------|
| `opt_build_and_diag_H_3rdnn.m` | Serial k-loop replacement (no `parfor`/`spmd`) |
| `parpool.m` | Stub — returns `[]`, no-op |
| `parcluster.m` | Stub — returns dummy struct |
| `gcp.m` | Stub — returns `[]`, no-op |

## Using with PCT (parallel execution)

If your cluster has a PCT license, do **not** copy these files to your working
directory (or delete them if already copied).  `run_twistb.sh` will not
overwrite existing files, so you can remove them selectively.

With PCT you also need `num_workers > 1` in `input.dat` and a SLURM header
such as:

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
```

Then set `num_workers: 16` in `input.dat`.
