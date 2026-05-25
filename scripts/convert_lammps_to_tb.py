#!/usr/bin/env python3
"""
Convert a LAMMPS molecular-style data file for a twisted TMD bilayer
into the TwisTB geometry file format expected by read_structure.m.

Atom positions are preserved exactly — no rotation is applied.
The theta values written to the file are the crystallographic orientation
angles of each layer's lattice frame, used by TwisTB to orient the
Slater-Koster hopping tables.

Usage:
  python3 convert_lammps_to_tb.py input.lammps [output.dat] [--n N --m M]

  --n N  --m M   Supercell indices from generate_twisted_2L_WSe2_new.py.
                 When provided, layer thetas are computed analytically
                 (exact and instant). Otherwise they are estimated from
                 W-W bond directions via ASE NeighborList (slower).
"""

import sys
import os
import time
import argparse
from collections import Counter

import numpy as np
from ase.io import read
from ase.neighborlist import NeighborList


def thetas_from_nm(N, M):
    """
    Compute the two TwisTB layer orientation angles exactly from (N, M).

    The generator rotates layer 1 by +twist/2 and layer 2 by -twist/2
    from a reference lattice whose W-W NN direction is at 30° from x.
    Angles are folded into [0°, 60°) as TwisTB expects.
    """
    cos_twist = (N**2 + 4*N*M + M**2) / (2 * (N**2 + N*M + M**2))
    twist_deg = np.degrees(np.arccos(cos_twist))
    theta1 = (30.0 + twist_deg / 2.0) % 60.0
    theta2 = (30.0 - twist_deg / 2.0) % 60.0
    return theta1, theta2, twist_deg


def _layer_theta_from_positions(W_atoms, cutoff=4.5):
    """
    Fallback: crystallographic orientation angle from W-W NN bond directions.
    Returns the median angle folded to [0°, 60°).
    ASE NeighborList handles PBC efficiently.
    """
    n = len(W_atoms)
    pos = W_atoms.get_positions()
    cell = W_atoms.get_cell()

    nl = NeighborList([cutoff / 2] * n, self_interaction=False, bothways=False)
    nl.update(W_atoms)

    angles = []
    for i in range(n):
        indices, offsets = nl.get_neighbors(i)
        for j, offset in zip(indices, offsets):
            dr = pos[j] + offset @ cell - pos[i]
            angles.append(np.degrees(np.arctan2(dr[1], dr[0])) % 60.0)

    if not angles:
        raise RuntimeError("No W-W neighbours found — check cutoff")
    return float(np.median(angles))


def write_tb_geometry(atoms, out_fname, nm=None):
    """
    nm : tuple (N, M) or None.
         If given, layer thetas are computed analytically from the supercell
         indices. Otherwise estimated from W-W bond directions.
    """
    mol_ids = atoms.get_array('mol-id')
    layers = sorted(set(mol_ids))
    nlayer = len(layers)

    cell = atoms.get_cell()
    a1_M = np.array(cell[0])
    a2_M = np.array(cell[1])
    a3   = np.array([0., 0., 40.])

    syms = atoms.get_chemical_symbols()
    layer_idx = {lyr: [i for i, m in enumerate(mol_ids) if m == lyr] for lyr in layers}
    nat_per_layer = len(layer_idx[layers[0]])
    n_prim = nat_per_layer // 3
    a_mono = np.linalg.norm(a1_M) / np.sqrt(n_prim)

    a1_mono = np.array([a_mono, 0., 0.])
    a2_mono = np.array([a_mono / 2., a_mono * np.sqrt(3) / 2., 0.])

    if nm is not None:
        N, M = nm
        theta1, theta2, twist_deg = thetas_from_nm(N, M)
        thetas = [theta1, theta2]
        print(f"  Thetas from (N={N}, M={M}): {theta1:.6f}°  {theta2:.6f}°"
              f"  (twist = {twist_deg:.4f}°)", flush=True)
    else:
        thetas = []
        for lyr in layers:
            W_idx = [i for i in layer_idx[lyr] if syms[i] == 'W']
            print(f"  Computing theta for layer {lyr} ({len(W_idx)} W atoms)...", flush=True)
            t0 = time.time()
            theta = _layer_theta_from_positions(atoms[W_idx])
            thetas.append(theta)
            print(f"  Layer {lyr} theta = {theta:.4f}° [{time.time()-t0:.1f}s]", flush=True)

    pos = atoms.get_positions()
    with open(out_fname, 'w') as f:
        f.write('  '.join(str(nat_per_layer) for _ in layers) + '\n')
        f.write('  '.join(f'{a_mono:.6f}' for _ in layers) + '\n')
        f.write('  '.join(f'{t:.6f}' for t in thetas) + '\n')
        f.write('  '.join('1.00' for _ in layers) + '\n')
        for v in [a1_mono, a2_mono, a3]:
            f.write(f'{v[0]:.8f}   {v[1]:.8f}   {v[2]:.8f}\n')
        for v in [a1_M, a2_M, a3]:
            f.write(f'{v[0]:.8f}   {v[1]:.8f}   {v[2]:.8f}\n')
        for lyr in layers:
            for i in layer_idx[lyr]:
                f.write(f"{syms[i]} {lyr} {pos[i,0]:.6f} {pos[i,1]:.6f} {pos[i,2]:.6f}\n")

    print(f"Written: {out_fname}")
    print(f"  {nlayer} layers, {nat_per_layer} atoms/layer, a_mono = {a_mono:.4f} Å")
    for lyr, t in zip(layers, thetas):
        print(f"  Layer {lyr}: theta = {t:.4f}°")
    if nlayer == 2:
        print(f"  Twist angle = {abs(thetas[1] - thetas[0]):.4f}°")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input',  help='LAMMPS data file')
    parser.add_argument('output', nargs='?', help='Output TwisTB file (default: input_TB.dat)')
    parser.add_argument('--n', type=int, default=None, metavar='N',
                        help='Supercell index N from generate_twisted_2L_WSe2_new.py')
    parser.add_argument('--m', type=int, default=None, metavar='M',
                        help='Supercell index M from generate_twisted_2L_WSe2_new.py')
    args = parser.parse_args()

    in_fname  = args.input
    out_fname = args.output or os.path.splitext(in_fname)[0] + '_TB.dat'
    nm = (args.n, args.m) if (args.n is not None and args.m is not None) else None

    t_start = time.time()

    print(f"==> Reading LAMMPS file: {in_fname}", flush=True)
    atoms = read(in_fname, format='lammps-data', atom_style='molecular')
    counts = Counter(atoms.get_chemical_symbols())
    print(f"    {len(atoms)} atoms: " + ", ".join(f"{n} {s}" for s, n in sorted(counts.items())),
          flush=True)

    if nm:
        print(f"==> Using analytical thetas for N={nm[0]}, M={nm[1]}", flush=True)
    else:
        print(f"==> Estimating thetas from W-W bond directions (no --n/--m given)", flush=True)

    print(f"==> Writing: {out_fname}", flush=True)
    write_tb_geometry(atoms, out_fname, nm=nm)

    print(f"==> Done in {time.time() - t_start:.1f}s", flush=True)


if __name__ == '__main__':
    main()
