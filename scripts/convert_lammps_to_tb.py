#!/usr/bin/env python3
"""
Convert a LAMMPS molecular-style data file for a twisted TMD bilayer
into the TwisTB geometry file format expected by read_structure.m.

Atom positions are preserved exactly — no rotation is applied.
The theta values written to the file are the crystallographic orientation
angles of each layer's lattice frame, used by TwisTB to orient the
Slater-Koster hopping tables.

Usage: python3 convert_lammps_to_tb.py input.lammps [output.dat]
"""

import numpy as np
import sys
import os
import time

# ── species identification from atomic mass ──────────────────────────────────
MASS_TO_SPECIES = {
    'Mo':  (95.94,  3.0),   # (nominal mass, tolerance)
    'W':   (183.85, 3.0),
    'S':   (32.06,  1.0),
    'Se':  (78.96,  1.0),
}

def identify_species(mass):
    for name, (ref, tol) in MASS_TO_SPECIES.items():
        if abs(mass - ref) < tol:
            return name
    raise ValueError(f"Cannot identify species for mass {mass:.3f}")

def is_metal(species):
    return species in ('Mo', 'W')

# ── LAMMPS parser ─────────────────────────────────────────────────────────────
def parse_lammps(fname):
    with open(fname) as f:
        lines = f.readlines()

    box = {}
    masses = {}
    atoms = []
    section = None

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue

        if 'xlo xhi' in stripped:
            p = stripped.split(); box['xlo'], box['xhi'] = float(p[0]), float(p[1])
        elif 'ylo yhi' in stripped:
            p = stripped.split(); box['ylo'], box['yhi'] = float(p[0]), float(p[1])
        elif 'zlo zhi' in stripped:
            p = stripped.split(); box['zlo'], box['zhi'] = float(p[0]), float(p[1])
        elif 'xy xz yz' in stripped:
            p = stripped.split(); box['xy'] = float(p[0])
        elif stripped == 'Masses':
            section = 'masses'
        elif stripped.startswith('Atoms'):
            section = 'atoms'
        elif section == 'masses':
            p = stripped.split()
            if len(p) >= 2:
                try:
                    masses[int(p[0])] = float(p[1])
                except ValueError:
                    pass
        elif section == 'atoms':
            p = stripped.split()
            if len(p) >= 6:
                try:
                    atom_id  = int(p[0])
                    mol_id   = int(p[1])   # layer index (1-based)
                    atype    = int(p[2])
                    x, y, z  = float(p[3]), float(p[4]), float(p[5])
                    ix = int(p[6]) if len(p) > 6 else 0
                    iy = int(p[7]) if len(p) > 7 else 0
                    atoms.append(dict(id=atom_id, layer=mol_id, atype=atype,
                                      x=x, y=y, z=z, ix=ix, iy=iy))
                except (ValueError, IndexError):
                    section = None   # hit "Velocities" or similar section

    Lx = box['xhi'] - box['xlo']
    Ly = box['yhi'] - box['ylo']
    xy = box.get('xy', 0.0)

    # apply image-flag unwrapping for periodic boundary correction
    for a in atoms:
        a['xu'] = a['x'] + a['ix'] * Lx + a['iy'] * xy
        a['yu'] = a['y'] + a['iy'] * Ly

    # map atom types to species names
    type_to_species = {t: identify_species(m) for t, m in masses.items()}
    for a in atoms:
        a['species'] = type_to_species[a['atype']]

    return atoms, Lx, Ly, xy, masses, type_to_species


# ── geometry computation ──────────────────────────────────────────────────────
def compute_geometry(atoms, Lx, Ly, xy):
    """Compute moiré cell, monolayer lattice parameter, and layer theta values."""
    a1_M = np.array([Lx, 0.0, 0.0])
    a2_M = np.array([xy, Ly, 0.0])

    layers = sorted(set(a['layer'] for a in atoms))
    nlayer = len(layers)
    nat_per_layer = sum(1 for a in atoms if a['layer'] == layers[0])

    # number of primitive TMD units per layer (3 atoms each: M + 2X)
    n_prim = nat_per_layer // 3

    # monolayer lattice parameter from moiré cell magnitude
    a_mono = np.linalg.norm(a1_M) / np.sqrt(n_prim)

    # monolayer unit cell (convention 2, 60°)
    a1_mono = np.array([a_mono, 0.0, 0.0])
    a2_mono = np.array([a_mono / 2.0, a_mono * np.sqrt(3) / 2.0, 0.0])

    # theta for each layer: angle of the W/Mo NN direction relative to x-axis
    thetas = []
    for lyr in layers:
        metals = [a for a in atoms if a['layer'] == lyr and is_metal(a['species'])]
        print(f"  Computing theta for layer {lyr} ({len(metals)} metal atoms)...", flush=True)
        theta = _layer_theta(metals, Lx, Ly, xy, a_mono, layer_label=lyr)
        thetas.append(theta)
        print(f"  Layer {lyr} theta = {theta:.4f}°", flush=True)

    return a1_M, a2_M, a1_mono, a2_mono, a_mono, thetas, layers, nat_per_layer


def _nn_vector(x0, y0, metals, Lx, Ly, xy, a_mono):
    """Find the nearest-neighbour metal-metal vector using periodic images."""
    best_d, best_v = np.inf, None
    for m in metals:
        dx0 = m['xu'] - x0
        dy0 = m['yu'] - y0
        for nx in range(-2, 3):
            for ny in range(-2, 3):
                dx = dx0 + nx * Lx + ny * xy
                dy = dy0 + ny * Ly
                d = np.hypot(dx, dy)
                if 0.5 * a_mono < d < 1.5 * a_mono and d < best_d:
                    best_d, best_v = d, (dx, dy)
    return best_v


def _layer_theta(metals, Lx, Ly, xy, a_mono, layer_label=None):
    """Robust theta estimate: median of all NN metal-metal bond angles folded to [0°, 60°)."""
    angles = []
    n = len(metals)
    report_every = max(1, n // 10)
    t0 = time.time()
    for i, ref in enumerate(metals):
        v = _nn_vector(ref['xu'], ref['yu'],
                       [m for m in metals if m['id'] != ref['id']],
                       Lx, Ly, xy, a_mono)
        if v is not None:
            # fold to [0°, 60°) — fundamental domain of the hexagonal lattice
            ang = np.degrees(np.arctan2(v[1], v[0])) % 60.0
            angles.append(ang)
        if (i + 1) % report_every == 0 or (i + 1) == n:
            elapsed = time.time() - t0
            tag = f" (layer {layer_label})" if layer_label is not None else ""
            print(f"    theta{tag}: {i+1}/{n} metal atoms processed  [{elapsed:.1f}s]",
                  flush=True)
    if not angles:
        raise RuntimeError("Could not find NN metal-metal vector")
    return float(np.median(angles))


# ── writer ────────────────────────────────────────────────────────────────────
def write_tb_geometry(atoms, a1_M, a2_M, a1_mono, a2_mono, a_mono,
                      thetas, layers, nat_per_layer, out_fname):
    nlayer = len(layers)
    a3 = np.array([0.0, 0.0, 40.0])  # vacuum in z

    with open(out_fname, 'w') as f:
        # line 1: nat per layer
        f.write('  '.join(str(nat_per_layer) for _ in layers) + '\n')
        # line 2: monolayer lattice parameter per layer
        f.write('  '.join(f'{a_mono:.6f}' for _ in layers) + '\n')
        # line 3: theta per layer (degrees)
        f.write('  '.join(f'{t:.6f}' for t in thetas) + '\n')
        # line 4: strain per layer (no strain)
        f.write('  '.join('1.00' for _ in layers) + '\n')
        # lines 5-7: monolayer unit cell vectors
        for v in [a1_mono, a2_mono, a3]:
            f.write(f'{v[0]:.8f}   {v[1]:.8f}   {v[2]:.8f}\n')
        # lines 8-10: moiré supercell vectors
        for v in [a1_M, a2_M, a3]:
            f.write(f'{v[0]:.8f}   {v[1]:.8f}   {v[2]:.8f}\n')
        # atoms: all layer-1 atoms first, then layer-2, etc.
        for lyr in layers:
            for a in atoms:
                if a['layer'] == lyr:
                    f.write(f"{a['species']} {lyr} "
                            f"{a['xu']:.6f} {a['yu']:.6f} {a['z']:.6f}\n")

    print(f"Written TwisTB geometry: {out_fname}")
    print(f"  {nlayer} layers, {nat_per_layer} atoms/layer")
    print(f"  a_mono = {a_mono:.4f} Å")
    for i, (lyr, t) in enumerate(zip(layers, thetas)):
        print(f"  Layer {lyr}: theta = {t:.4f}°")
    if nlayer == 2:
        print(f"  Twist angle = {abs(thetas[1]-thetas[0]):.4f}°")


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    in_fname = sys.argv[1]
    if len(sys.argv) >= 3:
        out_fname = sys.argv[2]
    else:
        base = os.path.splitext(in_fname)[0]
        out_fname = base + '_TB.dat'

    t_start = time.time()

    print(f"==> Reading LAMMPS file: {in_fname}", flush=True)
    atoms, Lx, Ly, xy, masses, type_to_species = parse_lammps(in_fname)
    species_counts = {}
    for a in atoms:
        species_counts[a['species']] = species_counts.get(a['species'], 0) + 1
    print(f"    {len(atoms)} atoms total: " +
          ", ".join(f"{n} {s}" for s, n in sorted(species_counts.items())),
          flush=True)

    print(f"==> Computing geometry (moiré cell, lattice parameter, layer thetas)...", flush=True)
    a1_M, a2_M, a1_mono, a2_mono, a_mono, thetas, layers, nat_per_layer = \
        compute_geometry(atoms, Lx, Ly, xy)

    print(f"==> Writing TwisTB geometry: {out_fname}", flush=True)
    write_tb_geometry(atoms, a1_M, a2_M, a1_mono, a2_mono, a_mono,
                      thetas, layers, nat_per_layer, out_fname)

    print(f"==> Done in {time.time() - t_start:.1f}s", flush=True)
    return out_fname


if __name__ == '__main__':
    main()
