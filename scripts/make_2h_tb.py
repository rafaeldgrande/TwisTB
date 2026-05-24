#!/usr/bin/env python3
"""
Generate a TwisTB geometry file for the 2H-stacked bilayer WSe2
from a Quantum ESPRESSO input file.
"""
import numpy as np

# ── QE cell parameters (Angstrom) ─────────────────────────────────────────────
a1_QE = np.array([2.81229396100000, 1.62367147400000, 0.0])
a2_QE = np.array([2.81229396100000, -1.62367147400000, 0.0])

# ── QE atomic positions (Angstrom) ────────────────────────────────────────────
atoms_qe = [
    ('W',  -0.0132720709,  0.0000000000, -0.7956758462),
    ('Se',  1.8615935375,  0.0000000000,  0.8710863829),
    ('Se',  1.8615808294,  0.0000000000, -2.4587612713),
    ('W',   1.8628291270,  0.0000000000,  5.5956647853),
    ('Se', -0.0120237744,  0.0000000000,  7.2587407104),
    ('Se', -0.0120359370,  0.0000000000,  3.9289452389),
]

# ── derived quantities ─────────────────────────────────────────────────────────
a_mono = float(np.linalg.norm(a1_QE[:2]))
theta  = float(np.degrees(np.arctan2(a1_QE[1], a1_QE[0])) % 60.0)

print(f"a_mono = {a_mono:.6f} Å")
print(f"theta  = {theta:.6f}°  (both layers, no twist)")

# monolayer unit cell: TwisTB convention 2 (a1 along x)
a1_mono = np.array([a_mono, 0.0, 0.0])
a2_mono = np.array([a_mono / 2.0, a_mono * np.sqrt(3) / 2.0, 0.0])
a3      = np.array([0.0, 0.0, 40.0])

# moiré cell = primitive cell; use actual QE vectors
a1_M = np.array([a1_QE[0], a1_QE[1], 0.0])
a2_M = np.array([a2_QE[0], a2_QE[1], 0.0])

# ── assign layers by z (split midway between the two W atoms) ─────────────────
z_W = sorted(a[3] for a in atoms_qe if a[0] == 'W')
z_split = (z_W[0] + z_W[1]) / 2.0   # midpoint between z≈-0.80 and z≈5.60
print(f"Layer split at z = {z_split:.3f} Å  (W1 z={z_W[0]:.3f}, W2 z={z_W[1]:.3f})")

atoms = [(sp, 1 if z < z_split else 2, x, y, z) for (sp, x, y, z) in atoms_qe]

# ── write TwisTB geometry file ─────────────────────────────────────────────────
out = "bilayer_2H_WSe2_TB.dat"
nlayer = 2
nat_per_layer = 3
thetas = [theta, theta]

with open(out, 'w') as f:
    f.write('  '.join(str(nat_per_layer) for _ in range(nlayer)) + '\n')
    f.write('  '.join(f'{a_mono:.6f}' for _ in range(nlayer)) + '\n')
    f.write('  '.join(f'{t:.6f}' for t in thetas) + '\n')
    f.write('  '.join('1.00' for _ in range(nlayer)) + '\n')
    for v in [a1_mono, a2_mono, a3]:
        f.write(f'{v[0]:.8f}   {v[1]:.8f}   {v[2]:.8f}\n')
    for v in [a1_M, a2_M, a3]:
        f.write(f'{v[0]:.8f}   {v[1]:.8f}   {v[2]:.8f}\n')
    for lyr in [1, 2]:
        for sp, layer, x, y, z in atoms:
            if layer == lyr:
                f.write(f"{sp} {lyr} {x:.6f} {y:.6f} {z:.6f}\n")

print(f"Written: {out}")
