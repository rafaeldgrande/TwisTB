#!/usr/bin/env python3
"""
Convert a TwisTB ASG geometry file (old 9-line format) to the
TwisTB v2.0 geometry file format (10-line format).

ASG format (9 header lines):
  nat1  nat2
  theta_twist          ← single value
  strain               ← single value
  a1_mono x y z
  a2_mono x y z
  a3 x y z
  a1_moire x y z
  a2_moire x y z
  a3 x y z
  atoms...

TwisTB v2.0 format (10 header lines):
  nat1  nat2
  a_mono1  a_mono2     ← bondlength per layer
  theta1  theta2       ← per-layer orientation angle folded to [0,60)
  strain1  strain2
  a1_mono x y z
  ...

Usage: python3 convert_asg_to_tb.py input.dat [output.dat]
"""
import numpy as np
import sys, os

def _nn_vector(x0, y0, metals, Lx, Ly, a_mono):
    best_d, best_v = np.inf, None
    for mx, my in metals:
        dx0, dy0 = mx - x0, my - y0
        for nx in range(-3, 4):
            for ny in range(-3, 4):
                dx = dx0 + nx * Lx
                dy = dy0 + ny * Ly
                d = np.hypot(dx, dy)
                if 0.5 * a_mono < d < 1.5 * a_mono and d < best_d:
                    best_d, best_v = d, (dx, dy)
    return best_v

def layer_theta(metals_xy, Lx, Ly, a_mono):
    angles = []
    for i, (x0, y0) in enumerate(metals_xy):
        others = [p for j, p in enumerate(metals_xy) if j != i]
        v = _nn_vector(x0, y0, others, Lx, Ly, a_mono)
        if v is not None:
            angles.append(np.degrees(np.arctan2(v[1], v[0])) % 60.0)
    if not angles:
        raise RuntimeError("No NN metal-metal bond found")
    return float(np.median(angles))

def convert(in_fname, out_fname):
    with open(in_fname) as f:
        lines = f.readlines()

    # ── parse header ──────────────────────────────────────────────────────────
    idx = 0
    nat = list(map(int, lines[idx].split())); idx += 1
    nlayer = len(nat)

    theta_twist = float(lines[idx].strip()); idx += 1
    strain_val  = float(lines[idx].strip()); idx += 1

    a1_mono = np.array(list(map(float, lines[idx].split()))); idx += 1
    a2_mono = np.array(list(map(float, lines[idx].split()))); idx += 1
    a3_mono = np.array(list(map(float, lines[idx].split()))); idx += 1

    a1_M = np.array(list(map(float, lines[idx].split()))); idx += 1
    a2_M = np.array(list(map(float, lines[idx].split()))); idx += 1
    a3_M = np.array(list(map(float, lines[idx].split()))); idx += 1

    # ── parse atoms ───────────────────────────────────────────────────────────
    atoms = []
    for line in lines[idx:]:
        s = line.strip()
        if not s:
            continue
        parts = s.split()
        sp, lyr = parts[0], int(parts[1])
        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
        atoms.append((sp, lyr, x, y, z))

    # ── compute a_mono from cell vector ───────────────────────────────────────
    a_mono_val = float(np.linalg.norm(a1_mono[:2]))

    # ── per-layer theta from ASG rotation convention ──────────────────────────
    # ASG starts from convention-2 (a1 along x), then rotates
    # layer 1 by +theta_twist/2 and layer 2 by -theta_twist/2.
    # NN bond angle folded to [0°, 60°):
    #   layer 1: +theta_twist/2 % 60
    #   layer 2: -theta_twist/2 % 60  =  60 - theta_twist/2  (for 0 < twist < 120)
    half = theta_twist / 2.0
    thetas = [half % 60.0, (60.0 - half) % 60.0]

    # ── write v2.0 format ─────────────────────────────────────────────────────
    with open(out_fname, 'w') as f:
        f.write('  '.join(str(n) for n in nat) + '\n')
        f.write('  '.join(f'{a_mono_val:.6f}' for _ in nat) + '\n')
        f.write('  '.join(f'{t:.6f}' for t in thetas) + '\n')
        f.write('  '.join(f'{strain_val:.2f}' for _ in nat) + '\n')
        for v in [a1_mono, a2_mono, a3_mono]:
            f.write(f'{v[0]:.8f}   {v[1]:.8f}   {v[2]:.8f}\n')
        for v in [a1_M, a2_M, a3_M]:
            f.write(f'{v[0]:.8f}   {v[1]:.8f}   {v[2]:.8f}\n')
        for sp, lyr, x, y, z in atoms:
            f.write(f'{sp} {lyr} {x:.6f} {y:.6f} {z:.6f}\n')

    print(f"  a_mono = {a_mono_val:.4f} Å")
    print(f"  theta (twist) = {theta_twist:.4f}°")
    for i, t in enumerate(thetas):
        print(f"  Layer {i+1}: theta = {t:.4f}°")
    if nlayer == 2:
        print(f"  Twist angle (from thetas) = {abs(thetas[1]-thetas[0]):.4f}°")
    print(f"Written: {out_fname}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    in_f = sys.argv[1]
    out_f = sys.argv[2] if len(sys.argv) >= 3 else os.path.splitext(in_f)[0] + '_v2.dat'
    print(f"Converting {in_f} → {out_f}")
    convert(in_f, out_f)
