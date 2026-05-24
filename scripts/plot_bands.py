#!/usr/bin/env python3
"""
Plot TwisTB band structure from a *_BS.dat file.

Usage:
    python3 plot_bands.py <BS.dat> [--ef EF] [--emin EMIN] [--emax EMAX]
                                   [--out OUTPUT.png]

Energy reference:
    By default, Ef = max energy of the last valence band (band 7*x/11,
    where x = total number of bands), so E - Ef = 0 is the VBM.
    Pass --ef to override with an absolute value (eV, raw TB scale).

K-path labels are read from the filename convention used by TwisTB
(Gamma-M-K-Gamma of the mini-BZ).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
import argparse
import os
import sys

STYLE = "/Users/rdelgrande/work/Codes/utilities/presentation.mplstyle"

def parse_bs(fname):
    """Return (k_coords, energies) where energies.shape = (nbands, nkpts)."""
    bands, block = [], []
    with open(fname) as f:
        for line in f:
            s = line.strip()
            if s:
                k, e = map(float, s.split())
                block.append((k, e))
            else:
                if block:
                    bands.append(block)
                    block = []
    if block:
        bands.append(block)

    k_coords = np.array([pt[0] for pt in bands[0]])
    energies  = np.array([[pt[1] for pt in b] for b in bands])
    return k_coords, energies


def kpath_ticks(k_coords):
    """
    TwisTB writes Gamma-M-K-Gamma for the mini-BZ.
    High-symmetry points are at equal thirds of the total k-path length.
    """
    kmin, kmax = k_coords[0], k_coords[-1]
    span = kmax - kmin
    ticks  = [kmin, kmin + span/3, kmin + 2*span/3, kmax]
    labels = [r'$\Gamma$', r'$M$', r'$K$', r'$\Gamma$']
    return ticks, labels


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('bs_file', help='*_BS.dat file from TwisTB')
    ap.add_argument('--ef',   type=float, default=None,
                    help='Fermi level in raw TB energy units (eV). '
                         'Default: max of band 7*nbands/11 (last valence band = VBM).')
    ap.add_argument('--emin', type=float, default=-3.0,
                    help='Lower energy limit relative to Ef (default: -3)')
    ap.add_argument('--emax', type=float, default=0.5,
                    help='Upper energy limit relative to Ef (default: +0.5)')
    ap.add_argument('--out',  default=None,
                    help='Output image filename (default: <bs_file_stem>_plot.png)')
    ap.add_argument('--noccs', type=int, default=None,
                    help='Number of occupied (valence) bands in the output. '
                         'Use when reduced_workspace=true: set to num_eigs value. '
                         'Default: 7*nbands/11 (full-spectrum formula).')
    args = ap.parse_args()

    k_coords, energies = parse_bs(args.bs_file)

    nbands = energies.shape[0]
    if args.noccs is not None:
        noccs = args.noccs
        print(f"Bands: {nbands} total = {noccs} valence + {nbands - noccs} conduction "
              f"(reduced_workspace mode, num_eigs={noccs})")
    else:
        # Derive band counts from the 11-orbital Fang model:
        #   total bands x = 11 * n_unitcells  →  valence = 7*x/11, conduction = 4*x/11
        if nbands % 11 != 0:
            print(f"WARNING: number of bands ({nbands}) is not a multiple of 11 — "
                  "noccs estimate may be wrong")
        noccs = 7 * nbands // 11
        print(f"Bands: {nbands} total = {noccs} valence + {nbands - noccs} conduction "
              f"({nbands // 11} formula units)")

    # Fermi level: default = max of the last valence band
    if args.ef is None:
        ef = float(energies[noccs - 1].max())
        print(f"Ef = {ef:.4f} eV  (max of band {noccs}, the last valence band)")
    else:
        ef = args.ef
        print(f"Ef = {ef:.4f} eV  (user-supplied)")

    energies_shifted = energies - ef

    # Output filename
    if args.out is None:
        stem = os.path.splitext(args.bs_file)[0]
        out_fname = stem + '_plot.png'
    else:
        out_fname = args.out

    # ── plot ─────────────────────────────────────────────────────────────────
    if os.path.isfile(STYLE):
        plt.style.use(STYLE)

    fig, ax = plt.subplots(figsize=(5, 5))

    for band in energies_shifted:
        ax.plot(k_coords, band, color='tab:red', lw=0.6)

    # Fermi level line
    ax.axhline(0, color='white', lw=0.8, ls='--', alpha=0.7)

    # High-symmetry ticks and vertical lines
    ticks, labels = kpath_ticks(k_coords)
    for t in ticks[1:-1]:
        ax.axvline(t, color='gray', lw=0.5, ls='--', alpha=0.5)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)

    ax.set_xlim(k_coords[0], k_coords[-1])
    ax.set_ylim(args.emin, args.emax)
    ax.set_ylabel(r'$E - E_{\rm F}$ (eV)')

    # Title: extract theta from filename
    base = os.path.basename(args.bs_file)
    ax.set_title(base.replace('_BS.dat', '').replace('_', ' '), fontsize=8)

    plt.tight_layout()
    fig.savefig(out_fname, dpi=150)
    print(f"Saved: {out_fname}")


if __name__ == '__main__':
    main()
