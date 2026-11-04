#!/usr/bin/env python3
"""Plot GridFire vs MESA (bbq) for the tables written by mesa_comparison_test.

    plot_comparison.py [--dir DIR] [--out DIR] [--species H-1 He-4 ...]

Reads mesa_comparison_<case>.csv from DIR (default: current directory, where the test writes them) and writes
mesa_comparison_<case>.png to OUT (default: DIR). Requires numpy and matplotlib.
"""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    with open(path) as f:
        rows = list(csv.reader(line for line in f if not line.startswith("#")))
    header, data = rows[0], np.array(rows[1:], dtype=float)
    return {name: data[:, i] for i, name in enumerate(header)}


def plot_case(path, out_dir, wanted):
    case = path.stem.removeprefix("mesa_comparison_")
    t = load(path)
    age = t["age"]
    species = [c.removeprefix("X_mesa_") for c in t if c.startswith("X_mesa_")]
    if wanted:
        species = [s for s in species if s in wanted]
    else:  # The species that matter: those reaching 1e-4 in either code
        species = [s for s in species if max(t[f"X_mesa_{s}"].max(), t[f"X_gridfire_{s}"].max()) >= 1e-4]

    fig, (ax_x, ax_xerr, ax_e) = plt.subplots(3, 1, figsize=(8, 10), sharex=True,
                                               gridspec_kw={"height_ratios": [3, 1.5, 1.5]})
    colors = plt.cm.tab20(np.linspace(0, 1, max(len(species), 1)))
    for s, color in zip(species, colors):
        xm, xg = t[f"X_mesa_{s}"], t[f"X_gridfire_{s}"]
        ax_x.loglog(age, xm, color=color, lw=2, alpha=0.5)
        ax_x.loglog(age, xg, color=color, lw=1, ls="--", label=s)
        mask = np.maximum(xm, xg) >= 1e-4
        rel = np.abs(xg - xm) / np.maximum(np.abs(xm), 1e-300)
        ax_xerr.loglog(age[mask], np.maximum(rel[mask], 1e-12), color=color, lw=1)

    em, eg = t["energy_mesa"], t["energy_gridfire"]
    ax_e.loglog(age, np.maximum(np.abs(eg - em) / np.maximum(np.abs(em), 1e-300), 1e-12), color="k", lw=1)

    ax_x.set_ylabel("mass fraction")
    ax_x.set_ylim(1e-8, 1.5)
    ax_x.legend(ncol=3, fontsize=7, title="solid: MESA, dashed: GridFire", title_fontsize=8)
    ax_xerr.set_ylabel("|ΔX| / X  (X ≥ 1e-4)")
    ax_e.set_ylabel("|ΔE_step| / E_step")
    ax_e.set_xlabel("burn time [s]")
    for ax in (ax_xerr, ax_e):
        ax.axhline(0.05, color="r", lw=0.8, ls=":")
    fig.suptitle(f"GridFire vs MESA (bbq): {case}")
    fig.tight_layout()
    out = out_dir / f"mesa_comparison_{case}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", type=Path, default=Path("."))
    parser.add_argument("--out", type=Path)
    parser.add_argument("--species", nargs="*")
    args = parser.parse_args()
    out_dir = args.out or args.dir
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = sorted(args.dir.glob("mesa_comparison_*.csv"))
    if not paths:
        raise SystemExit(f"no mesa_comparison_*.csv in {args.dir}; run mesa_comparison_test first")
    for path in paths:
        plot_case(path, out_dir, args.species)


if __name__ == "__main__":
    main()
