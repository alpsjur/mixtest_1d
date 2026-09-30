#!/usr/bin/env python3
# analysis/plot_c2c4_phase.py
"""
Phase-diagram plot of the empirical dimensionless mixing-completion time
t*_mix over the (structure.c4, GLS.C2) plane, for a single sweep manifest
(e.g. sweeps/gen_c2c4/manifest.yaml or sweeps/k-omega_c2c4/manifest.yaml).

Companion to analysis/mixing_timescale.py's plot_sensitivity(), which only
varies c4 at fixed C2. Here both closure coefficients are varied together
to see whether they interact (i.e. whether the sensitivity to c4 depends
on C2, and vice versa) rather than acting as separable, additive factors.

Usage:
    python analysis/plot_c2c4_phase.py sweeps/gen_c2c4/manifest.yaml \
        --title "GEN closure" --save figures/gen_c2c4_phase.png

    # Overlay c4 as a fraction of GLS.C1 (so gen and k-omega share an axis):
    python analysis/plot_c2c4_phase.py sweeps/k-omega_c2c4/manifest.yaml \
        --normalize-c4 --title "k-omega closure"
"""
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

import numpy as np

from utils.utils import load_yaml
from analysis.mixing_timescale import mixing_timescale


def collect_grid(manifest_path: str, plateau_frac: float = 0.5, length_scale: str = "H0",
                  normalize_c4: bool = False):
    """
    Run mixing_timescale() for every completed run in the manifest and
    return (c2_vals, c4_vals, T) where T[i, j] = t_star_mix for
    (c2_vals[i], c4_vals[j]), plus a dict of per-run diagnostics for
    bookkeeping (u_inf, Pstr_diag, NaNs, etc).
    """
    manifest = load_yaml(manifest_path)
    rows = []
    for r in manifest.get("runs", []):
        resolved_config = r["resolved_config"]
        if not os.path.isfile(resolved_config):
            continue
        try:
            result = mixing_timescale(resolved_config, plateau_frac=plateau_frac, length_scale=length_scale)
        except Exception as e:
            print(f"Skipping {r.get('run_name')}: {e}")
            continue
        p = r.get("params", {})
        c2 = float(p.get("GLS.C2"))
        c4 = float(p.get("structure.c4"))
        if normalize_c4:
            c1 = float(result["params"]["GLS"]["C1"])
            c4 = c4 / c1
        rows.append({"c2": c2, "c4": c4, "t_star_mix": result["t_star_mix"],
                      "Pstr_diag": result["Pstr_diag"], "run_name": r.get("run_name")})

    if not rows:
        raise RuntimeError(f"No completed runs found in {manifest_path}")

    c2_vals = sorted({row["c2"] for row in rows})
    c4_vals = sorted({row["c4"] for row in rows})
    T = np.full((len(c2_vals), len(c4_vals)), np.nan)
    for row in rows:
        i = c2_vals.index(row["c2"])
        j = c4_vals.index(row["c4"])
        T[i, j] = row["t_star_mix"]

    return np.array(c2_vals), np.array(c4_vals), T, rows


def plot_phase(manifest_path: str, ax=None, plateau_frac: float = 0.5, length_scale: str = "H0",
                normalize_c4: bool = False, title: str = None, vmax: float = None):
    import matplotlib.pyplot as plt
    try:
        from cmcrameri import cm as cmc
        cmap = cmc.batlow
    except ModuleNotFoundError:
        cmap = plt.get_cmap("viridis")

    c2_vals, c4_vals, T, rows = collect_grid(manifest_path, plateau_frac, length_scale, normalize_c4)

    if ax is None:
        _, ax = plt.subplots(figsize=(6.5, 5.5))

    n_nan = int(np.isnan(T).sum())
    if n_nan:
        print(f"{manifest_path}: {n_nan}/{T.size} combinations did not reach "
              f"phi_star<0.05 within the run (t_star_mix=NaN) -- shown blank.")

    im = ax.pcolormesh(c4_vals, c2_vals, T, cmap=cmap, shading="nearest", vmax=vmax)
    for row in rows:
        if np.isnan(row["t_star_mix"]):
            ax.plot(row["c4"], row["c2"], marker="x", color="red", markersize=8, mew=2)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r"$t^*_{mix}$ (dimensionless)")

    xlabel = r"structure.c4 / GLS.C1" if normalize_c4 else r"structure.c4"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"GLS.C2 ($c_{\psi2}$)")
    ax.set_title(title or os.path.basename(os.path.dirname(manifest_path)))
    return ax, (c2_vals, c4_vals, T)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", help="Path to a sweep manifest.yaml")
    parser.add_argument("--plateau-frac", type=float, default=0.5)
    parser.add_argument("--length-scale", choices=["H0", "pycnocline"], default="H0")
    parser.add_argument("--normalize-c4", action="store_true",
                         help="Plot c4/GLS.C1 on the x-axis instead of raw c4 "
                              "(needed to overlay closures with different C1)")
    parser.add_argument("--title", type=str, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--save", type=str, default=None, help="Output figure path")
    args = parser.parse_args()

    import matplotlib.pyplot as plt
    ax, (c2_vals, c4_vals, T) = plot_phase(
        args.manifest, plateau_frac=args.plateau_frac, length_scale=args.length_scale,
        normalize_c4=args.normalize_c4, title=args.title, vmax=args.vmax,
    )
    plt.tight_layout()
    if args.save:
        os.makedirs(os.path.dirname(args.save), exist_ok=True)
        plt.savefig(args.save, dpi=150)
        print(f"Saved {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
