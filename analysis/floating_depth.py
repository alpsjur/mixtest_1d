#!/usr/bin/env python3
# analysis/floating_depth.py
"""
Sensitivity analysis for floating (non-bottom-fixed) structures:
how does the structured-zone thickness (structure.depth_zero_below,
here tracked via the synthetic sweep parameter structure.depth_frac =
depth_zero_below/H0) affect mixing, restricted to the bfrc_cb == CD
(no-shear) case?

Primary metric: tau_x (default x_frac=0.10, "time to 10% mixing"), not
the full-completion t_star_mix/tau_mix used for bottom-fixed structures
-- see notes/floating_structure_sensitivity_analysis.md sec 2 for why.

Two companion sweeps are analyzed here:
  - templates/floating_depth_sweep.yaml (main sweep, NTIMES sized off
    tau_x_theory): tau_x sensitivity to depth_frac, c4, H0/zt.
  - templates/floating_depth_longrun_sweep.yaml (long-run subset, NTIMES
    sized off tau_mix_theory): phi(t) plateau / asymptotic mixed-fraction
    diagnostic vs depth_frac.

Usage:
    python analysis/floating_depth.py --sweep sweeps/floating_depth/manifest.yaml
    python analysis/floating_depth.py --longrun sweeps/floating_depth_longrun/manifest.yaml
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


def summarize_floating_sweep(manifest_path: str, x_frac: float = 0.10,
                              plateau_frac: float = 0.5, length_scale: str = "H0",
                              tail_frac: float = 0.1) -> list:
    """
    Compute mixing_timescale() (with the tau_x/plateau diagnostics) for
    every completed run in a floating-depth sweep manifest.

    Returns a list of dicts, one per run, with the swept parameters
    (structure.depth_frac, structure.c4, grid.H0, initial.temp_zt) plus
    tau_x_theory/tau_x_diagnostic (days), tau_mix_theory/tau_mix_diagnostic
    (days, reference scale only -- may not be reached), and the plateau
    diagnostics (plateaued, mixed_fraction_inf).
    """
    manifest = load_yaml(manifest_path)
    rows = []
    for r in manifest.get("runs", []):
        resolved_config = r["resolved_config"]
        if not os.path.isfile(resolved_config):
            continue
        try:
            result = mixing_timescale(
                resolved_config, plateau_frac=plateau_frac,
                length_scale=length_scale, x_frac=x_frac, tail_frac=tail_frac,
            )
        except Exception as e:
            print(f"Skipping {r.get('run_name')}: {e}")
            continue
        p = r.get("params", {})
        rows.append({
            "run_name": r.get("run_name"),
            "depth_frac": p.get("structure.depth_frac"),
            "c4": p.get("structure.c4"),
            "H0": p.get("grid.H0"),
            "temp_zt": p.get("initial.temp_zt"),
            "tau_x_theory_days": result["tau_x_theory"] / 86400.0,
            "tau_x_diagnostic_days": (
                result["tau_x_diagnostic"] / 86400.0
                if result["tau_x_diagnostic"] is not None and not np.isnan(result["tau_x_diagnostic"])
                else np.nan
            ),
            "tau_mix_theory_days": result["tau_mix_theory"] / 86400.0,
            "tau_mix_diagnostic_days": result["tau_mix_diagnostic"] / 86400.0,
            "plateaued": result["plateaued"],
            "mixed_fraction_inf": result["mixed_fraction_inf"],
        })
    return rows


def plot_tau_x_vs_depth_frac(manifest_path: str, ax=None, x_frac: float = 0.10,
                              plateau_frac: float = 0.5, length_scale: str = "H0"):
    """
    Plot diagnostic tau_x (days) vs structure.depth_frac, grouped by
    (structure.c4, grid.H0), with the theory prediction (tau_x_theory,
    which scales as 1/depth_frac via d_struct = depth_frac*H0 -- see
    utils.analytic_mixing_timescale) overlaid as dashed lines for
    comparison. This is the key plot for testing the P_str ~ d_struct
    physics hypothesis (notes/floating_structure_sensitivity_analysis.md).
    """
    import matplotlib.pyplot as plt
    from collections import defaultdict

    rows = summarize_floating_sweep(manifest_path, x_frac=x_frac,
                                     plateau_frac=plateau_frac, length_scale=length_scale)

    by_group = defaultdict(list)
    for row in rows:
        key = (row["c4"], row["H0"], row["temp_zt"])
        by_group[key].append(row)

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))

    cmap = plt.get_cmap("viridis")
    groups = sorted(by_group.keys())
    for i, key in enumerate(groups):
        c4, H0, zt = key
        group_rows = sorted(by_group[key], key=lambda r: r["depth_frac"])
        depth_fracs = [r["depth_frac"] for r in group_rows]
        tau_x_diag = [r["tau_x_diagnostic_days"] for r in group_rows]
        tau_x_th = [r["tau_x_theory_days"] for r in group_rows]
        color = cmap(i / max(len(groups) - 1, 1))
        ax.plot(depth_fracs, tau_x_diag, marker="o", color=color,
                label=f"c4={c4}, H0={H0}, diagnostic")
        ax.plot(depth_fracs, tau_x_th, marker="", linestyle="--", color=color,
                label=f"c4={c4}, H0={H0}, theory")

    ax.set_xlabel(r"depth_frac = depth_zero_below / H0")
    ax.set_ylabel(rf"$\tau_{{x={x_frac}}}$ (days)")
    ax.set_title(r"Time to $x$-fraction mixing vs structured-zone thickness")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="best")
    return ax


def summarize_plateau(manifest_path: str, plateau_frac: float = 0.5,
                       length_scale: str = "H0", tail_frac: float = 0.1) -> list:
    """
    From the long-run subset (templates/floating_depth_longrun_sweep.yaml),
    report the plateau diagnostic (mixed_fraction_inf, phi_inf) vs
    structure.depth_frac -- i.e. how much of the water column's initial
    stratification can a shallow structured zone ever mix, regardless of
    run duration.
    """
    rows = summarize_floating_sweep(manifest_path, plateau_frac=plateau_frac,
                                     length_scale=length_scale, tail_frac=tail_frac)
    return sorted(rows, key=lambda r: r["depth_frac"])


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep", type=str, help="Path to floating_depth manifest.yaml")
    parser.add_argument("--longrun", type=str, help="Path to floating_depth_longrun manifest.yaml")
    parser.add_argument("--x-frac", type=float, default=0.10)
    parser.add_argument("--plateau-frac", type=float, default=0.5)
    parser.add_argument("--length-scale", choices=["H0", "pycnocline"], default="H0")
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()

    if args.sweep:
        rows = summarize_floating_sweep(args.sweep, x_frac=args.x_frac,
                                         plateau_frac=args.plateau_frac,
                                         length_scale=args.length_scale)
        print(f"{'run_name':45s} {'depth_frac':>10s} {'c4':>6s} {'H0':>6s} "
              f"{'tau_x_th(d)':>12s} {'tau_x_diag(d)':>14s} {'mixed_inf':>10s}")
        for r in rows:
            print(f"{r['run_name']:45s} {r['depth_frac']:>10} {r['c4']:>6} {r['H0']:>6} "
                  f"{r['tau_x_theory_days']:>12.3f} {r['tau_x_diagnostic_days']:>14.3f} "
                  f"{r['mixed_fraction_inf']:>10.3f}")

        ax = plot_tau_x_vs_depth_frac(args.sweep, x_frac=args.x_frac,
                                       plateau_frac=args.plateau_frac,
                                       length_scale=args.length_scale)
        import matplotlib.pyplot as plt
        if args.save:
            filename = f"figures/floating_depth_tau_x_{args.length_scale}.png"
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            plt.savefig(filename)
            print(f"Plot saved to {filename}")
        else:
            plt.show()

    if args.longrun:
        rows = summarize_plateau(args.longrun, plateau_frac=args.plateau_frac,
                                  length_scale=args.length_scale)
        print(f"\n{'run_name':45s} {'depth_frac':>10s} {'plateaued':>10s} {'mixed_fraction_inf':>18s}")
        for r in rows:
            print(f"{r['run_name']:45s} {r['depth_frac']:>10} {str(r['plateaued']):>10s} "
                  f"{r['mixed_fraction_inf']:>18.3f}")

    if not args.sweep and not args.longrun:
        parser.error("Provide --sweep and/or --longrun manifest.yaml")


if __name__ == "__main__":
    main()
