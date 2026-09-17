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


def pycnocline_edge(params: dict, k: float = 1.0) -> float:
    """
    Approximate depth of the pycnocline's upper edge: z_t - k*temp_ht
    (the tanh thermocline profile's transition is ~complete within
    +-temp_ht of its centre z_t; k=1.0 is a reasonable first estimate of
    where the transition "begins", not an exact physical threshold --
    see notes/floating_structure_sensitivity_analysis.md sec 4.3).
    """
    zt = float(params["initial"]["temp_zt"])
    ht = float(params["initial"]["temp_ht"])
    return zt - k * ht


def summarize_floating_sweep(manifest_path, x_frac: float = 0.10,
                              plateau_frac: float = 0.5, length_scale: str = "H0",
                              tail_frac: float = 0.1, k_edge: float = 1.0) -> list:
    """
    Compute mixing_timescale() (with the tau_x/plateau diagnostics) for
    every completed run in a floating-depth sweep manifest (or a list of
    manifest paths, which are concatenated).

    Returns a list of dicts, one per run, with the swept parameters
    (structure.depth_frac, structure.c4, grid.H0, initial.temp_zt) plus
    tau_x_theory/tau_x_diagnostic (days), tau_mix_theory/tau_mix_diagnostic
    (days, reference scale only -- may not be reached), the plateau
    diagnostics (plateaued, mixed_fraction_inf), and the pycnocline-
    penetration diagnostics d_struct, pycnocline_edge, penetration_margin,
    penetrating (see pycnocline_edge()).
    """
    manifest_paths = [manifest_path] if isinstance(manifest_path, str) else list(manifest_path)
    rows = []
    for mp in manifest_paths:
        manifest = load_yaml(mp)
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
            edge = pycnocline_edge(result["params"], k=k_edge)
            penetration_margin = result["d_struct"] - edge
            rows.append({
                "run_name": r.get("run_name"),
                "depth_frac": p.get("structure.depth_frac"),
                "c4": p.get("structure.c4"),
                "H0": p.get("grid.H0"),
                "temp_zt": p.get("initial.temp_zt"),
                "d_struct": result["d_struct"],
                "pycnocline_edge": edge,
                "penetration_margin": penetration_margin,
                "penetrating": bool(penetration_margin >= 0),
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


def plot_tau_x_ratio_vs_penetration(manifest_paths, ax=None, x_frac: float = 0.10,
                                     plateau_frac: float = 0.5, k_edge: float = 1.0):
    """
    Plot the theory/diagnostic tau_x ratio against penetration_margin
    (d_struct - pycnocline_edge), across one or more sweep manifests,
    colored by grid.H0. This is the key diagnostic for testing whether
    "does the structured zone reach the pycnocline" (rather than
    depth_frac alone) controls the theory's over-/under-prediction
    crossover -- see
    notes/floating_structure_sensitivity_analysis.md sec 4.3.
    """
    import matplotlib.pyplot as plt

    rows = summarize_floating_sweep(manifest_paths, x_frac=x_frac,
                                     plateau_frac=plateau_frac, k_edge=k_edge)
    rows = [r for r in rows if not np.isnan(r["tau_x_diagnostic_days"])]

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))

    H0_vals = sorted({r["H0"] for r in rows})
    cmap = plt.get_cmap("viridis")
    for i, H0 in enumerate(H0_vals):
        group = sorted([r for r in rows if r["H0"] == H0], key=lambda r: r["penetration_margin"])
        margins = [r["penetration_margin"] for r in group]
        ratios = [r["tau_x_theory_days"] / r["tau_x_diagnostic_days"] for r in group]
        ax.plot(margins, ratios, marker="o", color=cmap(i / max(len(H0_vals) - 1, 1)),
                label=f"H0={H0}")

    ax.axvline(0.0, color="k", linestyle="--", linewidth=1, label="pycnocline edge")
    ax.axhline(1.0, color="k", linestyle=":", linewidth=1)
    ax.set_xlabel(r"penetration margin $= d_{struct} - (z_t - k\,h_t)$ (m)")
    ax.set_ylabel(r"$\tau_{x}$ theory / diagnostic ratio")
    ax.set_title("Theory/diagnostic ratio vs. pycnocline penetration")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    return ax


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


def plot_phi_longrun(manifest_path: str, ax=None, plateau_frac: float = 0.5,
                      length_scale: str = "H0"):
    """
    Plot phi_star(t) (in days, not the dimensionless t_star) for every run
    in the long-run/plateau-detection subset, one curve per
    structure.depth_frac. Illustrates whether -- and how slowly -- phi(t)
    approaches zero once the structured zone only occupies part of the
    water column (see notes/floating_structure_sensitivity_analysis.md
    sec 4.2): thinner structured zones take much longer to reach a given
    mixed fraction, but (at least for depth_frac >= 0.25) do eventually
    reach ~full mixing rather than plateauing at a permanent floor.
    """
    import matplotlib.pyplot as plt

    manifest = load_yaml(manifest_path)
    runs = sorted(manifest.get("runs", []),
                  key=lambda r: r.get("params", {}).get("structure.depth_frac", 0))

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))

    cmap = plt.get_cmap("viridis")
    n_plotted = 0
    n_total = len(runs)
    for r in runs:
        resolved_config = r["resolved_config"]
        if not os.path.isfile(resolved_config):
            continue
        try:
            result = mixing_timescale(resolved_config, plateau_frac=plateau_frac,
                                       length_scale=length_scale)
        except Exception as e:
            print(f"Skipping {r.get('run_name')}: {e}")
            continue
        depth_frac = r.get("params", {}).get("structure.depth_frac")
        ax.plot(
            result["days"], result["phi_star"],
            label=f"depth_frac={depth_frac}",
            color=cmap(n_plotted / max(n_total - 1, 1)),
        )
        n_plotted += 1

    ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel(r"$\phi(t)/\phi(0)$")
    ax.set_title("Long-run mixing: does $\\phi(t)$ plateau at a nonzero residual?")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    return ax


def gather_penetrating_curves(manifest_paths, length_scale: str = "H0",
                               x_frac: float = 0.10, plateau_frac: float = 0.5,
                               k_edge: float = 1.0):
    """
    Collect (run_name, H0, t_star, phi_star) for every completed run in
    one or more sweep manifests that is classified as "penetrating"
    (structured zone reaches into the pycnocline, per pycnocline_edge()),
    computing t_star with the given length_scale ("H0", "pycnocline", or
    "d_struct" -- see mixing_timescale()). Used for phase-A collapse
    exploration (notes/floating_structure_sensitivity_analysis.md sec
    4.4): only pre-threshold (t_star <= 1) portions are trustworthy for
    the "how far did it get" curve shape, but the full curve is returned
    here for inspection.
    """
    manifest_paths = [manifest_paths] if isinstance(manifest_paths, str) else list(manifest_paths)
    curves = []
    for mp in manifest_paths:
        manifest = load_yaml(mp)
        for r in manifest.get("runs", []):
            resolved_config = r["resolved_config"]
            if not os.path.isfile(resolved_config):
                continue
            try:
                result = mixing_timescale(resolved_config, plateau_frac=plateau_frac,
                                           length_scale=length_scale, x_frac=x_frac)
            except Exception as e:
                print(f"Skipping {r.get('run_name')}: {e}")
                continue
            edge = pycnocline_edge(result["params"], k=k_edge)
            if result["d_struct"] - edge < 0:
                continue  # non-penetrating -- deferred to phase B
            p = r.get("params", {})
            curves.append({
                "run_name": r.get("run_name"),
                "H0": p.get("grid.H0"),
                "depth_frac": p.get("structure.depth_frac"),
                "t_star": result["t_star"],
                "phi_star": result["phi_star"],
            })
    return curves


def collapse_spread_score(curves, t_star_max: float = 1.0, n_grid: int = 50) -> float:
    """
    Quantitative collapse-quality metric: interpolate all curves'
    phi_star(t_star) onto a common grid over [0, t_star_max] (the
    pre-threshold/early-time portion where curves overlap for essentially
    all runs), then compute the standard deviation of phi_star across
    curves at each grid point, and return the mean of that std-dev over
    the grid (lower = better collapse). Grid points where fewer than 2
    curves have data are skipped.

    A perfect collapse (all curves identical) gives a score of 0.
    """
    if not curves:
        return np.nan
    grid = np.linspace(0.0, t_star_max, n_grid)
    values_per_point = [[] for _ in grid]
    for c in curves:
        t_star = c["t_star"]
        phi_star = c["phi_star"]
        order = np.argsort(t_star)
        t_sorted = t_star[order]
        phi_sorted = phi_star[order]
        t_min, t_max = t_sorted[0], t_sorted[-1]
        for i, tg in enumerate(grid):
            if t_min <= tg <= t_max:
                values_per_point[i].append(float(np.interp(tg, t_sorted, phi_sorted)))
    stds = [np.std(v) for v in values_per_point if len(v) >= 2]
    return float(np.mean(stds)) if stds else np.nan


def plot_penetrating_collapse(manifest_paths, ax=None, length_scale: str = "H0",
                               x_frac: float = 0.10, plateau_frac: float = 0.5,
                               k_edge: float = 1.0, t_star_max: float = 1.2):
    """
    Plot phi_star(t_star) for only the "penetrating" runs (see
    gather_penetrating_curves), colored by grid.H0, for a given
    length_scale candidate. Prints the collapse_spread_score alongside.
    Phase-A exploration -- see
    notes/floating_structure_sensitivity_analysis.md sec 4.4.
    """
    import matplotlib.pyplot as plt

    curves = gather_penetrating_curves(manifest_paths, length_scale=length_scale,
                                        x_frac=x_frac, plateau_frac=plateau_frac,
                                        k_edge=k_edge)
    score = collapse_spread_score(curves, t_star_max=min(t_star_max, 1.0))

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))

    H0_vals = sorted({c["H0"] for c in curves})
    cmap = plt.get_cmap("viridis")
    for c in curves:
        i = H0_vals.index(c["H0"])
        ax.plot(c["t_star"], c["phi_star"], color=cmap(i / max(len(H0_vals) - 1, 1)),
                alpha=0.7, linewidth=1)

    ax.axvline(1.0, color="k", linestyle="--", linewidth=1)
    ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
    ax.set_xlim(0, t_star_max)
    ax.set_xlabel(rf"$t^*$ (length_scale={length_scale})")
    ax.set_ylabel(r"$\phi(t)/\phi(0)$")
    ax.set_title(f"Penetrating-only collapse (L={length_scale}), score={score:.4f}")
    ax.grid(True, alpha=0.3)
    # Manual legend by H0 color only (many overlapping runs otherwise).
    handles = [plt.Line2D([0], [0], color=cmap(i / max(len(H0_vals) - 1, 1)), label=f"H0={H0}")
               for i, H0 in enumerate(H0_vals)]
    ax.legend(handles=handles, fontsize=8, loc="best")
    return ax, score


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep", type=str, nargs="+", help="Path(s) to floating_depth manifest.yaml")
    parser.add_argument("--longrun", type=str, help="Path to floating_depth_longrun manifest.yaml")
    parser.add_argument("--penetration", type=str, nargs="+",
                         help="Path(s) to sweep manifest.yaml(s) to combine for the "
                              "penetration-margin diagnostic and/or phase-A collapse "
                              "exploration (e.g. --sweep + the two straddle sweeps)")
    parser.add_argument("--collapse", action="store_true",
                         help="Run phase-A penetrating-only collapse exploration "
                              "(requires --penetration)")
    parser.add_argument("--x-frac", type=float, default=0.10)
    parser.add_argument("--plateau-frac", type=float, default=0.5)
    parser.add_argument("--length-scale", choices=["H0", "pycnocline", "d_struct"], default="H0")
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()

    if args.sweep:
        sweep_paths = args.sweep if len(args.sweep) > 1 else args.sweep[0]
        rows = summarize_floating_sweep(sweep_paths, x_frac=args.x_frac,
                                         plateau_frac=args.plateau_frac,
                                         length_scale=args.length_scale)
        print(f"{'run_name':45s} {'depth_frac':>10s} {'c4':>6s} {'H0':>6s} "
              f"{'tau_x_th(d)':>12s} {'tau_x_diag(d)':>14s} {'mixed_inf':>10s}")
        for r in rows:
            print(f"{r['run_name']:45s} {r['depth_frac']:>10} {r['c4']:>6} {r['H0']:>6} "
                  f"{r['tau_x_theory_days']:>12.3f} {r['tau_x_diagnostic_days']:>14.3f} "
                  f"{r['mixed_fraction_inf']:>10.3f}")

        if len(args.sweep) == 1:
            ax = plot_tau_x_vs_depth_frac(args.sweep[0], x_frac=args.x_frac,
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

        ax3 = plot_phi_longrun(args.longrun, plateau_frac=args.plateau_frac,
                                length_scale=args.length_scale)
        import matplotlib.pyplot as plt
        if args.save:
            filename3 = f"figures/floating_depth_longrun_phi_{args.length_scale}.png"
            os.makedirs(os.path.dirname(filename3), exist_ok=True)
            plt.savefig(filename3)
            print(f"Plot saved to {filename3}")
        else:
            plt.show()

    if args.penetration:
        import matplotlib.pyplot as plt

        ax4 = plot_tau_x_ratio_vs_penetration(args.penetration, x_frac=args.x_frac,
                                               plateau_frac=args.plateau_frac)
        if args.save:
            filename4 = "figures/floating_depth_penetration_margin.png"
            plt.savefig(filename4)
            print(f"Plot saved to {filename4}")
        else:
            plt.show()

        if args.collapse:
            for ls in ["H0", "pycnocline", "d_struct"]:
                ax5, score = plot_penetrating_collapse(args.penetration, length_scale=ls,
                                                         x_frac=args.x_frac,
                                                         plateau_frac=args.plateau_frac)
                print(f"length_scale={ls}: collapse_spread_score={score:.4f}")
                if args.save:
                    filename5 = f"figures/floating_depth_penetrating_collapse_{ls}.png"
                    plt.savefig(filename5)
                    print(f"Plot saved to {filename5}")
                else:
                    plt.show()

    if not args.sweep and not args.longrun and not args.penetration:
        parser.error("Provide --sweep, --longrun, and/or --penetration manifest.yaml")


if __name__ == "__main__":
    main()

