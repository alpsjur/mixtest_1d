#!/usr/bin/env python3
# analysis/floating_structure_mixing.py
"""
Analysis/plotting for the floating-structure (spar / semisubmersible) draft,
geometry, and seasonality experiments -- see
notes/floating_structure_depth_distribution.md.

Unlike analysis/mixing_timescale.py, this module does NOT use
utils.compute_Pstr / utils.analytic_mixing_timescale for normalization,
because both assume a depth-uniform structure.str_a spanning the full water
column (structure.depth_zero_below >= grid.H0), which is false for every
run analysed here (partial-depth, "floating" structures). Everything below
is computed directly from phi(t) = utils.compute_phi(...) and compared in
real (dimensional) time instead of the dimensionless t_star used elsewhere
in this repo.

Usage (regenerate all three figures from the prepared runs under runs/):
    python analysis/floating_structure_mixing.py
"""
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.utils import open_roms_dataset, compute_phi, compute_time_vector

RUNS_DIR = os.path.join(ROOT_DIR, "runs")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

# draft (m) for every run prepared by tools/prep_structure_geometry_sweep.py
DRAFTS = [15, 20, 30, 40, 50, 65, 90, 120, 150]


def load_phi_star(run_name: str):
    """
    Load phi_star(t) = phi(t)/phi(0) and the matching real-time axis (days)
    for a single run, given its name under runs/.

    Returns
    -------
    t_days : np.ndarray
    phi_star : np.ndarray
    """
    cfg_path = os.path.join(RUNS_DIR, run_name, "resolved_config.yaml")
    ds, grid, params = open_roms_dataset(cfg_path)
    phi, _ = compute_phi(ds, grid, params)
    phi_star = (phi / phi.isel(ocean_time=0)).values
    t_days = compute_time_vector(params)
    ds.close()
    return t_days, phi_star


def crossing_time(t_days: np.ndarray, phi_star: np.ndarray, threshold: float) -> float:
    """
    First time (linearly interpolated) at which phi_star drops to or below
    `threshold`. Returns NaN if the threshold is never reached within the
    run.
    """
    below = np.where(phi_star <= threshold)[0]
    if len(below) == 0:
        return np.nan
    i1 = below[0]
    if i1 == 0:
        return float(t_days[0])
    t0, t1 = t_days[i1 - 1], t_days[i1]
    p0, p1 = phi_star[i1 - 1], phi_star[i1]
    if p0 == p1:
        return float(t1)
    return float(t0 + (threshold - p0) * (t1 - t0) / (p1 - p0))


def plot_draft_sweep(drafts=DRAFTS, out_path=None):
    """
    Figure 1 (Result A): phi_star(t) for every draft in the fixed-total_a
    draft sweep, plus t(phi_star=0.5) / t(phi_star=0.1) vs draft.
    """
    results = {}
    for d in drafts:
        t_days, phi_star = load_phi_star(f"geom_matched_draft{d}")
        results[d] = (t_days, phi_star)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    cmap = plt.cm.viridis
    for i, d in enumerate(drafts):
        t_days, phi_star = results[d]
        ax.plot(t_days, phi_star, color=cmap(i / (len(drafts) - 1)), label=f"{d} m")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("time (days)")
    ax.set_ylabel(r"$\phi^*(t) = \phi(t)/\phi(0)$")
    ax.set_title("Stratification decay vs. structure draft\n"
                 "(fixed total frontal area, pycnocline at 40 m)")
    ax.legend(title="draft", fontsize=7, ncol=2)

    t50 = [crossing_time(*results[d], 0.5) for d in drafts]
    t10 = [crossing_time(*results[d], 0.1) for d in drafts]

    ax2 = axes[1]
    ax2.plot(drafts, t50, "o-", label=r"$t(\phi^*=0.5)$")
    ax2.plot(drafts, t10, "s-", label=r"$t(\phi^*=0.1)$")
    ax2.axvline(40, color="gray", ls="--", lw=1, label="pycnocline depth (40 m)")
    ax2.set_xlabel("structure draft (m)")
    ax2.set_ylabel("time (days)")
    ax2.set_title("Mixing timescale vs. draft\n(fixed total frontal area = 1.5)")
    ax2.legend(fontsize=8)

    plt.tight_layout()
    if out_path is None:
        out_path = os.path.join(FIGURES_DIR, "draft_depth_distribution_mixing.png")
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"saved {out_path}")
    return results, t50, t10


def plot_realistic_vs_matched(out_path=None):
    """
    Figure 2 (Result B): matched (total_a=1.5) vs. realistic
    (geometry-derived) total frontal area, for spar (90 m) and
    semisubmersible (20 m) drafts.
    """
    styles = {
        "geom_matched_draft90":   dict(color="tab:blue", ls="--", label="spar (90 m), matched total_a=1.5"),
        "geom_realistic_spar":    dict(color="tab:blue", ls="-", label="spar (90 m), realistic total_a=1.65"),
        "geom_matched_draft20":   dict(color="tab:orange", ls="--", label="semisub (20 m), matched total_a=1.5"),
        "geom_realistic_semisub": dict(color="tab:orange", ls="-", label="semisub (20 m), realistic total_a=0.63"),
    }

    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    for name, style in styles.items():
        t_days, phi_star = load_phi_star(name)
        ax.plot(t_days, phi_star, **style, lw=1.8)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("time (days)")
    ax.set_ylabel(r"$\phi^*(t)=\phi(t)/\phi(0)$")
    ax.set_title("Effect of realistic (vs. matched) total frontal area\n"
                 "spar vs. semisubmersible")
    ax.legend(fontsize=8)
    plt.tight_layout()
    if out_path is None:
        out_path = os.path.join(FIGURES_DIR, "realistic_total_area_mixing.png")
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"saved {out_path}")


def plot_seasonality(out_path=None):
    """
    Figure 3 (Result C): phi_star(t) for spar (90 m) and semisubmersible
    (20 m) drafts under three stratification regimes (summer / default /
    winter).
    """
    styles = {
        "season_summer_semisub": dict(color="tab:orange", ls="-",
                                       label="semisub (20 m), summer (zt=15, dT=8)"),
        "geom_matched_draft20":  dict(color="tab:orange", ls="--",
                                       label="semisub (20 m), default (zt=40, dT=10)"),
        "season_winter_semisub": dict(color="tab:orange", ls=":",
                                       label="semisub (20 m), winter (zt=60, dT=3)"),
        "season_summer_spar":    dict(color="tab:blue", ls="-", label="spar (90 m), summer"),
        "geom_matched_draft90":  dict(color="tab:blue", ls="--", label="spar (90 m), default"),
        "season_winter_spar":    dict(color="tab:blue", ls=":", label="spar (90 m), winter"),
    }

    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    for name, style in styles.items():
        t_days, phi_star = load_phi_star(name)
        ax.plot(t_days, phi_star, **style, lw=1.8)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("time (days)")
    ax.set_ylabel(r"$\phi^*(t)=\phi(t)/\phi(0)$")
    ax.set_title("Seasonal sensitivity of structure-induced mixing\n"
                 "(matched total frontal area, fixed structure draft)")
    ax.legend(fontsize=7.5, loc="upper right")
    plt.tight_layout()
    if out_path is None:
        out_path = os.path.join(FIGURES_DIR, "seasonal_sensitivity_mixing.png")
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"saved {out_path}")


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    results, t50, t10 = plot_draft_sweep()
    print("\ndraft (m)  t(phi*=0.5)  t(phi*=0.1)")
    for d, a, b in zip(DRAFTS, t50, t10):
        print(f"{d:>9}  {a:>11.2f}  {b if np.isnan(b) else round(b, 2)}")
    plot_realistic_vs_matched()
    plot_seasonality()


if __name__ == "__main__":
    main()
