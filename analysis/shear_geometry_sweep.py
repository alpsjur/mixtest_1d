#!/usr/bin/env python3
# analysis/shear_geometry_sweep.py
"""
Analysis/plotting for the shear-profile-vs-geometry comparison -- see
notes/shear_geometry_interaction.md.

Question: does a surface-intensified (NCC-like) sheared background current
change the spar-vs-semisubmersible mixing comparison, relative to the
uniform current used in notes/floating_structure_depth_distribution.md and
notes/current_strength_sensitivity.md? The two geometries were built to be
"matched" (same depth-integrated structure.str_a*draft = total_a = 1.5)
under a *uniform* current -- this checks whether that equivalence survives
once the current is sheared, surface-intensified, and depth-mean-matched
instead.

Reuses:
  - runs/shear_spar, runs/shear_semisub (new, bodyforce.mode="profile",
    profiles/ncc_shear_surface_intensified.txt)
  - runs/current_spar_u0p15, runs/current_semisub_u0p15 (already-run
    uniform-current baselines, bodyforce.U0=0.15 m/s, matching the shear
    profile's thickness-weighted column mean of ~0.15 m/s)

Usage:
    python analysis/shear_geometry_sweep.py
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

from current_strength_sensitivity import load_phi_star, crossing_time

RUNS_DIR = os.path.join(ROOT_DIR, "runs")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

# (run_name, structure_type, draft, forcing label)
RUNS = [
    ("current_spar_u0p15", "spar", 90.0, "uniform"),
    ("shear_spar", "spar", 90.0, "sheared"),
    ("current_semisub_u0p15", "semisub", 20.0, "uniform"),
    ("shear_semisub", "semisub", 20.0, "sheared"),
]

COLORS = {"spar": "tab:blue", "semisub": "tab:orange"}
STYLES = {"uniform": "--", "sheared": "-"}


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    results = {}
    for run_name, struct, draft, forcing in RUNS:
        t_days, phi_star = load_phi_star(run_name)
        results[run_name] = (t_days, phi_star)
        t50 = crossing_time(t_days, phi_star, 0.5)
        print(f"{run_name:24s} struct={struct:8s} forcing={forcing:8s} "
              f"t(phi*=0.5)={t50:.2f} d")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left panel: phi*(t) for all 4 runs
    ax = axes[0]
    for run_name, struct, draft, forcing in RUNS:
        t_days, phi_star = results[run_name]
        ax.plot(t_days, phi_star, color=COLORS[struct], linestyle=STYLES[forcing],
                label=f"{struct} ({forcing})")
    ax.set_xlabel("time (days)")
    ax.set_ylabel(r"$\phi^*(t) = \phi(t)/\phi(0)$")
    ax.set_title("Stratification decay: uniform vs. sheared current")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # Right panel: mixing timescale (t at phi*=0.5) comparison, uniform vs.
    # sheared, per geometry -- a clean summary of how much the shear
    # profile speeds up mixing relative to a depth-mean-matched uniform
    # current, and whether that speedup differs by geometry.
    ax = axes[1]
    x = np.arange(2)
    width = 0.35
    for i, struct in enumerate(["spar", "semisub"]):
        t_u = crossing_time(*results[f"current_{struct}_u0p15"], 0.5)
        t_s = crossing_time(*results[f"shear_{struct}"], 0.5)
        xpos = x + (i - 0.5) * width
        ax.bar(xpos[0], t_u, width, color=COLORS[struct], hatch="///", edgecolor="white")
        ax.bar(xpos[1], t_s, width, color=COLORS[struct])
    ax.set_xticks(x)
    ax.set_xticklabels(["uniform\n(U0=0.15 m/s)", "sheared\n(same depth-mean)"])
    ax.set_ylabel(r"$t(\phi^*=0.5)$ (days)")
    ax.set_title("Time to half-mix: uniform vs. sheared current")
    # de-dup legend (bar() called multiple times per color)
    handles = [plt.Rectangle((0, 0), 1, 1, color=COLORS[s]) for s in ["spar", "semisub"]]
    ax.legend(handles, ["spar", "semisub"], fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    fig.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "shear_geometry_sweep.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    main()
