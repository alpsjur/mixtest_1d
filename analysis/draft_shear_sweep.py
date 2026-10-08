#!/usr/bin/env python3
# analysis/draft_shear_sweep.py
"""
Analysis/plotting for the "Goldilocks draft" sweep -- see
notes/goldilocks_draft_depth.md.

Fixed structure.str_a (realistic, physically-anchored -- see
tools/prep_draft_shear_sweep.py), sweeping structure.depth_zero_below
("draft") continuously under the surface-intensified NCC-like shear
profile. Plots mixing timescale (t(phi*=0.5)) vs. draft to look for a
non-monotonic "worst-case"/"best-case" draft depth, rather than assuming
mixing simply gets stronger or weaker monotonically with draft.

Usage:
    python analysis/draft_shear_sweep.py
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

DRAFTS = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]

PYCNOCLINE_DEPTH = 40.0  # temp_zt default, m


def run_name(draft: int) -> str:
    return f"draft_shear_d{draft:03d}"


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    t50 = []
    t10 = []
    for draft in DRAFTS:
        t_days, phi_star = load_phi_star(run_name(draft))
        t50.append(crossing_time(t_days, phi_star, 0.5))
        t10.append(crossing_time(t_days, phi_star, 0.1))
        print(f"draft={draft:4d} m  t(phi*=0.5)={t50[-1]:.3f} d  "
              f"t(phi*=0.1)={t10[-1]:.3f} d")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(DRAFTS, t50, "o-", color="tab:blue", label=r"$t(\phi^*=0.5)$")
    ax.plot(DRAFTS, t10, "s--", color="tab:blue", alpha=0.5, label=r"$t(\phi^*=0.1)$")
    ax.axvline(PYCNOCLINE_DEPTH, color="gray", ls=":", lw=1.2,
               label=f"pycnocline depth ({PYCNOCLINE_DEPTH:.0f} m)")

    i_min = int(np.nanargmin(t50))
    ax.annotate(f"fastest mixing:\ndraft={DRAFTS[i_min]} m",
                xy=(DRAFTS[i_min], t50[i_min]),
                xytext=(DRAFTS[i_min] + 15, t50[i_min] + 0.15 * max(t50)),
                arrowprops=dict(arrowstyle="->", color="black"), fontsize=9)

    ax.set_xlabel("structure draft (m)")
    ax.set_ylabel("time (days)")
    ax.set_title("Mixing timescale vs. structure draft\n"
                  "(fixed, realistic str_a; surface-intensified shear profile)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()

    out_path = os.path.join(FIGURES_DIR, "goldilocks_draft_depth.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    main()
