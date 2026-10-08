#!/usr/bin/env python3
# analysis/shear_residence_time_tradeoff.py
"""
Residence-time vs. mixing-timescale tradeoff, repeated with a sheared
(surface-intensified, NCC-like) current instead of a depth-uniform one --
see notes/shear_residence_time_tradeoff.md. Companion to
analysis/residence_time_tradeoff.py (the original, uniform-current
version).

Uses the shear x current-strength sweep (tools/prep_shear_current_strength_
sweep.py, runs/shearcur_*) plus the already-completed runs/shear_spar and
runs/shear_semisub for the 0.15 m/s point, instead of the uniform-current
current_strength_sensitivity sweep.

Usage:
    python analysis/shear_residence_time_tradeoff.py
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

from analysis.current_strength_sensitivity import load_phi_star, diagnosed_speed
from analysis.residence_time_tradeoff import FARM_LENGTHS_M, phi_star_at_residence_time

FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

# (run_name, structure_type, draft, nominal u_target [m/s])
RUNS = [
    ("shearcur_spar_u0p05", "spar", 90.0, 0.05),
    ("shearcur_spar_u0p10", "spar", 90.0, 0.10),
    ("shear_spar", "spar", 90.0, 0.15),
    ("shearcur_spar_u0p20", "spar", 90.0, 0.20),
    ("shearcur_spar_u0p30", "spar", 90.0, 0.30),
    ("shearcur_semisub_u0p05", "semisub", 20.0, 0.05),
    ("shearcur_semisub_u0p10", "semisub", 20.0, 0.10),
    ("shear_semisub", "semisub", 20.0, 0.15),
    ("shearcur_semisub_u0p20", "semisub", 20.0, 0.20),
    ("shearcur_semisub_u0p30", "semisub", 20.0, 0.30),
]


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    by_struct = {}
    for run_name, struct, draft, u_target in RUNS:
        t_days, phi_star = load_phi_star(run_name)
        u = diagnosed_speed(run_name, draft)
        by_struct.setdefault(struct, []).append((u, t_days, phi_star, u_target))
        print(f"{run_name:24s} struct={struct:8s} u_target={u_target:.2f}  "
              f"u_diagnosed(structure-layer)={u:.4f} m/s")

    for struct in by_struct:
        by_struct[struct].sort(key=lambda row: row[0])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)

    for ax, struct in zip(axes, ["spar", "semisub"]):
        rows = by_struct[struct]
        us = [r[0] for r in rows]
        cmap = plt.cm.viridis
        for i, L in enumerate(FARM_LENGTHS_M):
            phi_exit = [phi_star_at_residence_time(t, p, u, L) for u, t, p, _ in rows]
            ax.plot(us, phi_exit, "o-", color=cmap(i / (len(FARM_LENGTHS_M) - 1)),
                    label=f"L={L/1000:.0f} km")
        ax.set_xlabel("diagnosed background current speed (m/s)")
        ax.set_title(f"{struct}: stratification remaining at farm exit (sheared current)\n"
                     r"$\phi^*(T_{res}),\ T_{res}=L/u$")
        ax.legend(fontsize=8)
        ax.axhline(0, color="gray", lw=0.5)

    axes[0].set_ylabel(r"$\phi^*(T_{res})$ (0 = fully mixed)")
    fig.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "shear_residence_time_tradeoff.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")

    for struct in ["spar", "semisub"]:
        rows = by_struct[struct]
        print(f"\n{struct}:")
        for u, t, p, u_target in rows:
            vals = [phi_star_at_residence_time(t, p, u, L) for L in FARM_LENGTHS_M]
            print(f"  u_target={u_target:.2f} u_diag={u:.3f} m/s  " +
                  "  ".join(f"L={L/1000:.0f}km:phi*={v:.3f}" if not np.isnan(v) else
                            f"L={L/1000:.0f}km:phi*=N/A(T_res>run)" for L, v in zip(FARM_LENGTHS_M, vals)))


if __name__ == "__main__":
    main()
