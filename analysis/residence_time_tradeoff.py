#!/usr/bin/env python3
# analysis/residence_time_tradeoff.py
"""
Residence-time vs. mixing-timescale tradeoff for a finite-extent floating
wind farm -- see notes/residence_time_tradeoff.md.

Motivating idea (user-suggested): a real wind farm has a finite along-flow
extent L. A water column only experiences the farm's structure drag for a
residence time T_res(u) = L/u while crossing it -- so even though faster
currents mix a *stationary* water column faster (shorter t_mix, see
notes/current_strength_sensitivity.md), they also sweep each water column
through the farm faster, leaving less *time* (not necessarily less
*progress*) to exploit that faster mixing. This script evaluates
phi*(T_res(u)) using the already-completed current-strength sweep's own
phi*(t) curves (one curve per diagnosed current speed, already run -- no
new ROMS runs needed), for a few representative farm lengths L.

Usage:
    python analysis/residence_time_tradeoff.py
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

from analysis.current_strength_sensitivity import RUNS, load_phi_star, diagnosed_speed

FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

# Representative along-flow farm lengths (m). Hywind Tampen spans roughly
# 7-10 km across its 11-turbine array; Utsira Nord and other planned
# Norwegian floating-wind areas are planned at O(10 km) scale. 2, 5, 10,
# 20 km bracket "a handful of turbine spacings" to "a large array".
FARM_LENGTHS_M = [2000.0, 5000.0, 10000.0, 20000.0]


def phi_star_at_residence_time(t_days: np.ndarray, phi_star: np.ndarray,
                                u: float, L_m: float):
    """
    phi*(T_res), T_res = L_m / u (seconds), linearly interpolated from the
    run's own phi*(t) time series. Returns NaN if T_res exceeds the run's
    simulated duration (residence time longer than the whole run -- the
    farm-crossing mixing estimate is then only a lower bound).
    """
    t_res_days = (L_m / u) / 86400.0
    if t_res_days > t_days[-1]:
        return np.nan
    return float(np.interp(t_res_days, t_days, phi_star))


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    # Group runs by structure type, collecting (u_diagnosed, t_days, phi_star).
    by_struct = {}
    for run_name, struct, draft, u_target in RUNS:
        t_days, phi_star = load_phi_star(run_name)
        u = diagnosed_speed(run_name, draft)
        by_struct.setdefault(struct, []).append((u, t_days, phi_star))

    for struct in by_struct:
        by_struct[struct].sort(key=lambda row: row[0])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)

    for ax, struct in zip(axes, ["spar", "semisub"]):
        rows = by_struct[struct]
        us = [r[0] for r in rows]
        cmap = plt.cm.viridis
        for i, L in enumerate(FARM_LENGTHS_M):
            phi_exit = [phi_star_at_residence_time(t, p, u, L) for u, t, p in rows]
            ax.plot(us, phi_exit, "o-", color=cmap(i / (len(FARM_LENGTHS_M) - 1)),
                    label=f"L={L/1000:.0f} km")
        ax.set_xlabel("diagnosed background current speed (m/s)")
        ax.set_title(f"{struct}: stratification remaining at farm exit\n"
                     r"$\phi^*(T_{res}),\ T_{res}=L/u$")
        ax.legend(fontsize=8)
        ax.axhline(0, color="gray", lw=0.5)

    axes[0].set_ylabel(r"$\phi^*(T_{res})$ (0 = fully mixed)")
    fig.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "residence_time_tradeoff.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")

    # Print a small table for the note.
    for struct in ["spar", "semisub"]:
        rows = by_struct[struct]
        print(f"\n{struct}:")
        for u, t, p in rows:
            vals = [phi_star_at_residence_time(t, p, u, L) for L in FARM_LENGTHS_M]
            print(f"  u={u:.3f} m/s  " +
                  "  ".join(f"L={L/1000:.0f}km:phi*={v:.3f}" if not np.isnan(v) else
                            f"L={L/1000:.0f}km:phi*=N/A(T_res>run)" for L, v in zip(FARM_LENGTHS_M, vals)))


if __name__ == "__main__":
    main()
