#!/usr/bin/env python3
# analysis/current_strength_sensitivity.py
"""
Analysis/plotting for the background current-strength (BFRC_U) sensitivity
sweep -- see notes/current_strength_sensitivity.md.

Follows the conventions of analysis/floating_structure_mixing.py (real,
dimensional time in days; phi_star(t) = phi(t)/phi(0)), since the structures
analysed here are partial-depth ("floating") and so
utils.analytic_mixing_timescale's dimensionless-time normalization does not
apply (it assumes a full-depth structure).

In addition to phi_star(t), this script diagnoses the *actual* achieved
background current speed from each run's "u" field (late-time, structure-
layer average) and compares it against the nominal target used to choose
BFRC_U -- the point being that the quadratic-drag scaling used to pick
BFRC_U (see tools/prep_current_strength_sweep.py) is only approximate for
these partial-depth structures, so the diagnosed speed is what should
actually be trusted/plotted against the mixing timescale.

Usage (regenerate the figure from the prepared/run sweep):
    python analysis/current_strength_sensitivity.py
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

# (run_name, structure_type, draft, nominal u_target [m/s])
RUNS = [
    ("current_spar_u0p05", "spar", 90.0, 0.05),
    ("geom_matched_draft90", "spar", 90.0, 0.10),
    ("current_spar_u0p15", "spar", 90.0, 0.15),
    ("current_spar_u0p20", "spar", 90.0, 0.20),
    ("current_spar_u0p30", "spar", 90.0, 0.30),
    ("current_semisub_u0p05", "semisub", 20.0, 0.05),
    ("geom_matched_draft20", "semisub", 20.0, 0.10),
    ("current_semisub_u0p15", "semisub", 20.0, 0.15),
    ("current_semisub_u0p20", "semisub", 20.0, 0.20),
    ("current_semisub_u0p30", "semisub", 20.0, 0.30),
]


def load_phi_star(run_name: str):
    """Same as analysis/floating_structure_mixing.py:load_phi_star."""
    cfg_path = os.path.join(RUNS_DIR, run_name, "resolved_config.yaml")
    ds, grid, params = open_roms_dataset(cfg_path)
    phi, _ = compute_phi(ds, grid, params)
    phi_star = (phi / phi.isel(ocean_time=0)).values
    t_days = compute_time_vector(params)
    ds.close()
    return t_days, phi_star


def crossing_time(t_days, phi_star, threshold: float) -> float:
    """First (linearly-interpolated) time phi_star drops to `threshold`."""
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


def diagnosed_speed(run_name: str, draft: float, frac_late: float = 0.2) -> float:
    """
    Late-time, structure-layer-averaged |u| (m/s), diagnosed directly from
    the model output rather than assumed from the BFRC_U scaling formula.

    Averages over the last `frac_late` fraction of the run (quasi-steady
    plateau) and over depths above `draft` (i.e. within the structure's
    depth range, z measured as depth below the surface).
    """
    cfg_path = os.path.join(RUNS_DIR, run_name, "resolved_config.yaml")
    ds, grid, params = open_roms_dataset(cfg_path)
    # z_rho: 0 at surface, negative down -> structure layer is z_rho >= -draft.
    # Grid is flat-bottom/horizontally uniform, so a single column's z_rho(s_rho)
    # is representative and can be broadcast against u's s_rho dimension.
    z_rho_col = ds["z_rho"].isel(eta_rho=0, xi_rho=0)
    u = ds["u"].where(z_rho_col >= -draft)
    n_t = ds.sizes["ocean_time"]
    i0 = int(n_t * (1.0 - frac_late))
    u_late = u.isel(ocean_time=slice(i0, None))
    speed = float(np.abs(u_late).mean())
    ds.close()
    return speed


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    results = {}
    diagnosed = {}
    for run_name, struct, draft, u_target in RUNS:
        t_days, phi_star = load_phi_star(run_name)
        results[run_name] = (t_days, phi_star)
        diagnosed[run_name] = diagnosed_speed(run_name, draft)
        print(f"{run_name:24s} struct={struct:8s} draft={draft:5.0f}m  "
              f"u_target={u_target:.2f}  u_diagnosed={diagnosed[run_name]:.4f} m/s")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    cmap_spar = plt.cm.Blues
    cmap_semi = plt.cm.Oranges
    u_targets_sorted = sorted(set(r[3] for r in RUNS))
    for run_name, struct, draft, u_target in RUNS:
        t_days, phi_star = results[run_name]
        i = u_targets_sorted.index(u_target)
        frac = 0.3 + 0.7 * i / (len(u_targets_sorted) - 1)
        cmap = cmap_spar if struct == "spar" else cmap_semi
        ls = "-" if struct == "spar" else "--"
        ax.plot(t_days, phi_star, color=cmap(frac), ls=ls,
                label=f"{struct} u={u_target:.2f}")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("time (days)")
    ax.set_ylabel(r"$\phi^*(t) = \phi(t)/\phi(0)$")
    ax.set_title("Stratification decay vs. background current strength")
    ax.legend(fontsize=6, ncol=2)

    ax2 = axes[1]
    for struct, marker in [("spar", "o"), ("semisub", "s")]:
        rows = [r for r in RUNS if r[1] == struct]
        u_diag = [diagnosed[r[0]] for r in rows]
        t50 = [crossing_time(*results[r[0]], 0.5) for r in rows]
        ax2.plot(u_diag, t50, marker + "-", label=f"{struct}: t(phi*=0.5)")
    ax2.set_xlabel("diagnosed background current speed (m/s)")
    ax2.set_ylabel("time (days)")
    ax2.set_title("Mixing timescale vs. diagnosed current speed")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "current_strength_sensitivity.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    main()
