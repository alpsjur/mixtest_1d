#!/usr/bin/env python3
# analysis/gamma_observational_check.py
"""
Analysis/plotting for a first observational sanity check of the
structure-mixing parametrization -- see notes/gamma_observational_check.md.

Diagnoses the bulk "mixing coefficient" / mixing efficiency

    Gamma(z,t) = -B(z,t) / epsilon(z,t)
    B(z,t)     = (g/rho0) * AKt(z,t) * d(rho)/dz   (buoyancy destruction of TKE)
    epsilon(z,t) = utils.utils.compute_epsilon(ds, params)   (GLS-generic)

directly from already-completed model output (no new simulations), and
compares it against the field/LES-based bulk mixing-efficiency estimate of
8-14% reported by Schultze et al. (2020, JGR Oceans) for monopile wind
turbine foundation wakes -- the one directly relevant *observational*
benchmark found for this class of parametrization (see
notes/gamma_observational_check.md for the full discussion/caveats).

An epsilon-weighted average (sum(Gamma*eps)/sum(eps) over the structure's
depth range and an early quasi-steady time window, before the water column
is close to fully mixed) is used rather than a simple unweighted mean,
since epsilon can be locally very small (dividing B by it would otherwise
produce unphysical outliers dominating a plain average) -- the
epsilon-weighted average is the physically meaningful "how much of the
*actual* dissipated turbulence went into mixing" bulk number.

Usage:
    python analysis/gamma_observational_check.py
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

from utils.utils import open_roms_dataset, compute_epsilon, compute_time_vector, G, RHO0

RUNS_DIR = os.path.join(ROOT_DIR, "runs")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

# Schultze et al. (2020) field/LES bulk mixing-efficiency estimate for
# monopile wakes (fraction of structure-generated TKE that goes into
# mixing rather than being dissipated).
SCHULTZE_GAMMA_RANGE = (0.08, 0.14)

# (run_name, structure_type, draft, label)
RUNS = [
    ("current_spar_u0p15", "spar", 90.0, "u=0.15 m/s"),
    ("geom_matched_draft90", "spar", 90.0, "u=0.10 m/s"),
    ("current_semisub_u0p15", "semisub", 20.0, "u=0.15 m/s"),
    ("geom_matched_draft20", "semisub", 20.0, "u=0.10 m/s"),
]

COLORS = {"spar": "tab:blue", "semisub": "tab:orange"}


def diagnosed_gamma(run_name: str, draft: float, t_window=(0.1, 0.4)) -> float:
    """
    Epsilon-weighted bulk Gamma over the structure's depth range (z_w >=
    -draft) and an early time window (days) chosen to be well before the
    column is close to fully mixed (phi* collapses to ~0), so there is
    still meaningful stratification/buoyancy flux to diagnose from.
    """
    cfg_path = os.path.join(RUNS_DIR, run_name, "resolved_config.yaml")
    ds, grid, params = open_roms_dataset(cfg_path)
    t_days = compute_time_vector(params)

    drho_dz = grid.derivative(ds["rho"], "Z")
    B = (G / RHO0) * ds["AKt"] * drho_dz
    eps = compute_epsilon(ds, params)
    Gamma = -B / eps

    z_w_col = ds["z_w"].isel(eta_rho=0, xi_rho=0)
    mask = z_w_col >= -draft
    idx = np.where((t_days >= t_window[0]) & (t_days <= t_window[1]))[0]

    Gamma_sel = Gamma.where(mask).isel(ocean_time=idx)
    eps_sel = eps.where(mask).isel(ocean_time=idx)
    val = float((Gamma_sel * eps_sel).sum() / eps_sel.sum())
    ds.close()
    return val


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    results = []
    for run_name, struct, draft, label in RUNS:
        g = diagnosed_gamma(run_name, draft)
        results.append((run_name, struct, draft, label, g))
        print(f"{run_name:24s} struct={struct:8s} {label:12s} Gamma_diagnosed={g:.4f}")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    x = np.arange(len(results))
    colors = [COLORS[r[1]] for r in results]
    bars = ax.bar(x, [r[4] for r in results], color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r[1]}\n{r[3]}" for r in results])
    ax.axhspan(*SCHULTZE_GAMMA_RANGE, color="green", alpha=0.2,
               label="Schultze et al. (2020)\nfield/LES estimate (8-14%)")
    ax.set_ylabel(r"diagnosed bulk $\Gamma = -B/\epsilon$")
    ax.set_title("Structure-mixing efficiency vs. an observational benchmark\n"
                  "(same str_c4=0.44 for all runs)")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "gamma_observational_check.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    main()
