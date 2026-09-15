#!/usr/bin/env python3
# analysis/mixing_timescale.py
"""
Mixing-timescale sensitivity analysis (Carpenter et al. 2016 style).

For a single run this computes:
  - phi(t): stratification potential energy anomaly (utils.compute_phi)
  - Pstr(t): diagnostic (model-derived) depth-integrated power extraction
    (utils.compute_Pstr), averaged over the last `plateau_frac` of the run
    to get a single quasi-steady Pstr_diag (the initial transient before
    the body-force-driven flow reaches its steady balance is excluded).
  - tau_mix_theory: analytic prediction from utils.analytic_mixing_timescale
    (pre-run estimate, uses the *assumed* steady u_inf).
  - tau_mix_diagnostic: same formula but using Pstr_diag instead of the
    analytic Pstr -- i.e. "what tau_mix actually was, given the model's own
    power extraction", which can differ from tau_mix_theory if the flow
    doesn't reach the assumed balance, or if str_a/CD/H0 assumptions in the
    analytic formula don't quite hold.
  - t_star = t * Pstr_diag / (g * delta_rho * H0^2)   (dimensionless time)
  - phi_star = phi(t) / phi(0)                         (normalized stratification)

Across a sweep, plotting phi_star vs t_star for every run and checking
whether the curves collapse (and cross ~0 near t_star=1) is the direct
analogue of Carpenter et al.'s Fig 4b.


Usage (single run):
    python analysis/mixing_timescale.py runs/<name>/resolved_config.yaml

Usage (whole sweep, collapse plot + sensitivity plot):
    python analysis/mixing_timescale.py --sweep sweeps/mixing_timescale/manifest.yaml
"""
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

import numpy as np

from utils.utils import (
    open_roms_dataset,
    compute_time_vector,
    compute_phi,
    compute_Pstr,
    analytic_mixing_timescale,
    load_yaml,
    G,
)


def pycnocline_length_scale(params: dict) -> float:
    """
    Geometric mean of the pycnocline centre's distances to the two
    boundaries (surface and bed): L = sqrt(zt * (H0 - zt)).

    """
    H0 = float(params["grid"]["H0"])
    zt = float(params["initial"]["temp_zt"])
    return np.sqrt(zt * (H0 - zt))


def mixing_timescale(resolved_config_path: str, plateau_frac: float = 0.5,
                      length_scale: str = "H0") -> dict:
    """
    Compute phi(t), Pstr(t), and the theoretical/diagnostic mixing time
    scales for a single completed run.

    Parameters
    ----------
    resolved_config_path : str
        Path to a resolved_config.yaml produced by prep_experiment.
    plateau_frac : float
        Fraction of the run (from the end) over which Pstr(t) is averaged
        to get the quasi-steady Pstr_diag, excluding the initial spin-up
        transient. Default 0.5 (last half of the run).
    length_scale : str
        Which length scale to use for nondimensionalizing time and tau_mix:
        "H0" (Carpenter et al.'s original choice, the full water column
        depth) or "pycnocline" (see pycnocline_length_scale() -- the
        geometric mean of the pycnocline's distances to the two
        boundaries, found empirically to collapse the sweep much better
        when grid.H0 is varied at fixed initial.temp_zt).

    Returns
    -------
    dict with keys:
        days, phi, phi_star, Pstr, Pstr_diag,
        tau_mix_theory, tau_mix_diagnostic, t_star, params
    """
    ds, grid, params = open_roms_dataset(resolved_config_path)

    phi, _ = compute_phi(ds, grid, params)
    Pstr = compute_Pstr(ds, grid, params)
    days = compute_time_vector(params)

    phi_vals = phi.values
    phi0 = phi_vals[0]
    phi_star = phi_vals / phi0 if phi0 != 0 else np.full_like(phi_vals, np.nan)

    n = len(phi_vals)
    i0 = max(0, int(np.floor(n * (1.0 - plateau_frac))))
    Pstr_diag = float(np.mean(Pstr.values[i0:]))

    theory = analytic_mixing_timescale(params)
    tau_mix_theory = theory["tau_mix"]
    delta_rho = theory["delta_rho"]
    H0 = float(params["grid"]["H0"])

    if length_scale == "H0":
        L = H0
    elif length_scale == "pycnocline":
        L = pycnocline_length_scale(params)
    else:
        raise ValueError(f"Unknown length_scale {length_scale!r}: use 'H0' or 'pycnocline'.")

    tau_mix_diagnostic = G * delta_rho * L ** 2 / Pstr_diag if Pstr_diag > 0 else np.nan

    t_seconds = days * 86400.0
    t_star = t_seconds * Pstr_diag / (G * delta_rho * L ** 2)

    # Empirical "mixing complete" dimensionless time: first t_star at which
    # phi_star drops below 5% of its initial value. Unlike tau_mix_theory
    # (which only encodes the power *input* scale g*delta_rho*H^2/Pstr),
    # this captures how efficiently the GLS closure actually converts that
    # power into eddy diffusivity -- found empirically to depend strongly on
    # structure.c4 and grid.H0 (pycnocline position), but negligibly on
    # structure.CD or initial.temp_dT (which only set the power/buoyancy
    # scales already absorbed into t_star itself).
    below = np.where(phi_star < 0.05)[0]
    if len(below):
        i1 = below[0]
        if i1 > 0:
            # Linear interpolation between the last sample >= 0.05 and the
            # first sample < 0.05 for sub-output-step resolution -- with
            # hourly history output (NHIS=90, DT=40s), the raw t_star grid
            # spacing (~0.01 in t_star for the sweeps in this repo) can be
            # too coarse to resolve small A(c4) differences at low c4 (see
            # notes/mixing_timescale_analysis.md, c4 fine-sweep section),
            # where several adjacent c4 values would otherwise land on the
            # same discrete output timestep and appear identical.
            i0 = i1 - 1
            x0, x1 = phi_star[i0], phi_star[i1]
            t0, t1 = t_star[i0], t_star[i1]
            frac = (0.05 - x0) / (x1 - x0) if x1 != x0 else 0.0
            t_star_mix = float(t0 + frac * (t1 - t0))
        else:
            t_star_mix = float(t_star[i1])
    else:
        t_star_mix = np.nan

    return {
        "days": days,
        "phi": phi_vals,
        "phi_star": phi_star,
        "Pstr": Pstr.values,
        "Pstr_diag": Pstr_diag,
        "tau_mix_theory": tau_mix_theory,
        "tau_mix_diagnostic": tau_mix_diagnostic,
        "t_star": t_star,
        "t_star_mix": t_star_mix,
        "params": params,
    }


def summarize_sweep(manifest_path: str, plateau_frac: float = 0.5, length_scale: str = "H0") -> list:
    """
    Compute mixing_timescale() for every completed run in a sweep manifest.

    Returns a list of dicts (one per run) with the swept parameters plus
    tau_mix_theory/tau_mix_diagnostic/Pstr_diag, for tabulating how the
    actual (diagnostic) mixing time scale depends on each parameter -- in
    particular whether it depends on structure.c4 even though
    tau_mix_theory does not.
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
        rows.append({
            "run_name": r.get("run_name"),
            "params": r.get("params", {}),
            "tau_mix_theory_days": result["tau_mix_theory"] / 86400.0,
            "tau_mix_diagnostic_days": result["tau_mix_diagnostic"] / 86400.0,
            "Pstr_diag": result["Pstr_diag"],
            "t_star_mix": result["t_star_mix"],
        })
    return rows


def plot_sensitivity(manifest_path: str, ax=None, plateau_frac: float = 0.5, length_scale: str = "H0"):
    """
    Plot the empirical dimensionless mixing-completion time t_star_mix
    (first t_star at which phi_star < 0.05) against structure.c4, grouped
    by grid.H0, averaging over structure.CD/initial.temp_dT (which were
    found to have negligible effect once time is nondimensionalized by
    Pstr_diag). This isolates the closure-efficiency sensitivity that the
    simple Carpenter et al. power-based tau_mix scaling does not capture.

    With length_scale="pycnocline" the residual grouping by grid.H0 mostly
    disappears (see pycnocline_length_scale() and
    notes/mixing_timescale_analysis.md) -- kept as an argument here so this
    plot can be used to visually confirm that collapse.
    """
    import matplotlib.pyplot as plt
    from collections import defaultdict

    manifest = load_yaml(manifest_path)
    by_c4_H0 = defaultdict(list)
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
        key = (p.get("structure.c4"), p.get("grid.H0"))
        by_c4_H0[key].append(result["t_star_mix"])

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))

    H0_vals = sorted({k[1] for k in by_c4_H0})
    cmap = plt.get_cmap("viridis")
    for i, H0 in enumerate(H0_vals):
        c4_vals = sorted(k[0] for k in by_c4_H0 if k[1] == H0)
        means = [np.mean(by_c4_H0[(c4, H0)]) for c4 in c4_vals]
        stds = [np.std(by_c4_H0[(c4, H0)]) for c4 in c4_vals]
        ax.errorbar(
            c4_vals, means, yerr=stds, marker="o",
            label=f"H0={H0} m",
            color=cmap(i / max(len(H0_vals) - 1, 1)),
        )

    ax.set_xlabel(r"structure.c4")
    ax.set_ylabel(r"$t^*_{mix}$ (dimensionless time for $\phi^*<0.05$)")
    ax.set_title(f"Mixing-completion sensitivity to closure coefficient c4 (L={length_scale})")
    ax.grid(True, alpha=0.3)
    ax.legend()
    return ax


def plot_collapse(manifest_path: str, ax=None, plateau_frac: float = 0.5, length_scale: str = "H0"):
    """
    Plot phi_star(t_star) for every completed run in a sweep manifest onto
    one axis, to check for the Carpenter-et-al.-style collapse. Adds a
    vertical reference line at t_star=1 (the theoretical "fully mixed"
    point) and a horizontal line at phi_star=0.

    length_scale="H0" reproduces Carpenter et al.'s original choice (the
    full water column depth). length_scale="pycnocline" uses
    pycnocline_length_scale() instead -- found empirically to give a much
    tighter collapse when grid.H0 varies at fixed initial.temp_zt (see
    notes/mixing_timescale_analysis.md).
    """
    import matplotlib.pyplot as plt
    try:
        from cmcrameri import cm as cmc
        cmap = cmc.batlow
    except ModuleNotFoundError:
        cmap = plt.get_cmap("viridis")

    manifest = load_yaml(manifest_path)
    runs = manifest.get("runs", [])

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))

    n_plotted = 0
    for r in runs:
        resolved_config = r["resolved_config"]
        if not os.path.isfile(resolved_config):
            continue
        try:
            result = mixing_timescale(resolved_config, plateau_frac=plateau_frac, length_scale=length_scale)
        except Exception as e:
            print(f"Skipping {r.get('run_name')}: {e}")
            continue
        ax.plot(
            result["t_star"], result["phi_star"],
            label=r.get("run_name"),
            color=cmap(n_plotted / max(len(runs) - 1, 1)),
        )
        n_plotted += 1

    ax.axvline(1.0, color="k", linestyle="--", linewidth=1, label=r"$t^*=1$")
    ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
    if length_scale == "H0":
        ax.set_xlabel(r"$t^* = t\,P_{str}/(g\,\Delta\rho\,H_0^2)$")
    else:
        ax.set_xlabel(r"$t^* = t\,P_{str}/(g\,\Delta\rho\,L^2)$, $L=\sqrt{z_t(H_0-z_t)}$")
    ax.set_ylabel(r"$\phi(t)/\phi(0)$")
    ax.set_title(f"Mixing timescale collapse (L={length_scale})")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="best")
    return ax


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("resolved_config", nargs="?", help="Path to a single run's resolved_config.yaml")
    parser.add_argument("--sweep", type=str, help="Path to a sweep manifest.yaml (collapse plot)")
    parser.add_argument("--plateau-frac", type=float, default=0.5)
    parser.add_argument("--length-scale", choices=["H0", "pycnocline"], default="H0",
                         help="Length scale for nondimensionalizing t_star/tau_mix (default: H0, "
                              "Carpenter et al.'s original choice; 'pycnocline' uses "
                              "sqrt(zt*(H0-zt)), which collapses better when H0 is varied)")
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()

    if args.sweep:
        ax = plot_collapse(args.sweep, plateau_frac=args.plateau_frac, length_scale=args.length_scale)
        import matplotlib.pyplot as plt
        if args.save:
            filename = f"figures/mixing_timescale_collapse_{os.path.basename(os.path.dirname(args.sweep))}_{args.length_scale}.png"
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            plt.savefig(filename)
            print(f"Plot saved to {filename}")
        else:
            plt.show()

        ax2 = plot_sensitivity(args.sweep, plateau_frac=args.plateau_frac, length_scale=args.length_scale)
        if args.save:
            filename2 = f"figures/mixing_timescale_sensitivity_{os.path.basename(os.path.dirname(args.sweep))}_{args.length_scale}.png"
            plt.savefig(filename2)
            print(f"Plot saved to {filename2}")
        else:
            plt.show()
    elif args.resolved_config:
        result = mixing_timescale(args.resolved_config, plateau_frac=args.plateau_frac, length_scale=args.length_scale)
        print(f"tau_mix_theory:      {result['tau_mix_theory']/86400.0:.3f} days")
        print(f"tau_mix_diagnostic:  {result['tau_mix_diagnostic']/86400.0:.3f} days")
        print(f"Pstr_diag:           {result['Pstr_diag']:.6g} W/m2")

        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
        axes[0].plot(result["days"], result["phi"])
        axes[0].set_xlabel("Time (days)")
        axes[0].set_ylabel(r"$\phi$ (J m$^{-2}$)")
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(result["t_star"], result["phi_star"])
        axes[1].axvline(1.0, color="k", linestyle="--", linewidth=1)
        axes[1].axhline(0.0, color="k", linestyle=":", linewidth=1)
        axes[1].set_xlabel(r"$t^*$")
        axes[1].set_ylabel(r"$\phi/\phi(0)$")
        axes[1].grid(True, alpha=0.3)

        fig.tight_layout()
        if args.save:
            filename = f"figures/mixing_timescale_{args.resolved_config.split('/')[-2]}.png"
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            plt.savefig(filename)
            print(f"Plot saved to {filename}")
        else:
            plt.show()
    else:
        parser.error("Provide either a resolved_config path or --sweep manifest.yaml")


if __name__ == "__main__":
    main()
