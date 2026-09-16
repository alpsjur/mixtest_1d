#!/usr/bin/env python3
# analysis/predict_pea.py
"""
Predictive function for the potential energy anomaly PEA(t) = phi(t), given
STRUCTURE_MIXING/GLS parameters, combining:

  1. tau_mix_theory: the pure power-balance analytic time scale (Carpenter
     et al. 2016), computed BEFORE any run from utils.analytic_mixing_timescale,
     using the pycnocline length scale L = sqrt(z_t*(H0-z_t)) (found in
     notes/mixing_timescale_analysis.md sections 4.5/4.6 to collapse the
     H0/z_t geometry dependence far better than the raw water-column depth
     H0 used in Carpenter et al.'s original formula).
  2. A(c4): an empirical closure-efficiency correction factor, fit to the
     combined c4_fine_sweep (7 points, c4 in [0.10, 1.00]) and
     c4_near_c1_sweep (7 points, c4 in [0.875, 1.00], x 2 H0 x 2 z_t) data.
     A(c4) is well described by a smooth quadratic base plus a steeply
     rising term that only matters as c4 -> GLS.C1:

         A(c4) = a0 + a1*c4 + a2*c4^2 + a3*c4^n

     (no true singularity at c4 = C1 -- A(c4) stays finite there; it is
     simply a rapidly-but-smoothly accelerating function of c4). This
     single functional form fits smoothly across the whole valid range
     0 < c4 <= C1, so -- contrary to the original working hypothesis in
     notes/mixing_timescale_analysis.md section 6 -- *no* separate
     "low-c4"/"near-c1" piecewise functions are needed; one continuous fit
     suffices given fine enough sampling near c4 = C1.
  3. S(x): a universal dimensionless shape function for phi_star(t_star),
     fit once across ALL 74 regenerated runs (4 sweeps): phi_star, as a
     function of x = t_star / t_star_mix(c4), collapses onto a single
     curve to ~1-2% (the per-run scatter at fixed x), confirming a genuine
     *shape* collapse -- not just a completion-*time* collapse. S(x) is
     modeled as the complementary regularized incomplete Beta function
     (the survival function of a Beta(a,b) distribution on [0,1]):

         S(x) = 1 - I_x(a, b),   x in [0,1];   S(x) = 0 for x > 1

     which is bounded, monotonically decreasing, exactly 1 at x=0 and 0 at
     x=1, and only needs 2 free parameters. This fits the ensemble-mean
     shape to RMSE ~1.3% (vs ~2.8% for a plain linear ramp), essentially at
     the level of the inter-run scatter.

Putting these together:

    phi(t) ~= phi(0) * S( t / (tau_mix_theory_pycnocline * A(c4)) )

Usage:
    # Refit A(c4) and S(x) from the regenerated sweeps and print parameters:
    python analysis/predict_pea.py --fit

    # Validate predict_phi() against every completed run across all 4 sweeps
    # (reports per-run and overall RMSE/R^2, saves a comparison figure):
    python analysis/predict_pea.py --validate
"""
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import betainc

from utils.utils import (
    load_yaml,
    analytic_mixing_timescale,
    G,
)
from mixing_timescale import (
    mixing_timescale,
    summarize_sweep,
    pycnocline_length_scale,
)

# ---------------------------------------------------------------------------
# Fitted defaults (from the 74-run sweep set regenerated for this analysis:
# mixing_timescale, pycnocline_zt, c4_fine, c4_near_c1 -- see
# notes/mixing_timescale_analysis.md section 4.9). Re-derive with --fit.
# ---------------------------------------------------------------------------
DEFAULT_A_C4_PARAMS = dict(a0=2.16959109, a1=0.07336219, a2=0.11181842, a3=1.01193866, n=16.41007755)
DEFAULT_SHAPE_PARAMS = dict(a=0.95238207, b=0.86847507)

C4_SWEEP_MANIFESTS = [
    "sweeps/c4_fine/manifest.yaml",
    "sweeps/c4_near_c1/manifest.yaml",
]
ALL_SWEEP_MANIFESTS = [
    "sweeps/mixing_timescale/manifest.yaml",
    "sweeps/pycnocline_zt/manifest.yaml",
    "sweeps/c4_fine/manifest.yaml",
    "sweeps/c4_near_c1/manifest.yaml",
]


def A_c4_model(c4, a0, a1, a2, a3, n):
    """Blended quadratic + steep-power-law correction model for A(c4)."""
    return a0 + a1 * c4 + a2 * c4 ** 2 + a3 * c4 ** n


def shape_function(x, a, b):
    """
    Universal normalized phi_star(t_star/t_star_mix) shape: complementary
    regularized incomplete Beta function, S(x) = 1 - I_x(a, b) for x in
    [0, 1], clipped to 0 for x > 1 (mixing complete/plateaued).
    """
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    inside = (x >= 0) & (x < 1)
    xin = np.clip(x[inside], 1e-12, 1.0 - 1e-12)
    out[inside] = 1.0 - betainc(a, b, xin)
    out[x < 0] = 1.0
    return out


def fit_A_c4(manifests=C4_SWEEP_MANIFESTS, length_scale: str = "pycnocline") -> dict:
    """
    Fit A(c4) = a0 + a1*c4 + a2*c4^2 + a3*c4^n to t_star_mix(c4), averaged
    over grid.H0/initial.temp_zt (found to be a negligible residual
    dependence once t_star_mix is computed with the pycnocline length
    scale -- see notes/mixing_timescale_analysis.md sections 4.5/4.6/4.9).

    Returns dict of fitted parameters (a0, a1, a2, a3, n) plus the raw
    (c4, mean t_star_mix, std) data used, for diagnostics/plotting.
    """
    from collections import defaultdict
    by_c4 = defaultdict(list)
    for manifest in manifests:
        path = os.path.join(ROOT_DIR, manifest)
        rows = summarize_sweep(path, length_scale=length_scale)
        for r in rows:
            c4 = r["params"]["structure.c4"]
            by_c4[c4].append(r["t_star_mix"])

    c4_vals = np.array(sorted(by_c4.keys()))
    A_vals = np.array([np.mean(by_c4[c4]) for c4 in c4_vals])
    A_std = np.array([np.std(by_c4[c4]) for c4 in c4_vals])

    popt, _ = curve_fit(
        A_c4_model, c4_vals, A_vals,
        p0=[2.1, 0.2, 0.3, 1.0, 8.0], maxfev=20000,
    )
    pred = A_c4_model(c4_vals, *popt)
    rmse = float(np.sqrt(np.mean((pred - A_vals) ** 2)))

    return {
        "a0": float(popt[0]), "a1": float(popt[1]), "a2": float(popt[2]),
        "a3": float(popt[3]), "n": float(popt[4]),
        "rmse": rmse, "rel_rmse": rmse / float(np.mean(A_vals)),
        "c4_vals": c4_vals, "A_vals": A_vals, "A_std": A_std,
    }


def fit_shape_function(manifests=ALL_SWEEP_MANIFESTS, length_scale: str = "pycnocline",
                        n_bins: int = 40, x_max: float = 1.3) -> dict:
    """
    Aggregate phi_star(t_star/t_star_mix) across every completed run in the
    given sweeps, bin-average, and fit the Beta-CDF-complement shape
    function S(x; a, b) to the binned ensemble mean (x <= 1.05, i.e. up to
    and just past the nominal mixing-completion point).

    Returns dict with fitted (a, b), fit RMSE, and the binned (x, mean, std)
    diagnostic arrays used for the fit / for per-run goodness-of-fit checks
    against the piecewise-linear ramp hypothesis.
    """
    all_x, all_y = [], []
    n_used = 0
    for manifest in manifests:
        path = os.path.join(ROOT_DIR, manifest)
        m = load_yaml(path)
        for r in m.get("runs", []):
            resolved_config = r["resolved_config"]
            if not os.path.isabs(resolved_config):
                resolved_config = os.path.join(ROOT_DIR, resolved_config)
            if not os.path.isfile(resolved_config):
                continue
            try:
                result = mixing_timescale(resolved_config, length_scale=length_scale)
            except Exception:
                continue
            tmix = result["t_star_mix"]
            if not np.isfinite(tmix) or tmix <= 0:
                continue
            x = result["t_star"] / tmix
            y = result["phi_star"]
            mask = x <= x_max
            all_x.append(x[mask])
            all_y.append(y[mask])
            n_used += 1

    x = np.concatenate(all_x)
    y = np.concatenate(all_y)

    bins = np.linspace(0, x_max, n_bins)
    idx = np.digitize(x, bins)
    centers, means, stds = [], [], []
    for i in range(1, len(bins)):
        sel = idx == i
        if sel.sum() > 5:
            centers.append(0.5 * (bins[i - 1] + bins[i]))
            means.append(y[sel].mean())
            stds.append(y[sel].std())
    centers = np.array(centers)
    means = np.array(means)
    stds = np.array(stds)

    fit_mask = centers <= 1.05
    popt, _ = curve_fit(shape_function, centers[fit_mask], means[fit_mask],
                         p0=[1.0, 1.0], maxfev=10000)
    pred = shape_function(centers[fit_mask], *popt)
    rmse = float(np.sqrt(np.mean((pred - means[fit_mask]) ** 2)))
    lin_pred = np.clip(1 - centers[fit_mask], 0, None)
    lin_rmse = float(np.sqrt(np.mean((lin_pred - means[fit_mask]) ** 2)))

    return {
        "a": float(popt[0]), "b": float(popt[1]),
        "rmse": rmse, "linear_ramp_rmse": lin_rmse,
        "n_runs_used": n_used,
        "x_centers": centers, "y_means": means, "y_stds": stds,
    }


def predict_phi(t_seconds, params: dict, phi0: float = None,
                 A_c4_params: dict = None, shape_params: dict = None):
    """
    Predict phi(t) given a run's resolved config `params` and a time axis
    `t_seconds` (array or scalar, seconds since run start).

    phi(t) ~= phi(0) * S( t / (tau_mix_theory_pycnocline * A(c4)) )

    where tau_mix_theory_pycnocline is the analytic power-balance mixing
    time scale (utils.analytic_mixing_timescale), renormalized to use the
    pycnocline length scale L = sqrt(z_t*(H0-z_t)) instead of the raw H0,
    A(c4) is the closure-efficiency correction (fit_A_c4), and S(x) is the
    universal normalized shape function (fit_shape_function).

    If `phi0` (phi at t=0) is not supplied, it is computed analytically
    from the initial condition: phi(0) = (rho0/8) * g * TCOEF * temp_dT *
    temp_ht * H0^2 for the tanh thermocline profile used in this repo's
    initial conditions -- NOT implemented here (would require importing
    the initial-condition profile shape); callers should generally pass
    the actual phi(0) diagnosed from compute_phi()[0] instead.

    Returns
    -------
    phi_pred : np.ndarray (or float), same shape as t_seconds
    """
    if A_c4_params is None:
        A_c4_params = DEFAULT_A_C4_PARAMS
    if shape_params is None:
        shape_params = DEFAULT_SHAPE_PARAMS
    if phi0 is None:
        raise ValueError(
            "predict_phi requires phi0 (phi(t=0)); compute it from the "
            "initial condition via utils.compute_phi(ds, grid, params)[0][0] "
            "on the run's own history file, or from a matching run's "
            "diagnosed phi(0)."
        )

    theory = analytic_mixing_timescale(params)
    delta_rho = theory["delta_rho"]
    Pstr_theory = theory["Pstr"]

    L = pycnocline_length_scale(params)
    tau_mix_theory_pyc = G * delta_rho * L ** 2 / Pstr_theory

    c4 = float(params["structure"]["c4"])
    A = A_c4_model(c4, A_c4_params["a0"], A_c4_params["a1"], A_c4_params["a2"],
                    A_c4_params["a3"], A_c4_params["n"])

    x = np.asarray(t_seconds, dtype=float) / (tau_mix_theory_pyc * A)
    S = shape_function(x, shape_params["a"], shape_params["b"])
    return phi0 * S


def validate_predict_phi(manifests=ALL_SWEEP_MANIFESTS, A_c4_params=None, shape_params=None,
                          length_scale: str = "pycnocline") -> list:
    """
    For every completed run in `manifests`, compute predict_phi(t) (using
    the run's own diagnosed phi(0) as the anchor) and compare against the
    actual phi(t) time series. Returns a list of per-run dicts with
    run_name, params, rmse (relative to phi(0)), and R^2.
    """
    rows = []
    for manifest in manifests:
        path = os.path.join(ROOT_DIR, manifest)
        m = load_yaml(path)
        for r in m.get("runs", []):
            resolved_config = r["resolved_config"]
            if not os.path.isabs(resolved_config):
                resolved_config = os.path.join(ROOT_DIR, resolved_config)
            if not os.path.isfile(resolved_config):
                continue
            try:
                result = mixing_timescale(resolved_config, length_scale=length_scale)
            except Exception as e:
                print(f"Skipping {r.get('run_name')}: {e}")
                continue

            phi_actual = result["phi"]
            phi0 = phi_actual[0]
            t_seconds = result["days"] * 86400.0
            phi_pred = predict_phi(t_seconds, result["params"], phi0=phi0,
                                    A_c4_params=A_c4_params, shape_params=shape_params)

            resid = (phi_pred - phi_actual) / phi0
            rmse = float(np.sqrt(np.mean(resid ** 2)))
            ss_res = np.sum((phi_actual - phi_pred) ** 2)
            ss_tot = np.sum((phi_actual - np.mean(phi_actual)) ** 2)
            r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan

            rows.append({
                "run_name": r.get("run_name"),
                "params": r.get("params", {}),
                "rmse_rel": rmse,
                "r2": r2,
            })
    return rows


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fit", action="store_true", help="Refit A(c4) and S(x) from the sweeps and print parameters")
    parser.add_argument("--validate", action="store_true", help="Validate predict_phi() against all completed runs")
    parser.add_argument("--save", action="store_true", help="Save validation figure to figures/")
    args = parser.parse_args()

    if args.fit:
        print("=== Fitting A(c4) (c4_fine + c4_near_c1 sweeps) ===")
        A_fit = fit_A_c4()
        print(f"a0={A_fit['a0']:.6f} a1={A_fit['a1']:.6f} a2={A_fit['a2']:.6f} "
              f"a3={A_fit['a3']:.6f} n={A_fit['n']:.6f}")
        print(f"RMSE={A_fit['rmse']:.5f}  relative RMSE={A_fit['rel_rmse']*100:.3f}%")
        for c4, a, s in zip(A_fit["c4_vals"], A_fit["A_vals"], A_fit["A_std"]):
            print(f"  c4={c4:.3f}  A_data={a:.4f}  std={s:.4f}  A_model={A_c4_model(c4, *[A_fit[k] for k in ('a0','a1','a2','a3','n')]):.4f}")

        print("\n=== Fitting universal shape function S(x) (all 4 sweeps) ===")
        S_fit = fit_shape_function()
        print(f"a={S_fit['a']:.6f} b={S_fit['b']:.6f}")
        print(f"RMSE={S_fit['rmse']:.5f} (linear-ramp RMSE for comparison: {S_fit['linear_ramp_rmse']:.5f})")
        print(f"n_runs_used={S_fit['n_runs_used']}")
        return

    if args.validate:
        rows = validate_predict_phi()
        rmses = np.array([r["rmse_rel"] for r in rows])
        r2s = np.array([r["r2"] for r in rows])
        print(f"n_runs validated: {len(rows)}")
        print(f"RMSE (relative to phi(0)): mean={rmses.mean():.4f} median={np.median(rmses):.4f} "
              f"max={rmses.max():.4f}")
        print(f"R^2: mean={r2s.mean():.4f} median={np.median(r2s):.4f} min={r2s.min():.4f}")
        worst = sorted(rows, key=lambda r: -r["rmse_rel"])[:5]
        print("\nWorst 5 runs by relative RMSE:")
        for r in worst:
            print(f"  {r['run_name']}: rmse={r['rmse_rel']:.4f} r2={r['r2']:.4f} params={r['params']}")

        if args.save:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(6, 5))
            ax.hist(rmses, bins=30)
            ax.set_xlabel("Relative RMSE of predict_phi(t) vs actual phi(t)")
            ax.set_ylabel("Number of runs")
            ax.set_title(f"predict_phi validation across {len(rows)} runs")
            os.makedirs(os.path.join(ROOT_DIR, "figures"), exist_ok=True)
            fname = os.path.join(ROOT_DIR, "figures", "predict_pea_validation_rmse_hist.png")
            plt.savefig(fname)
            print(f"\nSaved validation histogram to {fname}")
        return

    parser.error("Provide --fit and/or --validate")


if __name__ == "__main__":
    main()
