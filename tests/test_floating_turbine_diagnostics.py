#!/usr/bin/env python3
"""
tests/test_floating_turbine_diagnostics.py

Verification tests for floating turbine foundation modifications in:
- utils.utils (analytic_mixing_timescale, compute_Pd, compute_Pstr, compute_phi)
- analysis.mixing_timescale (pycnocline_length_scale)
"""

import os
import sys
import numpy as np
import xarray as xr
import xgcm

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import (
    analytic_mixing_timescale,
    compute_Pd,
    compute_Pstr,
    compute_phi,
    RHO0,
    R0,
    TCOEF,
    G,
)
from analysis.mixing_timescale import pycnocline_length_scale


def test_analytic_mixing_timescale():
    # Base bottom-fixed config
    cfg_bf = {
        "structure": {"CD": 0.63, "str_a": 0.01, "depth_zero_below": 1e9, "c4": 0.44},
        "bodyforce": {"BFRC_U": 315.0, "BFRC_V": 0.0},
        "grid": {"H0": 150.0},
        "initial": {"temp_dT": 10.0, "temp_zt": 30.0},
        "GLS": {"C1": 1.0},
    }
    res_bf = analytic_mixing_timescale(cfg_bf)
    assert res_bf["Hturb"] == 150.0
    u_inf = np.sqrt(2.0 * 3.15e-5 / (0.63 * 0.01))
    Pd = 0.5 * 0.63 * 0.01 * u_inf**3
    Pstr_bf = RHO0 * 150.0 * Pd
    assert np.isclose(res_bf["Pstr"], Pstr_bf)

    # Floating config with Hturb = 50 m
    cfg_fl = {
        "structure": {"CD": 0.63, "str_a": 0.01, "depth_zero_below": 50.0, "c4": 0.44},
        "bodyforce": {"BFRC_U": 315.0, "BFRC_V": 0.0},
        "grid": {"H0": 150.0},
        "initial": {"temp_dT": 10.0, "temp_zt": 30.0},
        "GLS": {"C1": 1.0},
    }
    res_fl = analytic_mixing_timescale(cfg_fl)
    assert res_fl["Hturb"] == 50.0
    assert np.isclose(res_fl["Pstr"], RHO0 * 50.0 * Pd)
    assert np.isclose(res_fl["Pstr"], res_bf["Pstr"] * (50.0 / 150.0))

    # Test tau_mix ratio: tau_mix = g * delta_rho * Hturb^2 / Pstr = g * delta_rho * Hturb / (RHO0 * Pd)
    # So tau_mix_fl / tau_mix_bf = Hturb / H0 = 50 / 150
    assert np.isclose(res_fl["tau_mix"] / res_bf["tau_mix"], 50.0 / 150.0)

    # Check validation error when temp_zt >= Hturb
    cfg_invalid = {
        "structure": {"CD": 0.63, "str_a": 0.01, "depth_zero_below": 30.0, "c4": 0.44},
        "bodyforce": {"BFRC_U": 315.0, "BFRC_V": 0.0},
        "grid": {"H0": 150.0},
        "initial": {"temp_dT": 10.0, "temp_zt": 30.0},
        "GLS": {"C1": 1.0},
    }
    try:
        analytic_mixing_timescale(cfg_invalid)
        assert False, "Should have raised ValueError for temp_zt >= Hturb"
    except ValueError:
        pass
    print("test_analytic_mixing_timescale: PASS")


def test_pycnocline_length_scale():
    # Bottom-fixed
    p_bf = {
        "grid": {"H0": 150.0},
        "structure": {"depth_zero_below": 1e9},
        "initial": {"temp_zt": 35.0},
    }
    L_bf = pycnocline_length_scale(p_bf)
    assert np.isclose(L_bf, np.sqrt(35.0 * (150.0 - 35.0)))

    # Floating: Hturb = 60 m, zt = 25 m
    p_fl = {
        "grid": {"H0": 150.0},
        "structure": {"depth_zero_below": 60.0},
        "initial": {"temp_zt": 25.0},
    }
    L_fl = pycnocline_length_scale(p_fl)
    assert np.isclose(L_fl, np.sqrt(25.0 * (60.0 - 25.0)))

    # Check error when zt >= Hturb
    p_err = {
        "grid": {"H0": 150.0},
        "structure": {"depth_zero_below": 40.0},
        "initial": {"temp_zt": 45.0},
    }
    try:
        pycnocline_length_scale(p_err)
        assert False, "Should have raised ValueError for zt >= Hturb"
    except ValueError:
        pass
    print("test_pycnocline_length_scale: PASS")


def _build_synthetic_dataset_and_grid(H0=150.0, N=30, u_val=0.1):
    """Create a 1D synthetic column dataset with xgcm Grid."""
    s_rho = np.arange(N)
    s_w = np.arange(N + 1)
    eta_rho = np.arange(3)
    xi_rho = np.arange(3)
    xi_u = np.arange(2)
    eta_v = np.arange(2)

    # Linear vertical coordinates
    z_w_1d = np.linspace(-H0, 0, N + 1)
    z_rho_1d = 0.5 * (z_w_1d[:-1] + z_w_1d[1:])
    dz_1d = np.diff(z_w_1d)

    # 3D arrays (s_rho, eta_rho, xi_rho)
    shape_rho = (N, len(eta_rho), len(xi_rho))
    z_rho_3d = np.broadcast_to(z_rho_1d[:, None, None], shape_rho).copy()
    dz_3d = np.broadcast_to(dz_1d[:, None, None], shape_rho).copy()

    # Metrics
    dx = 100.0
    dy = 100.0
    dA = dx * dy
    dV = dA * dz_3d

    # Synthetic density with pycnocline at z = -30 m
    # Stratification: delta_rho ~ 2 kg/m3 across z = -30 m
    rho_1d = 1025.0 - 1.0 * np.tanh((z_rho_1d + 30.0) / 10.0)
    rho_3d = np.broadcast_to(rho_1d[:, None, None], shape_rho).copy()

    ds = xr.Dataset(
        data_vars={
            "u": (("ocean_time", "s_rho", "eta_rho", "xi_u"), np.full((2, N, len(eta_rho), len(xi_u)), u_val)),
            "v": (("ocean_time", "s_rho", "eta_v", "xi_rho"), np.zeros((2, N, len(eta_v), len(xi_rho)))),
            "rho": (("ocean_time", "s_rho", "eta_rho", "xi_rho"), np.broadcast_to(rho_3d, (2, *shape_rho)).copy()),
            "z_rho": (("s_rho", "eta_rho", "xi_rho"), z_rho_3d),
            "dz": (("s_rho", "eta_rho", "xi_rho"), dz_3d),
            "dA": (("eta_rho", "xi_rho"), np.full((len(eta_rho), len(xi_rho)), dA)),
            "dV": (("s_rho", "eta_rho", "xi_rho"), dV),
            "pm": (("eta_rho", "xi_rho"), np.full((len(eta_rho), len(xi_rho)), 1.0 / dx)),
            "pn": (("eta_rho", "xi_rho"), np.full((len(eta_rho), len(xi_rho)), 1.0 / dy)),
        },
        coords={
            "ocean_time": [0.0, 3600.0],
            "s_rho": s_rho,
            "s_w": s_w,
            "eta_rho": eta_rho,
            "xi_rho": xi_rho,
            "xi_u": xi_u,
            "eta_v": eta_v,
        },
    )

    coords = {
        "X": {"center": "xi_rho", "inner": "xi_u"},
        "Y": {"center": "eta_rho", "inner": "eta_v"},
        "Z": {"center": "s_rho", "outer": "s_w"},
    }
    metrics = {
        ("Z",): ["dz"],
        ("X", "Y"): ["dA"],
    }
    grid = xgcm.Grid(ds, coords=coords, metrics=metrics, padding='periodic', autoparse_metadata=False)
    return ds, grid


def test_compute_Pd_and_Pstr():
    ds, grid = _build_synthetic_dataset_and_grid(H0=150.0, N=30, u_val=0.1)

    # 1. Bottom-fixed
    params_bf = {
        "structure": {"CD": 0.63, "str_a": 0.01, "depth_zero_below": 1e9},
        "grid": {"H0": 150.0},
    }
    Pd_bf = compute_Pd(ds, grid, params_bf)
    expected_Pd = 0.5 * 0.63 * 0.01 * (0.1**3)
    assert np.allclose(Pd_bf.values, expected_Pd)

    Pstr_bf = compute_Pstr(ds, grid, params_bf)
    expected_Pstr_bf = RHO0 * 150.0 * expected_Pd
    assert np.isclose(float(Pstr_bf.values[0]), expected_Pstr_bf, rtol=1e-3)

    # 2. Floating turbine with depth_zero_below = 50 m
    params_fl = {
        "structure": {"CD": 0.63, "str_a": 0.01, "depth_zero_below": 50.0},
        "grid": {"H0": 150.0},
    }
    Pd_fl = compute_Pd(ds, grid, params_fl)
    # Check that Pd_fl is 0 for depth > 50 m
    z_rho = ds["z_rho"].values[:, 0, 0]
    for k, z in enumerate(z_rho):
        depth_from_surface = -z
        if depth_from_surface > 50.0:
            assert np.allclose(Pd_fl.values[:, k, :, :], 0.0)
        else:
            assert np.allclose(Pd_fl.values[:, k, :, :], expected_Pd)

    Pstr_fl = compute_Pstr(ds, grid, params_fl)
    # Vertically integrated Pstr should be close to RHO0 * 50.0 * expected_Pd
    expected_Pstr_fl = RHO0 * 50.0 * expected_Pd
    assert np.isclose(float(Pstr_fl.values[0]), expected_Pstr_fl, rtol=0.05)
    print("test_compute_Pd_and_Pstr: PASS")


def test_compute_phi_floating():
    ds, grid = _build_synthetic_dataset_and_grid(H0=150.0, N=30, u_val=0.1)

    params_fl = {
        "structure": {"depth_zero_below": 60.0},
        "grid": {"H0": 150.0},
    }
    phi, rho_mix0 = compute_phi(ds, grid, params_fl)

    # Initial stratified state should have phi > 0
    assert float(phi.values[0]) > 0.0

    # If we artificially homogenize the upper 60 m in timestep 1:
    ds_homo = ds.copy(deep=True)
    mask = (-ds_homo["z_rho"] <= 60.0).values
    # Replace upper layer rho with its mean
    rho_mean = ds_homo["rho"].values[0, mask].mean()
    ds_homo["rho"].values[1, mask] = rho_mean

    phi_homo, _ = compute_phi(ds_homo, grid, params_fl)
    # phi at timestep 1 should be essentially 0
    assert np.isclose(float(phi_homo.values[1]), 0.0, atol=1e-3)
    print("test_compute_phi_floating: PASS")


if __name__ == "__main__":
    test_analytic_mixing_timescale()
    test_pycnocline_length_scale()
    test_compute_Pd_and_Pstr()
    test_compute_phi_floating()
    print("ALL FLOATING TURBINE DIAGNOSTIC TESTS PASSED!")
