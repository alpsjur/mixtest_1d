#!/usr/bin/env python3
"""
tools/make_ini.py

Create a ROMS-compatible initial conditions NetCDF using a resolved config dict.

Preferred usage (from orchestrator):
    from tools.make_ini import make_ini_from_config
    ini_path = make_ini_from_config(cfg_dict)

CLI fallback:
    python tools/make_ini.py path/to/config.yaml
"""

import os
import sys
import numpy as np
import netCDF4 as nc
import yaml
from utils.utils import compute_z_r

# ---------------------------------------------------------------------------
# Initialization Functions
# ---------------------------------------------------------------------------

def zeta_initial(x_rho, y_rho, cfg):
    """
    Define the initial free surface (`zeta`) as a function of `x_rho` and `y_rho`.
    
    Args:
        x_rho: X-coordinates at RHO points.
        y_rho: Y-coordinates at RHO points.
        cfg: A dictionary of configuration parameters (e.g., zeta_slope).

    Returns:
        A 2D array representing the initial `zeta`.
    """
    return np.zeros_like(x_rho, dtype=np.float64) 

def ubar_initial(eta_u, xi_u, cfg):
    """
    Define the initial vertically integrated u-momentum (`ubar`).

    Args:
        eta_u: Number of points in the eta direction for U points.
        xi_u: Number of points in the xi direction for U points.
        cfg: A dictionary of configuration parameters (e.g., ubar_constant).

    Returns:
        A 2D array representing the initial `ubar`.
    """
    return np.zeros((eta_u, xi_u), dtype=np.float64)

def vbar_initial(eta_v, xi_v, cfg):
    """
    Define the initial vertically integrated v-momentum (`vbar`).

    Args:
        eta_v: Number of points in the eta direction for V points.
        xi_v: Number of points in the xi direction for V points.
        cfg: A dictionary of configuration parameters (e.g., vbar_constant).

    Returns:
        A 2D array representing the initial `vbar`.
    """
    return np.zeros((eta_v, xi_v), dtype=np.float64)

def u_initial(N, eta_u, xi_u, cfg):
    """
    Define the initial 3D u-momentum component (`u`).

    Args:
        N: Number of vertical levels.
        eta_u: Number of points in the eta direction for U points.
        xi_u: Number of points in the xi direction for U points.
        cfg: A dictionary of configuration parameters (e.g., u_constant).

    Returns:
        A 3D array representing the initial `u`.
    """
    return np.zeros((N, eta_u, xi_u), dtype=np.float64)

def v_initial(N, eta_v, xi_v, cfg):
    """
    Define the initial 3D v-momentum component (`v`).

    Args:
        N: Number of vertical levels.
        eta_v: Number of points in the eta direction for V points.
        xi_v: Number of points in the xi direction for V points.
        cfg: A dictionary of configuration parameters (e.g., v_constant).

    Returns:
        A 3D array representing the initial `v`.
    """
    return np.zeros((N, eta_v, xi_v), dtype=np.float64)

def temp_initial(z_r, cfg):
    """
    Define the initial temperature profile (`temp`) as a function of depth
    using the hyperbolic tangent model.

    Args:
        z_r: 3D array of vertical coordinates at RHO points.
        cfg: A dictionary of configuration parameters (e.g., temp_T0, temp_dT, temp_zt, temp_ht).

    Returns:
        A 3D array representing the initial `temp`.
    """
    # Retrieve parameters with defaults
    temp_T0 = cfg["initial"]["temp_T0"]  # Surface temperature (T_s)
    temp_dT = cfg["initial"]["temp_dT"]  # Total temperature change across the thermocline (ΔT)
    temp_zt = cfg["initial"]["temp_zt"]  # Depth of the thermocline center (z_t)
    temp_ht = cfg["initial"]["temp_ht"]   # Thickness of the thermocline (h_t)

    # Hyperbolic tangent temperature profile
    return (temp_T0-temp_dT) + (temp_dT / 2.0) * (1 + np.tanh((z_r + temp_zt) / temp_ht))

def salt_initial(z_r, cfg):
    """
    Define the initial salinity profile (`salt`) as a function of depth.

    Args:
        z_r: 3D array of vertical coordinates at RHO points.
        cfg: A dictionary of configuration parameters (e.g., salt_S0).

    Returns:
        A 3D array representing the initial `salt`.
    """
    salt_S0 = cfg["initial"]["salt_S0"]
    return np.full_like(z_r, salt_S0, dtype=np.float64)

# ---------------------------------------------------------------------------
# Main Initial Conditions File Creation Function
# ---------------------------------------------------------------------------

def make_ini_from_config(cfg: dict) -> str:
    """
    Create the initial conditions file using values from a resolved config dict.
    """
    input_dir = cfg["io"]["input_dir"]
    grd_name  = cfg["files"]["grd"]
    ini_name  = cfg["files"]["ini"]

    grd_path = os.path.join(input_dir, grd_name)
    ini_path = os.path.join(input_dir, ini_name)
    os.makedirs(os.path.dirname(ini_path) or ".", exist_ok=True)

    # Read required config values
    N = int(cfg["grid"]["N"])

    Vtransform  = int(cfg["vertical"]["Vtransform"])
    Vstretching = int(cfg["vertical"]["Vstretching"])
    THETA_S     = float(cfg["vertical"]["THETA_S"])
    THETA_B     = float(cfg["vertical"]["THETA_B"])
    HC          = float(cfg["vertical"]["HC"])

    ocean_time_seconds = float(cfg["initial"]["ocean_time_seconds"])

    # Read grid geometry
    with nc.Dataset(grd_path, "r") as grd:
        h     = grd.variables["h"][:]
        x_rho = grd.variables["x_rho"][:]
        y_rho = grd.variables["y_rho"][:]
        x_u   = grd.variables["x_u"][:]
        y_u   = grd.variables["y_u"][:]
        x_v   = grd.variables["x_v"][:]
        y_v   = grd.variables["y_v"][:]

        xi_rho  = len(grd.dimensions["xi_rho"])
        eta_rho = len(grd.dimensions["eta_rho"])
        xi_u    = xi_rho - 1
        eta_u   = eta_rho
        xi_v    = xi_rho
        eta_v   = eta_rho - 1

    # Compute vertical coordinates at RHO-points
    z_r = compute_z_r(h, HC, THETA_S, THETA_B, N)  # shape: (N, eta_rho, xi_rho)

    # Interpolate z_r to staggered points
    z_r_u = 0.5 * (z_r[:, :, :-1] + z_r[:, :, 1:])  # (N, eta_u, xi_u)
    z_r_v = 0.5 * (z_r[:, :-1, :] + z_r[:, 1:, :])  # (N, eta_v, xi_v)

    # Allocate initial fields using parameterized functions
    zeta = zeta_initial(x_rho, y_rho, cfg)
    ubar = ubar_initial(eta_u, xi_u, cfg)
    vbar = vbar_initial(eta_v, xi_v, cfg)
    u_3d = u_initial(N, eta_u, xi_u, cfg)
    v_3d = v_initial(N, eta_v, xi_v, cfg)
    temp = temp_initial(z_r, cfg)
    salt = salt_initial(z_r, cfg)

    # Write initial conditions NetCDF
    with nc.Dataset(ini_path, "w", format="NETCDF4") as f:
        # Global attributes
        f.title = "ROMS Initial Conditions (parameterized)"
        f.history = "Created by tools/make_ini.py"
        f.description = "Initial conditions with parameterized functions for initialization"
        f.source = "Generated from grid file"

        # Dimensions
        f.createDimension("xi_rho", xi_rho)
        f.createDimension("eta_rho", eta_rho)
        f.createDimension("xi_u",   xi_u)
        f.createDimension("eta_u",  eta_u)
        f.createDimension("xi_v",   xi_v)
        f.createDimension("eta_v",  eta_v)
        f.createDimension("s_rho",  N)
        f.createDimension("ocean_time", None)

        # ocean_time
        ocean_time = f.createVariable("ocean_time", "f8", ("ocean_time",))
        ocean_time.long_name = "time since simulation start"
        ocean_time.units = "seconds since 0001-01-01 00:00:00"
        ocean_time.calendar = "360.0 days in every year"
        ocean_time[0] = ocean_time_seconds

        # Variables
        f.createVariable("zeta", "f8", ("ocean_time", "eta_rho", "xi_rho"))[0, :, :] = zeta
        f.createVariable("ubar", "f8", ("ocean_time", "eta_u", "xi_u"))[0, :, :] = ubar
        f.createVariable("vbar", "f8", ("ocean_time", "eta_v", "xi_v"))[0, :, :] = vbar
        f.createVariable("u", "f8", ("ocean_time", "s_rho", "eta_u", "xi_u"))[0, :, :, :] = u_3d
        f.createVariable("v", "f8", ("ocean_time", "s_rho", "eta_v", "xi_v"))[0, :, :, :] = v_3d
        f.createVariable("temp", "f8", ("ocean_time", "s_rho", "eta_rho", "xi_rho"))[0, :, :, :] = temp
        f.createVariable("salt", "f8", ("ocean_time", "s_rho", "eta_rho", "xi_rho"))[0, :, :, :] = salt

    #print(f"Initial conditions file written: {ini_path}")
    return ini_path


if __name__ == "__main__":
    # CLI fallback: accept a single config path
    if len(sys.argv) != 2:
        print("Usage: python tools/make_ini.py path/to/config.yaml", file=sys.stderr)
        sys.exit(1)
    with open(sys.argv[1], "r") as f:
        cfg = yaml.safe_load(f)
    make_ini_from_config(cfg)