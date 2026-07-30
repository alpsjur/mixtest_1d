#!/usr/bin/env python3
"""
tools/make_grd.py

Creates a ROMS-compatible grid NetCDF using a resolved configuration dict.

- Preferred usage (from orchestrator):
    from tools.make_grd import make_grid_from_config
    grid_path = make_grid_from_config(cfg_dict)

- CLI usage (for convenience):
    python tools/make_grd.py path/to/config.yaml
"""

import os
import sys
import numpy as np
import netCDF4 as nc
import yaml

# Ensure project root is importable when running from tools/
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import compute_z_r



def build_str_a(z_r, str_a_value, depth_zero_below, eta_limits=None, xi_limits=None):
    """
    Uniform str_a everywhere, zeroed below the given depth from the surface and outside specified limits.
    
    Parameters:
        z_r (np.ndarray): Vertical coordinate array.
        str_a_value (float): Structure area density value.
        depth_zero_below (float): Depth below which str_a is zero.
        eta_limits (list or tuple): [min_eta, max_eta] limits for eta dimension.
        xi_limits (list or tuple): [min_xi, max_xi] limits for xi dimension.
    """
    N_lev, eta_rho, xi_rho = z_r.shape
    str_a = np.full((N_lev, eta_rho, xi_rho), str_a_value, dtype=np.float64)

    # Apply depth condition
    dist_from_surface = -z_r  # positive; 0 at surface
    str_a[dist_from_surface > depth_zero_below] = 0.0

    # Apply eta and xi limits if provided
    if eta_limits:
        str_a[:, :eta_limits[0], :] = 0.0
        str_a[:, eta_limits[1] + 1:, :] = 0.0
    if xi_limits:
        str_a[:, :, :xi_limits[0]] = 0.0
        str_a[:, :, xi_limits[1] + 1:] = 0.0

    return str_a


def write_grid(output, Lm, Mm, N, Vtransform, Vstretching, THETA_S, THETA_B, HC,
               H0, DX, DY, F0, str_a_value, depth_zero_below, eta_limits=None, xi_limits=None):
    """
    Create a ROMS-compatible grid NetCDF at 'output'.
    """
    out_dir = os.path.dirname(output) or "."
    os.makedirs(out_dir, exist_ok=True)

    # ROMS staggered C-grid dimensions
    xi_rho  = Lm + 2
    eta_rho = Mm + 2
    xi_psi  = Lm + 1
    eta_psi = Mm + 1
    xi_u    = Lm + 1
    eta_u   = eta_rho
    xi_v    = xi_rho
    eta_v   = Mm + 1

    # Coordinates (uniform Cartesian grid)
    xi_idx  = np.arange(xi_rho,  dtype=float)
    eta_idx = np.arange(eta_rho, dtype=float)

    x_rho = (xi_idx[np.newaxis,  :] + 0.5) * DX
    y_rho = (eta_idx[:, np.newaxis] + 0.5) * DY

    x_psi = (np.arange(1, xi_rho,  dtype=float)) * DX
    y_psi = (np.arange(1, eta_rho, dtype=float)) * DY
    x_psi, y_psi = np.meshgrid(x_psi, y_psi)

    x_u   = (np.arange(1, xi_rho,  dtype=float)) * DX
    y_u   = y_rho[:, :xi_u].copy()
    x_u   = np.tile(x_u, (eta_u, 1))

    x_v   = x_rho[:eta_v, :].copy()
    y_v   = (np.arange(1, eta_rho, dtype=float)) * DY
    y_v   = np.tile(y_v[:, np.newaxis], (1, xi_v))

    # Scalars/metrics
    h      = np.full((eta_rho, xi_rho), H0,       dtype=np.float64)
    f      = np.full((eta_rho, xi_rho), F0,       dtype=np.float64)
    pm     = np.full((eta_rho, xi_rho), 1.0 / DX, dtype=np.float64)
    pn     = np.full((eta_rho, xi_rho), 1.0 / DY, dtype=np.float64)
    xl_val = float(Lm * DX)
    el_val = float(Mm * DY)

    # Vertical coordinates for str_a construction
    z_r = compute_z_r(h, HC, THETA_S, THETA_B, N)

    # Build str_a
    str_a = build_str_a(z_r, str_a_value, depth_zero_below, eta_limits, xi_limits)

    # Write NetCDF
    with nc.Dataset(output, "w", format="NETCDF4") as ds:
        ds.title       = "ROMS grid file (config-driven)"
        ds.type        = "ROMS grid file"
        ds.Conventions = "CF-1.6"

        ds.createDimension("xi_rho",  xi_rho)
        ds.createDimension("eta_rho", eta_rho)
        ds.createDimension("xi_psi",  xi_psi)
        ds.createDimension("eta_psi", eta_psi)
        ds.createDimension("xi_u",    xi_u)
        ds.createDimension("eta_u",   eta_u)
        ds.createDimension("xi_v",    xi_v)
        ds.createDimension("eta_v",   eta_v)
        ds.createDimension("s_rho",   N)
        ds.createDimension("one",     1)

        def scalar_var(name, value, long_name, units):
            v = ds.createVariable(name, "f8", ("one",))
            v.long_name = long_name
            v.units     = units
            v[:] = value

        scalar_var("xl", xl_val, "domain length in the XI-direction", "meter")
        scalar_var("el", el_val, "domain length in the ETA-direction", "meter")

        vs = ds.createVariable("spherical", "i4", ("one",))
        vs.long_name = "grid type logical switch"
        vs.flag_values   = [0, 1]
        vs.flag_meanings = "Cartesian spherical"
        vs[:] = 0

        def write2d(name, dims, data, long_name, units, dtype="f8"):
            v = ds.createVariable(name, dtype, dims)
            v.long_name = long_name
            v.units     = units
            v[:]        = data

        write2d("h",     ("eta_rho", "xi_rho"), h,    "bathymetry at RHO-points", "meter")
        write2d("f",     ("eta_rho", "xi_rho"), f,    "Coriolis parameter at RHO-points", "second-1")
        write2d("pm",    ("eta_rho", "xi_rho"), pm,   "metric in XI", "meter-1")
        write2d("pn",    ("eta_rho", "xi_rho"), pn,   "metric in ETA", "meter-1")
        write2d("x_rho", ("eta_rho", "xi_rho"), x_rho, "x-locations of RHO-points", "meter")
        write2d("y_rho", ("eta_rho", "xi_rho"), y_rho, "y-locations of RHO-points", "meter")
        write2d("x_psi", ("eta_psi", "xi_psi"), x_psi, "x-locations of PSI-points", "meter")
        write2d("y_psi", ("eta_psi", "xi_psi"), y_psi, "y-locations of PSI-points", "meter")
        write2d("x_u",   ("eta_u",   "xi_u"),   x_u,   "x-locations of U-points", "meter")
        write2d("y_u",   ("eta_u",   "xi_u"),   y_u,   "y-locations of U-points", "meter")
        write2d("x_v",   ("eta_v",   "xi_v"),   x_v,   "x-locations of V-points", "meter")
        write2d("y_v",   ("eta_v",   "xi_v"),   y_v,   "y-locations of V-points", "meter")

        v = ds.createVariable("str_a", "f8", ("s_rho", "eta_rho", "xi_rho"), fill_value=9.99e36)
        v.long_name = "structure area density"
        v.units     = "meter-1"
        v.valid_min = np.float64(0.0)
        v[:]        = str_a

    #print(f"Grid file written: {output}")


def make_grid_from_config(cfg: dict) -> str:
    """
    Create the grid using values from a resolved config dict.
    """
    input_dir = cfg["io"]["input_dir"]
    grd_name  = cfg["files"]["grd"]
    output = os.path.join(input_dir, grd_name)

    Lm = int(cfg["grid"]["Lm"])
    Mm = int(cfg["grid"]["Mm"])
    N  = int(cfg["grid"]["N"])
    H0 = float(cfg["grid"]["H0"])
    DX = float(cfg["grid"]["DX"])
    DY = float(cfg["grid"]["DY"])

    Vtransform  = int(cfg["vertical"]["Vtransform"])
    Vstretching = int(cfg["vertical"]["Vstretching"])
    THETA_S     = float(cfg["vertical"]["THETA_S"])
    THETA_B     = float(cfg["vertical"]["THETA_B"])
    HC          = float(cfg["vertical"]["HC"])

    F0 = float(cfg["phys"]["F0"])

    str_a_value      = float(cfg["structure"]["str_a"])
    depth_zero_below = float(cfg["structure"]["depth_zero_below"])

    # Get eta and xi limits from configuration
    eta_limits = cfg["structure"].get("eta_limits", None)
    xi_limits = cfg["structure"].get("xi_limits", None)

    write_grid(
        output=output,
        Lm=Lm, Mm=Mm, N=N,
        Vtransform=Vtransform, Vstretching=Vstretching,
        THETA_S=THETA_S, THETA_B=THETA_B, HC=HC,
        H0=H0, DX=DX, DY=DY, F0=F0,
        str_a_value=str_a_value,
        depth_zero_below=depth_zero_below,
        eta_limits=eta_limits,
        xi_limits=xi_limits,
    )

    return output


if __name__ == "__main__":
    # CLI fallback: accept a single config path, load YAML, and call the function
    if len(sys.argv) != 2:
        print("Usage: python tools/make_grd.py path/to/config.yaml", file=sys.stderr)
        sys.exit(1)

    config_path = sys.argv[1]
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)
    make_grid_from_config(cfg)