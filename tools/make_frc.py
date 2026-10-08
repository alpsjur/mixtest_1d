#!/usr/bin/env python3
"""
tools/make_frc.py

Create a ROMS-compatible surface wind-stress forcing NetCDF for mixtest_1d,
adapted from mixtest_3d/tools/make_frc.py.

Used for the structure-mixing / surface-mixing synergy experiment (see
notes/surface_mixing_synergy.md): structure-induced mixing in this repo so
far has always been the *only* source of turbulence (the body force alone
drives a uniform, shear-free flow -- see
notes/floating_structure_depth_distribution.md Sec. 2). Adding a real
surface wind stress gives this idealised column a second, independent TKE
source so the GLS closure's nonlinear response to *combined* forcing can be
tested, which is what this file enables.

Preferred usage (from the orchestrator):
    from tools.make_frc import make_frc_from_config
    frc_path = make_frc_from_config(cfg_dict)

CLI fallback:
    python tools/make_frc.py path/to/resolved_config.yaml
"""

import os
import sys
import numpy as np
import netCDF4 as nc
import yaml

# ---------------------------------------------------------------------------
# Surface stress functions
# ---------------------------------------------------------------------------

def sustr_forcing(frc_time, eta_u, xi_u, cfg):
    """
    Surface zonal wind stress on u-points, in *kinematic* units (m2/s2,
    i.e. tau/rho0) to match ROMS' sustr/svstr convention (see
    ana_smflux.h) -- NOT N/m^2 as in mixtest_3d's version.

    cfg["forcing"]["sustr_tau0"] is the (constant in time and space)
    kinematic stress; a positive value drives an eastward surface current.
    """
    n_time = len(frc_time)
    tau0 = cfg.get("forcing", {}).get("sustr_tau0", 0.0)
    return np.full((n_time, eta_u, xi_u), tau0, dtype=np.float64)


def svstr_forcing(frc_time, eta_v, xi_v, cfg):
    """Surface meridional wind stress on v-points (kinematic, m2/s2)."""
    n_time = len(frc_time)
    tau0 = cfg.get("forcing", {}).get("svstr_tau0", 0.0)
    return np.full((n_time, eta_v, xi_v), tau0, dtype=np.float64)


# ---------------------------------------------------------------------------
# Main forcing file creation function
# ---------------------------------------------------------------------------

def make_frc_from_config(cfg: dict) -> str:
    """
    Create the surface forcing file using values from a resolved config dict.

    Always generates a file (with sustr=svstr=0.0 if cfg["forcing"] is
    absent or its tau0 values are unset) -- since ANA_SMFLUX is undefined
    in this application's compiled executable (roms-related/Include/
    mixtest_1d.h), ROMS reads sustr/svstr from FRCNAME unconditionally, so
    every run needs this file to exist, not just the ones that want a
    non-zero wind stress.

    The forcing has its own time axis defined by cfg["forcing"]["dt_frc"]
    (seconds); ROMS linearly interpolates between snapshots at run time.
    A constant-in-time stress only needs two snapshots (start and end).
    """
    input_dir = cfg["io"]["input_dir"]
    grd_name = cfg["files"]["grd"]
    frc_name = cfg["files"].get("frc", "mixtest_1d_frc.nc")

    grd_path = os.path.join(input_dir, grd_name)
    frc_path = os.path.join(input_dir, frc_name)
    os.makedirs(os.path.dirname(frc_path) or ".", exist_ok=True)

    DT = float(cfg["time_stepping"]["DT"])
    NTIMES = int(cfg["time_stepping"]["NTIMES"])
    t_start = 0.0
    t_end = NTIMES * DT
    dt_frc = float(cfg.get("forcing", {}).get("dt_frc", t_end))  # default: one interval spanning the whole run
    frc_time = np.arange(t_start, t_end + dt_frc, dt_frc)

    with nc.Dataset(grd_path, "r") as grd:
        xi_rho = len(grd.dimensions["xi_rho"])
        eta_rho = len(grd.dimensions["eta_rho"])

    xi_u, eta_u = xi_rho - 1, eta_rho
    xi_v, eta_v = xi_rho, eta_rho - 1

    sustr = sustr_forcing(frc_time, eta_u, xi_u, cfg)
    svstr = svstr_forcing(frc_time, eta_v, xi_v, cfg)

    with nc.Dataset(frc_path, "w", format="NETCDF4") as f:
        f.title = "ROMS Surface Forcing (mixtest_1d)"
        f.history = "Created by tools/make_frc.py"
        f.description = "Surface wind stress forcing (sustr, svstr), kinematic (m2/s2)"

        f.createDimension("xi_rho", xi_rho)
        f.createDimension("eta_rho", eta_rho)
        f.createDimension("xi_u", xi_u)
        f.createDimension("eta_u", eta_u)
        f.createDimension("xi_v", xi_v)
        f.createDimension("eta_v", eta_v)
        f.createDimension("sms_time", None)

        ot = f.createVariable("sms_time", "f8", ("sms_time",))
        ot.long_name = "time since simulation start"
        ot.units = "seconds since 0001-01-01 00:00:00"
        ot.calendar = "360.0 days in every year"
        ot[:] = frc_time

        su = f.createVariable("sustr", "f8", ("sms_time", "eta_u", "xi_u"))
        su.long_name = "surface u-momentum stress"
        su.units = "meter2 second-2"
        su[:] = sustr

        sv = f.createVariable("svstr", "f8", ("sms_time", "eta_v", "xi_v"))
        sv.long_name = "surface v-momentum stress"
        sv.units = "meter2 second-2"
        sv[:] = svstr

    return frc_path


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python tools/make_frc.py path/to/resolved_config.yaml", file=sys.stderr)
        sys.exit(1)
    with open(sys.argv[1], "r") as fh:
        cfg = yaml.safe_load(fh)
    path = make_frc_from_config(cfg)
    print(path)
