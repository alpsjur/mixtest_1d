#!/usr/bin/env python3
# prep_profiles.py
import os
import sys
import argparse
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Make project utils importable (adjust paths if needed)
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import prep_ds, load_yaml

def prep_profiles(resolved_config_path: str, var: str, timestep: int = -1, subtract_ini: bool = False):

    # Load configuration
    params = load_yaml(resolved_config_path)
    out_dir = params["io"]["output_dir"]
    his_file = params["files"]["his"]
    ds_path = os.path.join(out_dir, his_file)

    # Open dataset; use_cftime to avoid warnings for non-proleptic Gregorian calendars
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds = xr.open_dataset(ds_path, decode_times=time_coder)

    # Prepare dataset and grid (project-specific)
    ds, grid = prep_ds(ds, params)

    if var not in ds:
        ds.close()
        raise KeyError(f"Variable '{var}' not found in dataset.")

    da = ds[var]

    # Select the desired timestep
    da = da.isel(ocean_time=timestep)

    # average over horizontal dimensions (X, Y)
    da_avg = grid.average(da, axis=("X", "Y")).squeeze()
    if subtract_ini:
        da_ini = grid.average(ds[var].isel(ocean_time=0), axis=("X", "Y")).squeeze()
        da_avg = da_avg - da_ini

    # check which vertical dimension is present and return it along with the averaged data
    if "s_rho" in da_avg.dims:
        z = grid.average(ds.z_rho, axis=("X", "Y")).squeeze()
    elif "s_w" in da_avg.dims:
        z = grid.average(ds.z_w, axis=("X", "Y")).squeeze()
    else:
        raise ValueError("No recognized vertical dimension found in the data array.")
    return da_avg, z