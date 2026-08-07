#!/usr/bin/env python3
# prep_profiles.py
import os
import sys

# Make project utils importable
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import open_roms_dataset


def prep_profiles(resolved_config_path: str, var: str, timestep: int = -1, subtract_ini: bool = False):
    """
    Compute the horizontally averaged vertical profile of a variable.

    Parameters
    ----------
    resolved_config_path : str
        Path to a resolved_config.yaml produced by prep_experiment.
    var : str
        Variable name to read from the history file.
    timestep : int
        Ocean time index to use (default -1 = last timestep).
    subtract_ini : bool
        If True, subtract the initial (timestep=0) profile so the result
        shows the change relative to initial conditions.

    Returns
    -------
    da_avg : DataArray
        Horizontally averaged profile at the selected timestep.
    z : DataArray
        Corresponding depth coordinate (m, negative down).
    """
    ds, grid, params = open_roms_dataset(resolved_config_path)

    if var not in ds:
        ds.close()
        raise KeyError(f"Variable '{var}' not found in dataset.")

    da = grid.average(ds[var].isel(ocean_time=timestep), axis=("X", "Y")).squeeze()

    if subtract_ini:
        da_ini = grid.average(ds[var].isel(ocean_time=0), axis=("X", "Y")).squeeze()
        da = da - da_ini

    # Return the matching depth coordinate based on the vertical dimension present
    if "s_rho" in da.dims:
        z = grid.average(ds.z_rho, axis=("X", "Y")).squeeze()
    elif "s_w" in da.dims:
        z = grid.average(ds.z_w, axis=("X", "Y")).squeeze()
    else:
        raise ValueError("No recognized vertical dimension found in the data array.")

    return da, z


def main():
    import argparse
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description="Plot horizontally averaged vertical profile of a variable.")
    parser.add_argument("resolved_config", type=str, help="Path to resolved_config.yaml")
    parser.add_argument("variable", type=str, help="Variable name to plot")
    parser.add_argument("--timestep", type=int, default=-1, help="Ocean time index to use (default: last timestep)")
    parser.add_argument("--subtract_ini", action="store_true", help="Subtract initial profile (timestep=0)")

    args = parser.parse_args()

    da_avg, z = prep_profiles(args.resolved_config, args.variable, args.timestep, args.subtract_ini)

    plt.figure(figsize=(6, 8))
    plt.plot(da_avg, z)
    plt.gca().invert_yaxis()
    plt.xlabel(f"{args.variable} ({da_avg.units})")
    plt.ylabel("Depth (m)")
    plt.title(f"Horizontally averaged profile of {args.variable}")
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()