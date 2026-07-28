import xarray as xr
import yaml
import numpy as np
import xgcm
import os 
import sys

# Ensure project root is importable when running from tools/
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import compute_depths

def load_yaml(path: str) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)

def prep_ds(ds, params):
    # Rename dimensions to match expected naming conventions for xgcm
    ds = ds.rename({'eta_u': 'eta_rho', 'xi_v': 'xi_rho', 'xi_psi': 'xi_u', 'eta_psi': 'eta_v'})

    # Define the coordinates dictionary for xgcm grid object
    coords = {
        'X': {'center': 'xi_rho', 'inner': 'xi_u'}, 
        'Y': {'center': 'eta_rho', 'inner': 'eta_v'}, 
        'Z': {'center': 's_rho', 'outer': 's_w'}
    }

    # Create the grid object using xgcm
    grid = xgcm.Grid(ds, coords=coords, autoparse_metadata=False, padding='periodic')

    # Calculate the vertical coordinates (z-levels) at rho and w points
    z_r, z_w = compute_depths(ds.isel(ocean_time=0).h.squeeze().values, 
                              params['vertical']['HC'], 
                              params['vertical']['THETA_S'], 
                              params['vertical']['THETA_B'], 
                              params['grid']['N']
                              )

    # Add the calculated z-coordinates to the dataset
    ds.coords['z_w'] = (('s_w', 'eta_rho', 'xi_rho'), z_w)
    ds.coords['z_rho'] = (('s_rho', 'eta_rho', 'xi_rho'), z_r)

    # Add calculated z-coordinates to dataset
    #ds.coords['z_w'] = z_w.where(ds.mask_rho, 0).transpose('ocean_time', 's_w', 'eta_rho', 'xi_rho')
    #ds.coords['z_rho'] = z_r.where(ds.mask_rho, 0).transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho')

    # Interpolate grid metrics to u, v, and psi points
    ds['pm_v'] = grid.interp(ds.pm, 'Y')
    ds['pn_u'] = grid.interp(ds.pn, 'X')
    ds['pm_u'] = grid.interp(ds.pm, 'X')
    ds['pn_v'] = grid.interp(ds.pn, 'Y')
    ds['pm_psi'] = grid.interp(grid.interp(ds.pm, 'Y'), 'X')  # Interpolated to psi points
    ds['pn_psi'] = grid.interp(grid.interp(ds.pn, 'X'), 'Y')  # Interpolated to psi points

    # Calculate grid spacings (dx, dy) at various grid points
    ds['dx'] = 1 / ds.pm
    ds['dx_u'] = 1 / ds.pm_u
    ds['dx_v'] = 1 / ds.pm_v
    ds['dx_psi'] = 1 / ds.pm_psi

    ds['dy'] = 1 / ds.pn
    ds['dy_u'] = 1 / ds.pn_u
    ds['dy_v'] = 1 / ds.pn_v
    ds['dy_psi'] = 1 / ds.pn_psi

    # Calculate vertical grid spacing differences
    ds['dz'] = grid.diff(ds.z_w, 'Z', padding='fill')
    ds['dz_w'] = grid.diff(ds.z_rho, 'Z', padding='fill')
    ds['dz_u'] = grid.interp(ds.dz, 'X')
    ds['dz_w_u'] = grid.interp(ds.dz_w, 'X')
    ds['dz_v'] = grid.interp(ds.dz, 'Y')
    ds['dz_w_v'] = grid.interp(ds.dz_w, 'Y')

    # Calculate grid cell areas
    ds['dA'] = ds.dx * ds.dy

    # Define metrics for xgcm grid object
    metrics = {
        ('X',): ['dx', 'dx_u', 'dx_v', 'dx_psi'],  # X distances
        ('Y',): ['dy', 'dy_u', 'dy_v', 'dy_psi'],  # Y distances
        ('Z',): ['dz', 'dz_u', 'dz_v', 'dz_w', 'dz_w_u', 'dz_w_v'],  # Z distances
        ('X', 'Y'): ['dA']  # Areas
    }

    # Re-create the grid object with the new metrics
    grid = xgcm.Grid(ds, coords=coords, metrics=metrics, padding='periodic', autoparse_metadata=False)

    return ds, grid