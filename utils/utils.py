# python
import os
import numpy as np
import yaml
import xarray as xr
import xgcm

def compute_stretching(theta_s, theta_b, N):
    """
    Compute vertical stretching curves (s_rho, Cs_r) and (s_w, Cs_w)
    for Vstretching=5 (Souza et al., 2015).

    Parameters
    ----------
    theta_s : float
        Surface control parameter (0 <= theta_s <= 10 recommended)
    theta_b : float
        Bottom control parameter (0 <= theta_b <= 4 recommended)
    N : int
        Number of vertical rho levels (N >= 2)

    Returns
    -------
    s_rho : ndarray (N,), normalized S-coordinates at rho-points
    Cs_r : ndarray (N,), stretching curves at rho-points
    s_w : ndarray (N+1,), normalized S-coordinates at w-points
    Cs_w : ndarray (N+1,), stretching curves at w-points
    """
    if N < 2:
        raise ValueError("N must be >= 2 for the quadratic Legendre sigma definition.")

    rN = float(N)

    # Indices
    k_w = np.arange(0, N + 1, dtype=float)        # W-points: k = 0..N
    k_r = np.arange(1, N + 1, dtype=float) - 0.5  # Rho-points: k-0.5, with k = 1..N

    # Sigma at W-points 
    s_w = -(k_w**2 - 2.0 * k_w * rN + k_w + rN**2 - rN) / (rN**2 - rN) \
          - 0.01 * (k_w**2 - k_w * rN) / (1.0 - rN)

    # Sigma at Rho-points (note the 0.5 shift)
    s_rho = -(k_r**2 - 2.0 * k_r * rN + k_r + rN**2 - rN) / (rN**2 - rN) \
            - 0.01 * (k_r**2 - k_r * rN) / (1.0 - rN)

    # Surface refinement C_sur
    if theta_s > 0.0:
        Csur_rho = (1.0 - np.cosh(theta_s * s_rho)) / (np.cosh(theta_s) - 1.0)
        Csur_w   = (1.0 - np.cosh(theta_s * s_w))   / (np.cosh(theta_s) - 1.0)
    else:
        Csur_rho = -(s_rho ** 2)
        Csur_w   = -(s_w ** 2)

    # Bottom refinement (second stretching)
    if theta_b > 0.0:
        denom = 1.0 - np.exp(-theta_b)  # as in the provided description
        Cs_r = (np.exp(theta_b * Csur_rho) - 1.0) / denom
        Cs_w = (np.exp(theta_b * Csur_w)   - 1.0) / denom
    else:
        Cs_r = Csur_rho
        Cs_w = Csur_w

    return s_rho, Cs_r, s_w, Cs_w


def _vtransform2_depths(h, hc, s, Cs):
    """
    Compute depths for Vtransform=2 with zeta=0:
        z = h * (hc*s + h*Cs) / (hc + h)

    Parameters
    ----------
    h : ndarray (eta_rho, xi_rho), positive downward [m]
    hc : scalar
    s : ndarray (K,), normalized S-coordinates (rho or w)
    Cs : ndarray (K,), stretching function (rho or w)

    Returns
    -------
    z : ndarray (K, eta_rho, xi_rho), negative (sea level = 0)
        Index [0] is the bottom-most level, [K-1] the top-most.
    """
    h3d = h[np.newaxis, :, :]
    s3d = s[:, np.newaxis, np.newaxis]
    Cs3d = Cs[:, np.newaxis, np.newaxis]
    return h3d * (hc * s3d + h3d * Cs3d) / (hc + h3d)


def compute_depths(h, hc, theta_s, theta_b, N):
    """
    Compute depths at rho- and w-levels (z_r, z_w) for Vtransform=2, zeta=0.

    Parameters
    ----------
    h : ndarray (eta_rho, xi_rho), positive downward [m]
    hc, theta_s, theta_b, N : scalar ROMS parameters

    Returns
    -------
    z_r : ndarray (N, eta_rho, xi_rho)
        Negative (sea level = 0). Index [0] is bottom-most (k=1 in Fortran),
        index [N-1] is top-most (k=N in Fortran).
    z_w : ndarray (N+1, eta_rho, xi_rho)
        W-level depths with same ordering (bottom to surface).
    """
    s_rho, Cs_r, s_w, Cs_w = compute_stretching(theta_s, theta_b, N)
    z_r = _vtransform2_depths(h, hc, s_rho, Cs_r)
    z_w = _vtransform2_depths(h, hc, s_w, Cs_w)
    return z_r, z_w


# Thin wrapper for backward-compatibility
def compute_z_r(h, hc, theta_s, theta_b, N):
    z_r, _ = compute_depths(h, hc, theta_s, theta_b, N)
    return z_r

def compute_z_w(h, hc, theta_s, theta_b, N):
    _, z_w = compute_depths(h, hc, theta_s, theta_b, N)
    return z_w



def load_yaml(path: str) -> dict:
    """Load a YAML file and return its contents as a dict."""
    with open(path, "r") as f:
        return yaml.safe_load(f)


def save_yaml(path: str, data: dict) -> None:
    """Write a dict to a YAML file."""
    with open(path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def ensure_dir(path: str) -> None:
    """Create a directory (and any parents) if it does not already exist."""
    os.makedirs(path, exist_ok=True)


def open_roms_dataset(config_path: str):
    """
    Open a ROMS history file and return a prepared dataset, grid, and config.

    Reads the resolved config at `config_path`, opens the history NetCDF file
    referenced by that config, and calls `prep_ds` to attach grid metrics.

    Parameters
    ----------
    config_path : str
        Path to a resolved_config.yaml produced by prep_experiment.

    Returns
    -------
    ds : xarray.Dataset
        History dataset with vertical/horizontal metrics attached.
    grid : xgcm.Grid
        xgcm Grid with coordinate metrics for averaging/differencing.
    params : dict
        The loaded config dict (useful for accessing run parameters).
    """
    params = load_yaml(config_path)
    his_path = os.path.join(params["io"]["output_dir"], params["files"]["his"])
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds = xr.open_dataset(his_path, decode_times=time_coder)
    ds, grid = prep_ds(ds, params)
    return ds, grid, params

def prep_ds(ds, params):
    """
    Prepare dataset and xgcm Grid with correct vertical and horizontal metrics.

    Adds:
      - z_rho, z_w
      - dx/dy at rho, u, v, psi
      - dz at rho (layer thickness) and dual dz_w at w
      - dA at rho, and dA_u/dA_v/dA_psi at staggered points
      - 3D volume metrics dV, dV_u, dV_v, dV_w

    Returns
    -------
    ds : xarray.Dataset
    grid : xgcm.Grid
    """
    # 1) Rename dims if necessary (only those present)
    rename_map = {
        'eta_u': 'eta_rho',
        'xi_v': 'xi_rho',
        'xi_psi': 'xi_u',
        'eta_psi': 'eta_v',
    }
    ds = ds.rename({k: v for k, v in rename_map.items() if k in ds.dims})

    # 2) Axis mapping for xgcm
    coords = {
        'X': {'center': 'xi_rho', 'inner': 'xi_u'},
        'Y': {'center': 'eta_rho', 'inner': 'eta_v'},
        'Z': {'center': 's_rho',  'outer': 's_w'},
    }

    # 3) Provisional grid (no metrics yet)
    grid = xgcm.Grid(ds, coords=coords, autoparse_metadata=False, padding='periodic')

    # 4) Vertical coordinates
    z_r, z_w = compute_depths(
        ds.isel(ocean_time=0).h.squeeze().values,
        params['vertical']['HC'],
        params['vertical']['THETA_S'],
        params['vertical']['THETA_B'],
        params['grid']['N'],
    )

    ds = ds.assign_coords(
        z_w=(('s_w', 'eta_rho', 'xi_rho'), z_w),
        z_rho=(('s_rho', 'eta_rho', 'xi_rho'), z_r),
    )

    # 5) Horizontal metrics at staggered points
    ds['pm_v']   = grid.interp(ds.pm, 'Y')
    ds['pn_u']   = grid.interp(ds.pn, 'X')
    ds['pm_u']   = grid.interp(ds.pm, 'X')
    ds['pn_v']   = grid.interp(ds.pn, 'Y')
    ds['pm_psi'] = grid.interp(grid.interp(ds.pm, 'Y'), 'X')
    ds['pn_psi'] = grid.interp(grid.interp(ds.pn, 'X'), 'Y')

    ds['dx']      = 1.0 / ds.pm
    ds['dx_u']    = 1.0 / ds.pm_u
    ds['dx_v']    = 1.0 / ds.pm_v
    ds['dx_psi']  = 1.0 / ds.pm_psi

    ds['dy']      = 1.0 / ds.pn
    ds['dy_u']    = 1.0 / ds.pn_u
    ds['dy_v']    = 1.0 / ds.pn_v
    ds['dy_psi']  = 1.0 / ds.pn_psi

    # 6) Vertical metrics
    # 6a) rho-layer thickness (positive)
    dz_rho = grid.diff(ds.z_w, 'Z')  # (s_rho, eta_rho, xi_rho)
    dz_rho = abs(dz_rho)
    ds['dz'] = dz_rho

    # 6b) dual thickness at w (positive): average dz to w and half endpoints
    dz_w = grid.interp(ds['dz'], 'Z', padding='extend')  # (s_w, eta_rho, xi_rho)
    w_s = xr.DataArray(np.ones(dz_w.sizes['s_w']), dims=['s_w'])
    w_s[0]  = 0.5
    w_s[-1] = 0.5
    dz_w = dz_w * w_s
    ds['dz_w'] = dz_w

    # u/v variants
    ds['dz_u']   = grid.interp(ds['dz'],   'X')
    ds['dz_v']   = grid.interp(ds['dz'],   'Y')
    ds['dz_w_u'] = grid.interp(ds['dz_w'], 'X')
    ds['dz_w_v'] = grid.interp(ds['dz_w'], 'Y')

    # 7) 2D areas at all staggered points
    ds['dA']     = ds['dx']     * ds['dy']
    ds['dA_u']   = ds['dx_u']   * ds['dy_u']
    ds['dA_v']   = ds['dx_v']   * ds['dy_v']
    ds['dA_psi'] = ds['dx_psi'] * ds['dy_psi']

    # 8) Optional: 3D volumes (helpful for explicit weighting)
    ds['dV']   = ds['dA']   * ds['dz']
    ds['dV_u'] = ds['dA_u'] * ds['dz_u']
    ds['dV_v'] = ds['dA_v'] * ds['dz_v']
    ds['dV_w'] = ds['dA']   * ds['dz_w']  # area at rho, thickness at w

    # 9) Metrics dictionary (include u/v/psi areas to avoid interpolation)
    metrics = {
        ('X',): ['dx', 'dx_u', 'dx_v', 'dx_psi'],
        ('Y',): ['dy', 'dy_u', 'dy_v', 'dy_psi'],
        ('Z',): ['dz', 'dz_u', 'dz_v', 'dz_w', 'dz_w_u', 'dz_w_v'],
        ('X', 'Y'): ['dA', 'dA_u', 'dA_v', 'dA_psi'],
    }

    grid = xgcm.Grid(ds, coords=coords, metrics=metrics, padding='periodic', autoparse_metadata=False)
    return ds, grid