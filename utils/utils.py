# python
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
    with open(path, "r") as f:
        return yaml.safe_load(f)

def prep_ds(ds, params):
    """
    Prepare dataset and xgcm Grid with correct vertical metrics for volume integrals.

    Steps:
    - Rename dims to xgcm-friendly names (rho/u/v/psi alignment)
    - Build xgcm Grid
    - Compute z_rho and z_w from bathymetry and vertical params
    - Add z-coordinates to ds
    - Interpolate horizontal metrics to u, v, psi points and compute dx/dy
    - Compute dz (rho thickness) and dz_w (dual thickness at w points) correctly
    - Build complete metrics dict and re-create the Grid

    Parameters
    ----------
    ds : xarray.Dataset
        ROMS dataset with h, pm, pn, and standard ROMS grid dims/coords.
    params : dict
        {
          'vertical': {
              'HC': float, 'THETA_S': float, 'THETA_B': float
          },
          'grid': {'N': int}
        }

    Returns
    -------
    ds : xarray.Dataset
        Dataset augmented with z_rho, z_w, dx/dy, dz, dz_w and their u/v/psi variants.
    grid : xgcm.Grid
        Grid object configured with the correct metrics.
    """
    # 1) Rename dims that match
    rename_map = {
        'eta_u': 'eta_rho',
        'xi_v': 'xi_rho',
        'xi_psi': 'xi_u',
        'eta_psi': 'eta_v',
    }
    existing = {k: v for k, v in rename_map.items() if k in ds.dims}
    if existing:
        ds = ds.rename(existing)

    # 2) Define coords mapping for xgcm
    coords = {
        'X': {'center': 'xi_rho', 'inner': 'xi_u'},
        'Y': {'center': 'eta_rho', 'inner': 'eta_v'},
        'Z': {'center': 's_rho',  'outer': 's_w'},
    }

    # 3) Create a provisional Grid (no metrics yet)
    grid = xgcm.Grid(ds, coords=coords, autoparse_metadata=False, padding='periodic')

    # 4) Compute vertical coordinates (z-levels) at rho and w points
    # Assumes compute_depths returns z_r (s_rho, eta_rho, xi_rho) and z_w (s_w, eta_rho, xi_rho)
    z_r, z_w = compute_depths(
        ds.isel(ocean_time=0).h.squeeze().values,
        params['vertical']['HC'],
        params['vertical']['THETA_S'],
        params['vertical']['THETA_B'],
        params['grid']['N'],
    )

    # 5) Add z-coordinates to dataset
    ds = ds.assign_coords(
        z_w=(('s_w', 'eta_rho', 'xi_rho'), z_w),
        z_rho=(('s_rho', 'eta_rho', 'xi_rho'), z_r),
    )

    # 6) Interpolate horizontal metrics to u, v, and psi points
    ds['pm_v']   = grid.interp(ds.pm, 'Y')
    ds['pn_u']   = grid.interp(ds.pn, 'X')
    ds['pm_u']   = grid.interp(ds.pm, 'X')
    ds['pn_v']   = grid.interp(ds.pn, 'Y')
    ds['pm_psi'] = grid.interp(grid.interp(ds.pm, 'Y'), 'X')
    ds['pn_psi'] = grid.interp(grid.interp(ds.pn, 'X'), 'Y')

    # 7) Horizontal grid spacings
    ds['dx']      = 1.0 / ds.pm
    ds['dx_u']    = 1.0 / ds.pm_u
    ds['dx_v']    = 1.0 / ds.pm_v
    ds['dx_psi']  = 1.0 / ds.pm_psi

    ds['dy']      = 1.0 / ds.pn
    ds['dy_u']    = 1.0 / ds.pn_u
    ds['dy_v']    = 1.0 / ds.pn_v
    ds['dy_psi']  = 1.0 / ds.pn_psi

    # 8) Vertical metrics
    # 8a) Layer thickness at rho points (N layers), computed from z_w
    # grid.diff(z_w,'Z') returns differences between outer (w) to center (rho), i.e. dz at s_rho
    dz_rho = grid.diff(ds.z_w, 'Z')  # (s_rho, eta_rho, xi_rho)
    # Make thickness positive regardless of z sign convention
    dz_rho = abs(dz_rho)
    ds['dz'] = dz_rho  # thickness for rho-centered fields

    # 8b) Dual-cell thickness at w points for w-centered fields
    # Interpolate dz (rho) to w, then halve the end points
    dz_w = grid.interp(ds['dz'], 'Z', padding='extend')  # (s_w, eta_rho, xi_rho)
    # Halve top and bottom control volumes using a 1D weight along s_w
    w_s = xr.DataArray(np.ones(dz_w.sizes['s_w']), dims=['s_w'])
    w_s[0]  = 0.5
    w_s[-1] = 0.5
    dz_w = dz_w * w_s
    ds['dz_w'] = dz_w

    # Also build u/v versions
    ds['dz_u']   = grid.interp(ds['dz'],   'X')
    ds['dz_v']   = grid.interp(ds['dz'],   'Y')
    ds['dz_w_u'] = grid.interp(ds['dz_w'], 'X')
    ds['dz_w_v'] = grid.interp(ds['dz_w'], 'Y')

    # 9) Cell area at rho points (and can interpolate as needed)
    ds['dA'] = ds['dx'] * ds['dy']

    # 10) Define metrics and re-create the Grid with metrics
    metrics = {
        ('X',): ['dx', 'dx_u', 'dx_v', 'dx_psi'],
        ('Y',): ['dy', 'dy_u', 'dy_v', 'dy_psi'],
        ('Z',): ['dz', 'dz_u', 'dz_v', 'dz_w', 'dz_w_u', 'dz_w_v'],
        ('X', 'Y'): ['dA'],
    }

    grid = xgcm.Grid(ds, coords=coords, metrics=metrics, padding='periodic', autoparse_metadata=False)

    return ds, grid