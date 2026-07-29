# python
import numpy as np
import yaml 
import xarray as xr
import xgcm

def compute_stretching(theta_s, theta_b, N):
    """
    Compute vertical stretching curves (s_rho, Cs_r) and (s_w, Cs_w)
    assuming zeta=0 (mean sea level) for Vtransform=2 and Vstretching=5.

    Parameters
    ----------
    h : ndarray (eta_rho, xi_rho), positive downward [m]
        (not used directly; included for API symmetry with other ROMS code)
    hc, theta_s, theta_b, N : scalar ROMS parameters

    Returns
    -------
    s_rho : ndarray (N,), normalized S-coordinates at rho-points
    Cs_r : ndarray (N,), stretching curves at rho-points
    s_w : ndarray (N+1,), normalized S-coordinates at w-points
    Cs_w : ndarray (N+1,), stretching curves at w-points
    """
    # Compute s_rho (N levels) and s_w (N+1 levels)
    k_rho = np.arange(1, N + 1, dtype=float)  # Rho-points
    k_w = np.arange(0, N + 1, dtype=float)    # W-points
    rN = float(N)

    # Normalized S-coordinates for rho-points
    s_rho = -(k_rho**2 - 2.0 * k_rho * rN + k_rho + rN**2 - rN) / (rN**2 - rN) \
            - 0.01 * (k_rho**2 - k_rho * rN) / (1.0 - rN)
    # Normalized S-coordinates for w-points
    s_w = -(k_w**2 - 2.0 * k_w * rN + k_w + rN**2 - rN) / (rN**2 - rN) \
          - 0.01 * (k_w**2 - k_w * rN) / (1.0 - rN)

    # Stretching function for rho-points
    if theta_s > 0:
        Csur_rho = (1.0 - np.cosh(theta_s * s_rho)) / (np.cosh(theta_s) - 1.0)
    else:
        Csur_rho = -(s_rho ** 2)

    if theta_b > 0:
        Cs_r = (np.exp(theta_b * Csur_rho) - 1.0) / (1.0 - np.exp(-theta_b))
    else:
        Cs_r = Csur_rho

    # Stretching function for w-points
    if theta_s > 0:
        Csur_w = (1.0 - np.cosh(theta_s * s_w)) / (np.cosh(theta_s) - 1.0)
    else:
        Csur_w = -(s_w ** 2)

    if theta_b > 0:
        Cs_w = (np.exp(theta_b * Csur_w) - 1.0) / (1.0 - np.exp(-theta_b))
    else:
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