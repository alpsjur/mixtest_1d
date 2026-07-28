# python
import numpy as np

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