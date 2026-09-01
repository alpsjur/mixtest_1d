# python
import os
import numpy as np
import yaml
import xarray as xr
import xgcm

# ---------------------------------------------------------------------------
# Physical constants
#
# These must stay in sync with the hardcoded EOS/Boussinesq values in
# templates/mixtest_1d.in.j2 (RHO0, R0, TCOEF are not currently exposed as
# config parameters there).
# ---------------------------------------------------------------------------
G = 9.81          # gravitational acceleration [m s-2]
RHO0 = 1025.0     # Boussinesq reference density [kg m-3] (template RHO0)
R0 = 1027.0       # linear EOS reference density [kg m-3]  (template R0)
TCOEF = 1.7e-4    # thermal expansion coefficient [degC-1] (template TCOEF)

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


def compute_phi(ds, grid, params):
    """
    Time evolution of the stratification potential energy anomaly phi(t),
    following Carpenter et al. (2016), Eq. 7:

        phi(t) = integral_0^H  g * z * (rho_mix - rho(z,t))  dz

    with z measured *upward from the seabed* (z=0 at the bed, z=H at the
    surface) -- the opposite convention to ROMS' z_rho (0 at the surface,
    negative downward). rho_mix is the density the water column would have
    if instantaneously and completely mixed; since this idealised setup has
    no surface/bottom buoyancy fluxes and no tracer sources or nudging
    (LtracerSrc/LtracerCLM == F in the .in template), the volume-averaged
    density is conserved and rho_mix is just the (constant-in-time)
    volume average of the initial density field.

    phi > 0 for stable stratification and phi -> 0 as the water column
    is mixed. Throughout this function "rho" refers to whatever density
    field is stored in the history file (ROMS' idDano "density anomaly");
    since phi only involves differences, the constant offset used by that
    convention is irrelevant.

    Parameters
    ----------
    ds : xarray.Dataset
        Prepared history dataset (see open_roms_dataset / prep_ds).
    grid : xgcm.Grid
        Grid with metrics attached (see prep_ds).
    params : dict
        Resolved run config (used for grid.H0).

    Returns
    -------
    phi : xarray.DataArray
        Time series of phi(t), dims (ocean_time,), units J m-2 (per unit
        horizontal area).
    rho_mix0 : xarray.DataArray
        Scalar volume-averaged initial density (the "rho_mix" datum).
    """
    if "rho" not in ds:
        raise KeyError(
            "Variable 'rho' (density anomaly, idDano) not found in dataset. "
            "Enable Hout(idDano) in the ROMS input file."
        )

    H0 = float(params["grid"]["H0"])

    # rho_mix: volume-averaged density at t=0 (conserved for all t in this
    # closed-column setup -- see docstring). Kept as a DataArray so it
    # broadcasts cleanly against the full field below.
    rho_mix0 = grid.average(ds["rho"].isel(ocean_time=0), axis=("X", "Y", "Z"))

    # z measured upward from the seabed (Carpenter et al. convention).
    z_from_bed = ds["z_rho"] + H0

    integrand = G * (rho_mix0 - ds["rho"]) * z_from_bed

    # Integrate vertically first (matches metric dims exactly), then average
    # horizontally -- equivalent to averaging first since the domain is
    # horizontally homogeneous, but avoids metric/dim mismatches.
    phi_xy = grid.integrate(integrand, "Z")
    phi = grid.average(phi_xy, axis=("X", "Y")).squeeze()
    phi.name = "phi"
    phi.attrs["long_name"] = "stratification potential energy anomaly"
    phi.attrs["units"] = "J m-2"

    return phi, rho_mix0


def compute_Pd(ds, grid, params):
    """
    Structure-drag turbulence production rate P_d(z,t), following the
    STRUCTURE_MIXING parametrization (Rennau, Schimmels & Burchard 2012):

        P_d = 0.5 * CD * str_a * (u^2 + v^2)^(3/2)     [m2 s-3]

    str_a and CD are taken from the resolved config (params["structure"]),
    i.e. this assumes str_a is spatially uniform over the water column, as
    is the case for the baseline/sensitivity-sweep configs (depth_zero_below
    set far below H0). If a depth-varying str_a is used, this needs to be
    read from the grid file instead.

    Returns
    -------
    Pd : xarray.DataArray
        P_d(z,t) at rho-points, dims (ocean_time, s_rho, eta_rho, xi_rho).
    """
    CD = float(params["structure"]["CD"])
    str_a = float(params["structure"]["str_a"])

    u_r = grid.interp(ds["u"], "X")
    v_r = grid.interp(ds["v"], "Y")
    spd2 = u_r ** 2 + v_r ** 2

    Pd = 0.5 * CD * str_a * spd2 * np.sqrt(spd2)
    Pd.name = "Pd"
    Pd.attrs["long_name"] = "structure-drag turbulence production"
    Pd.attrs["units"] = "m2 s-3"
    return Pd


def compute_Pstr(ds, grid, params):
    """
    Time series of P_str(t), the depth-integrated (per unit horizontal area)
    power extracted from the mean flow by structure drag:

        P_str(t) = rho0 * integral_0^H  P_d(z,t)  dz     [W m-2]

    This is the diagnostic (model-derived) counterpart to the analytic
    P_str used in analytic_mixing_timescale, and is meant to be evaluated
    from actual u,v time series (accounts for spin-up transients etc.).

    Returns
    -------
    Pstr : xarray.DataArray
        Time series, dims (ocean_time,), units W m-2.
    """
    Pd = compute_Pd(ds, grid, params)
    Pstr_xy = grid.integrate(RHO0 * Pd, "Z")
    Pstr = grid.average(Pstr_xy, axis=("X", "Y")).squeeze()
    Pstr.name = "Pstr"
    Pstr.attrs["long_name"] = "depth-integrated structure-drag power extraction"
    Pstr.attrs["units"] = "W m-2"
    return Pstr


def analytic_mixing_timescale(cfg: dict) -> dict:
    """
    Analytic (pre-run) estimate of the mixing time scale tau_mix, following
    Carpenter et al. (2016).

    Note: the body force in this idealised setup represents *any* steady
    background current (e.g. a mean coastal current, or a tidally-averaged
    residual current) -- the drag/mixing balance below only depends on the
    resulting flow speed, not on what physically drives it.

    Assumes the body-force-driven flow reaches a quasi-steady balance with
    structure drag (same balance as tests/test_STRUCTURE_DRAG.py):

        u_inf = sqrt(2 * BFRC_U / (CD * str_a))

    This requires BFRC_V == 0 (else the 2D force/drag balance is not a
    simple closed-form expression) and str_a active over the full depth
    (structure.depth_zero_below >= grid.H0), since P_str = rho0 * H * P_d
    assumes P_d is depth-uniform.

    tau_mix = g * delta_rho * H^2 / P_str
            = g * delta_rho * H / (rho0 * P_d)

    Use this to size NTIMES for a run *before* it has been executed (e.g. in
    a sensitivity sweep where each parameter combination implies a different
    tau_mix). To check whether the actual model run mixes on this predicted
    time scale, compare against the diagnostic P_str from compute_Pstr()
    and the phi(t) curve from compute_phi() after the fact.

    Parameters
    ----------
    cfg : dict
        Resolved (or about-to-be-resolved) run config, with at least
        structure.{CD,str_a,depth_zero_below}, bodyforce.{BFRC_U,BFRC_V},
        grid.H0, initial.temp_dT.

    Returns
    -------
    dict with keys:
        u_inf         quasi-steady speed [m/s]
        Pd            structure-drag production [m2/s3]
        Pstr          depth-integrated power extraction [W/m2]
        delta_rho     top-to-bottom density difference [kg/m3]
        tau_mix       predicted mixing time scale [s]
    """
    CD = float(cfg["structure"]["CD"])
    str_a = float(cfg["structure"]["str_a"])
    depth_zero_below = float(cfg["structure"]["depth_zero_below"])
    H0 = float(cfg["grid"]["H0"])
    BFRC_U = float(cfg["bodyforce"]["BFRC_U"]) * 1e-7  # config stores e-7 m/s2
    BFRC_V = float(cfg["bodyforce"].get("BFRC_V", 0.0)) * 1e-7
    temp_dT = float(cfg["initial"]["temp_dT"])

    if abs(BFRC_V) > 0.0:
        raise ValueError(
            "analytic_mixing_timescale assumes BFRC_V == 0 (simple 1D force/"
            "drag balance). For BFRC_V != 0, estimate tau_mix from a model "
            "run instead (compute_Pstr + compute_phi)."
        )
    if depth_zero_below < H0:
        raise ValueError(
            "analytic_mixing_timescale assumes str_a is uniform over the "
            "full water column (structure.depth_zero_below >= grid.H0)."
        )
    if str_a <= 0.0 or CD <= 0.0:
        raise ValueError("analytic_mixing_timescale requires str_a > 0 and CD > 0.")

    c4 = cfg.get("structure", {}).get("c4")
    c1 = cfg.get("GLS", {}).get("C1")
    if c4 is not None and c1 is not None and float(c4) > float(c1):
        raise ValueError(
            f"structure.c4 ({c4}) > GLS.C1 ({c1}). In Carpenter et al., c4 > C1 "
            "corresponds to reduced structure-induced mixing efficiency "
            "(c4 = C1 is the neutral point). In this implementation that regime "
            "is also outside the validated/numerically robust range, so "
            "analytic_mixing_timescale() rejects it."
        )

    u_inf = np.sqrt(2.0 * BFRC_U / (CD * str_a))
    Pd = 0.5 * CD * str_a * u_inf ** 3
    Pstr = RHO0 * H0 * Pd
    delta_rho = R0 * TCOEF * temp_dT

    if delta_rho <= 0.0:
        raise ValueError("analytic_mixing_timescale requires initial.temp_dT > 0.")

    tau_mix = G * delta_rho * H0 ** 2 / Pstr

    return {
        "u_inf": u_inf,
        "Pd": Pd,
        "Pstr": Pstr,
        "delta_rho": delta_rho,
        "tau_mix": tau_mix,
    }


def compute_time_vector(params: dict) -> np.ndarray:
    """
    Build a time axis in days from run config parameters.

    Uses DT (timestep in seconds), NHIS (output interval in timesteps),
    and NTIMES (total timesteps) to produce one value per output record.
    """
    DT = params["time_stepping"]["DT"]
    NHIS = params["time_stepping"]["NHIS"]
    NTIMES = params["time_stepping"]["NTIMES"]
    N = NTIMES // NHIS
    dt = DT * NHIS  # seconds per output record
    t_seconds = np.arange(N+1, dtype=np.float64) * dt
    return t_seconds / (60 * 60 * 24)  # convert to days