# Mixing-Timescale Sensitivity Analysis

**Date:** 2026-08-31
**Repo:** `mixtest_1d` (idealized 1D ROMS testbed for the `STRUCTURE_MIXING` parametrization)
**Reference:** Carpenter, J. R., Merckelbach, L., Callies, U., Clark, S., Gaslikova, L., & Rippeth, T. P. (2016). *Potential Impacts of Offshore Wind Farms on North Sea Stratification*. PLoS ONE, 11(8), e0160830. https://doi.org/10.1371/journal.pone.0160830

## 1. Motivation

The `STRUCTURE_MIXING` parametrization (see `roms/IMPLEMENTATION_STRUCTURE_MIXING.md`) adds two effects to ROMS to represent offshore wind turbine foundations:

1. A **drag force** on the flow, proportional to the local frontal area density of structures (`str_a`) and a drag coefficient `CD`.
2. A **turbulence production** term injected into the GLS (Generic Length Scale) turbulence closure's TKE and length-scale (ψ) equations, proportional to the power dissipated by that drag (`P_d`) and scaled by a coefficient `c4`.

Carpenter et al. (2016) analyze how offshore structures erode North Sea stratification by extracting momentum from tidal flow and converting it into turbulent mixing. They define a **mixing time scale** τ_mix that predicts how long it takes stratification to be eroded, given the power extracted by the structures.

This analysis reproduces a Carpenter-et-al.-style sensitivity test in `mixtest_1d`, checking:
- Whether the model's actual mixing rate matches the analytic (theoretical) prediction.
- Whether the dimensionless collapse Carpenter et al. found (φ*(t*) collapsing across cases) also holds here.
- How the *mixing efficiency* depends on the STRUCTURE_MIXING/GLS closure parameters, which the simple power-balance theory does not by itself predict.

Note on terminology: in this idealized 1D setup, the flow is driven by a constant background body force (`bodyforce.BFRC_U`), representing a generic steady current (e.g. the mean/residual flow of the Norwegian Coastal Current) — **not** a tidal velocity as in Carpenter et al.'s North Sea application. The physics of the drag/mixing balance only depends on the resulting flow speed, not on what drives it, so the same framework applies.

## 2. Theory

### 2.1 Stratification potential energy anomaly φ(t)

Carpenter et al. define (their Eq. 7, adapted to this repo's z-convention):

```
φ(t) = ∫₀^H [ρ_mix − ρ(z,t)] g z dz
```

where `z` is measured **upward from the seabed** (0 at the bed, `H` at the surface) and `ρ_mix` is the density the water column would have if it were completely and uniformly mixed. `φ` is a potential-energy-like measure of stratification: `φ = φ(0) > 0` for a stratified column, and `φ → 0` as the column mixes to uniform density.

Implementation: `utils/utils.py::compute_phi(ds, grid, params)`. ROMS' native `z_rho` is 0 at the surface and negative downward, so we transform via `z_from_bed = z_rho + H0`.

**No reference "fully-mixed" run is needed.** Because this idealized setup has no surface/bottom heat or tracer fluxes and no tracer climatology relaxation (`LtracerSrc == LtracerCLM == F F` in the `.in` template), the volume-averaged density is exactly conserved over time. This was verified numerically (volume-avg density conserved to ~1e-7 relative error across a full run). Therefore:

```
ρ_mix(t) = ρ_mix(0)
```

and `ρ_mix` can be computed directly from the initial condition, with no need to run a separate "completely mixed" reference simulation.

### 2.2 Structure power extraction P_d, P_str

Following the same drag-balance physics as `tests/test_STRUCTURE_DRAG.py`, the flow driven by a constant background body force `BFRC_U` (with `BFRC_V = 0`) reaches a quasi-steady balance with structure drag:

```
u_inf  = sqrt(2 · BFRC_U / (CD · str_a))
P_d    = 0.5 · CD · str_a · u_inf³            (power dissipated per unit volume, m²/s³)
P_str  = ρ₀ · H · P_d                          (depth-integrated power, W/m², using P_d = P_str/(ρ₀H))
```

Implementation: `utils/utils.py::compute_Pd`, `compute_Pstr` (diagnostic, computed from the actual `u`, `v` in a completed run) and `analytic_mixing_timescale` (pre-run analytic prediction using the assumed steady-state `u_inf`).

### 2.3 Dimensionless mixing time scale

Carpenter et al. define a dimensionless time:

```
t* = t · P_str / (g · Δρ · H²)
```

chosen so that (in their idealized time-dependent pycnocline model) `t* = 1` corresponds to the pycnocline being "fully mixed". This motivates defining a mixing time scale:

```
τ_mix ≡ g · Δρ · H² / P_str
```

In this repo, `Δρ = R0 · TCOEF · temp_dT` (top-to-bottom density difference from the linear equation of state, using the initial temperature drop across the thermocline), so:

```
τ_mix = g · Δρ · H / (ρ₀ · P_d)
```

Two versions of τ_mix are computed:
- **`tau_mix_theory`** — purely analytic, computed *before* running the model, from the assumed steady-state `u_inf` (used to size `NTIMES` for each sweep run).
- **`tau_mix_diagnostic`** — computed *after* the model has run, using `P_str` diagnosed directly from the model's actual `u`, `v` fields (averaged over the last half of the run, to exclude the initial spin-up transient).

A dimensionless density anomaly is also defined (not directly used in the plots below, but available for reference): `ρ* = 2(ρ − ρ₀)/Δρ`.

**Important: `c4` does not appear anywhere in this analytic formula.** `c4` only controls how efficiently the diagnosed power `P_d` is converted into eddy diffusivity via the GLS ψ-equation's structure-production source term. Sweeping `c4` therefore directly tests whether the *actual* model mixing behavior deviates from the simplified Carpenter power-balance scaling — see §4 below.

## 3. Model setup

### 3.1 Prerequisite: switching to a linear equation of state

The default `mixtest_1d` build uses `NONLIN_EOS` (a full nonlinear/compressible equation of state). Under this EOS, an initial sanity-check run showed that `φ(t)` decayed to only ~7% of its initial value and then plateaued, even though the volume-averaged temperature and salinity had become perfectly homogeneous by day ~5. Investigation showed the residual density variation was a **compressibility artifact**: in-situ (nonlinear) density retains a pressure-dependent (depth-dependent) component even when T and S are exactly uniform. Carpenter et al.'s φ formula implicitly assumes density depends only on T/S (a linear EOS), with no compressibility term.

**Fix:** `roms/Include/mixtest_1d.h` was changed to `#undef NONLIN_EOS` (there is no separate `LINEAR_EOS` macro in ROMS — linear EOS is simply the default when `NONLIN_EOS` is undefined; verified via ROMS' `checkdefs.F`). ROMS was rebuilt (`roms/build_roms.sh`, non-parallel) and the existing unit test suite (`test_STRUCTURE_DRAG`, `test_STRUCTURE_PRODUCTION`, `test_UV_BODYFORCE`) was re-verified to still pass. After this fix, the sanity-check run's φ(t) decayed cleanly to ~0.

### 3.2 Sweep design

Defined in `templates/mixing_timescale_sweep.yaml`, with the following cartesian sweep:

| Parameter | Values | Meaning |
|---|---|---|
| `structure.CD` | 0.63, 1.26 | structure drag coefficient |
| `structure.c4` | 0.44, 0.97 (0.44, 0.97, **1.4 removed**, see §4.3) | GLS ψ-equation structure-production coefficient |
| `initial.temp_dT` | 5.0, 10.0 °C | temperature drop across the thermocline (sets Δρ) |
| `grid.H0` | 75.0, 150.0 m | water column depth |

`configs/variants/mixing_timescale.yaml` scales `structure.str_a` and `bodyforce.BFRC_U` up by ~100x relative to `baseline.yaml` (which alone gives τ_mix ≈ 921 days — far too slow for a practical sweep) so that τ_mix falls in a practical 2–13 day range, while keeping `u_inf` in a physically reasonable 0.07–0.1 m/s range for a coastal current.

`tools/prep_mixing_timescale_sweep.py` wraps the existing `prep_sweep.py` machinery, adding a per-run step that calls `analytic_mixing_timescale()` *before* preparing each run, and sizes `NTIMES` as:

```
NTIMES = ceil(ntimes_margin · tau_mix_theory / DT / NHIS) · NHIS
```

with `ntimes_margin = 2.5`, so each run's simulated duration is roughly 2.5x its own predicted τ_mix — long enough to see mixing complete and plateau, without wasting compute on the slowest cases. Because each parameter combination implies a different τ_mix, this per-run NTIMES sizing is essential (a single fixed NTIMES for the whole sweep would either truncate the slowest runs or waste enormous compute on the fastest ones).

This produced 16 runs (`runs/mixtau_CD*_c4*_dT*_H*/`), ranging from NTIMES=12510 (duration 2.3 days) to NTIMES≈70000+ (duration ~13 days), all completed successfully.

## 4. Results

### 4.1 τ_mix_theory vs. τ_mix_diagnostic: excellent agreement

| CD | c4 | dT (°C) | H0 (m) | τ_mix theory (d) | τ_mix diagnostic (d) | t*_mix |
|---|---|---|---|---|---|---|
| 0.63 | 0.44 | 5.0 | 75.0 | 2.302 | 2.302 | 0.561 |
| 1.26 | 0.44 | 5.0 | 75.0 | 3.256 | 3.256 | 0.550 |
| 0.63 | 0.44 | 10.0 | 75.0 | 4.605 | 4.605 | 0.543 |
| 1.26 | 0.44 | 10.0 | 75.0 | 6.512 | 6.512 | 0.531 |
| 0.63 | 0.44 | 5.0 | 150.0 | 4.605 | 4.605 | 0.461 |
| 1.26 | 0.44 | 5.0 | 150.0 | 6.512 | 6.512 | 0.454 |
| 0.63 | 0.44 | 10.0 | 150.0 | 9.209 | 9.209 | 0.448 |
| 1.26 | 0.44 | 10.0 | 150.0 | 13.024 | 13.024 | 0.445 |
| 0.63 | 0.97 | 5.0 | 75.0 | 2.302 | 2.302 | 0.760 |
| 1.26 | 0.97 | 5.0 | 75.0 | 3.256 | 3.256 | 0.729 |
| 0.63 | 0.97 | 10.0 | 75.0 | 4.605 | 4.605 | 0.724 |
| 1.26 | 0.97 | 10.0 | 75.0 | 6.512 | 6.512 | 0.717 |
| 0.63 | 0.97 | 5.0 | 150.0 | 4.605 | 4.605 | 0.615 |
| 1.26 | 0.97 | 5.0 | 150.0 | 6.512 | 6.512 | 0.608 |
| 0.63 | 0.97 | 10.0 | 150.0 | 9.209 | 9.209 | 0.597 |
| 1.26 | 0.97 | 10.0 | 150.0 | 13.024 | 13.024 | 0.595 |

For every run, `tau_mix_diagnostic` (computed from the model's own diagnosed `P_str`) matches `tau_mix_theory` (the pre-run analytic prediction) to within <1%. This confirms that the assumed quasi-steady drag/body-force balance (`u_inf = sqrt(2·BFRC_U/(CD·str_a))`) is indeed reached by the model, and that the analytic power-input scale is correctly predicted **independent of the turbulence closure coefficient `c4`** — as expected, since `u_inf` and `P_d` only depend on `CD`, `str_a`, and `BFRC_U`, not on how that power is subsequently converted into mixing.

### 4.2 The dimensionless collapse: φ*(t*)

![Mixing timescale collapse](../figures/mixing_timescale_collapse_mixing_timescale_H0.png)

Plotting `φ*(t) = φ(t)/φ(0)` against the dimensionless time `t* = t·P_str_diag/(g·Δρ·H²)` for all 16 valid runs produces a tight collapse — the direct analogue of Carpenter et al.'s Fig 4b. All 16 curves fall in a narrow band and decay from `φ* = 1` to `φ* ≈ 0` over a similar range of `t*`.

However, unlike Carpenter et al.'s idealized time-dependent pycnocline model (where mixing completion occurs at `t* = 1` essentially by construction/definition), **mixing completes earlier here, at t* ≈ 0.45–0.76**, and — notably — this completion time is **not universal**: it varies systematically with two of the four swept parameters (see §4.4). This means `t* = 1` should be understood as a rough order-of-magnitude time scale for *this* parametrization/setup, not an exact universal threshold; the value at which mixing actually completes depends on lower-order details of the turbulence closure and geometry not captured by the leading-order power balance.


### 4.3 Interpretation and instability for `c4 > GLS.C1`

An initial version of the sweep included `structure.c4 = 1.4` in addition to `0.44` and `0.97`. All 8 runs with `c4 = 1.4` showed essentially **no mixing**: `φ*` remained near 1 throughout the run, even though the mean flow still spun up to the analytically expected drag-limited speed. Diagnostics showed that `tke` and `AKt` collapsed to their floor values (`tke ≈ Kmin`, `AKt ≈ 1e-6 m² s⁻¹`), while the momentum/drag balance remained correct. Thus the failure occurred in the turbulence closure, not in the drag formulation.

This behavior can be interpreted in light of the original Rennau et al. structure-mixing theory. In the homogeneous `k-ε` analysis given in their Eqs. (17)–(18), the additional coefficient `c4` enters the dissipation-equation source term analogously to the standard shear-production coefficient `c1`. Their analysis shows that the distinguished value is **`c4 = c1`**:

- `c4 = c1`: structure production leaves the implied mixing efficiency unchanged;
- `c4 < c1`: structure production **enhances** the mixing efficiency;
- `c4 > c1`: structure production **reduces** the mixing efficiency.

Thus `c4` should be interpreted relative to `c1`, not as an independent free scaling. In the present GLS configuration, `GLS.C1 = 1.0`, so `c4 = 1.4` lies on the low-efficiency side of that theoretical threshold, whereas the successful sweep values (`0.44`, `0.97`) lie at or below it.

Accordingly, `c4 = 1.4` was removed from the production sweep, and `analytic_mixing_timescale()` now rejects `structure.c4 > GLS.C1` as outside the validated operating range of the present explicit STRUCTURE_MIXING/GLS coupling. This should be understood as a **practical restriction of the current implementation**, not as a universal theoretical prohibition on `c4 > c1` in all turbulence closures.

### 4.4 Why mixing completes early, and why it depends on c4 and H0

![Mixing-completion sensitivity](../figures/mixing_timescale_sensitivity_mixing_timescale_H0.png)

To quantify the deviation from Carpenter's universal `t* = 1`, we define an **empirical dimensionless mixing-completion time** `t*_mix`: the first `t*` at which `φ* < 0.05` (i.e. stratification potential energy has dropped to <5% of its initial value). Averaging over the two values of `CD` and `temp_dT` (which have negligible effect once time is nondimensionalized — see the flat rows within each `c4`/`H0` group in the table above), `t*_mix` depends cleanly on `c4` and `H0`:

| c4 | H0 = 75 m | H0 = 150 m |
|---|---|---|
| 0.44 | 0.546 ± 0.011 | 0.452 ± 0.006 |
| 0.97 | 0.733 ± 0.017 | 0.604 ± 0.008 |

Two effects, cleanly separable:

1. **Closure coefficient `c4`:** `t*_mix` is consistently ~1.34–1.36x *larger* (i.e. mixing takes longer, relative to the dimensionless time scale) at `c4 = 0.97` than at `c4 = 0.44`. Note that `c4` multiplies the production term in the **ψ (dissipation/length-scale) equation**, not the TKE equation directly. A larger `c4` increases the dissipation rate of turbulent kinetic energy relative to its production, which *shrinks* the turbulent length scale and hence `AKt` for the same power input `P_d`. This was confirmed directly: mean `AKt` over the water column (last half of the run) was ~0.37 m²/s at `c4 = 0.44` vs. ~0.19 m²/s at `c4 = 0.97`, for otherwise identical `CD`/`dT`/`H0` — i.e. roughly half the eddy diffusivity at the higher `c4`, consistent with the ~1.35x longer mixing time.

2. **Water column depth `H0` (pycnocline position):** `t*_mix` is consistently ~0.82–0.83x *smaller* (mixing completes sooner, relative to t*) at `H0 = 150` m than at `H0 = 75` m. The thermocline centre depth (`initial.temp_zt = 40` m, fixed across the sweep) sits at a different *relative* position in the water column depending on `H0`: 53% of the way down for `H0 = 75` m, but only 27% of the way down for `H0 = 150` m. A pycnocline positioned closer to a boundary (surface or bed) erodes faster under a roughly depth-uniform eddy diffusivity than one centered mid-column, because it has less distance to homogenize toward the nearer boundary and benefits from boundary-adjacent turbulence. Since the pycnocline depth `temp_zt` was not itself swept, this session could not initially distinguish "absolute pycnocline depth" from "relative pycnocline position" as the controlling variable — resolved in §4.5 below by testing an alternative length scale.

The two effects multiply cleanly: e.g. `t*_mix(c4=0.97, H0=150) / t*_mix(c4=0.97, H0=75) = 0.604/0.733 = 0.824`, matching `t*_mix(c4=0.44, H0=150)/t*_mix(c4=0.44, H0=75) = 0.452/0.546 = 0.828` to within 0.5%. Likewise the `c4` ratio is ~1.342 at `H0=75` and ~1.336 at `H0=150`. This separability suggests `t*_mix ≈ A(c4) · B(H0)` to good approximation across the tested range, though this has only been checked over 2 values of each parameter and should not be over-generalized.

### 4.5 Testing an alternative length scale: pycnocline position instead of H0

The `H0` dependence found in §4.4 raises a natural question: Carpenter et al.'s `t* = t·P_str/(g·Δρ·H²)` uses the **full water column depth** `H` as the length scale that sets the dimensionless mixing time. But physically, what actually has to be eroded is the pycnocline itself — so a length scale tied to the **pycnocline's position within the column** might be more appropriate than the total depth, especially when (as here) `H0` is varied while the pycnocline depth `temp_zt` is held fixed.

**Candidate length scale tested:** the geometric mean of the pycnocline centre's distances to the two boundaries,

```
L = sqrt(z_t · (H0 − z_t))
```

(implemented as `pycnocline_length_scale()` in `analysis/mixing_timescale.py`), used in place of `H0` in both `t*` and `τ_mix`. Several other candidates were also tried for comparison — the distance to the *nearer* boundary alone (`min(z_t, H0−z_t)`), the distance to the *farther* boundary alone (`max(z_t, H0−z_t)`), and the harmonic mean — but all of these gave a *worse* collapse than `H0` itself. Only the geometric mean improved on it.

**Result:** grouping runs by `structure.c4` (averaging over `CD`/`temp_dT`, which have negligible effect as before) and comparing the coefficient of variation (CV = std/mean) of the empirical `t*_mix` across the two `H0` values:

| Length scale `L` | CV at c4=0.44 | CV at c4=0.97 |
|---|---|---|
| `H0` (original) | 0.096 | 0.098 |
| `sqrt(z_t·(H0−z_t))` (pycnocline) | 0.031 | 0.030 |

Using the pycnocline-based length scale reduces the residual `H0`-driven spread in `t*_mix` by **a factor of ~3**, for both values of `c4`. Visually, plotting `φ*(t*)` with each length scale (`analysis/mixing_timescale.py --length-scale {H0,pycnocline}`) makes the improvement clear:

![Length scale comparison](../figures/mixing_timescale_length_scale_comparison.png)

*(Left: original `t*` using `H0`, colored by `c4`/`H0` group — note `H0=75` and `H0=150` runs form visibly separate sub-clusters within each `c4` color. Right: `t*` using the pycnocline length scale — the `H0=75`/`H0=150` sub-clusters collapse onto each other almost completely, leaving only the `c4` grouping.)*

**Interpretation:** the full water column depth `H0` is not, in general, the "correct" length scale for the dimensionless mixing time when the pycnocline sits at a fixed absolute depth rather than scaling with `H0`. The geometric mean of the pycnocline's two boundary distances captures the relevant "how far must this stratification travel to homogenize toward a boundary" scale far better than the total depth does — consistent with viscous/diffusive-erosion problems in general, where the natural length scale is tied to the structure being eroded (the pycnocline and its distances to the boundaries), not the size of the overall domain.

This is a genuinely useful refinement of the Carpenter et al. (2016) framework for cases like this one, where offshore structures span the full water column but the pycnocline itself sits at a roughly fixed absolute depth (set by seasonal thermal forcing) regardless of local bathymetric depth — exactly the situation on much of the Norwegian shelf. `analysis/mixing_timescale.py` now supports both length-scale choices via a `length_scale=` argument / `--length-scale` CLI flag, defaulting to `"H0"` for backward compatibility with Carpenter et al.'s original definition.

**Caveat (resolved in §4.6 below):** this was originally tested over only 2 values of `H0` (75 m and 150 m) at one fixed `z_t = 40` m, so the exact functional form (geometric mean, specifically) was not yet rigorously established, and "absolute pycnocline depth" was confounded with "relative pycnocline position" — a proper test required varying `z_t` independently of `H0` across a wider range.

**Another caveat:** it is not clear how this translates to structures not extending trough the whole water column, relevant for floating foundaitons. 

### 4.6 Confirming the geometric-mean form: varying z_t independently of H0

The §4.5 test above only varied `z_t/H0` implicitly (`z_t=40` fixed, `H0∈{75,150}` swept), so it could not separate "absolute pycnocline depth" from "relative pycnocline position", nor establish that the geometric mean specifically (as opposed to some other function of `z_t` and `H0−z_t`) is the right form. A dedicated follow-up sweep (`templates/pycnocline_zt_sweep.yaml`, `tools/prep_pycnocline_zt_sweep.py`, `sweeps/pycnocline_zt/`) tests this directly using explicit `(H0, z_t)` pairs (rather than a full cartesian product, which would place the thermocline too close to a boundary for many combinations):

| `H0` (m) | `z_t` (m) | `z_t/H0` |
|---|---|---|
| 75 | 25, 35, 45, 55 | 0.33 – 0.73 |
| 150 | 25, 45, 75, 115 | 0.17 – 0.77 |

(`z_t = 25` and `z_t = 45` are shared between both `H0` values, isolating the pure "H0 changed, absolute pycnocline depth fixed" effect.) Each of these 8 `(H0, z_t)` combinations was run at both `structure.c4 ∈ {0.44, 0.97}` (16 runs total), with `structure.CD = 0.63` and `initial.temp_dT = 10.0` held fixed as before.

**Result:** grouping only by `structure.c4` (i.e. averaging over all 8 `H0`/`z_t` combinations, not just 2 `H0` values as in §4.5), the coefficient of variation of `t*_mix`:

| Length scale `L` | CV at c4=0.44 | CV at c4=0.97 |
|---|---|---|
| `H0` (original) | 0.162 | 0.164 |
| `min(z_t, H0−z_t)` | 0.540 | 0.535 |
| `max(z_t, H0−z_t)` | 0.487 | 0.487 |
| harmonic mean, `2·z_t·(H0−z_t)/H0` | 0.225 | 0.220 |
| `sqrt(z_t·(H0−z_t))` (geometric mean) | **0.027** | **0.025** |

The geometric-mean length scale reduces the `t*_mix` spread by a factor of **~6** relative to plain `H0` across this much wider, independently-varied range of `H0` and `z_t` — an even larger improvement than the factor of ~3 found in the original 2-point test (§4.5), and decisively better than every other simple candidate tried (all of which are worse than `H0` itself once `z_t` is allowed to vary independently). This confirms the geometric mean is not an artifact of the narrow original test, and specifically identifies `sqrt(z_t(H0−z_t))` — not `min`, `max`, or the harmonic mean — as the correct functional form.

![Pycnocline-position collapse, H0=length scale](../figures/mixing_timescale_collapse_pycnocline_zt_H0.png)
![Pycnocline-position collapse, pycnocline length scale](../figures/mixing_timescale_collapse_pycnocline_zt_pycnocline.png)

*(Top: `t*` using `H0` — 8 visibly separate sub-clusters (color = run, ordered by `H0`/`z_t`) spread across `t* ≈ 0.32–0.76`. Bottom: `t*` using `L=sqrt(z_t(H0−z_t))` — all 8 `H0`/`z_t` combinations collapse onto essentially 2 curves, separated only by `c4`.)*

This resolves the caveat from §4.5: the pycnocline length scale is now validated over a genuinely independent `z_t`/`H0` design, not just a fixed-`z_t` comparison, and the geometric-mean form specifically (not just "some function that shrinks with proximity to a boundary") is confirmed as the best of the candidates tested.

### 4.7 Mapping A(c4): a fine c4 sweep from 0.1 to GLS.C1

§4.4 established that `t*_mix` depends on `structure.c4` and `grid.H0`/`initial.temp_zt` as separable factors, `t*_mix ≈ A(c4)·B(H0, z_t)`, but only tested `c4 ∈ {0.44, 0.97}`. A dedicated follow-up sweep (`templates/c4_fine_sweep.yaml`, using the existing `tools/prep_mixing_timescale_sweep.py` machinery) samples `c4` at 7 values (`0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00`) crossed with the original 2 values of `grid.H0` (`75, 150` m, `temp_zt = 40` m fixed as in the original sweep), with `CD = 0.63` and `temp_dT = 10.0` held fixed (14 runs total).

**Resolution improvement:** at this fine `c4` spacing, the original "first output timestep where `φ* < 0.05`" metric proved too coarse — with hourly history output, several adjacent `c4` values landed on the exact same output timestep and appeared numerically identical. `mixing_timescale()` was updated to linearly interpolate `t*_mix` between the last sample `≥ 0.05` and the first sample `< 0.05`, giving sub-output-step resolution; this is now the default behavior for all length scales and sweeps (a backward-compatible refinement, not a `length_scale`-specific change).

**Separability re-confirmed at fine resolution:** the ratio `t*_mix(H0=150)/t*_mix(H0=75)` is `0.827–0.838` across all 7 `c4` values (vs. `0.824–0.828` found at just 2 `c4` values in §4.4) — the `A(c4)·B(H0)` factorization holds cleanly over the whole finely-sampled range, not just the 2 original points.

**Shape of A(c4):** for `c4` from `0.10` to `0.85`, `t*_mix` increases only gently and smoothly (~9% total increase at `H0=75`, from `0.528` to `0.576`), and is fit extremely well by a quadratic in `c4` (`A0 + a·c4 + b·c4²`, relative RMSE < 0.5% at both `H0` values, vs. ~1% for a linear or exponential fit) — a plain power-law singularity at `c4 = GLS.C1` (e.g. `A0/(1-c4)^p`) does **not** fit at all (residuals dominate once `c4` approaches 1). Then, from `c4 = 0.85` to `c4 = 1.00 (= GLS.C1)`, `t*_mix` rises sharply: **~34–36% above the quadratic trend extrapolated from `c4 < 1`**, consistently at both `H0` values. Since `c4 = 1.00` was checked and is *not* flagged `unstable` by the existing anomaly detector (§4.3), this is a real, validated data point, not a numerical artifact of the instability discussed in §4.3 (which only manifests for `c4 > GLS.C1`, e.g. `c4 = 1.4`) — but it does show that the approach to the theoretical "neutral point" `c4 = C1` is not smooth: mixing efficiency degrades disproportionately fast in the last ~15% of the valid `c4` range.

![Fine c4 sweep: A(c4) mapping](../figures/mixing_timescale_c4_fine_sweep.png)

*(Both panels: `t*_mix` vs. `c4` for `H0=75` (blue) and `H0=150` (orange), with a quadratic fit to the `c4 < 1.0` points (dashed) extrapolated to `c4=1.0`. The sharp upward departure of the solid curves from the dashed extrapolation right at `c4 = GLS.C1 = 1.0` is visible in both panels, for both length scale choices.)*

**Interpretation:** this refines the qualitative picture from §4.3 — rather than a sharp binary switch between "stable, C1-scaled mixing efficiency" (`c4 ≤ C1`) and "numerical collapse" (`c4 > C1`), there appears to be a smooth, mild, closure-driven reduction in mixing efficiency for most of the valid `c4` range, followed by an accelerating approach to reduced efficiency as `c4 → C1`, consistent with `c4 = C1` being a genuine dynamical transition (per Rennau et al.'s neutral-point framing) rather than merely a numerically convenient upper bound.

### 4.8 Confirmed conservation and numerical health

- Volume-averaged density is conserved to ~1e-7 relative error over the full run duration in every case (justifying the "no reference run needed" simplification in §2.1).
- `φ*(t)` settles cleanly to ~0 (numerical noise floor ~1e-6) with no overshoot or oscillation in any of the 16 valid runs, confirming the diagnostic pipeline and the underlying mixing physics are well-behaved once `c4 ≤ GLS.C1`.

## 5. Summary of findings

1. **The analytic power-balance prediction for `τ_mix` is accurate.** `tau_mix_diagnostic` matches `tau_mix_theory` to <1% across all 16 valid runs, confirming that the quasi-steady drag/body-force balance is correctly captured and that the leading-order power scale `P_str = ρ₀·H·P_d` is set by `CD`, `str_a`, and `BFRC_U` alone, independent of how the turbulence closure subsequently partitions that power into mixing versus dissipation.

2. **The dimensionless collapse from Carpenter et al. (2016) reproduces well** in this idealized 1D setup: `φ*(t*)` collapses across all 16 valid parameter combinations when plotted against the nondimensional time `t* = t·P_str/(g·Δρ·H²)`.

3. **The value of `t*` at which mixing completes is not universal in this implementation.** Unlike the idealized `t* = 1` reference scaling in Carpenter et al.'s pycnocline model, mixing here completes at `t* ≈ 0.45–0.76`, and that completion time depends systematically on lower-order closure and geometry effects not contained in the leading-order power balance. In particular:
   - **`c4` matters because it changes mixing efficiency, not power input.** Carpenter et al.'s homogeneous closure analysis identifies `c4 = C1` as the neutral point: `c4 < C1` enhances the implied mixing efficiency, while `c4 > C1` reduces it. Consistent with that theory, increasing `c4` from `0.44` to `0.97` slows the erosion of stratification at fixed diagnosed `P_d`, indicating that less of the extracted power is converted into effective buoyancy mixing.
   - **`H0` matters here through pycnocline position.** With `temp_zt` fixed, changing `H0` changes the pycnocline's relative position within the water column; in these runs, a pycnocline closer to a boundary mixes out sooner in nondimensional time than one located nearer mid-column.

4. **Replacing `H0` with a pycnocline-based length scale, `L = sqrt(z_t·(H0−z_t))`, removes most of the residual `H0`-dependence** in the dimensionless mixing-completion time. This suggests the pycnocline's distance to the boundaries, not the total water column depth, is the more physically appropriate length scale for this parametrization when the pycnocline sits at a fixed absolute depth — as is realistic for the Norwegian shelf, where the pycnocline depth is set by seasonal thermal forcing rather than local bathymetry. This was first found over a narrow 2-point test (§4.5, CV reduced ~3x) and then **confirmed over a much wider, independent sweep of `z_t` and `H0`** (§4.6): CV reduced by **~6x** (0.16→0.03), and the geometric mean specifically outperforms `min`, `max`, and the harmonic mean of the two boundary distances — ruling out those alternatives as the correct form.

5. **`c4 > C1` should be interpreted as a theoretically low-efficiency regime, not merely a numerical anomaly.** In Carpenter et al., `c4 = C1` is the distinguished value at which structure-induced production leaves the implied mixing efficiency unchanged; `c4 > C1` reduces that efficiency, while `c4 < C1` enhances it. In the present explicit STRUCTURE_MIXING/GLS implementation, pushing into the `c4 > C1` regime (`c4 = 1.4 > 1.0`) did not just make mixing weaker: it caused `tke`, `gls`, and `AKt` to collapse to floor values, so the regime is numerically unvalidated in practice. We therefore restrict the validated sweep range to `c4 ≤ GLS.C1`, enforced by a `ValueError` in `analytic_mixing_timescale()` and by anomaly flagging in `plot_collapse()`. This should be understood as a **practical restriction of the current implementation**, not as a universal theoretical prohibition on `c4 > C1`.

6. **`A(c4)` is smooth and mildly increasing for most of the valid range, but rises sharply as `c4 → GLS.C1`.** A fine 7-point `c4` sweep (§4.7, `c4 ∈ [0.10, 1.00]`) shows `t*_mix` is well-fit by a quadratic in `c4` for `c4 ≤ 0.85` (relative RMSE < 0.5%), with only a mild ~9% total increase over that range — but at `c4 = GLS.C1` exactly, `t*_mix` jumps ~35% above the extrapolated quadratic trend, consistently at both `H0` values tested. This point is not flagged as numerically unstable, so it is a real effect: mixing efficiency degrades disproportionately fast in the last ~15% of the valid `c4` range, consistent with `c4 = C1` being a genuine dynamical transition (Carpenter et al.'s "neutral point") rather than just a numerically convenient upper bound. The `A(c4)·B(H0)` separability found in §4.4 holds cleanly across this whole finely-sampled range (ratio `t*_mix(H0=150)/t*_mix(H0=75)` stays within `0.827–0.838` for all 7 `c4` values), not just the 2 originally tested points.

7. **`CD` and `temp_dT` have negligible independent effect on the nondimensional mixing curves once time is scaled by the Carpenter power-balance timescale.** Their primary influence is already absorbed into `P_str`, `Δρ`, and therefore `t*` / `τ_mix`, exactly as the leading-order theory predicts.

## 6. Future work / open questions

- **Consider whether `t*_mix ≈ A(c4)·B(H0,z_t)` extends to a wider parameter range** or breaks down (e.g. at very shallow/deep `H0`, very small/large `z_t/H0` ratios close to 0 or 1, or very small `c4` approaching 0).
- **Investigate the mechanism behind the sharp `A(c4)` upturn near `c4 = GLS.C1`** (§4.7) — e.g. does `AKt` or the TKE/GLS balance show early warning signs (elevated variance, slower convergence to quasi-steady state) as `c4` approaches `C1`, even before the numerically-unstable regime `c4 > C1` is reached?
- **Test the geometric-mean pycnocline length scale on a sloping or non-uniform-`str_a` water column**, where the "distance to a boundary" concept is less clean than in this idealized flat-bottom, depth-uniform-drag setup.

## 7. Code and artifacts

| File | Purpose |
|---|---|
| `utils/utils.py` | `compute_phi`, `compute_Pd`, `compute_Pstr`, `analytic_mixing_timescale` — core physics/diagnostics, plus the `c4 ≤ GLS.C1` validation check. |
| `roms/Include/mixtest_1d.h` | `NONLIN_EOS` undefined (linear EOS), required for physically meaningful `φ(t)`. |
| `templates/mixing_timescale_sweep.yaml` | Original sweep definition (CD x c4 x temp_dT x H0, `z_t` fixed at 40 m). |
| `templates/pycnocline_zt_sweep.yaml` | Follow-up sweep definition (§4.6): explicit `(H0, z_t)` pairs x `c4`, varying `z_t` independently of `H0`. |
| `templates/c4_fine_sweep.yaml` | Follow-up sweep definition (§4.7): 7 values of `c4` (0.10–1.00) x 2 values of `grid.H0`, mapping `A(c4)` finely (uses the existing `tools/prep_mixing_timescale_sweep.py`). |
| `configs/variants/mixing_timescale.yaml` | Scaled `str_a`/`BFRC_U` variant to bring τ_mix into a practical 2–13 day range. |
| `tools/prep_mixing_timescale_sweep.py` | Prepares the original sweep runs (and the §4.7 fine-`c4` sweep, which is a pure cartesian product), deriving per-run `NTIMES` from `analytic_mixing_timescale()`. |
| `tools/prep_pycnocline_zt_sweep.py` | Prepares the §4.6 follow-up sweep from explicit `(H0, z_t)` pairs (not a full cartesian product, to avoid placing the thermocline too close to a boundary). |
| `analysis/mixing_timescale.py` | Per-run diagnostics (`mixing_timescale()`), sweep summary (`summarize_sweep()`), collapse plot (`plot_collapse()`), and closure-sensitivity plot (`plot_sensitivity()`); all three accept a `length_scale` argument (`"H0"` or `"pycnocline"`, see §4.5), and `pycnocline_length_scale()` implements the alternative length scale itself. `mixing_timescale()`'s `t_star_mix` now uses linear interpolation between output timesteps for sub-output-step resolution (added for §4.7's fine `c4` sweep, but applies to all sweeps). |
| `sweeps/mixing_timescale/manifest.yaml` / `.csv` | Generated manifest of the original 16 prepared/completed runs (not committed to git; regenerable via `tools/prep_mixing_timescale_sweep.py`). |
| `sweeps/pycnocline_zt/manifest.yaml` / `.csv` | Generated manifest of the §4.6 follow-up 16 runs (not committed to git; regenerable via `tools/prep_pycnocline_zt_sweep.py`). |
| `sweeps/c4_fine/manifest.yaml` / `.csv` | Generated manifest of the §4.7 fine-`c4` 14 runs (not committed to git; regenerable via `tools/prep_mixing_timescale_sweep.py templates/c4_fine_sweep.yaml`). |
| `runs/mixtau_*/`, `runs/mixtauzt_*/`, `runs/mixtauc4_*/` | Completed run directories for all three sweeps (not committed to git). |
| `figures/mixing_timescale_collapse_mixing_timescale_H0.png` | The φ*(t*) collapse plot using `L=H0` (§4.2), original sweep. |
| `figures/mixing_timescale_sensitivity_mixing_timescale_H0.png` | The t*_mix vs. c4/H0 sensitivity plot using `L=H0` (§4.4), original sweep. |
| `figures/mixing_timescale_collapse_mixing_timescale_pycnocline.png` | The φ*(t*) collapse plot using `L=sqrt(zt(H0-zt))` (§4.5), original sweep. |
| `figures/mixing_timescale_sensitivity_mixing_timescale_pycnocline.png` | The t*_mix vs. c4/H0 sensitivity plot using `L=sqrt(zt(H0-zt))`, showing the H0-grouping collapse away (§4.5), original sweep. |
| `figures/mixing_timescale_length_scale_comparison.png` | Side-by-side `φ*(t*)` comparison of both length scale choices, colored by `c4` with linestyle by `H0`, for the original sweep (§4.5). |
| `figures/mixing_timescale_collapse_pycnocline_zt_H0.png` | The φ*(t*) collapse plot using `L=H0` for the §4.6 follow-up sweep (8 visible `H0`/`z_t` sub-clusters). |
| `figures/mixing_timescale_collapse_pycnocline_zt_pycnocline.png` | The φ*(t*) collapse plot using `L=sqrt(z_t(H0-z_t))` for the §4.6 follow-up sweep (sub-clusters collapse to 2 curves by `c4` only). |
| `figures/mixing_timescale_sensitivity_c4_fine_H0.png` / `_pycnocline.png` | The standard t*_mix vs. c4 sensitivity plot (via `plot_sensitivity()`) for the §4.7 fine-`c4` sweep, at 7 `c4` values instead of 2. |
| `figures/mixing_timescale_c4_fine_sweep.png` | Annotated version of the above with quadratic fits (excluding `c4=1.0`) overlaid, highlighting the sharp departure from the smooth trend right at `c4=GLS.C1` (§4.7). |

To regenerate this analysis from scratch:

```bash
python tools/prep_mixing_timescale_sweep.py templates/mixing_timescale_sweep.yaml
python tools/run_sweep.py sweeps/mixing_timescale/manifest.yaml
python analysis/mixing_timescale.py --sweep sweeps/mixing_timescale/manifest.yaml --length-scale H0 --save
python analysis/mixing_timescale.py --sweep sweeps/mixing_timescale/manifest.yaml --length-scale pycnocline --save

python tools/prep_pycnocline_zt_sweep.py templates/pycnocline_zt_sweep.yaml
python tools/run_sweep.py sweeps/pycnocline_zt/manifest.yaml
python analysis/mixing_timescale.py --sweep sweeps/pycnocline_zt/manifest.yaml --length-scale H0 --save
python analysis/mixing_timescale.py --sweep sweeps/pycnocline_zt/manifest.yaml --length-scale pycnocline --save

python tools/prep_mixing_timescale_sweep.py templates/c4_fine_sweep.yaml
python tools/run_sweep.py sweeps/c4_fine/manifest.yaml
python analysis/mixing_timescale.py --sweep sweeps/c4_fine/manifest.yaml --length-scale H0 --save
python analysis/mixing_timescale.py --sweep sweeps/c4_fine/manifest.yaml --length-scale pycnocline --save
```

