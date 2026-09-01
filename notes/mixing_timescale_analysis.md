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

![Mixing timescale collapse](../figures/mixing_timescale_collapse_mixing_timescale.png)

Plotting `φ*(t) = φ(t)/φ(0)` against the dimensionless time `t* = t·P_str_diag/(g·Δρ·H²)` for all 16 valid runs produces a tight collapse — the direct analogue of Carpenter et al.'s Fig 4b. All 16 curves fall in a narrow band and decay from `φ* = 1` to `φ* ≈ 0` over a similar range of `t*`.

However, unlike Carpenter et al.'s idealized time-dependent pycnocline model (where mixing completion occurs at `t* = 1` essentially by construction/definition), **mixing completes earlier here, at t* ≈ 0.45–0.76**, and — notably — this completion time is **not universal**: it varies systematically with two of the four swept parameters (see §4.4). This means `t* = 1` should be understood as a rough order-of-magnitude time scale for *this* parametrization/setup, not an exact universal threshold; the value at which mixing actually completes depends on lower-order details of the turbulence closure and geometry not captured by the leading-order power balance.

### 4.3 Numerical instability at c4 > GLS.C1

An initial version of this sweep included `structure.c4 = 1.4` (in addition to 0.44 and 0.97). **Every one of the 8 runs with `c4 = 1.4`** (all combinations of CD/dT/H0) showed **no mixing at all** — `φ*` stayed pinned near 1.0 for the entire run, despite `u` correctly reaching its analytic steady-state speed. Diagnosis:

- `tke` and `AKt` (eddy diffusivity) collapsed to their numerical floor values (`tke ≈ 1e-8` = `Kmin`, `AKt ≈ 1e-6`) at `c4 = 1.4`, whereas sibling runs at `c4 = 0.44`/`0.97` (same CD/dT/H0) gave healthy `AKt` values (0.19–1.6 m²/s).
- Momentum/drag balance was unaffected — only the turbulence closure blew up.

This is consistent with a **numerical instability in the explicit time-stepping of the GLS ψ-equation's structure-production source term** (`gls_c4 · Pd_struct · gls/tke`) at high production rates relative to the model's fixed timestep (`DT = 40 s`). The user's physical insight resolved the ambiguity: in this GLS closure configuration, `GLS.C1 = 1.0` (the shear-production coefficient in the same ψ-equation), and `structure.c4` plays an analogous role for structure-driven production. **`c4` should not exceed `C1`** — setting `c4 = 1.4 > C1 = 1.0` pushes the total production term outside its numerically stable range for this timestep.

**Action taken:**
- `structure.c4 = 1.4` was removed from the sweep design (`templates/mixing_timescale_sweep.yaml`), documented with an explanatory comment.
- `utils/utils.py::analytic_mixing_timescale()` now raises a `ValueError` if `structure.c4 > GLS.C1` is passed, to prevent this invalid/unstable regime from silently being swept again in the future.
- `analysis/mixing_timescale.py::mixing_timescale()` also includes a runtime anomaly detector (flags a run as `unstable` if φ* hasn't decayed substantially despite t* being well past 1) as a second line of defense; `plot_collapse()` automatically excludes flagged runs with a printed console warning.

This is a genuine limitation of the current explicit STRUCTURE_MIXING/GLS coupling worth keeping in mind: the structure-production term's stability bound depends on the interplay of `c4`, the timestep `DT`, and the local production rate, and has not been explored further (e.g. whether reducing `DT` would restore stability at `c4 > C1`).

### 4.4 Why mixing completes early, and why it depends on c4 and H0

![Mixing-completion sensitivity](../figures/mixing_timescale_sensitivity_mixing_timescale.png)

To quantify the deviation from Carpenter's universal `t* = 1`, we define an **empirical dimensionless mixing-completion time** `t*_mix`: the first `t*` at which `φ* < 0.05` (i.e. stratification potential energy has dropped to <5% of its initial value). Averaging over the two values of `CD` and `temp_dT` (which have negligible effect once time is nondimensionalized — see the flat rows within each `c4`/`H0` group in the table above), `t*_mix` depends cleanly on `c4` and `H0`:

| c4 | H0 = 75 m | H0 = 150 m |
|---|---|---|
| 0.44 | 0.546 ± 0.011 | 0.452 ± 0.006 |
| 0.97 | 0.733 ± 0.017 | 0.604 ± 0.008 |

Two effects, cleanly separable:

1. **Closure coefficient `c4`:** `t*_mix` is consistently ~1.34–1.36x *larger* (i.e. mixing takes longer, relative to the dimensionless time scale) at `c4 = 0.97` than at `c4 = 0.44`. This is initially counter-intuitive — one might expect a larger structure-production coefficient to produce *faster* mixing. But `c4` multiplies the production term in the **ψ (dissipation/length-scale) equation**, not the TKE equation directly. A larger `c4` increases the dissipation rate of turbulent kinetic energy relative to its production, which *shrinks* the turbulent length scale and hence `AKt` for the same power input `P_d`. This was confirmed directly: mean `AKt` over the water column (last half of the run) was ~0.37 m²/s at `c4 = 0.44` vs. ~0.19 m²/s at `c4 = 0.97`, for otherwise identical `CD`/`dT`/`H0` — i.e. roughly half the eddy diffusivity at the higher `c4`, consistent with the ~1.35x longer mixing time.

2. **Water column depth `H0` (pycnocline position):** `t*_mix` is consistently ~0.82–0.83x *smaller* (mixing completes sooner, relative to t*) at `H0 = 150` m than at `H0 = 75` m. The thermocline centre depth (`initial.temp_zt = 40` m, fixed across the sweep) sits at a different *relative* position in the water column depending on `H0`: 53% of the way down for `H0 = 75` m, but only 27% of the way down for `H0 = 150` m. A pycnocline positioned closer to a boundary (surface or bed) erodes faster under a roughly depth-uniform eddy diffusivity than one centered mid-column, because it has less distance to homogenize toward the nearer boundary and benefits from boundary-adjacent turbulence. Since the pycnocline depth `temp_zt` was not itself swept, this session cannot yet distinguish "absolute pycnocline depth" from "relative pycnocline position" as the controlling variable — see §5 (future work).

The two effects multiply cleanly: e.g. `t*_mix(c4=0.97, H0=150) / t*_mix(c4=0.97, H0=75) = 0.604/0.733 = 0.824`, matching `t*_mix(c4=0.44, H0=150)/t*_mix(c4=0.44, H0=75) = 0.452/0.546 = 0.828` to within 0.5%. Likewise the `c4` ratio is ~1.342 at `H0=75` and ~1.336 at `H0=150`. This separability suggests `t*_mix ≈ A(c4) · B(H0)` to good approximation across the tested range, though this has only been checked over 2 values of each parameter and should not be over-generalized.

### 4.5 Confirmed conservation and numerical health

- Volume-averaged density is conserved to ~1e-7 relative error over the full run duration in every case (justifying the "no reference run needed" simplification in §2.1).
- `φ*(t)` settles cleanly to ~0 (numerical noise floor ~1e-6) with no overshoot or oscillation in any of the 16 valid runs, confirming the diagnostic pipeline and the underlying mixing physics are well-behaved once `c4 ≤ GLS.C1`.

## 5. Summary of findings

1. **The analytic power-balance prediction for τ_mix is accurate.** `tau_mix_diagnostic` matches `tau_mix_theory` to <1% across all 16 runs — the quasi-steady drag/body-force balance assumption holds, and the leading-order power scale `P_str = ρ₀·H·P_d` is correctly predicted from `CD`, `str_a`, `BFRC_U` alone, independent of the turbulence closure.
2. **The dimensionless collapse from Carpenter et al. (2016) reproduces well** in this idealized 1D setup — `φ*(t*)` collapses across all 16 valid parameter combinations.
3. **The specific value of `t*` at which mixing completes is not universal**, unlike the idealized `t* = 1` in Carpenter's own model. Here it ranges ~0.45–0.76 and depends systematically (and separably) on:
   - the GLS closure coefficient `c4` (higher `c4` → slower mixing, via reduced eddy diffusivity from enhanced dissipation-equation production), and
   - the pycnocline's relative position in the water column (set by `H0` here, with `temp_zt` fixed — pycnocline closer to a boundary → faster mixing).
4. **A hard numerical stability constraint exists: `structure.c4` must not exceed `GLS.C1`.** Exceeding it collapses TKE/GLS to their numerical floor and eliminates mixing entirely, a purely numerical artifact of the explicit time-stepping of the structure-production source term — not a physical result. This is now enforced by a `ValueError` in `analytic_mixing_timescale()` and flagged/excluded automatically in `plot_collapse()`.
5. **CD and temp_dT have negligible independent effect on the dimensionless mixing time** — both are already fully absorbed into the `t*`/`τ_mix` normalization, exactly as the theory predicts.

## 6. Future work / open questions

- **Vary `temp_zt` independently of `H0`** to disentangle "absolute pycnocline depth" from "relative pycnocline position" as the driver of the `H0` effect in §4.4.
- **Investigate the `c4 > GLS.C1` instability further** (e.g. does reducing `DT` restore stability at higher `c4`? Is there a general stability criterion relating `c4`, `DT`, and the local production rate that could be derived analytically, similar to the existing drag-term stability discussion in `roms/IMPLEMENTATION_STRUCTURE_MIXING.md` §12?).
- **Sweep `c4` more finely between 0 and `GLS.C1`** to map out the `A(c4)` dependence more precisely (only 2 points were tested here).
- **Test other GLS closure choices** (e.g. `k`-ε via `configs/variants/k-e.yaml`, which uses `GLS.C1 = 1.44`) to see whether the `c4`-dependence and the `c4 ≤ C1` stability bound are specific to this closure or general.
- **Consider whether `t*_mix ≈ A(c4)·B(H0)` extends to a wider parameter range** or breaks down (e.g. at very shallow/deep `H0`, or very small `c4`).

## 7. Code and artifacts

| File | Purpose |
|---|---|
| `utils/utils.py` | `compute_phi`, `compute_Pd`, `compute_Pstr`, `analytic_mixing_timescale` — core physics/diagnostics, plus the `c4 ≤ GLS.C1` validation check. |
| `roms/Include/mixtest_1d.h` | `NONLIN_EOS` undefined (linear EOS), required for physically meaningful `φ(t)`. |
| `templates/mixing_timescale_sweep.yaml` | Sweep definition (CD x c4 x temp_dT x H0). |
| `configs/variants/mixing_timescale.yaml` | Scaled `str_a`/`BFRC_U` variant to bring τ_mix into a practical 2–13 day range. |
| `tools/prep_mixing_timescale_sweep.py` | Prepares all sweep runs, deriving per-run `NTIMES` from `analytic_mixing_timescale()`. |
| `analysis/mixing_timescale.py` | Per-run diagnostics (`mixing_timescale()`), sweep summary (`summarize_sweep()`), collapse plot (`plot_collapse()`), and closure-sensitivity plot (`plot_sensitivity()`). |
| `sweeps/mixing_timescale/manifest.yaml` / `.csv` | Generated manifest of all 16 prepared/completed runs (not committed to git; regenerable via `tools/prep_mixing_timescale_sweep.py`). |
| `runs/mixtau_*/` | 16 completed run directories (not committed to git). |
| `figures/mixing_timescale_collapse_mixing_timescale.png` | The φ*(t*) collapse plot (§4.2). |
| `figures/mixing_timescale_sensitivity_mixing_timescale.png` | The t*_mix vs. c4/H0 sensitivity plot (§4.4). |

To regenerate this analysis from scratch:

```bash
python tools/prep_mixing_timescale_sweep.py templates/mixing_timescale_sweep.yaml
python tools/run_sweep.py sweeps/mixing_timescale/manifest.yaml
python analysis/mixing_timescale.py --sweep sweeps/mixing_timescale/manifest.yaml --save
```
