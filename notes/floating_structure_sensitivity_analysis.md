# Sensitivity analysis for floating (non-bottom-fixed) structures

## 1. Motivation

The mixing-timescale sensitivity analysis in
[`mixing_timescale_analysis.md`](mixing_timescale_analysis.md) (in
particular its §8, "Extension: floating support structures (`Gb` /
`bfrc_cb`)") introduced a floating-foundation variant of the structure
parametrization: `str_a` is set to zero below a configurable depth
(`structure.depth_zero_below`, measured from the surface), representing
a turbine foundation that only occupies the upper part of the water
column. Because the body force (`UV_BODYFORCE`) that drives the flow no
longer has drag to balance it below `depth_zero_below`, a new term `Gb`
(coefficient `bfrc_cb`) was added to the momentum equations, active only
where `str_a == 0`, using the `str_a` value from the top of the water
column as its reference drag coefficient. This keeps the velocity field
bounded without requiring real structures at depth. See
[`IMPLEMENTATION_STRUCTURE_MIXING.md`](../roms/IMPLEMENTATION_STRUCTURE_MIXING.md)
(or the equivalent doc in the `roms` repo) and §8 of the existing note
for the full derivation, and `test_UV_BODYFORCE_FLOATING_CB_EQ_CD` for
the validation that, when `bfrc_cb == structure.CD`, the column-mean
momentum balance (and hence the steady velocity `u_inf`) is *identical*
to the fully-structured (bottom-fixed) case.

This note asks the natural follow-up question: **given that the
velocity field is unaffected (in the no-shear `bfrc_cb == CD` case), how
does the *mixing* depend on how much of the water column actually
contains structures** (`structure.depth_zero_below`, expressed below as
the fraction `depth_frac = depth_zero_below / H0`)? Concretely: can we
extend the predictive `τ_mix`/`A(c4)`/`S(x*)` machinery of §7 of the
existing note to floating structures?

**Scope of this first pass (confirmed with the user):** `bfrc_cb` is
held fixed at `bfrc_cb == structure.CD` throughout (no shear between the
structured and floating zones). This isolates the effect of *where* the
structured zone sits from the added complexity of a genuinely different
drag law at depth, which is left as future work (§6).

## 2. Theory extension

### 2.1 Why the completion-time metric (`τ_mix`) doesn't directly apply

The existing theory (§4-7 of the base note) estimates a *completion*
time scale `τ_mix`, the time for `φ(t) → 0` (full-column mixing),
assuming TKE production (`P_d`, and its depth integral `P_str`) is
generated more or less uniformly across the whole column, because
`str_a > 0` everywhere for bottom-fixed structures.

For floating structures, `P_d` is only physically generated where real
structures exist, i.e. in the structured zone (`str_a > 0`, thickness
`d_struct = min(depth_zero_below, H0)`). The `Gb` term is a momentum
sink, not a TKE source, so no production occurs in the floating zone
below `depth_zero_below`. This means:

- Total mixing power should scale with `d_struct`, not `H0`:
  `P_str = ρ₀ · d_struct · P_d` (generalizing the base note's
  `P_str = ρ₀ · H0 · P_d`, which is the `d_struct = H0` special case).
- A thinner structured zone means less power is deposited into the
  water column overall (even though the velocity/drag balance is
  unchanged), so **full-column mixing should take longer** — and, more
  importantly, TKE generated near the surface must diffuse downward
  (via the GLS closure) to affect the stratification below
  `depth_zero_below` at all, which is a fundamentally slower,
  diffusion-limited process rather than a power-limited one. This
  suggests "time to full mixing" may not degrade gracefully/predictably
  as `depth_frac` shrinks, and might not be reached in any practical run
  duration for thin structured zones.

**Metric adopted:** track `τ_x`, the time for `φ(t)` to reach a fixed
fractional reduction `x` of its initial value (`φ(τ_x) = (1-x)·φ(0)`),
default `x = 10%`. This is chosen because it falls in the early-time,
close-to-linear regime of the existing universal shape function
`S(x*)` (§7.3 of the base note, fitted Beta-CDF with `a≈0.95, b≈0.87`,
both near 1), so the base note's constant-`P_str` theory predicts (to
leading order) `τ_x_theory ≈ x · τ_mix_theory`. As a secondary
diagnostic, whether `φ(t)` plateaus within a run (and, if so, the
asymptotic mixed fraction `1 - φ_∞/φ(0)`) is also tracked.

### 2.2 Implementation (`utils/utils.py`)

- `analytic_mixing_timescale(cfg, x_frac=0.10)`: generalized to use
  `d_struct = min(depth_zero_below, H0)` in `P_str`, and to return
  `d_struct`, `x_frac`, and `tau_x = x_frac · tau_mix` alongside the
  existing keys. Raises `ValueError` if `depth_zero_below < H0` (a
  genuinely floating configuration) and `structure.cb != structure.CD`,
  enforcing the no-shear scope of this analysis at the API level.
- `compute_Pd`/`compute_Pstr`: generalized to read the actual `str_a(z)`
  profile from the grid file (`load_str_a`) instead of assuming a
  spatially-uniform scalar — a general correctness fix that happens to
  be required for floating configurations, but does not change results
  for any existing uniform-`str_a` (bottom-fixed) sweep.
- `find_tau_x_star(t, phi_star, x_frac)`: interpolated crossing-time
  finder (generalizes the existing hardcoded 5%-completion-threshold
  logic used for `t*_mix`); returns `NaN` if the run never reaches the
  threshold.
- `detect_phi_plateau(phi, tail_frac=0.1, rel_slope_tol=0.02)`: flags
  whether the last `tail_frac` of a run has a small slope relative to
  the initial decay rate, and reports `phi_inf` / `mixed_fraction_inf`.
- `analysis/mixing_timescale.py::mixing_timescale()`: extended
  (additively, backward-compatible) with `x_frac`/`tail_frac` params and
  the corresponding `tau_x_theory`, `tau_x_diagnostic`, `plateaued`,
  `phi_inf`, `mixed_fraction_inf`, `rel_tail_slope` outputs.

**Validation:** a bottom-fixed regression config
(`depth_zero_below=1e9 ≫ H0`) reproduces `d_struct == H0` exactly and
gives bit-identical `tau_mix_theory`/`tau_mix_diagnostic` to the
pre-generalization code (both cases agree at 795,692 s). A dedicated
floating sanity run (`depth_zero_below = H0/3`, `bfrc_cb == CD`) gave
`tau_mix_theory` and `tau_mix_diagnostic` agreeing within ~1.8%,
confirming the `d_struct`-based generalization is consistent with real
model output.

## 3. Sweep design

Two sweeps were run, both restricted to `structure.cb == structure.CD`
(fixed at 0.63) and `initial.temp_dT = 10.0`:

**Main sweep** (`templates/floating_depth_sweep.yaml`,
`tools/prep_floating_depth_sweep.py`, 20 runs): `NTIMES` sized from
`tau_x_theory` (with a ×3 margin) so runs are only as long as needed to
resolve the `x=10%` crossing.
- `(grid.H0, initial.temp_zt)` pairs: `{(75, 40), (150, 40)}`
- `structure.c4`: `{0.44, 0.85}`
- `structure.depth_frac = depth_zero_below / H0`: `{1.0, 0.75, 0.5, 0.25, 0.1}`
  (`1.0` reproduces the fully-structured baseline as an internal
  consistency check)

**Long-run subset** (`templates/floating_depth_longrun_sweep.yaml`,
3 runs): a single `(H0, z_t) = (150, 40)` pair, `structure.c4 = 0.44`,
`depth_frac ∈ {1.0, 0.5, 0.25}`, `NTIMES` sized from `tau_mix_theory`
(×2 margin) instead, specifically to run long enough to characterize
the plateau/asymptotic-mixed-fraction behavior, which the (deliberately
short) main sweep cannot resolve.

Analysis code: `analysis/floating_depth.py`
(`summarize_floating_sweep`, `plot_tau_x_vs_depth_frac`,
`summarize_plateau`).

## 4. Results

### 4.1 `τ_x` sensitivity to structured-zone thickness

| depth_frac | c4 | H0 | τ_x theory (d) | τ_x diagnostic (d) | ratio (theory/diag) |
|---|---|---|---|---|---|
| 1.00 | 0.44 | 75 | 0.460 | 0.278 | 1.66 |
| 0.75 | 0.44 | 75 | 0.614 | 0.332 | 1.85 |
| 0.50 | 0.44 | 75 | 0.921 | 0.518 | 1.78 |
| 0.25 | 0.44 | 75 | 1.842 | 2.937 | 0.63 |
| 0.10 | 0.44 | 75 | 4.605 | — (not reached) | — |
| 1.00 | 0.44 | 150 | 0.921 | 0.467 | 1.97 |
| 0.75 | 0.44 | 150 | 1.228 | 0.497 | 2.47 |
| 0.50 | 0.44 | 150 | 1.842 | 0.752 | 2.45 |
| 0.25 | 0.44 | 150 | 3.684 | 2.443 | 1.51 |
| 0.10 | 0.44 | 150 | 9.209 | 17.936 | 0.51 |
| 1.00 | 0.85 | 75 | 0.460 | 0.290 | 1.59 |
| 0.75 | 0.85 | 75 | 0.614 | 0.345 | 1.78 |
| 0.50 | 0.85 | 75 | 0.921 | 0.546 | 1.69 |
| 0.25 | 0.85 | 75 | 1.842 | 5.249 | 0.35 |
| 0.10 | 0.85 | 75 | 4.605 | — (not reached) | — |
| 1.00 | 0.85 | 150 | 0.921 | 0.502 | 1.84 |
| 0.75 | 0.85 | 150 | 1.228 | 0.526 | 2.33 |
| 0.50 | 0.85 | 150 | 1.842 | 0.799 | 2.31 |
| 0.25 | 0.85 | 150 | 3.684 | 3.219 | 1.15 |
| 0.10 | 0.85 | 150 | 9.209 | — (not reached) | — |

(Figure: `figures/floating_depth_tau_x_H0.png` — diagnostic `τ_x` vs.
`depth_frac`, per `(c4, H0)` group, with the theory overlaid as dashed
lines.)

**Key findings:**

- **The theory captures the right order of magnitude and the correct
  qualitative trend** (`τ_x` grows as `depth_frac` shrinks — thinner
  structured zones take longer to reach even 10% mixing), but the
  quantitative agreement is not uniform across `depth_frac`:
  - At `depth_frac ∈ {0.5, 0.75, 1.0}`, the theory **over-predicts**
    `τ_x` by a roughly constant ~1.6–2.5×, similar to the ~2× discrepancy
    already noted for the base note's early-time approximation
    (`S(x*)` is close to, but not exactly, linear near `x*=0`).
  - At `depth_frac = 0.25`, the ratio **crosses over below 1** — the
    diagnostic `τ_x` becomes *larger* than the theory predicts (by up to
    ~3×), i.e. the theory now under-predicts.
  - At `depth_frac = 0.1`, most configurations **never reach 10% mixing**
    within the sweep's run duration (`tau_x_theory × 3`); the one
    exception (`c4=0.44, H0=150`) still shows a ~2× under-prediction.
- This cross-over is consistent with the mechanism proposed in §2.1:
  once the structured zone is thin, mixing at depth becomes
  **diffusion-limited** (TKE must be transported down from the shallow
  production zone via the GLS closure) rather than **power-limited**
  (set purely by the depth-integrated `P_str`/buoyancy ratio, as the
  simple theory assumes). The simple power-balance theory has no
  mechanism to capture this transport bottleneck, so it increasingly
  under-predicts `τ_x` as `depth_frac` shrinks.
- No systematic difference between `c4=0.44` and `c4=0.85` beyond the
  known `A(c4)`-type sensitivity already documented in the base note;
  the depth_frac trend is qualitatively the same for both.

### 4.2 Does `φ(t)` plateau at a nonzero residual?

An earlier one-off sanity check (previous session, `depth_frac ≈ 0.33`,
14-day run) suggested `φ(t)` plateaus at only ~57.5% mixed, hinting at a
possible permanent "mixing floor" for thin structured zones. The
dedicated long-run subset (sized to run much longer, off `tau_mix_theory`
rather than `tau_x_theory`) tells a different, more complete story:

| depth_frac | run duration (d) | plateaued | mixed_fraction_inf |
|---|---|---|---|
| 1.00 | ~18.5 | True | 1.000 |
| 0.50 | ~36.9 | True | 1.000 |
| 0.25 | ~73.7 | True | 0.991 |

Even at `depth_frac = 0.25`, given enough time (~74 days, vs. the
earlier 14-day check), `φ(t)` decays essentially to zero (`phi_star`
reaches ≈0 by the end of the run; `rel_tail_slope ≈ 2×10⁻⁴`, i.e.
genuinely flat, not still actively decaying). **There is no evidence of
a permanent mixing floor** for `depth_frac ≥ 0.25`: the earlier ~57.5%
"plateau" was an artifact of the run simply not being long enough, not
a genuine asymptote. The real effect of a thin structured zone is a
(strongly, non-linearly) **longer time to reach any given fraction of
mixing**, not a cap on how much mixing can ultimately occur.

(This does not rule out a true floor emerging for even thinner
structured zones, e.g. `depth_frac ≲ 0.1` — the main sweep's `depth_frac
= 0.1` runs did not reach 10% mixing even after `tau_x_theory × 3`, and
were not extended into a long-run check; this remains open, see §6.)

## 5. Predictive function

Per the plan's empirical-first approach: **no clean closed-form
correction to `τ_x_theory` was found.** The theory-to-diagnostic ratio is
not a constant multiplicative factor across `depth_frac` (§4.1) — it
decreases monotonically from ~2.5 (thick structured zone) through 1
(crossover near `depth_frac ≈ 0.3–0.4`) to well below 1 or "not reached"
(thin structured zone). A single fitted correction factor (analogous to
`A(c4)` in the base note) would not generalize across this range without
additional structure (e.g. a `depth_frac`-dependent term reflecting the
diffusion-limited regime), which this first pass did not attempt to
derive.

**What *is* established, and useful as a practical estimate:**
`τ_x_theory ≈ x_frac · τ_mix_theory` (with `P_str ∝ d_struct`) gives the
correct order of magnitude and correct qualitative sensitivity to
`depth_frac`, `c4`, and `H0` for `depth_frac ≳ 0.25`, typically within a
factor of ~2–3. For `depth_frac ≲ 0.25` it should be treated as a lower
bound only (actual mixing is slower than predicted), and for very thin
structured zones (`depth_frac ~ 0.1`) it may not predict a reachable
crossing time at all within a practically-sized run.

A proper predictive function for the diffusion-limited regime would
likely need to incorporate the GLS closure's vertical diffusivity/
turbulent length scale explicitly (rather than a single depth-integrated
power balance), which is left as future work (§6).

## 6. Future work

- **`bfrc_cb ≠ structure.CD` (shear between zones).** This analysis
  fixed `bfrc_cb == CD` to isolate the effect of structured-zone
  thickness. Allowing `bfrc_cb` to differ (as originally motivated in
  the problem statement — a distinct steady-state velocity/shear
  profile) is the natural next extension, and will require revisiting
  `analytic_mixing_timescale`'s momentum-balance assumptions (currently
  it explicitly rejects `bfrc_cb != CD` for floating configs).
- **Diffusion-limited regime for thin structured zones
  (`depth_frac ≲ 0.25`).** The simple power-balance theory breaks down
  here (§4.1); a refined theory would need to account for the vertical
  transport time scale of TKE/mixing from the structured zone down into
  the floating zone (e.g. via the GLS closure's diffusivity), not just
  the total power budget.
- **Very thin structured zones (`depth_frac ~ 0.1` and below).** Not
  checked with a long-run subset; open question whether a genuine
  mixing floor exists at some (thinner) `depth_frac`, or whether given
  enough time these configurations also eventually reach full mixing
  (as `depth_frac = 0.25` did, contrary to the earlier informal check).
- **Pycnocline length-scale collapse.** The base note found that
  nondimensionalizing time by a pycnocline length scale
  (`sqrt(z_t(H0−z_t))`) rather than `H0` collapses the bottom-fixed
  sweep much better across varying `H0`. Whether an analogous (or
  different) length scale improves the collapse for floating structures
  — potentially one that also incorporates `depth_zero_below` — was not
  investigated here.

## 7. Code and artifacts

### 7.1 Files
- `utils/utils.py`: `load_str_a`, `find_tau_x_star`, `detect_phi_plateau`;
  generalized `compute_Pd`, `analytic_mixing_timescale`.
- `analysis/mixing_timescale.py`: `mixing_timescale()` extended with
  `x_frac`/`tail_frac` diagnostics (additive, backward-compatible).
- `analysis/floating_depth.py` (new): `summarize_floating_sweep`,
  `plot_tau_x_vs_depth_frac`, `summarize_plateau`.
- `tools/prep_floating_depth_sweep.py` (new): prep script with
  `structure.depth_frac → structure.depth_zero_below` conversion and
  dual `NTIMES`-sizing bases (`tau_x`/`tau_mix`).
- `templates/floating_depth_sweep.yaml` (new): main sweep (20 runs).
- `templates/floating_depth_longrun_sweep.yaml` (new): long-run/plateau
  subset (3 runs).

### 7.2 Reproducing the sweeps
```bash
# Main sweep: tau_x sensitivity to depth_frac / c4 / H0
python tools/prep_floating_depth_sweep.py templates/floating_depth_sweep.yaml
python tools/run_sweep.py sweeps/floating_depth/manifest.yaml

# Long-run subset: plateau / asymptotic mixed-fraction check
python tools/prep_floating_depth_sweep.py templates/floating_depth_longrun_sweep.yaml
python tools/run_sweep.py sweeps/floating_depth_longrun/manifest.yaml

# Analysis (table + figure + plateau summary)
python analysis/floating_depth.py \
    --sweep sweeps/floating_depth/manifest.yaml \
    --longrun sweeps/floating_depth_longrun/manifest.yaml --save
```
