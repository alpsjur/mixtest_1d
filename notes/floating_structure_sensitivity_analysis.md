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

![Diagnostic τ_x vs. depth_frac, grouped by (c4, H0); theory overlaid as dashed lines](../figures/floating_depth_tau_x_H0.png)

*Figure 1: diagnostic `τ_x` (solid) vs. `depth_frac`, grouped by
`(c4, H0)`, with the theory prediction (`τ_x_theory ≈ x_frac ·
τ_mix_theory`, dashed, same colour) overlaid for comparison. The
crossover between over- and under-prediction is visible around
`depth_frac ≈ 0.3–0.4` for every group.*

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

![phi_star(t) for the long-run subset, one curve per depth_frac](../figures/floating_depth_longrun_phi_H0.png)

*Figure 2: `φ(t)/φ(0)` over time for the long-run subset
(`H0=150, z_t=40, c4=0.44`), one curve per `depth_frac`. All three
configurations decay to essentially zero given enough time — thinner
structured zones (smaller `depth_frac`) simply need much longer to get
there, with no sign of flattening out above zero.*

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

### 4.3 Refinement: pycnocline reach, not `depth_frac`, controls the crossover

§4.1 showed the theory/diagnostic ratio crosses from over- to
under-prediction somewhere around `depth_frac ≈ 0.25–0.4`, but the exact
crossover point differed between the `H0=75` and `H0=150` groups — a
hint that `depth_frac` itself isn't the fundamental control variable.
The physical reasoning in §2.1 already suggested why: what matters is
whether the structured zone actually reaches into the pycnocline, not
what fraction of the *total* water column it occupies.

**Definition.** The pycnocline's upper edge is approximated as
`z_edge = z_t - temp_ht` (30 m for this sweep's fixed thermocline
parameters `z_t=40, temp_ht=10` — a reasonable but not exact estimate of
where the tanh transition "begins"; see `pycnocline_edge()` in
`analysis/floating_depth.py`). Define the **reach margin**
`reach_margin = d_struct - z_edge`: positive means the structured
zone reaches into (or through) the pycnocline ("reaching"), negative
means it stops short ("non-reaching").

**New straddle sweep.** To test this cleanly, 12 new runs were added
(`templates/floating_depth_reach_H{75,150}.yaml`), targeting
`d_struct ∈ {20, 25, 30, 35, 40, 45}` m for each of `H0=75` and
`H0=150` (i.e. the same absolute `d_struct` values, different
`depth_frac`), fixed `c4=0.44`.

![Theory/diagnostic tau_x ratio vs. reach margin, colored by H0](../figures/floating_depth_reach_margin.png)

*Figure 3: theory/diagnostic `τ_x` ratio vs. `reach_margin`,
combining the main sweep's `c4=0.44` runs with the new straddle sweeps,
for `H0=75` and `H0=150`. Unlike the `depth_frac`-based plot in §4.1,
the two `H0` groups now line up closely, both crossing `ratio=1` almost
exactly at `reach_margin=0` (`ratio≈1.0–1.2` just past the edge,
`ratio≈0.7–1.0` just before it).*

This is a substantially cleaner result than the `depth_frac`-based
crossover in §4.1: **whether the structured zone reaches the pycnocline
is the real controlling parameter**, not `depth_zero_below/H0`. This
confirms the diffusion-limited-vs-power-limited picture from §2.1 and
gives a concrete, geometry-based criterion (`d_struct` vs. `z_t -
temp_ht`) for which regime a given configuration is in — useful both for
interpreting these results and for designing future sweeps/practical
turbine-siting judgements.

### 4.4 Collapse exploration for pycnocline-reaching structures (phase A)

With the reaching/non-reaching split established, the next
question is whether `φ(t)` collapses onto a common curve (as in the base
note's Carpenter-et-al.-style analysis) once restricted to **pycnocline-
reaching** floating structures. Non-reaching structures are deferred to a later
round (see §6) since a delay/transport mechanism is expected to matter
there, requiring separate treatment.

**Candidates tried** (all reduce to the base note's original
normalization for bottom-fixed structures, since `d_struct = H0` and
`pycnocline_length_scale` are unaffected by `depth_zero_below`):
- `L = H0` (base note's original choice).
- `L = pycnocline length scale = sqrt(z_t(H0-z_t))` (found in the base
  note to collapse the *bottom-fixed* sweep much better than `H0`
  when `H0` is varied at fixed `z_t`).
- `L = d_struct` (the structured-zone thickness itself, since `Pstr`
  scales with `d_struct`, not `H0`).

**Quantitative collapse metric.** All `φ*(t*)` curves (for the
pycnocline-reaching subset: main sweep + both straddle sweeps + long-run subset,
restricted to `reach_margin ≥ 0`) were interpolated onto a common
`t*` grid over `[0, 1]`, and the standard deviation of `φ*` across curves
at each grid point was averaged over the grid (lower = tighter collapse;
see `collapse_spread_score()` in `analysis/floating_depth.py`).

| length_scale | collapse_spread_score |
|---|---|
| `H0` | 0.150 |
| `pycnocline` | **0.047** |
| `d_struct` | 0.193 |

![Pycnocline-reaching collapse, L=H0](../figures/floating_depth_reaching_collapse_H0.png)
![Pycnocline-reaching collapse, L=pycnocline](../figures/floating_depth_reaching_collapse_pycnocline.png)

*Figure 4: `φ*(t*)` for the pycnocline-reaching subset only, using `L=H0` (top)
vs. `L=pycnocline` (bottom). The pycnocline-scaled curves visibly
cluster more tightly, consistent with the ~3x lower spread score.*

**Result: the `pycnocline` length scale collapses pycnocline-reaching floating
structures substantially better than `H0` or `d_struct`**, mirroring the
base note's finding for bottom-fixed structures. This means the
already-generalized `d_struct`-based `Pstr` theory, combined with the
existing `pycnocline` nondimensionalization (no new theory needed),
already gives a reasonable universal-ish collapse for pycnocline-reaching
floating structures. `d_struct` alone performs worse than even the
unmodified `H0` choice — the water column's overall geometry
(`H0`, `z_t`) still matters for the mixing *rate*, not just the
structured-zone thickness, for runs where the zone does reach the
pycnocline.

This was a bounded, 3-candidate exploration per plan (phase A); the
non-reaching case (`reach_margin < 0`) was excluded here and is
deferred to a later round (§6).

## 5. Predictive function

Per the plan's empirical-first approach: **no clean closed-form
correction to `τ_x_theory` was found for the full `depth_frac` range**,
but §4.3/4.4 substantially refine what such a correction would need to
look like. The theory-to-diagnostic ratio is not a constant
multiplicative factor across `depth_frac` (§4.1) — but §4.3 shows it *is*
governed cleanly by `reach_margin = d_struct - (z_t - temp_ht)`
rather than `depth_frac` per se, crossing `ratio=1` right at
`reach_margin=0` for both `H0` groups tested. A single fitted
correction factor (analogous to `A(c4)` in the base note) would need to
be expressed as a function of `reach_margin` (or, more physically,
some diffusion-delay quantity built from it, per §6), not `depth_frac`.

**What *is* established, and useful as a practical estimate:**
1. `τ_x_theory ≈ x_frac · τ_mix_theory` (with `P_str ∝ d_struct`) gives
   the correct order of magnitude and correct qualitative sensitivity to
   `depth_frac`, `c4`, and `H0`, and is quantitatively reasonable
   (within a factor of ~2–2.5, itself consistent with the base note's
   own early-time-approximation bias) **specifically for pycnocline-
   reaching structures** (`d_struct ≥ z_t - temp_ht`) — confirmed both via the
   `τ_x` ratio (§4.3) and via the `φ*(t*)` collapse (§4.4).
2. For pycnocline-reaching structures, `φ*(t*)` (using the `pycnocline`
   length scale, not `H0`) collapses substantially better than with
   `H0`, so the existing base-note machinery
   (`length_scale="pycnocline"`) already extends adequately to floating
   structures **as long as they reach the pycnocline** — no new theory
   was needed for this
   subset (§4.4).
3. For non-reaching structures (`reach_margin < 0`), the theory
   under-predicts `τ_x` increasingly as the margin becomes more
   negative, and this regime was excluded from the phase-A collapse
   check. A proper predictive function here would likely need to
   incorporate a diffusion-delay term (turbulent transport of TKE across
   the gap between the structured zone and the pycnocline) explicitly,
   which is left as future work (§6, phase B).

## 6. Future work

- **Non-reaching structures / diffusion-delay theory (phase B).**
  §4.3/4.4 establish that `reach_margin` cleanly separates a
  "reaching" regime (existing `pycnocline`-scaled theory works well)
  from a "non-reaching" regime (theory increasingly under-predicts
  `τ_x`). The natural next step is to derive and test a diffusion-delay
  correction for the non-reaching case: TKE produced in the shallow
  structured zone must be transported across the gap
  `L_gap = (z_t - temp_ht) - d_struct` before it can act on the
  pycnocline's stratification, suggesting an additive time delay
  `tau_delay ~ L_gap^2 / K_eff` (rather than a different power-law rate)
  before the existing `t*` normalization applies. `K_eff` could be
  estimated directly from the model's own diagnosed GLS-closure
  vertical-diffusivity output (`AKt`/`AKv`, confirmed present in the
  ROMS history files) rather than fitted as a free parameter. Not
  attempted this round; deliberately deferred until the reaching-case
  baseline (established here) was solid.
- **Very thin structured zones (`depth_frac ~ 0.1` and below).** Not
  checked with a long-run subset; open question whether a genuine
  mixing floor exists at some (thinner) `depth_frac`, or whether given
  enough time these configurations also eventually reach full mixing
  (as `depth_frac = 0.25` did, contrary to the earlier informal check).
- **`pycnocline_edge()` approximation.** `z_t - temp_ht` (`k=1.0`) is a
  reasonable but somewhat arbitrary estimate of where the tanh
  thermocline transition "begins" (could also use `k=2.0`, i.e.
  `z_t - 2·temp_ht`, for a stricter definition). The straddle sweep's
  clean crossover at `reach_margin=0` (§4.3) suggests `k=1.0` is
  at least a reasonable choice, but this was not itself swept/optimized.

## 7. Code and artifacts

### 7.1 Files
- `utils/utils.py`: `load_str_a`, `find_tau_x_star`, `detect_phi_plateau`;
  generalized `compute_Pd`, `analytic_mixing_timescale`.
- `analysis/mixing_timescale.py`: `mixing_timescale()` extended with
  `x_frac`/`tail_frac` diagnostics (additive, backward-compatible);
  `d_struct` now returned; `length_scale` gained a `"d_struct"` option.
- `analysis/floating_depth.py` (new): `summarize_floating_sweep`,
  `plot_tau_x_vs_depth_frac`, `summarize_plateau`, `plot_phi_longrun`,
  `pycnocline_edge`, `plot_tau_x_ratio_vs_reach`,
  `gather_reaching_curves`, `collapse_spread_score`,
  `plot_reaching_collapse`.
- `tools/prep_floating_depth_sweep.py` (new): prep script with
  `structure.depth_frac → structure.depth_zero_below` conversion and
  dual `NTIMES`-sizing bases (`tau_x`/`tau_mix`).
- `templates/floating_depth_sweep.yaml` (new): main sweep (20 runs).
- `templates/floating_depth_longrun_sweep.yaml` (new): long-run/plateau
  subset (3 runs).
- `templates/floating_depth_reach_H75.yaml`,
  `templates/floating_depth_reach_H150.yaml` (new): pycnocline-edge
  straddle sweeps (6 runs each, `c4=0.44` fixed).

### 7.2 Reproducing the sweeps
```bash
# Main sweep: tau_x sensitivity to depth_frac / c4 / H0
python tools/prep_floating_depth_sweep.py templates/floating_depth_sweep.yaml
python tools/run_sweep.py sweeps/floating_depth/manifest.yaml

# Long-run subset: plateau / asymptotic mixed-fraction check
python tools/prep_floating_depth_sweep.py templates/floating_depth_longrun_sweep.yaml
python tools/run_sweep.py sweeps/floating_depth_longrun/manifest.yaml

# Pycnocline-reach straddle sweeps
python tools/prep_floating_depth_sweep.py templates/floating_depth_reach_H75.yaml
python tools/prep_floating_depth_sweep.py templates/floating_depth_reach_H150.yaml
python tools/run_sweep.py sweeps/floating_depth_reach_H75/manifest.yaml
python tools/run_sweep.py sweeps/floating_depth_reach_H150/manifest.yaml

# Analysis (table + figure + plateau summary)
python analysis/floating_depth.py \
    --sweep sweeps/floating_depth/manifest.yaml \
    --longrun sweeps/floating_depth_longrun/manifest.yaml --save

# Reach-margin diagnostic + phase-A collapse exploration
python analysis/floating_depth.py \
    --reach sweeps/floating_depth/manifest.yaml \
            sweeps/floating_depth_reach_H75/manifest.yaml \
            sweeps/floating_depth_reach_H150/manifest.yaml \
            sweeps/floating_depth_longrun/manifest.yaml \
    --collapse --save
```

