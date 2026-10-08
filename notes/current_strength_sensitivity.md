# How strongly does the background current speed control structure-induced mixing?

**Context.** Every prior floating-structure experiment in this repo
(`notes/floating_structure_depth_distribution.md`) used a single background
current speed (`u ≈ 0.1 m/s`, chosen only to make runs tractable, not
because it is a particularly representative value). But the Norwegian
Coastal Current (NCC) -- the current regime most relevant to the eventual
regional-model application along the Norwegian coast -- has a wide natural
range, from a few cm/s to several tens of cm/s depending on location,
season and wind forcing. This note asks the obvious next question: how
sensitive is the structure-induced mixing timescale to the background
current speed, for the matched spar-like and semisubmersible-like
geometries already used in the depth-distribution note?

This is an exploratory extension, run autonomously; the experimental
choices below are defensible but not yet reviewed by a human.

## 1. A methodological pitfall worth documenting

The first version of this sweep varied `bodyforce.BFRC_U` (the config key
that also feeds `utils.analytic_mixing_timescale()`'s pre-run tau_mix
estimate) under the assumption that it was the actual simulated forcing
knob, scaled via the quasi-steady quadratic drag balance
`u_inf = sqrt(2*BFRC_U/(CD*str_a))`. It is not: `tools/make_grd.py`'s
`build_bfrc()` shows that for `bodyforce.mode: "uniform"` (used throughout
this repo), the simulated body-force acceleration is

```
bfrc_u = BFRC_CD * U0
```

i.e. it depends on `bodyforce.U0` (the literal target velocity), not on
`BFRC_U` at all. `BFRC_U` is only consumed by `analytic_mixing_timescale()`,
and only for full-depth structures -- for the partial-depth "floating"
structures studied in this repo, it is inert bookkeeping copied from an
earlier full-depth sweep and has no effect on the simulation. The first
8-run attempt at this sweep consequently produced 8 runs with *identical*
forcing (`U0=0.1 m/s` throughout, the baseline default) despite 8 different
nominal "target speeds" (`BFRC_U` values spanning a 36x range).

This was only caught because the diagnosed/measured background speed
(averaged over the structure's depth range and the last 20% of the run,
directly from the model's `u` field) was checked against the nominal label
rather than trusted -- it came back numerically identical across all 8 runs.
**The practical lesson: when driving this body-force setup, always vary
`bodyforce.U0` directly, and always diagnose the actual resulting speed
from model output rather than from input-parameter arithmetic -- `BFRC_U`
is not a safe proxy for anything except the one full-depth analytic
formula it was designed for.** The sweep was corrected
(`tools/prep_current_strength_sweep.py` now overrides `bodyforce.U0`
directly) and re-run; all results below are from the corrected runs.

## 2. Experiment design

Same matched-`total_a=1.5` geometries as the depth-distribution note:

- **spar-like**: draft = 90 m, `str_a = 1.5/90 ≈ 0.0167 m⁻¹`
- **semisub-like**: draft = 20 m, `str_a = 1.5/20 = 0.075 m⁻¹`

Five target current speeds per geometry: `U0 ∈ {0.05, 0.10, 0.15, 0.20,
0.30} m/s` (the `0.10` point reuses the existing
`geom_matched_draft{90,20}` runs). All other settings (stratification,
`CD=0.63`, fixed `NTIMES=99540` i.e. 46.08 days, `BFRC_CD=0.005`) are
unchanged from the depth-distribution note, so results are directly
comparable to it.

Because the quadratic structure drag is itself a function of the flow
speed it also damps, the *achieved* quasi-steady speed is somewhat lower
than the nominal target `U0` (the structure's own drag is an additional
sink the linear `BFRC_CD` damping has to overcome). The diagnosed speed is
what is actually plotted against the mixing timescale below.

| geometry | nominal `U0` (m/s) | diagnosed speed (m/s) |
|---|---|---|
| spar   | 0.05 | 0.048 |
| spar   | 0.10 | 0.092 |
| spar   | 0.15 | 0.133 |
| spar   | 0.20 | 0.171 |
| spar   | 0.30 | 0.244 |
| semisub | 0.05 | 0.043 |
| semisub | 0.10 | 0.077 |
| semisub | 0.15 | 0.108 |
| semisub | 0.20 | 0.135 |
| semisub | 0.30 | 0.185 |

(Diagnosed speed is the mean of `|u|` over the structure's depth range and
the last 20% of the 46-day run; see `analysis/current_strength_sensitivity.py`.)

## 3. Results

![Stratification decay and mixing timescale vs. current speed](../figures/current_strength_sensitivity.png)

**Left panel:** `phi*(t) = phi(t)/phi(0)` for every run. Faster currents
mix the column dramatically faster, and -- as in the depth-distribution
note -- the spar-like (full-column-reaching, deep-draft) geometry mixes
faster than the semisub-like (shallow-draft) geometry at every matched
current speed, for the same total frontal area.

**Right panel:** `t(phi*=0.5)` vs. diagnosed current speed, for both
geometries (log-scale would show it even more clearly, but the number of
points is small). The mixing timescale falls extremely steeply with
current speed -- roughly consistent with the `t_mix ∝ u⁻³` scaling implied
by Carpenter et al.'s `P_d ∝ CD·str_a·u³`: going from the slowest to the
fastest spar run (0.048 → 0.244 m/s, a 5x increase in speed) cuts the
mixing time from 16.1 days to 0.13 days (~120x), and for the semisub run
that reached `phi*=0.5` at both ends (0.077 → 0.185 m/s, a 2.4x increase)
the mixing time drops by ~16x (24.97 → 1.54 days) -- both broadly
consistent with a cubic-ish scaling, though the comparison is too noisy
(5 points, non-matching speed ratios between geometries) to claim a
precise exponent.

The weakest semisub run (`u ≈ 0.043 m/s`) does not complete mixing within
the fixed 46-day window at all (`phi*(46 d) = 0.776`) -- itself a
meaningful result, not a sizing artefact: a slow-moving semisubmersible
farm in a weak part of the NCC may simply not mix the local water column
within a residence-time-relevant period.

## 4. Practical implication for the Norwegian Coastal Current

The NCC's climatological speed is commonly in the few-cm/s to ~0.3 m/s
range depending on location/season/wind forcing, i.e. it spans essentially
the whole range tested here. Given how steeply `t_mix` falls with current
speed, this means structure-induced mixing from a floating wind farm is
likely to be a strongly intermittent, event-like process in a real
(unsteady) current regime -- dominated by the current's faster episodes
(e.g. wind-driven pulses or the along-coast jet core) rather than behaving
like a slow, steady background diffusivity increase. A regional model
with a realistic, time-varying NCC should therefore be expected to show
patchy/pulsed mixing signatures rather than a spatially smooth, constant
enhancement -- consistent with how nonlinearly (and quickly) `t_mix`
responds to `u` in this idealised test.

## 5. Caveats

- Only 5 discrete speeds per geometry, and the two geometries don't share
  exactly the same diagnosed-speed grid (since the diagnosed speed is an
  emergent, not directly imposed, quantity) -- not enough to fit a reliable
  power-law exponent.
- No rotation (`F0=0`) and no tidal/oscillatory component to the current;
  real NCC flow is unsteady on multiple timescales this idealised setup
  cannot represent (see `notes/current_strength_sensitivity.md`'s
  deprioritised sibling idea, oscillatory forcing, which would need a
  Fortran/recompile change not attempted here).
- `BFRC_CD=0.005` (fixed throughout) sets how tightly the flow is nudged
  toward `U0`; this was not itself varied, so the diagnosed-speed-vs-`U0`
  relationship documented in Section 2 is specific to this damping rate.

## Reproducing this note

```
python tools/prep_current_strength_sweep.py
python tools/run_sweep.py sweeps/current_strength_sensitivity/manifest.yaml
python analysis/current_strength_sensitivity.py
```
