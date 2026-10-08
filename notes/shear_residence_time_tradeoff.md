# Residence time vs. mixing timescale under a sheared current: does the earlier conclusion survive?

**Context.** `notes/residence_time_tradeoff.md` found that, for a
depth-uniform current, faster currents still deliver *more* mixing by farm
exit despite a shorter residence time inside a finite-length farm
(`phi*(T_res)` decreases monotonically with current speed, for every farm
length tested). That used the *uniform*-current `current_strength_sensitivity`
sweep. Since `notes/shear_geometry_interaction.md` showed that a realistic,
surface-intensified sheared current materially changes the spar-vs-semisub
mixing comparison, it's natural to ask: **does the residence-time tradeoff
conclusion still hold once the current is realistically sheared, or could
shear reverse it?**

## Method

New sweep: `tools/prep_shear_current_strength_sweep.py`, 8 new runs (+2
reused from `shear_geometry_interaction.md`). Same matched-`total_a=1.5`
spar/semisub geometries as every other geometry comparison in this project.
Same *shape* of sheared current profile as
`profiles/ncc_shear_surface_intensified.txt`, uniformly rescaled to 5
depth-mean targets (0.05, 0.10, 0.15, 0.20, 0.30 m/s —
`profiles/ncc_shear_scale_u{0p05,0p10,0p20,0p30}.txt`, plus the existing
0.15 m/s point from `runs/shear_spar`/`runs/shear_semisub`), matching the
uniform-current sweep's targets exactly for a fair comparison.

As before, "diagnosed" current speed (not the nominal target) is used for
`T_res = L/u`: the late-time, structure-layer-averaged `|u|`, following
`analysis/current_strength_sensitivity.py:diagnosed_speed` — for a sheared
profile this is the speed averaged over *only* the depths the structure
actually occupies, so a shallow semisubmersible's diagnosed speed samples
the fast near-surface part of the profile, while the deep spar's averages
over most of the profile including its decayed, slower deep-water part.

## Result: the tradeoff conclusion survives shear

![Stratification remaining at farm exit vs. current speed, sheared current, both geometries](../figures/shear_residence_time_tradeoff.png)

For every farm length tested (2, 5, 10, 20 km) and both geometries,
`phi*(T_res)` still **decreases monotonically with current speed** under
the sheared current — the same qualitative conclusion as the uniform-current
case. Faster currents (even though shear means "faster" now also implies a
disproportionately faster *near-surface* current) still mix a transiting
water column more by the time it exits the farm, not less.

A few shear-specific differences worth noting:

- **Semisub is less farm-exit-sensitive to speed than spar, under shear.**
  At the highest tested speed, the spar's `phi*(T_res)` for a 20 km farm has
  already collapsed to ≈0 (fully mixed well before/at exit) for diagnosed
  `u` ≳ 0.24 m/s, while the semisub only reaches `phi*≈0.39` at its highest
  diagnosed speed (0.42 m/s) for the same farm length — i.e. the curves are
  noticeably flatter for the semisub than for the spar. This is consistent
  with `notes/shear_geometry_interaction.md` and
  `notes/draft_depth_mixing_sweep.md`: once a structure's draft reaches the
  pycnocline, additional current speed keeps buying disproportionately more
  mixing (saturating-but-still-steep response), while a shallow structure's
  response is comparatively damped.
- **Diagnosed speeds are noticeably non-linear in the nominal scale factor**
  for the semisub: the `u_target=0.30` point diagnoses at **0.42 m/s**
  (40% higher than target), substantially more than the spar's corresponding
  overshoot (0.33 vs. 0.30 target). This is because the semisub's
  structure-layer average samples only the shallow, most strongly
  surface-intensified part of the profile, which is scaled by the same
  linear factor as the whole profile but starts from a higher baseline —
  not a sign of any instability, but a reminder that "current strength" is
  not a single number once the profile is sheared and geometry-dependent.

## Caveats

- Same caveats as the uniform-current version apply here (see
  `notes/residence_time_tradeoff.md`): `T_res = L/u` ignores farm-drag
  feedback on the along-flow current itself, assumes a single steady
  transit speed, and only addresses one along-flow pass through an
  idealised, horizontally homogeneous 1D column.
- Only one shear *shape* was tested (uniformly rescaled in amplitude, same
  surface-intensification length scale throughout) — a current that becomes
  more surface-trapped (not just stronger) at higher speeds was not tested
  and could behave differently.

## Practical implication

The regional-model takeaway from `notes/residence_time_tradeoff.md` stands
even once realistic current shear is included: **don't expect
slower-current regions/episodes of the NCC to see comparable net
mixing-per-farm-transit to faster ones** — if anything, the shear case
sharpens the geometry-dependence of this effect, with semisubmersible-type
floaters showing a flatter (less speed-sensitive) farm-exit response than
spar-type floaters.

## Reproduction

```bash
mamba activate roms   # velocity-shear-bodyforce branch, roms + mixtest_1d
python tools/prep_shear_current_strength_sweep.py
python tools/run_sweep.py sweeps/shear_current_strength_sweep/manifest.yaml
python analysis/shear_residence_time_tradeoff.py   # -> figures/shear_residence_time_tradeoff.png
```
