# Does it matter *where* in the water column a floating wind farm sits? Draft, geometry, and season

**Context.** All prior notes in this repo used a structure that occupies the
*full* water column (`depth_zero_below >= H0`) -- a reasonable idealisation
for a bottom-fixed monopile farm. But the long-term motivation for this
project is floating wind (spar and semisubmersible substructures) along the
Norwegian coast, and floating substructures only occupy the *upper* part of
the water column. This note asks: does the vertical placement of the
structure's drag (its draft) change how effectively it breaks down
stratification, independent of how much total drag it has? And does using
realistic substructure geometry (rather than an idealised, "matched" drag
budget) change the comparison between a spar and a semisubmersible?

This is an exploratory extension, run autonomously overnight; the
experimental choices below are defensible but not yet reviewed by a human.

## 1. Real-world geometry basis

Two Norwegian/public-domain floating substructure references were used as
concrete anchors (rather than inventing arbitrary numbers):

- **Spar**: Equinor's **Hywind Tampen** (operating, Norwegian North Sea,
  8.6 MW turbines): draft ≈ 90 m, hull diameter ≈ 14.7 m.
- **Semisubmersible**: the **UMaine VolturnUS-S** reference design for the
  **IEA 15 MW** reference turbine (NREL/TP-5000-76773): draft 20 m, one
  central column (10 m diameter) and three offset columns (12.5 m diameter)
  at 51.75 m radius from the tower axis, connected by pontoons.

**Caveat (stated once, applies throughout):** the two real structures are
sized for different turbine ratings (8.6 MW vs 15 MW), and the
semisubmersible's three-legged plan geometry is approximated here as a
single "effective" frontal width (~2 column diameters ≈ 25 m, a rough
average over current heading) and the spar as a uniform cylinder -- both are
simplifications of the real (tapered, multi-column) structures. The model's
`str_a(z)` is also a strict step function (uniform drag density down to
`depth_zero_below`, zero below); it cannot represent the taper of a real
spar or the multi-column plan form of a semisubmersible. All results below
should be read as "what the step-function approximation of each type
predicts," not a substitute for a geometry-resolving simulation.

## 2. Experiment design

Everything here reuses the inflated (but internally consistent) drag/forcing
scaling from `configs/variants/mixing_timescale.yaml` (`u_inf ≈ 0.1 m/s`
background current, `str_a` and `BFRC_U` scaled together so stratification
actually breaks down within a tractable run length). The default
stratification (`temp_zt=40 m`, `temp_dT=10` unless noted) and `CD=0.63`
are unchanged from the baseline config.

Because `utils.utils.analytic_mixing_timescale()` (used in earlier sweeps to
auto-size `NTIMES`) requires a full-depth structure, it cannot be used here.
Instead, every run below uses the same fixed `NTIMES=99540` (`DT=40 s`,
46.08 days of simulated time) -- 5x the analytic mixing time of the
full-depth reference case -- so that all runs are directly comparable in
real time, and "did not finish mixing within the window" is itself a
meaningful (negative) result rather than an artefact of mis-sized `NTIMES`.

Diagnostic: `phi*(t) = phi(t)/phi(0)`, where `phi` is the stratification
potential-energy anomaly (`utils.utils.compute_phi`, Carpenter et al. 2016
Eq. 7). `phi* -> 0` means fully mixed; `phi*` close to 1 means little
change. (`compute_Pd`/`compute_Pstr`/`analytic_mixing_timescale`, which
assume a depth-uniform `str_a`, are **not** used here since they silently
mis-handle partial-depth structures; `phi*` depends only on the density
field and is structure-agnostic.)

Three sub-experiments:

**(A) Draft sweep at fixed total frontal area.** Define
`total_a = str_a * depth_zero_below` (the depth-integrated drag density).
Holding `total_a = 1.5` fixed (equal to the existing validated full-depth
case, `str_a=0.01, H0=150 m`) and varying the draft
`d ∈ {15,20,30,40,50,65,90,120,150} m` (with `str_a = total_a/d` adjusted to
compensate) isolates the effect of *where* a fixed drag budget is placed in
the water column, independent of how much drag there is in total. The
pycnocline is fixed at 40 m depth.

**(B) Realistic (non-matched) total frontal area.** Using the real
geometry above at a common hypothetical 500x500 m turbine spacing (the
worked example already used in the repo's README):

| type | effective diameter | draft | `a = d/spacing^2` | `total_a = a x draft` |
|---|---|---|---|---|
| monopile (bottom-fixed reference) | 8 m | 150 m (full) | 3.2e-5 | 4.8e-3 |
| spar (Hywind Tampen) | 14.7 m | 90 m | 5.88e-5 | 5.29e-3 |
| semisubmersible (VolturnUS-S) | ~25 m (2-column effective) | 20 m | 1.0e-4 | 2.0e-3 |

Rescaled by `(total_a_type / total_a_monopile) x 1.5` to keep the same
absolute magnitude as experiment (A) (for comparable run lengths) while
preserving the realistic *ratio* between types: `total_a_spar=1.65`,
`total_a_semisub=0.63` (vs. the monopile/matched value of 1.5).

**(C) Seasonality.** The same spar (90 m) and semisub (20 m) geometries
from (A), each run against two additional idealised Norwegian-coastal
stratification regimes on top of the default: a "summer" case (shallow,
strong pycnocline: `temp_zt=15 m, temp_dT=8`) and a "winter" case (deep,
weak pycnocline: `temp_zt=60 m, temp_dT=3`).

All 15 new runs completed successfully (`returncode=0`), run sequentially.

## 3. Result A: draft sweep (fixed total frontal area)

![Draft-depth distribution of mixing](../figures/draft_depth_distribution_mixing.png)

| draft (m) | `str_a` | t(φ*=0.5) (days) | t(φ*=0.1) (days) | fully mixed by day 46? |
|---|---|---|---|---|
| 15  | 0.100   | 39.5 | never | no (φ*=0.46) |
| 20  | 0.075   | 25.0 | never | no (φ*=0.30) |
| 30  | 0.050   | 11.3 | 31.3  | yes |
| 40  | 0.0375  | 6.9  | 18.0  | yes |
| 50  | 0.030   | 5.2  | 13.3  | yes |
| 65  | 0.0231  | 3.1  | 8.2   | yes |
| 90  | 0.0167  | 2.3  | 5.4   | yes |
| 120 | 0.0125  | 2.2  | 3.9   | yes (fastest) |
| 150 (full) | 0.010 | 2.5 | 4.6 | yes |

Two findings stand out:

1. **A shallow-draft structure (less than the pycnocline depth) is
   dramatically less effective, even carrying the exact same total drag
   budget.** At draft=15-20 m (shallower than the 40 m pycnocline), the
   water column has not even reached 50% mixed by the end of the 46-day
   window; at draft≥30 m it mixes completely, and the time to do so falls
   off steeply as the draft approaches and then exceeds the pycnocline
   depth. This is consistent with the structure needing to reach down to
   (or past) the stratified interface to disrupt it directly, rather than
   only stirring the surface mixed layer above it.
2. **Full-depth placement is not actually optimal.** Mixing keeps getting
   faster as draft increases up to ~120 m, then turns around slightly at
   150 m (full depth). The effect is small (2.16 vs 2.52 days to φ*=0.5)
   but consistent across both the 50%- and 10%-mixed thresholds, so it is
   unlikely to be noise. A plausible mechanism: leaving a thin
   structure-free gap just above the seabed lets a return flow develop
   below the drag layer, adding an extra shear interface (between the
   decelerated, structure-affected water above and the freely-accelerating
   water in the gap) on top of the one at the pycnocline -- two shear
   sources instead of one. This is a second-order effect, not a dominant
   one, but worth keeping in mind if someone is tempted to assume "more
   draft is always better."

## 4. Result B: realistic vs. matched total frontal area

![Realistic total frontal area](../figures/realistic_total_area_mixing.png)

| case | draft | `total_a` | fully mixed by day 46? |
|---|---|---|---|
| spar, matched | 90 m | 1.50 | yes |
| spar, realistic | 90 m | 1.65 | yes (marginally faster) |
| semisub, matched | 20 m | 1.50 | no (φ*=0.30) |
| semisub, realistic | 20 m | 0.63 | no, and *less* mixed (φ*=0.49) |

When the comparison is forced to use real geometry-derived frontal areas
instead of an artificially equalised drag budget, the gap between spar and
semisubmersible *widens*, it does not narrow. This is because the two
disadvantages compound in the same direction for the semisubmersible: it is
already handicapped by its shallow draft (Result A), and on top of that its
real frontal area (relative to a spar, at the same farm footprint) is
roughly 3x smaller in this geometry approximation. A spar farm is
essentially unaffected by the correction (it was already fast, and if
anything the realistic spar has marginally *more* frontal area than the
"matched" baseline, so it mixes a little faster still).

## 5. Result C: seasonal sensitivity

![Seasonal sensitivity](../figures/seasonal_sensitivity_mixing.png)

| case | t(φ*=0.5) (days) | fully mixed by day 46? |
|---|---|---|
| spar, summer (zt=15, dT=8)   | 0.78 | yes |
| spar, default (zt=40, dT=10) | 2.28 | yes |
| spar, winter (zt=60, dT=3)   | 0.98 | yes |
| semisub, summer (zt=15, dT=8) | 6.35  | yes |
| semisub, default (zt=40, dT=10) | 25.0 | **no** (φ*=0.30) |
| semisub, winter (zt=60, dT=3)   | 12.5 | yes |

The spar is effectively insensitive to season: it fully mixes the column
within a few days regardless of stratification regime, because its draft
(90 m) exceeds all three pycnocline depths tested, so it always acts
directly on the stratified interface.

The semisubmersible (20 m draft) is, by contrast, strongly season-dependent,
and the result is not simply "stronger stratification = slower mixing."
The semisub fully mixes in *both* the summer case (pycnocline shallower
than its draft -- direct disruption, despite strong stratification) and the
winter case (much weaker stratification, despite the pycnocline being far
below its draft -- less potential energy to overcome). It is the "default"
regime specifically -- a combination of *relatively strong* stratification
*and* a pycnocline well below the structure's reach -- that is hardest for
a shallow structure to break down. In other words, a semisubmersible farm's
regional mixing footprint should be expected to depend on the local
seasonal cycle in a non-trivial, two-sided way: it is weakest in exactly
the conditions where the pycnocline sits below the structure's draft *and*
stratification is strong -- plausibly a late-summer/autumn situation in
parts of the Norwegian coastal current.

## 6. Practical implications for the regional model

- **Substructure type is a first-order control on the mixing footprint
  predicted by this parametrization, and the control is not well summarised
  by total installed drag alone.** A spar farm (deep draft) is predicted to
  mix the water column almost completely within days in this idealised
  setting, essentially independent of season. A semisubmersible farm
  (shallow draft) is predicted to have a much weaker, strongly
  season/stratification-dependent effect, and can fail to fully mix the
  column at all within a realistic simulation window when the pycnocline
  sits below its draft.
- **If the regional model can only afford to resolve one substructure type
  as a worst case / upper bound, use a spar-type (or other deep-draft)
  geometry** -- it gives a more pronounced and more robust (less
  season-sensitive) signal. A semisubmersible-only farm may be a poor
  candidate for *detecting* structure-induced mixing impacts in a regional
  simulation unless the local pycnocline is unusually shallow.
- **Draft relative to the local pycnocline depth matters more than draft in
  absolute terms.** Before committing resolution/compute budget to a
  region, it would help to know the typical seasonal pycnocline depth range
  along the relevant stretch of coast and compare it to the planned
  farm's substructure draft(s) -- this analysis suggests the structure
  needs to reach to or past the pycnocline to meaningfully affect it.
- These are all step-function-drag, single-column, no-rotation,
  one-realization idealisations. They are best read as hypotheses to check
  against the regional run (e.g.: "does a spar farm show a bigger and more
  season-independent SST/mixed-layer-depth signature than a semisub farm of
  comparable capacity, in the regional model?") rather than as quantitative
  predictions.

## 7. Reproducibility

- Config: `configs/variants/floating_structure_mixing.yaml` (fixed
  `NTIMES=99540`, built on `configs/variants/mixing_timescale.yaml`).
- Prep script (bespoke, not a cartesian `prep_sweep.py` sweep, since draft
  and `str_a` are paired, not independent, axes):
  `python tools/prep_structure_geometry_sweep.py`
  (writes `sweeps/floating_structure_geometry/manifest.yaml`, 15 runs).
- Run: `python tools/run_sweep.py sweeps/floating_structure_geometry/manifest.yaml`
  (sequential; do not parallelise, see repo-wide note on `romsS` races).
- Figures: `python analysis/floating_structure_mixing.py` (reads the 15 run
  directories under `runs/` by name and regenerates all three figures in
  `figures/`). Does **not** reuse `analysis/mixing_timescale.py`'s
  `mixing_timescale()`/`t_star` machinery, since that depends on
  `utils.compute_Pstr`/`utils.analytic_mixing_timescale`, both of which
  assume a full-depth, depth-uniform `str_a` -- false for every run here.
  Instead it computes `phi_star(t)` directly and compares runs in real
  (dimensional) time.
- Repo state: `mixtest_1d` and `roms` both on branch
  `velocity-shear-bodyforce` (`roms` at commit `75ac6191`).
- Run outputs live under `runs/<name>/` (gitignored); this note, the
  analysis script, and the figures are the persistent record.
