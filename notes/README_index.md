# Index: structure-induced mixing exploration notes

A navigation guide to all notes in this directory, grouped by theme,
newest work last within each group. See each note for full methodology,
figures, and reproduction commands.

## 1. Foundational sensitivity analysis

- **`mixing_timescale_analysis.md`** — the original sensitivity analysis:
  how the structure-mixing parametrization's dimensionless mixing timescale
  responds to `structure.c4`, closure choice, and other base parameters.
  Establishes the Carpenter et al. (2016) dimensionless-time framework used
  throughout later notes.
- **`c2_c4_closure_dependence.md`** — follow-up looking at the GLS closure
  coefficients `c2`/`c4` together, motivated by the Rennau et al. (2012)
  theoretical result that `c4 = c1` should make the structure contribution
  to the mixing coefficient Γ neutral.
- **`c2_dependence_of_structure_mixing.md`** — why `GLS.C2` alone can
  determine whether structure-induced mixing "works" at all (sign flips in
  the background term).
- **`c2_background_term_gamma_decomposition.md`** — resolves an apparent
  contradiction from the `c2` sensitivity sweep: the "danger zone" where
  mixing collapses for `c2 < c1` comes from Γ's *background* term (not the
  structure/`c4` term), i.e. it is a structure-independent effect of the
  closure's stationary-solution background balance, not a failure of the
  `c4=c1` neutrality argument itself. **This was the direct motivation for
  the present exploration session.**
- **`kmin_pmin_sensitivity.md`** — sensitivity to the GLS `Kmin`/`Pmin`
  floors (numerical safety parameters, not physical closure coefficients).

## 2. Geometry: where in the water column does a structure sit?

- **`floating_structure_depth_distribution.md`** — establishes the
  realistic spar (Hywind Tampen, draft≈90 m) vs. semisubmersible (VolturnUS-S,
  draft≈20 m) geometry anchors used in every later geometry comparison, and
  the project-wide `str_a` rescaling convention (×312.5) needed to keep
  idealised-run mixing timescales tractable.
- **`draft_depth_mixing_sweep.md`** — a continuous draft sweep (10–150 m,
  fixed realistic `str_a`, sheared current): no "Goldilocks" non-monotonic
  draft exists; the mixing timescale decreases monotonically with draft but
  saturates sharply beyond the pycnocline depth — draft differences matter
  most exactly in the upper ~40 m where semisubmersibles sit.

## 3. Background current: strength, shear, and exposure time

- **`current_strength_sensitivity.md`** — how strongly background current
  speed controls mixing (steep, roughly cubic `Pd ∝ u³`-consistent
  sensitivity); also documents a real bug caught by diagnosing the model's
  actual `u` field rather than trusting nominal sweep labels.
- **`residence_time_tradeoff.md`** — does a finite farm length reverse the
  current-speed trend, since faster currents give a water column less time
  inside the farm? No: faster currents still deliver *more* mixing by farm
  exit despite the shorter residence time (`T_res/t_mix ∝ u²` grows with
  speed).
- **`shear_geometry_interaction.md`** — does a realistic, surface-intensified
  (NCC-like) sheared current change the spar-vs-semisub comparison, since
  semisubmersibles sit where real currents are fastest? Yes: shear closes
  most (not all) of the gap — semisub `t(phi*=0.5)` drops 8.73→2.55 days,
  spar barely changes.

## 4. Observational grounding

- **`gamma_observational_check.md`** — a first sanity check of the
  parametrization's free parameters (`str_Cd`, `str_c4`) against an
  independent observational/LES benchmark (Schultze et al. 2020, 8–14%
  bulk mixing efficiency for monopile wakes): the same `c4=0.44` diagnoses
  a bulk Γ above this range for the spar and below it for the semisub — one
  `c4` cannot simultaneously match both geometries to this benchmark.

## 5. Visualization

- **`mixing_mechanism_visualizations.md`** — animations (`mixing_race.gif`,
  `mixing_mechanism.gif`) and static Hovmöller diagrams
  (`dissipation_hovmoller.png`, `gamma_hovmoller.png`) built from the
  shear-geometry runs, showing density erosion, dissipation rate,
  buoyancy destruction of TKE, and the local mixing coefficient Γ(z,t)
  evolving in depth and time. Finds that (a) dissipation declines smoothly
  with depth rather than stepping at the structure draft — background
  shear production dominates the upper water column almost everywhere —
  and (b) diagnosed Γ is persistently well above the Schultze et al.
  observational band along the active mixing front, reinforcing the
  bulk-Γ overshoot found in Section 4 as a spatially/temporally
  widespread feature, not a quirk of how the bulk average was computed.

## 6. Residence time revisited under shear

- **`shear_residence_time_tradeoff.md`** — repeats the residence-time
  tradeoff (Section 3) using the sheared current instead of a uniform one:
  the original conclusion survives (`phi*` at farm exit still decreases
  monotonically with current speed for both geometries), but the
  semisubmersible's farm-exit response is noticeably flatter/less
  speed-sensitive than the spar's once shear is included.

---

**Overall practical takeaways for the regional, realistic model** (repeated
across several notes, collected here for convenience):
- Mixing impact is dominated by the current's faster episodes/regions, at
  both the single-column and farm-footprint (residence-time) scale.
- Draft/geometry differences between floater types matter most in the
  upper ~40 m — exactly where semisubmersibles, not spars, sit — and this
  sensitivity is amplified, not reduced, by realistic current shear.
- The `c4` parameter cannot be chosen once to match an observational
  mixing-efficiency benchmark for both spar- and semisub-like geometries
  simultaneously with the current rescaling convention; this is a known
  open calibration gap, not a validated result.
