# Visualizing the mixing mechanism: animations and Hovmöller diagrams

A set of visualization-focused companions to
`notes/shear_geometry_interaction.md`, built purely from already-completed
simulation output (no new runs). All three use the same pair of runs:

| run | geometry | draft | `str_a` | `total_a` (= str_a × draft) |
|---|---|---:|---:|---:|
| `runs/shear_spar`    | spar (Hywind Tampen-like)      | 90 m | 0.01667 | 1.50 |
| `runs/shear_semisub` | semisub (VolturnUS-S-like)     | 20 m | 0.075   | 1.50 |

Both use the **same, matched `total_a = 1.5`** (same total depth-integrated
drag budget, just redistributed over a shallower or deeper layer — this is
the "matched-budget" convention, *not* the fixed-`str_a` convention used in
the separate draft sweep, see `notes/draft_depth_mixing_sweep.md`), the same
closure (`gen`, `str_c4 = 0.44`, `str_Cd = 0.63`), and the same **sheared**
current forcing: `bodyforce.mode = "profile"`, reading
`profiles/ncc_shear_surface_intensified.txt` (the surface-intensified,
NCC-like shear profile with `u_floor=0.08`, `u_surf=0.35 m/s`, `scale=40 m`,
depth-mean 0.15 m/s — see `notes/shear_geometry_interaction.md` for how this
profile was built). This is **not** the uniform-current baseline used
elsewhere in the project — shear is the whole point of these runs, since it
is what lets the near-surface semisubmersible access energetic surface
currents that a depth-uniform current would not distinguish from the spar's.

## 1. `figures/mixing_race.gif` — the original "mixing race"

`analysis/mixing_race_animation.py`. Two panels, animated over the first 6
days (144 frames at 5 fps):
- **left**: the density profile eroding over time, spar (blue) vs. semisub
  (orange), with each structure's draft shaded.
- **right**: `phi*(t)` (fraction of initial stratification potential energy
  remaining) for both, with a moving marker — literally a race to `phi*→0`.

Confirms visually what `notes/shear_geometry_interaction.md` found
quantitatively: under shear, the semisub's `t(phi*=0.5)` (2.55 days) is much
closer to the spar's (0.62 days) than under a uniform current (8.73 vs. 0.76
days) — still slower, but no longer in a different league.

## 2. `figures/dissipation_hovmoller.png` — where/when does mixing happen?

`analysis/dissipation_hovmoller.py`. A static 2×2 depth-time heatmap: rows
are spar/semisub, columns are density anomaly and dissipation rate
`epsilon(z,t)` (derived via `utils.utils.compute_epsilon`, log colour scale).
The structure's draft is marked with a dashed line.

**Finding**: `epsilon` does **not** show a sharp step at the structure
draft — it declines smoothly with depth in both runs, with the e-folding
scale set by the shear profile's surface intensification, not by the
structure's footprint. This means, in this parametrization/setup, the
near-surface *background* shear production dominates the dissipation
budget over the *structure-added* production almost everywhere except deep
within the stratified interior — a useful sanity check on how "localized"
the structure's signature actually is in the dissipation field (it is more
visible in the *density* panel, where the erosion front clearly stalls at
the structure's draft, than in the `epsilon` panel itself).

## 3. `figures/mixing_mechanism.gif` — Γ = -B/ε in motion

`analysis/mixing_mechanism_animation.py`. A 2×2 animated version extending
the race: density profile, dissipation rate `epsilon(z)`, buoyancy
destruction of TKE `B(z) = (g/ρ₀) A_Kt ∂ρ/∂z` (the numerator of Γ), and the
`phi*(t)` race track — all four evolving together, frame by frame, for
spar vs. semisub. This is the most direct visual of the Rennau et al.
(2012) theory from `notes/c2_c4_closure_dependence.md`: you can watch `B`
(strongly negative where stratification is actively suppressing turbulence)
and `epsilon` evolve together, and see where/when their ratio (Γ) would be
largest.

Note: the very first output timestep (t ≈ 0.04 days) contains a transient
adjustment spike in `B` about 2 orders of magnitude larger than anything
afterward (initial GLS spin-up artifact); axis limits are set from
percentiles excluding that frame so the meaningful later-time structure
remains visible.

## 4. `figures/gamma_hovmoller.png` — where/when is the local mixing coefficient largest?

`analysis/gamma_hovmoller.py`. A static depth-time heatmap of the local
diagnosed mixing coefficient `Gamma(z,t) = -B(z,t)/epsilon(z,t)` itself
(not just its two ingredients separately, as in panel 2 above, nor the
single bulk/epsilon-weighted number from
`notes/gamma_observational_check.md`). Cells where `epsilon` is below the
20th percentile are masked grey (the ratio is numerically meaningless
there — both numerator and denominator are near zero). The colour scale
is annotated with the Schultze et al. (2020) 8-14% observational band for
direct visual comparison.

**Finding**: diagnosed Γ is **well above** the observational band (often
saturating the ±0.3 colour limit) along the active mixing front for both
geometries, consistent with the bulk-Γ numbers in
`notes/gamma_observational_check.md` — this reinforces that visually, not
just as a single summary statistic: the overshoot is a persistent feature
across the whole active mixing region and most of the simulated period,
not an artifact of how the bulk average was computed.

## Reproduction

```bash
mamba activate roms   # velocity-shear-bodyforce branch, roms + mixtest_1d
# runs/shear_spar and runs/shear_semisub already exist (from
# tools/prep_shear_geometry_sweep.py) -- no new simulations needed.
python analysis/mixing_race_animation.py        # -> figures/mixing_race.gif
python analysis/dissipation_hovmoller.py        # -> figures/dissipation_hovmoller.png
python analysis/mixing_mechanism_animation.py   # -> figures/mixing_mechanism.gif
python analysis/gamma_hovmoller.py              # -> figures/gamma_hovmoller.png
```
