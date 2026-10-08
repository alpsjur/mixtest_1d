# Does a surface-intensified sheared current change the spar-vs-semisubmersible mixing comparison?

**Motivation.** `notes/floating_structure_depth_distribution.md` compared
a spar-like structure (drag spread over a 90 m draft) against a
semisubmersible-like structure (drag concentrated in a shallow 20 m
draft), matched to have the same depth-integrated drag
(`structure.str_a * draft = total_a = 1.5`). Under a *uniform* background
current, that comparison found the spar mixes the water column much
faster than the semisubmersible, because its drag acts throughout the
full stratified column rather than only near the surface.

All of that used a spatially uniform current (`bodyforce.mode:
"uniform"`). But the project owner pointed out that real currents --
including the Norwegian Coastal Current -- are commonly surface-
intensified (strongest near the surface, decaying with depth, driven by
wind/buoyancy/Ekman dynamics), and a semisubmersible's shallow draft sits
right in that near-surface high-velocity layer, while a spar's draft
reaches down into comparatively slack water. Since structure-drag
production scales roughly as `str_a * u^3`, this raised the question:
**does the semisubmersible catch up on the spar once the current
driving it is realistically sheared, rather than depth-uniform?**

## Method

The idealised setup already supports `bodyforce.mode: "profile"` (an
arbitrary depth-varying target current, read from a text file, see
`utils.resolve_velocity_profile`), used here for the first time.

A representative surface-intensified profile was built
(`profiles/ncc_shear_surface_intensified.txt`):

```
u(z) = u_floor + (u_surf - u_floor) * exp(-depth/scale)
u_floor = 0.08 m/s, u_surf = 0.35 m/s, scale = 40 m
```

(an e-folding depth of 40 m is consistent with the Ekman/plume-type
decay scale reported for the NCC's surface-intensified shear). The
floor/surface values were tuned so that this profile's depth-thickness-
weighted column mean, on this model's actual vertical grid, is 0.150
m/s -- matching the already-completed **uniform**-current runs
`current_spar_u0p15` / `current_semisub_u0p15` (`bodyforce.U0=0.15
m/s`) exactly. This makes the comparison "fair": both forcings put the
same total momentum into the water column, only its *vertical
distribution* differs.

Two new runs were prepared and executed (`tools/prep_shear_geometry_sweep.py`):
`shear_spar` and `shear_semisub`, identical to the matched-total_a
geometries used throughout this project, but with `bodyforce.mode:
"profile"` pointing at the shear profile above. The `ANA_SMFLUX`
surface-forcing infrastructure built for the wind-stress experiment was
reused/verified as unaffected (`bodyforce.mode` is independent of
`forcing`/wind stress).

## Result: the semisubmersible catches up dramatically; the spar barely changes

![Stratification decay and time-to-half-mix, uniform vs. sheared current, for both geometries](../figures/shear_geometry_sweep.png)

| geometry | forcing | t(phi*=0.5) |
|---|---|---|
| spar | uniform (U0=0.15 m/s) | 0.76 d |
| spar | sheared (same depth-mean) | 0.62 d |
| semisub | uniform (U0=0.15 m/s) | 8.73 d |
| semisub | sheared (same depth-mean) | **2.55 d** |

Going from a uniform to a depth-mean-matched, surface-intensified
current:
- the **spar's** mixing timescale barely changes (0.76 -> 0.62 days,
  ~18% faster) -- unsurprising, since its drag already spans the whole
  stratified column and so already experiences a broad range of current
  speeds (deep nearly-floor-value, shallow nearly-surface-value) even
  under the uniform case's one speed, in contrast the shear profile
  changes the balance only modestly for an already-distributed drag;
- the **semisubmersible's** mixing timescale collapses by a factor of
  **~3.4x** (8.73 -> 2.55 days) -- its shallow, concentrated drag now
  sits almost entirely within the fast part of the current profile
  (near-surface u up to 0.35 m/s, more than double the uniform case's
  0.15 m/s), so the cubic `u^3` production scaling gives it a large
  boost that a uniform current of the same bulk/mean strength could
  never deliver.

This confirms the project owner's hypothesis directly: **the
geometry comparison is not robust to how the background current's
momentum is distributed in the vertical.** A semisubmersible that looks
far less effective than a spar at distributed/bulk mixing under a
uniform current can become much more competitive (though still not
fully equivalent here -- 2.55 d vs. 0.62 d, roughly 4x slower) once a
realistic, surface-intensified current profile is used instead.

## Caveats

- Only one shear-profile shape/strength was tested (one e-folding
  scale, one surface/floor pair, chosen to match a plausible NCC
  description and one depth-mean value). The sensitivity to the shear
  profile's shape (e-folding depth relative to each structure's draft)
  is not explored -- a shallower or sharper profile would presumably
  favour the semisubmersible even more; a deeper/weaker one would
  favour it less. This is a single data point demonstrating the
  *existence* and *size* of the effect, not a full sensitivity sweep.
- The comparison is "fair" only in the sense of matching the
  depth-mean/column-integrated momentum input; it does not match
  (nor attempt to match) any other bulk quantity (e.g. depth-integrated
  kinetic energy, which is *not* conserved by this rescaling and is in
  fact higher for the sheared case by construction -- part of why it
  mixes faster even for the spar).
- As with the rest of this project, this is a 1D/column, horizontally
  homogeneous idealisation with no horizontal transport, no wave
  effects, and a fixed draft/geometry; the real structures' wake
  dynamics and 3D interactions are not represented.

## Practical implication

For the regional model: **do not assume the uniform-current
spar-vs-semisubmersible mixing comparison transfers directly to a
real, vertically sheared NCC.** The semisubmersible's disadvantage
relative to the spar, while real, is likely substantially smaller in
practice than the uniform-current idealisation suggests -- because
semisubmersibles are specifically positioned in the part of the water
column where real currents tend to be strongest. This is directly
relevant to comparing the regional mixing footprint of spar-based vs.
semisubmersible-based floating wind farm designs along the Norwegian
coast.

## Reproducing this note

```
python tools/prep_shear_geometry_sweep.py   # prepares shear_spar, shear_semisub
python tools/run_sweep.py sweeps/shear_geometry_sweep/manifest.yaml
python analysis/shear_geometry_sweep.py     # requires current_spar_u0p15 /
                                             # current_semisub_u0p15 already run
                                             # (see notes/current_strength_sensitivity.md)
```
