# No "Goldilocks" draft: mixing timescale vs. structure draft (fixed drag density, sheared current)

## Motivation

The shear × geometry experiment (`notes/shear_geometry_interaction.md`) compared only
two drafts (spar, 90 m; semisub, 20 m) at a *matched total drag budget*
(`total_a = 1.5` for both). That comparison can't tell us whether draft itself
has a non-monotonic effect — e.g. a hypothetical worst-case (or best-case) draft
for mixing, for a structure of some fixed, realistic drag density that simply
sits deeper or shallower in the water column.

This note sweeps draft continuously (10–150 m, 15 runs) under the same
surface-intensified shear profile, this time holding `structure.str_a` (drag
per unit depth) **fixed** at a realistic value instead of matching total drag.
That means `total_a = str_a × draft` **increases** with draft in this sweep —
a deeper structure of the same cross-sectional drag density simply has more
total drag. This is a deliberately different convention from the matched-budget
comparisons used elsewhere in this project, and the two should not be confused.

`str_a = 0.018375` is anchored to Hywind Tampen's real spar geometry (hull
diameter 14.7 m, 500×500 m spacing → physical frontal-area density
`8/500²·14.7/8 ≈ 3.2e-5 m⁻¹`... using the project's established scaling
convention from `notes/floating_structure_depth_distribution.md`:
`str_a = (14.7/500²) × 312.5 = 0.018375`, so draft=90 m reproduces the real
Hywind Tampen spar's total drag, and the other 14 drafts are hypothetical
"same structure, different draft" extensions of it.

| draft (m) | total_a (= str_a × draft) |
|---:|---:|
| 10  | 0.18 |
| 20  | 0.37 |
| 40  | 0.74 |
| 90  | 1.65 (real Hywind Tampen spar) |
| 150 | 2.76 |

## Result

![Mixing timescale vs. structure draft](../figures/goldilocks_draft_depth.png)

There is **no Goldilocks minimum or maximum** — `t(phi*=0.5)` and
`t(phi*=0.1)` both decrease **monotonically** as draft increases, with a
strongly saturating (diminishing-returns) shape:

- Shallow structures (draft ≤ 40 m, i.e. draft ≲ pycnocline depth) show the
  steepest sensitivity: `t(phi*=0.5)` falls from 14.5 days (draft=10 m) to
  1.4 days (draft=40 m) — a >10x change over just 30 m of draft.
- Beyond the pycnocline depth (40 m), the curve flattens sharply: from
  draft=90 m (0.58 days) to draft=150 m (0.50 days) is only a 15% change
  despite a 60 m increase in draft (and nearly doubling `total_a`, from 1.65
  to 2.76).

In other words: **once a structure's draft reaches down to (and past) the
pycnocline, adding more draft keeps buying faster mixing, but with rapidly
diminishing returns.** The interesting part of the curve — where draft
differences matter most for mixing speed — is entirely in the upper ~40 m,
which is exactly where semisubmersibles (draft ≈ 20 m) sit and spars
(draft ≈ 90 m) do not.

## Interpretation

This result is consistent with (not contradictory to) the earlier
shear-geometry finding. There, we found that shear closes most of the
spar/semisub mixing-timescale *gap* (semisub speeds up 3.4x under shear,
spar barely changes) but doesn't fully erase it. This sweep shows why: the
semisub's draft (20 m) sits on the steep part of this curve, where small
differences in effective draft (or in how much of the sheared, energetic
near-surface current it can tap) translate into large differences in mixing
speed. The spar's draft (90 m) sits on the flat part, where its total drag
budget is already large enough that it operates in a saturated,
input-insensitive regime.

Practically, this means: **for OWF regional-impact purposes, draft
differences between floater types matter most exactly in the upper-ocean
layer where semisubmersibles live** — predicted mixing impact is far more
sensitive to small differences in semisubmersible draft/geometry than to
equivalent differences for spars, which only amplify an already-large
starting effect.

### Caveats

- This sweep **confounds** two effects that both scale with draft here: (1)
  more total drag (`total_a`), and (2) more of the water column exposed to
  the surface-intensified shear current. It cannot separate "deeper draft
  increases total drag" from "deeper draft taps a different part of the
  current profile" — disentangling these would require an additional sweep
  holding `total_a` fixed while varying draft under the *same* shear profile
  (as opposed to the matched-budget, 2-point `shear_spar`/`shear_semisub`
  comparison already done).
- Only one current shear shape (`profiles/ncc_shear_surface_intensified.txt`)
  and one pycnocline depth (40 m, the project default) were tested; the kink
  location is almost certainly profile-dependent.

## Reproduction

```bash
mamba activate roms   # velocity-shear-bodyforce branch, roms + mixtest_1d
python tools/prep_draft_shear_sweep.py     # prepares 15 runs, runs/draft_shear_d*
python tools/run_sweep.py sweeps/draft_shear_sweep/manifest.yaml
python analysis/draft_shear_sweep.py       # -> figures/goldilocks_draft_depth.png
```

## Bonus: "mixing race" animation

As a companion visualization (not quantitative, just intuition-building),
`analysis/mixing_race_animation.py` animates the spar vs. semisub density
profiles (from the matched-budget `shear_spar`/`shear_semisub` runs) eroding
side by side, together with a `phi*(t)` "race track":

```bash
python analysis/mixing_race_animation.py   # -> figures/mixing_race.gif
```

![mixing race animation (static preview)](../figures/mixing_race.gif)
