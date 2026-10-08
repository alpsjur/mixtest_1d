# Does surface wind stress combine with structure-induced mixing, or does it get swamped?

**Context.** Every mixing experiment in this repo so far uses structure
drag as the *only* source of turbulence -- the body force itself drives a
uniform, shear-free background flow with no TKE production of its own.
Real floating wind farms sit in a water column that is *also* being mixed
by wind-driven surface shear, and the long-term goal of this project is to
understand regional, realistic impacts -- so a natural next question is
whether structure-induced and wind-induced mixing simply add, or interact
nonlinearly (plausible, since the GLS TKE/dissipation equations are
nonlinear in their production terms). This note builds the infrastructure
to apply a real surface wind stress in `mixtest_1d`, and reports what the
first attempt at this additivity/synergy test found.

This is an exploratory extension, run autonomously; the experimental
choices below are defensible but not yet reviewed by a human.

## 1. Infrastructure built

`mixtest_1d` previously had no way to apply a non-zero surface momentum
stress: `roms-related/Include/mixtest_1d.h` had `ANA_SMFLUX` defined, and
the compiled analytic routine (`ana_smflux.h`) has no custom case for this
application, so it silently fell through to `sustr=svstr=0`. The following
pieces were added (modelled on `mixtest_3d/tools/make_frc.py`, used only as
a design reference -- that implementation there is itself an unfinished
stub, not something that could be imported directly):

- **`tools/make_frc.py`**: writes a ROMS forcing NetCDF (`sustr`/`svstr` on
  the u/v grids, kinematic units `m2/s2` -- note this is a genuine unit
  difference from `mixtest_3d`'s version, which labels the same fields
  "N/m2"; ROMS' own `ana_smflux.h` comments confirm the kinematic
  convention is what the model expects). Always generates a file (constant
  in time by default, zero stress unless `cfg["forcing"]` overrides it),
  since the compiled model now has no analytic fallback.
- **`configs/baseline.yaml`**: new `files.frc` key and a commented-out
  `forcing:` section documenting `sustr_tau0`/`svstr_tau0`/`dt_frc`.
- **`templates/mixtest_1d.in.j2`**: `FRCNAME` now points at the generated
  file instead of a hardcoded, non-functional placeholder.
- **`tools/prep_experiment.py`**: wired in the new `make_frc_from_config()` call.
- **`roms-related/Include/mixtest_1d.h`**: `ANA_SMFLUX` undefined; ROMS
  rebuilt successfully (`roms-related/build_roms.sh`).
- **Regression check**: re-ran the existing `geom_matched_draft90` config
  through the new pipeline (zero-stress forcing file) and compared the
  full history output against the original (pre-rebuild) run -- **every
  field (`temp`, `salt`, `rho`, `u`, `v`, `tke`) is bit-identical**
  (max abs diff = 0.0). The switch from analytic to file-based zero stress
  changed nothing; every previous result in this repo remains valid.

## 2. Experiment design

A 2x2 matrix at the matched spar-like geometry (draft=90 m,
`str_a=1.5/90`), fixed background current (`U0=0.1 m/s`), same
stratification and 46.08-day window as the other floating-structure
experiments:

- structure drag: **off** (`str_a=0`) or **on**
- surface wind stress: **off** or **on**, a single "moderate" 10 m/s wind
  via the standard bulk formula `tau = rho_air * Cd10 * W10^2` (`Cd10 =
  1.3e-3`, `rho_air = 1.225 kg/m3`) -> `tau = 0.159 Pa` -> kinematic
  stress `1.55e-4 m2/s2`, constant in time, directed along the background
  current (+x).

## 3. Result: wind stress has essentially no effect in this setup

![Stratification decay and near-surface TKE, structure vs. wind](../figures/wind_structure_synergy.png)

**Left panel:** the structure-off/wind-off and structure-off/wind-on
curves are visually (and numerically, to 3+ significant figures)
identical, as are the structure-on/wind-off and structure-on/wind-on
curves. **Adding wind stress makes no detectable difference to bulk
mixing, with or without the structure.** This was checked again with a
10x stronger ("storm-force", ~35 m/s equivalent) wind stress
(`synergy_windcheck_strong`) -- still no detectable change in `phi*(t)`.

**Right panel:** the wind stress is not simply failing to couple into the
model -- near-surface TKE does respond (control ~1e-8 m2/s2 floor value;
moderate wind ~1e-7–1e-6; storm wind higher still), confirming
`sustr`/`svstr` are read and used. But the response stays confined to the
top ~1-2 m and never builds up enough to measurably affect the
bulk/column-integrated stratification.

## 4. Why: the body force is a full-column Rayleigh damper

The reason is structural, not a units or sign-convention bug.
`UV_BODYFORCE` (used to sustain the background current in this idealised
setup) adds an explicit relaxation term `du/dt = bfrc_u - BFRC_CD*u`
**at every depth level, including the surface**, with `BFRC_CD = 0.005
s⁻¹` fixed throughout. This gives a relaxation timescale `1/BFRC_CD = 200
s` -- any wind-driven perturbation to the surface current is damped back
toward the prescribed target (`U0=0.1 m/s`) within a few minutes, long
before enough shear or TKE can be generated or mixed down into the
stratified interior. Structure drag, by contrast, is *not* damped this way
-- it acts through the (undamped) quadratic term in the TKE production
budget, which is why it alone produces dramatic mixing while a comparably
realistic wind stress produces essentially none.

**This is a property of the idealised body-force mechanism used throughout
this repo to sustain a background current, not a statement about real
wind/structure mixing interactions.** A free-running model (no artificial
volumetric relaxation at every depth -- i.e. any realistic regional
configuration, including the planned Norwegian-coast regional model) would
not have this suppression, since real currents are sustained by pressure
gradients and boundary/surface forcing rather than a uniform body force
nudging every grid cell back to a fixed target.

## 5. Practical implication

**This idealised 1D framework, as currently built, cannot be used to test
wind/structure mixing additivity or synergy** -- the artificial
full-column current-sustaining mechanism masks it. Two options for anyone
wanting to pursue this question further:

1. Redesign the body force so it does not directly damp the surface
   layer(s) (e.g. restrict `UV_BODYFORCE` to depths below some distance
   from the surface, or switch to a pressure-gradient-driven setup) --
   would require new Fortran/CPP work, not attempted here.
2. Test the additivity/synergy question directly in the regional,
   realistic model instead, where the background current is not
   artificially sustained by a volumetric relaxation term. Given the
   infrastructure built here is verified and reusable
   (`tools/make_frc.py`'s approach/unit convention), the main transferable
   result for the regional model is this: **don't assume wind-driven and
   structure-driven mixing contributions are independent/additive without
   checking** -- the mechanisms compete for the same near-surface momentum
   budget, and in a free-running model (unlike here) that competition
   could go either way (reinforcing or suppressing), not just vanish.

## Reproducing this note

```
python tools/prep_wind_structure_synergy.py
python tools/run_sweep.py sweeps/wind_structure_synergy/manifest.yaml
python analysis/wind_structure_synergy.py
```

(the storm-force diagnostic run `synergy_windcheck_strong` was prepared
and run as a one-off check, not via a sweep script; see this note's
history for the exact override used: `forcing.sustr_tau0 = 1.5e-3`, all
else identical to `synergy_struct_off_wind_off`.)
