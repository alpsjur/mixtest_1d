# Can the structure-mixing parametrization be checked against observations?

**Motivation.** The structure-mixing parametrization in this project has
two free parameters that were not derived from first principles for this
specific application: `structure.CD` (drag coefficient, taken from
general structure-hydrodynamics literature) and `structure.c4` (the GLS
psi-equation structure-production weight, currently `c4=0.44`, chosen as
"a reasonable guess for realistic mixing" based on the sensitivity
behaviour explored in `notes/c2_c4_closure_dependence.md`, rather than
fit to any specific dataset). This note asks: is there any observational
or field-validated benchmark we can check the parametrization's *output*
against, even approximately, given how idealised this whole 1D setup is?

## The benchmark: Schultze et al. (2020)

A literature search for field/observational studies of structure-induced
mixing in offshore wind contexts turned up a directly relevant paper by
(part of) the same research group behind the Rennau et al. (2012)
parametrization used here:

> Schultze, L.K.P., Merckelbach, L.M., Horstmann, J., Raasch, S., &
> Carpenter, J.R. (2020). *Increased mixing and turbulence in the wake
> of offshore wind farm foundations.* Journal of Geophysical Research:
> Oceans, 125, e2019JC015858.

Using glider/ADCP-based field observations combined with LES around a
real monopile foundation, they report a **bulk mixing efficiency of
8-14%** -- i.e. roughly a tenth of the turbulent kinetic energy
generated in the structure's wake goes into irreversibly mixing the
stratified water column, with the rest dissipated as heat. This
"mixing efficiency" is, by construction, exactly the same physical
quantity as this project's `Gamma = -B/epsilon` (buoyancy destruction of
TKE divided by dissipation) -- it is an **independent, field-constrained
number that this model's diagnosed `Gamma` can be checked against,
without needing to run anything new.**

(Caveat up front: Schultze et al.'s number is for a monopile -- a rigid,
near-full-depth, bottom-fixed cylinder in ~30-40 m water -- not a
floating spar/semisubmersible. It is the closest available field-
validated number for "a GLS-type structure-mixing parametrization in an
ocean-model-relevant regime," but should be read as an order-of-magnitude
sanity check, not a precise calibration target for this project's
specific geometries.)

## Method: diagnosing Gamma directly from already-completed runs

No new simulations were needed for this check -- `Gamma(z,t)` was
computed from the existing current-strength/geometry runs'
history-file output, following exactly the same derivation already
validated in `notes/c2_c4_closure_dependence.md`:

```
B(z,t)       = (g/rho0) * AKt(z,t) * d(rho)/dz        (buoyancy flux)
epsilon(z,t) = CMU0^(3+P/N) * k^(3/2+M/N) * psi^(-1/N) (GLS-generic, from tke/gls)
Gamma(z,t)   = -B(z,t) / epsilon(z,t)
```

`Gamma` is then averaged, **weighted by `epsilon`** (not a plain mean --
`epsilon` can be locally very small, which would otherwise let division
noise dominate a plain average), over each structure's depth range and
an early time window (0.1-0.4 days into the run -- well before the
column approaches full mixing, so there is still a meaningful buoyancy
flux to diagnose).

See `analysis/gamma_observational_check.py`.

## Result: the single fixed c4=0.44 does not land consistently inside the observed range

![Diagnosed bulk Gamma for spar and semisubmersible geometries vs. the Schultze et al. (2020) field/LES range](../figures/gamma_observational_check.png)

| run | geometry | current | diagnosed bulk Gamma |
|---|---|---|---|
| `current_spar_u0p15` | spar (draft=90 m) | 0.15 m/s | 0.217 |
| `geom_matched_draft90` | spar (draft=90 m) | 0.10 m/s | 0.237 |
| `current_semisub_u0p15` | semisub (draft=20 m) | 0.15 m/s | 0.018 |
| `geom_matched_draft20` | semisub (draft=20 m) | 0.10 m/s | 0.026 |

With the project's current `c4=0.44` choice, the two "matched" geometries
diagnose to **very different** mixing efficiencies, on opposite sides of
the Schultze et al. observed band:

- the **spar** (drag spread over the full 90 m column) diagnoses to
  Gamma ~ 0.22-0.24 -- roughly **1.5-3x above** the observed 8-14% range;
- the **semisubmersible** (drag concentrated in the shallow 20 m draft)
  diagnoses to Gamma ~ 0.02-0.03 -- roughly **3-7x below** it.

This is consistent (diagnosed Gamma is not very sensitive to current
speed within a geometry) but strongly geometry-dependent. In other
words: **a single, geometry-independent `c4` cannot simultaneously
reproduce the Schultze et al. benchmark for both structure types as
currently configured** -- at least not for the "matched total depth-
integrated drag" geometries used throughout this project.

## Interpretation

This is not necessarily a flaw in the parametrization itself -- `Gamma`
in the Rennau et al. (2012) derivation is a property of the *local*
flow/turbulence state (via the `Pd/epsilon` ratio and the closure
constants), not a universal constant, so some geometry dependence is
expected in principle. But the size and *direction* of the mismatch here
is informative:

- the spar's drag, spread over the full water column including the
  weakly-turbulent deep interior, generates comparatively more buoyancy
  flux per unit dissipation in this diagnostic window than the field
  benchmark suggests is realistic;
- the semisubmersible's drag, concentrated at the energetic, well-mixed
  near-surface layer where ambient turbulence/dissipation is already
  large, generates comparatively *less* mixing per unit dissipation than
  the benchmark suggests.

This gives a concrete, falsifiable path to at least a first-order
calibration exercise, rather than treating `c4` as "a reasonable guess"
indefinitely:

1. **Rerun a small `c4` sensitivity sweep per geometry** (the
   infrastructure already exists from `notes/c2_c4_closure_dependence.md`)
   and check, with the *same* diagnostic used here, what `c4` value
   brings each geometry's diagnosed `Gamma` into the observed 8-14%
   band. If very different `c4` values are needed for the two
   geometries, that is itself a finding worth reporting (either `c4`
   genuinely should depend on how a structure's drag is distributed in
   depth, or the matched-total_a framework used to compare geometries
   throughout this project needs revisiting).
2. Treat the Schultze et al. monopile number as a floor/sanity check
   only -- a literature search for any mixing-efficiency estimates
   specific to floating structures (spar/semisubmersible mooring lines,
   catenary/taut systems) would substantially strengthen this
   comparison if such data exists; none was found in this pass.
3. Beyond `Gamma`, Schultze et al. also report the *spatial* structure
   of the wake (narrow, energetic, confined to within ~100 m
   downstream) -- a comparison this 1D, horizontally homogeneous model
   cannot make at all (it has no concept of "downstream distance"). A
   3D or along-flow-resolving regional-model run would be needed to
   check that aspect; this remains squarely a limitation of the
   idealised approach used throughout this project, and is flagged here
   rather than attempted.

## Caveats

- The Schultze et al. benchmark is for a rigid monopile in ~30-40 m
  water depth under North Sea conditions -- not a floating
  spar/semisubmersible, not NCC conditions, and not derived from the
  same GLS-type closure used here (their number comes from LES +
  field glider/ADCP data, not a RANS/GLS ocean model). Treat the
  comparison as an order-of-magnitude sanity check only.
- The diagnosed `Gamma` here depends on the (somewhat arbitrary) choice
  of depth mask and time window; a brief robustness check across three
  time windows (0.05-0.2 d, 0.1-0.5 d, 0.2-0.6 d) showed the diagnosed
  values are stable to within ~15% of each other for a given run,
  which is small compared to the order-of-magnitude mismatch reported
  above, but this was not exhaustively tested.
- No attempt was made here to analytically back-solve for a "corrected"
  `c4` using the stationary Rennau et al. formula -- the diagnosed
  `Gamma` is a property of the full, transient, spatially-resolved
  model (not the idealised 0-D stationary system the formula assumes),
  so any such back-solved estimate would be unreliable without actually
  rerunning at the suggested value and re-diagnosing.

## Reproducing this note

```
# (requires current_spar_u0p15, current_semisub_u0p15, geom_matched_draft90,
#  geom_matched_draft20 already run -- see notes/current_strength_sensitivity.md
#  and notes/floating_structure_depth_distribution.md)
python analysis/gamma_observational_check.py
```
