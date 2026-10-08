# Residence time vs. mixing timescale: does a finite farm length reverse the current-speed trend?

**Context.** `notes/current_strength_sensitivity.md` found that faster
background currents mix a water column dramatically faster (the mixing
timescale falls steeply with current speed, roughly consistent with the
cubic `P_d ∝ u³` scaling in Carpenter et al.). But a real wind farm has a
finite along-flow extent `L`: a water column is only exposed to the farm's
structure drag for a residence time `T_res = L/u` while it crosses the
farm. As current speed increases, a column both (a) mixes faster *and* (b)
spends less time inside the farm -- two effects pulling in opposite
directions. This note (prompted by a direct question from the project
owner) asks which one wins: does the amount of mixing actually delivered
to a water column crossing a real, finite-length farm increase or decrease
with current speed?

This reuses the already-completed `current_strength_sensitivity` sweep
(no new ROMS runs) -- each run's own `phi*(t)` curve is read off at
`t = T_res(u) = L/u` instead of at a fixed absolute time, for a few
representative farm lengths.

## Method

For each of the 5 diagnosed current speeds per geometry (spar-like,
draft=90 m; semisub-like, draft=20 m -- same matched-`total_a=1.5` runs as
the current-strength note), and for farm lengths `L ∈ {2, 5, 10, 20} km`
(spanning "a handful of turbine spacings" to "a large array" -- Hywind
Tampen's 11-turbine array spans roughly this scale), the residence time
`T_res = L/u` is computed and `phi*(T_res)` is linearly interpolated from
that run's own `phi*(t)` time series (every `T_res` tested here is well
within the 46-day simulated window, so no extrapolation is needed).

## Result: faster currents still deliver *more* mixing by farm exit, not less

![Stratification remaining at farm exit vs. current speed, for several farm lengths](../figures/residence_time_tradeoff.png)

For every farm length tested and both geometries, `phi*(T_res)` **decreases
monotonically with current speed** -- i.e. a faster current still leaves a
water column *more* mixed by the time it exits the farm, despite spending
less absolute time inside it. The effect can be large: for the spar-like
geometry and a 10 km farm, a column crossing at 0.048 m/s retains 92% of
its original stratification (`phi*=0.92`) by exit, while one crossing at
0.244 m/s is already fully mixed (`phi*≈0`) well before it exits.

This is the expected outcome of the scaling argument once both timescales
are written out: the stratification-mixing timescale falls roughly as
`t_mix ∝ u⁻³` (cubic structure-drag production), while the residence time
only falls as `T_res ∝ u⁻¹`. Their ratio `T_res/t_mix ∝ u²` therefore
*grows* with speed -- a water column gets proportionally more "mixing
timescales" of exposure per farm transit as the current speeds up, not
fewer, even though the absolute transit time shrinks. The residence-time
penalty is real (it does reduce the *absolute* exposure time), but it is
outpaced by how much faster the mixing process itself becomes.

## Caveats / limits of this result

- This uses the already-run current-speed grid (5 points, `u` up to
  ~0.24-0.25 m/s); the monotonic decrease is demonstrated over that range,
  not proven as `u → ∞`. The underlying scaling argument (`T_res/t_mix ∝
  u²`) suggests the trend should continue, not reverse, at higher speeds,
  but this has not been tested with new runs.
- `T_res = L/u` assumes a column crosses the whole farm length at the
  (uniform, steady) background current speed -- it ignores any
  current-speed reduction *caused by* the farm's own drag (the structure
  drag that produces the mixing also removes momentum from the flow,
  which would tend to slow the column down inside the farm and so
  *increase* its real residence time beyond this simple estimate -- a
  second-order effect not captured here, and in the opposite direction
  needed to reverse this note's conclusion, i.e. ignoring it is
  conservative with respect to the result above).
- Real currents are unsteady (tidal, wind-driven pulses); a single
  steady-state transit speed is a simplification, consistent with the
  rest of this idealised 1D framework (see `current_strength_sensitivity.md`
  Section 5 for the same caveat).
- This only addresses a single along-flow pass. Spatial mixing patterns
  inside a real (2D/3D, heterogeneous-density) farm footprint, including
  wake superposition between turbines noted in the literature search
  earlier in this project, are not represented in this idealised
  column-only result.

## Practical implication

For the regional model: **do not assume that weaker/slower-current regions
of the Norwegian Coastal Current see proportionally more "exposure time"
to a given farm's mixing effect and therefore comparable net impact to
faster regions.** This result suggests the opposite -- the faster-current
segments of the NCC (or faster tidal/wind-driven episodes) should produce
disproportionately *more* net mixing per farm transit, not less, reinforcing
the intermittency conclusion from `current_strength_sensitivity.md`
(mixing dominated by the current's faster episodes) at the farm-footprint
scale as well as the single-column scale.

## Reproducing this note

```
# (requires the current_strength_sensitivity sweep already run, see that note)
python analysis/residence_time_tradeoff.py
```
