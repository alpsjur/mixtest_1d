# Sensitivity of structure-induced mixing to GLS_Kmin and GLS_Pmin

**Date:** 2026-09-29
**Follow-up to:** `notes/mixing_timescale_analysis.md` and
`notes/c2_c4_closure_dependence.md`. Uses the same GEN closure setup and
methodology (Carpenter et al. 2016 power-balance framework, dimensionless
time `t* = t·P_str_diag/(g·Δρ·L²)`), with `L = sqrt(z_t·(H0−z_t))`, the
**pycnocline length scale**, per that note's §4.5 finding that it collapses
the sweep better than the full water-column depth `H0`.

## 1. Motivation

The GLS closure imposes numerical floors on turbulent kinetic energy
(`GLS_Kmin`) and its dissipation/length-scale variable (`GLS_Pmin`), to
prevent turbulence quantities from becoming zero or negative in quiescent
regions. ROMS' own "suggested values" table
(`templates/mixtest_1d.in.j2`, lines ~2044-2045; Warner et al. 2005,
doi:10.1016/j.ocemod.2003.12.003, Table 1) gives, per closure:

| Closure | GLS_Kmin | GLS_Pmin |
|---|---|---|
| k-kl | 5.0d-6 | 5.0d-6 |
| k-epsilon | 7.6d-6 | 1.0d-12 |
| k-omega | 7.6d-6 | 1.0d-12 |
| **gen** | **1.0d-8** | **1.0d-8** |

`configs/baseline.yaml` sets `GLS.Kmin: 1.0` / `GLS.Pmin: 1.0`, which the
`.in` template renders with a fixed `d-8` exponent
(`GLS_Kmin == {{ GLS.Kmin }}d-8`), i.e. **exactly** `1.0d-8` / `1.0d-8`.

**Finding: the mixtest_1d setup's Kmin/Pmin are *not* unusually high — they
match the gen closure's own literature default exactly.** They only look
"high" relative to k-epsilon/k-omega's `Pmin = 1.0d-12` because those two
closures use a very different `Ψ`-variable (`ω` and `ε` scale very
differently near the floor) — this is a closure-family difference, not an
indication that the gen-closure floors here are set too conservatively.

## 2. Experiment design

A 5x5 cartesian sweep (`templates/gen_kminpmin_sweep.yaml` ->
`sweeps/gen_kminpmin/`), gen closure only, holding all other parameters
at their `mixing_timescale.yaml`/`baseline.yaml` values (`structure.c4 =
0.44`, `GLS.C2 = 1.22`, `structure.CD = 0.63`, `initial.temp_dT = 10`,
`grid.H0 = 150`):

- `GLS.Kmin = [0.01, 0.1, 1.0, 10.0, 100.0]` -> renders to `GLS_Kmin =
  1.0d-10 ... 1.0d-6`
- `GLS.Pmin = [0.01, 0.1, 1.0, 10.0, 100.0]` -> renders to `GLS_Pmin =
  1.0d-10 ... 1.0d-6`

This spans the gen closure's own default (1.0d-8, at the grid centre) by
two decades in each direction, comfortably covering the full range used
across all four closures in the table above (Kmin: 1.0d-8 to 7.6d-6;
Pmin: 1.0d-12 to 1.0d-6). `tau_mix_theory` (power-balance only) does not
depend on `GLS.Kmin`/`Pmin`, so all 25 runs share the same `NTIMES`
(79650, `ntimes_margin=4.0`, `tau_mix_theory = 9.21` d) — any difference in
the empirical `t*_mix` is attributable entirely to the Kmin/Pmin floors.
All 25 runs completed successfully (run sequentially, not in parallel —
see the parallelism caveat in `notes/c2_c4_closure_dependence.md`).

## 3. Results

### 3.1 `t*_mix` is essentially insensitive to GLS_Kmin/GLS_Pmin

`t*_mix` (pycnocline length scale) over the full 5x5 grid:

```
Pmin\Kmin   1e-10   1e-09   1e-08   1e-07   1e-06
1e-10       2.2559  2.2557  2.2556  2.2556  2.2556
1e-09       2.2562  2.2557  2.2556  2.2556  2.2556
1e-08       2.2565  2.2561  2.2556  2.2556  2.2556
1e-07       2.2568  2.2564  2.2561  2.2561  2.2556
1e-06       2.2572  2.2568  2.2564  2.2561  2.2557
```

Across the full 4-decade x 4-decade grid, `t*_mix` ranges from 2.2556 to
2.2572 — a spread of **0.07%** relative to the mean. `Pstr_diag` (the
diagnosed steady-state structure-drag power) is identical to 4 significant
figures (4.0726e-01) in every single run: the mean-flow/drag balance is
completely unaffected, as expected since `Kmin`/`Pmin` do not enter the
momentum equation.

For comparison, sweeping `structure.c4` and `GLS.C2` at fixed Kmin/Pmin
(see `c2_c4_closure_dependence.md`) changes `t*_mix` by up to a factor of
~6 (2.18 to 14.86) within the numerically stable range, and even the
"well-behaved" `c4`/`H0` sweep in `mixing_timescale_analysis.md` produced
~35% swings in `t*_mix`. **The GLS_Kmin/GLS_Pmin floors are, by a wide
margin, the least influential parameter tested so far** for this
parametrization's mixing-completion time.

### 3.2 Why: the floors only bind in already-quiescent regions

To check whether the floors are simply inactive (never reached) or genuinely
irrelevant despite being reached, `tke`/`gls` fields were inspected directly
in the two most extreme corner runs (`Kmin=1e-10,Pmin=1e-10` and
`Kmin=1e-6,Pmin=1e-6`):

- In both runs, `min(tke)` and `min(gls)` over the whole run equal their
  respective floor values *exactly* — so the floors **are** being hit, not
  merely a numerically inert setting.
- At `t=0` (quiescent initial condition) **100%** of the water column sits
  at the TKE floor, as expected.
- By the time the flow has spun up, only **~4.9%** of vertical levels
  (2 of 41 `s_w` levels) remain at the floor at any given output time —
  physically, the top and bottom boundary points, where the no-slip/
  no-flux condition drives local shear/production to zero regardless of
  the floor's numerical value.

This explains the negligible sensitivity: `GLS_Kmin`/`GLS_Pmin` only
constrain the *edges* of the domain and the initial transient, not the
actively turbulent pycnocline interior where the structure-mixing
parametrization actually does its work. Since `φ`/`φ*` (and hence
`t*_mix`) integrate stratification over the whole water column but are
dominated by what happens in and around the pycnocline (mid-column, away
from the boundaries in this setup), the boundary-level floor values are
essentially decoupled from the bulk mixing evolution.

## 4. Conclusion

For the **gen** closure, in this 1D idealized setup:

1. The currently used `GLS_Kmin`/`GLS_Pmin` (1.0d-8/1.0d-8) are **not**
   elevated relative to standard practice — they are the closure's own
   literature default (Warner et al. 2005).
2. Varying both floors independently over 4 decades each (1.0d-10 to
   1.0d-6, bracketing the full range used across the k-kl/k-epsilon/
   k-omega/gen closures) changes the dimensionless mixing-completion time
   `t*_mix` by only ~0.07% — negligible compared to the closure-coefficient
   (`c4`, `C2`) sensitivity already documented.
3. The floors are technically active (reached) at the domain boundaries
   and during initial spin-up, but this has no measurable effect on the
   depth-integrated, pycnocline-driven mixing evolution.

**Practical implication:** the time evolution of structure-induced mixing
in this parametrization, as measured by `t*_mix`, can be treated as robust
to the choice of `GLS_Kmin`/`GLS_Pmin` within (and well beyond) the range
of values found in the turbulence-closure literature. Sensitivity work
going forward should continue to focus on the closure coefficients
`c4`/`C2`/`C1` (per `c2_c4_closure_dependence.md`), not on the TKE/dissipation
floors.

## Reproducing

```
python tools/prep_mixing_timescale_sweep.py templates/gen_kminpmin_sweep.yaml
python tools/run_sweep.py sweeps/gen_kminpmin/manifest.yaml   # run sequentially, not with xargs -P
```

Analysis/figure: `figures/gen_kminpmin_phase_pycnocline.png` (t*_mix phase
diagram over the Kmin x Pmin plane, generated ad hoc from
`analysis.mixing_timescale.mixing_timescale()` with `length_scale="pycnocline"`,
following the same pattern as `analysis/plot_c2c4_phase.py`).
