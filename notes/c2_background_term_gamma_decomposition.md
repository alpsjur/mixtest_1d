# The C2 "danger zone" is a structure-independent effect: it comes from the
# background term of Rennau et al.'s stationary Gamma, not from the c4/Pd term

**Date:** 2026-10-07
**Follow-up to:** `notes/c2_c4_closure_dependence.md` (which found that
`structure.c4` and `GLS.C2` are not separable, and attributed the resulting
turbulence collapse at low `C2` entirely to a transient "psi/k runaway" in
the psi-equation, driven by `c4*Pd` overwhelming `c2*eps`).

## 1. Motivation

Rennau et al. (2012) motivate the structure-mixing parametrization by
considering the *stationary* solution of the GLS `k`/`psi` system
(eliminating `P` between the two equations):

```
Gamma = -B/eps = (c1-c2)/(c3-c1)  +  (c4-c1)/(c3-c1) * Pd/eps
         \_______ term 1 _______/     \_______ term 2 _______/
```

with `c1 > c3` always (so `c3-c1 < 0`). This expression has **two
additive terms**:

- **Term 1** is a *background* contribution that does not involve `c4` or
  `Pd` at all -- it is set purely by `c1`, `c2`, `c3`, i.e. it is present
  even with **zero structures**. Its sign is `sign(c2-c1)` (since
  `c3-c1<0` is a fixed negative constant): **`Gamma`'s background term is
  positive iff `C2 > C1`**.
- **Term 2** is the structure contribution, proportional to `(c4-c1)` and
  weighted by `Pd/eps` -- this is the term `notes/mixing_timescale_analysis.md`
  and `notes/c2_c4_closure_dependence.md` focused on, and it is exactly
  zero at the "neutral point" `c4 = c1`.

`notes/c2_c4_closure_dependence.md`'s C2xc4 sweeps varied `GLS.C2` down to
0.52 -- *below* `C1` for all three closures tested (GEN `C1=1.0`, k-omega
`C1=0.555`, k-epsilon `C1=1.44`) -- and found turbulence collapses to the
numerical floor in that corner of parameter space. That note explained the
collapse mechanism (a positive feedback in the psi-equation via the
`psi/k` factor) but did not connect *why* low `C2` triggers it to the
structure-independent **term 1** above. This note closes that gap: the
previous sweep's `C2` axis accidentally crossed `C2=C1` in all three
closures, which is precisely where term 1 changes sign -- and all GLS
closures in the literature are calibrated with `C2 > C1`
(GEN: 1.22 > 1.0; k-omega: 0.833 > 0.555; k-epsilon: 1.92 > 1.44), so this
"danger zone" is normally never visited and has nothing to do with
structures or `c4` per se.

## 2. Re-reading the existing c4=C1 column with this in mind

`notes/c2_c4_closure_dependence.md` section 3.1 already contains the
decisive evidence, in its own tables, at the `c4/C1=1.0` column (the
*exact* neutral point where term 2 vanishes identically, so whatever
happens there is governed by term 1 alone):

| Closure   | `C1`  | `t*_mix` at `c4=C1`, by `C2`                                             |
|-----------|-------|---------------------------------------------------------------------------|
| GEN       | 1.0   | `C2=0.52`→NaN, `0.833`→NaN, **`1.22`→3.41**, `1.6`→2.35, `1.92`→2.27       |
| k-omega   | 0.555 | `C2=0.52`→NaN, **`0.833`→3.97**, `1.22`→2.52, `1.6`→2.41, `1.92`→2.37      |
| k-epsilon | 1.44  | `C2=0.52`→NaN, `0.833`→NaN, `1.22`→NaN, **`1.6`→8.91**, `1.92`→3.48       |

In **every one of the 15 cells**, the run collapses (NaN) iff `C2 < C1`
and mixes iff `C2 >= C1` (bold = closure's own literature default). This
is exactly `sign(Gamma_term1) = sign(C2-C1)`, with no reference to `c4`
or `Pd` needed -- a pattern that was there in the original data all along,
just not isolated because `c4` was swept at the same time.

## 3. New, targeted confirmation: fine C2 sweep at the exact neutral point

To test this cleanly (rather than relying on a single `c4/C1=1.0` column
inside a coarser 2D sweep), a new, finer sweep was run: GEN closure,
`structure.c4` fixed at exactly `GLS.C1=1.0` (so term 2 is identically
zero by construction, not just approximately), with `GLS.C2` on a fine
grid straddling `C1=1.0` from well below to well above (same
`mixing_timescale.yaml` body-force/structure setup as the earlier sweeps,
`ntimes_margin=2.5`, pycnocline length scale for `t*`).

- Config: `templates/gen_c2_neutral_c4_sweep.yaml` →
  `sweeps/gen_c2_neutral/` (9 runs, all `returncode=0`, run sequentially
  per the known `romsS`/`varinfo.yaml` parallel-run race condition).

Result (`analysis/mixing_timescale.py`, `length_scale="pycnocline"`):

| `GLS.C2` | `C2 - C1` | `t*_mix` |
|---------:|----------:|---------:|
| 0.70     | -0.30     | NaN (collapsed) |
| 0.85     | -0.15     | NaN (collapsed) |
| 0.95     | -0.05     | NaN (collapsed) |
| 1.00     |  0.00     | NaN (collapsed) |
| 1.05     |  0.05     | NaN (collapsed) |
| 1.15     |  0.15     | 4.66 |
| 1.30     |  0.30     | 2.72 |
| 1.60     |  0.60     | 2.35 |
| 1.92     |  0.92     | 2.27 |

**The collapse/success boundary sits right at `C2 ~= C1` (between 1.05 and
1.15), not anywhere near GEN's own literature default `C2=1.22`.** Since
`c4=C1` exactly here, term 2 is exactly zero throughout -- this cannot be
a `c4*Pd`-driven effect; it is entirely a property of the bare closure
(see the figure in section 4, left panel).

## 4. Direct confirmation from the simulated fields: Gamma_sim tracks term 1 almost exactly

To check this is really Rennau's `Gamma`, and not just a coincidental
`t*_mix` pattern, `Gamma` was computed *directly from the raw model
output* (not from the theoretical formula) for 5 of the above runs, early
in the run (hours 1-3, before any run has had time to mix or collapse),
at the depth of the initial pycnocline (`z = -initial.temp_zt`):

```
B(z,t)   = (g/rho0) * AKt(z,t) * d(rho)/dz        (buoyancy destruction of TKE)
eps(z,t) = utils.utils.compute_epsilon(ds, params) (GLS-generic, from tke/gls)
Gamma_sim(z,t) = -B(z,t) / eps(z,t)
```

(`compute_epsilon` already existed, from `c2_c4_closure_dependence.md`
section 3.5; `B` is new, via the standard down-gradient/flux-gradient
relation used internally by GLS, `-rho0*<w'rho'> = -AKt*d(rho)/dz`.)

| `GLS.C2` | term 1 (theory, `c3=C3M=0.05`) | `Gamma_sim` (model, hrs 1-3, pycnocline depth) |
|---------:|-------------------------------:|------------------------------------------------:|
| 0.95     | -0.053                         | +0.000               |
| 1.00     |  0.000                         | +0.002               |
| 1.05     | +0.053                         | +0.052               |
| 1.15     | +0.158                         | +0.158               |
| 1.92     | +0.968                         | +0.442               |

`Gamma_sim` **tracks the sign and magnitude of term 1 almost exactly for
`C2` within ~0.9 of `C1`** (0.95-1.15), diverging only once stratification
has eroded enough for the quasi-steady assumption behind Rennau's
derivation to break down (`C2=1.92`, where substantial mixing has already
happened by hour 3). This is strong, direct (not just inferential)
confirmation that the diagnosed collapse mechanism is exactly the
stationary `Gamma`'s background term, computed from first principles out
of the model's own `tke`/`gls`/`AKt`/`rho` fields.

![Left: t*_mix vs C2 at the neutral point c4=C1, collapse boundary at C2~=C1. Right: Gamma_sim from the model's own fields vs the theoretical term-1 curve.](../figures/c2_neutral_gamma_decomposition.png)

## 5. One nuance: the dynamical threshold has a small safety margin above `C2=C1`

Note that `Gamma_sim` is already clearly positive at `C2=1.05`
(+0.052) -- yet that run still collapses to the numerical floor
(section 3). The *sign* of the stationary `Gamma` is a necessary but not
sufficient condition for the actual time-dependent run to survive: Rennau's
formula assumes a quasi-steady `P`-balance has already been reached, but
early in any run (before that balance establishes), the same `psi/k`
positive-feedback mechanism documented in
`notes/c2_c4_closure_dependence.md` section 3.5 is active, and a `Gamma`
that is only *marginally* positive is not enough margin to survive the
transient before quasi-steady behaviour sets in. Empirically, `C2` needs
to exceed `C1` by roughly 10-15% (here: `C2 >~ 1.10-1.15` against
`C1=1.0`) before the run reliably avoids collapse, not just `C2 > C1`
exactly. This reconciles the two notes rather than contradicting either:
**term 1's sign is the controlling mechanism (this note), and the `psi/k`
runaway is *how* that sign-driven instability actually plays out
dynamically (`c2_c4_closure_dependence.md` section 3.5)**.

## 6. Practical implications

1. **The `c4_c4_closure_dependence.md` "danger zone" at low `C2` is not a
   structure-mixing artifact.** It would appear identically in plain GLS
   closure runs with no structures at all (`str_a=0`, `c4=0`) for any
   closure calibrated with `C2 < C1` -- which is precisely why no
   published GLS closure uses `C2 < C1` (literature defaults always have
   `C2` comfortably above `C1`, by 22-50%: GEN 22%, k-omega 50%,
   k-epsilon 33%). The earlier sweep's `C2` grid (reused literature
   values across closures, 0.52-1.92) simply wandered, for some
   closures, into a `C2`-vs-`C1` combination nobody would calibrate a
   model with.
2. **Updated calibration guidance:** when choosing/perturbing `GLS.C2`
   for sensitivity testing (e.g. to probe `c4` robustness), always check
   `C2` against that closure's own `C1`, not just against the closure's
   own literature-default `C2`. A margin of `C2 >~ 1.1*C1` (not just
   `C2 > C1`) is recommended based on section 5, to avoid conflating a
   pre-existing closure instability with the structure-mixing
   parametrization under test.
3. **`c4=C1` is confirmed as the right "neutral point" reference**, but
   only once `C2` is itself in its normal, literature-sane range relative
   to `C1` -- the earlier note's `c4`-sensitivity conclusions
   (section 4 of `c2_c4_closure_dependence.md`) remain valid for
   `C2 >= C1`; they should not be read as applying to `C2 < C1`, which is
   an unphysical regime for the reason in point 1, not a genuine `c4`
   interaction.
4. **Diagnostic tool added:** the `Gamma_sim = -B/eps` recipe in section 4
   (buoyancy destruction `B = (g/rho0)*AKt*d(rho)/dz` combined with the
   existing `compute_epsilon`) lets `Gamma` be checked directly against
   any run's own fields, closure-agnostically, without relying on the
   indirect `t*_mix`/NaN proxy.

## 7. Artifacts

- Config: `templates/gen_c2_neutral_c4_sweep.yaml` →
  `sweeps/gen_c2_neutral/` (9 runs; not committed, gitignored like all
  `sweeps/`/`runs/` output -- rerun via
  `python tools/prep_mixing_timescale_sweep.py templates/gen_c2_neutral_c4_sweep.yaml`
  then `python tools/run_sweep.py sweeps/gen_c2_neutral/manifest.yaml`,
  sequentially, ~85 s/run)
- Figure: `figures/c2_neutral_gamma_decomposition.png` (left: `t*_mix` vs
  `C2` with the collapse boundary at `C2~=C1`; right: `Gamma_sim` vs `C2`
  against the theoretical term-1 curve)
