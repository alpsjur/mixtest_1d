# The role of GLS closure coefficients c2 and c4 in structure-induced mixing

**Date:** 2026-09-28 (updated with k-epsilon and pycnocline length scale)
**Follow-up to:** `notes/mixing_timescale_analysis.md` (which fixed `GLS.C2` at each
closure's literature default and varied `structure.c4` alone).

## 1. Motivation

The structure-mixing parametrization (see `roms/README_STRUCTURE_MIXING.md`) adds
one new GLS coefficient, `c4` (`c_psi4`), which plays the same role for
structure-drag production `P_d` that `c1` (`c_psi1`) plays for shear production
in the psi-equation:

```
d(psi)/dt = ... + (psi/k) * (c1*P + c4*Pd + c3*G - c2*eps)
```

Everything so far (`mixing_timescale_analysis.md`, `c4_fine_sweep.yaml`,
`c4_near_c1_sweep.yaml`) varied `c4` at a *fixed* `C2` (the closure's own
literature default) and found that the empirical mixing-completion time
`t*_mix` depends on `c4` mainly very close to `C1` (the "neutral point"). This
raised an obvious question: is that `c4`-sensitivity actually independent of
`C2`, or does `C2` control *how much* `c4` matters? `C2` is the destruction
coefficient for `psi` -- exactly the term that balances `c4*Pd` -- so a
coupling is physically plausible. We then extended the question to a third
closure, k-epsilon, to see whether any (c2,c4) interaction found is generic to
the GLS family or specific to one closure's exponents.

**Methodological note on the length scale used for `t*`:** the dimensionless
time `t* = t / tau_mix_theory` uses a length scale `L` inside
`tau_mix_theory ~ (rho0*L^2 / P)^(1/3)` (Carpenter et al. power-balance
estimate). The first pass of this analysis used `L = H0` (full water-column
depth, 150 m). That is not the physically relevant scale -- the parametrized
mixing acts on the initial pycnocline, whose thickness scale is
`L_pycnocline = sqrt(zt*(H0-zt))` with `zt = initial.temp_zt = 40 m`, giving
`L_pycnocline ~ 66.3 m` versus `H0 = 150 m`, i.e. `(H0/L_pycnocline)^2 ~ 5.11`.
All `t*_mix` values below therefore come out ~5x larger than an earlier,
now-superseded H0-based pass, but the sweep design, mechanism, and qualitative
conclusions are unchanged -- only the absolute numbers shifted. All figures and
numbers in this note use `--length-scale pycnocline`
(`analysis/mixing_timescale.py:mixing_timescale(..., length_scale="pycnocline")`).

## 2. Experiment design

Three 5x5 cartesian sweeps (`structure.c4` x `GLS.C2`), one per closure:

- `templates/gen_c2c4_sweep.yaml` -> `sweeps/gen_c2c4/` (GEN closure,
  `GLS.C1=1.0`, the baseline closure)
- `templates/k-omega_c2c4_sweep.yaml` -> `sweeps/k-omega_c2c4/` (k-omega
  closure, `GLS.C1=0.555`; coefficients in `configs/variants/k-omega.yaml`,
  following Warner et al. 2005 (doi:10.1357/002224003322005087) Table 1 /
  Umlauf & Burchard 2003: P=-1, M=0.5, N=-1, CMU0=0.5477, C3M=-0.6, C3P=1.0,
  SIGK=SIGP=2.0)
- `templates/k-e_c2c4_sweep.yaml` -> `sweeps/k-e_c2c4/` (k-epsilon closure,
  `GLS.C1=1.44`; coefficients in `configs/variants/k-e.yaml`: P=3.0, M=1.5,
  N=-1.0, C2=1.92, C3M=-0.63, C3P=1.0, SIGK=2, SIGP=1.3)

All three sweeps use:
- `GLS.C2 = {0.52, 0.833, 1.22, 1.6, 1.92}` -- these are the literature `C2`
  values for the k-kl, k-omega, gen, and k-epsilon closures respectively
  (1.6 fills the gen-k-epsilon gap), used identically across all three
  sweeps so the closures' `C2` sensitivity is directly comparable.
- `structure.c4 = {0.1, 0.3, 0.55, 0.8, 1.0} x GLS.C1` -- five fractions of
  each closure's own `C1`, so all three closures can be compared via `c4/C1`
  on a common axis (`c4` must stay `<= C1`, see section 4.3 of
  `mixing_timescale_analysis.md`).
- `structure.CD/str_a`, `initial.temp_dT`, `grid.H0` fixed at the
  `mixing_timescale.yaml` variant's values (same as `c4_fine_sweep.yaml`),
  so `tau_mix_theory` (the Carpenter et al. power-based estimate, which
  does **not** depend on GLS coefficients) is identical (9.21 d,
  `NTIMES=79650`) across every run in all three sweeps -- any *difference*
  in the empirical `t*_mix` is therefore attributable entirely to the
  closure coefficients, not to a differently sized problem.
- `ntimes_margin: 4.0` (vs 2.5-3.0 in the c4-only sweeps), anticipating
  that varying `C2` as well as `c4` might slow mixing further.

**Analysis tool added:** `analysis/plot_c2c4_phase.py` -- pcolormesh of
`t*_mix` over the (`c4`, `C2`) plane (optionally normalized to `c4/C1`, and
selectable `--length-scale {H0,pycnocline}`), marking non-converged
(`t*_mix = NaN`) cells with a red X.

**Infrastructure note:** running the 75 combined runs via
`xargs -P 4 tools/run_experiment.py` intermittently caused runs (about half
of k-omega's, none of GEN's, none of k-epsilon's) to fail immediately with a
spurious `Fortran runtime error: Bad real number in item 1 of list input` in
`inp_decode.f90`, on `.in` files that were verified byte-for-byte correct.
Rerunning the *same* `resolved_config.yaml` sequentially (no parallelism)
always succeeded. This looks like a race condition when several `romsS`
processes read a shared file at `roms-related/` (plausibly `varinfo.yaml`)
concurrently -- not a physical instability. **Takeaway: don't run more than
one romsS process at a time sharing the same `roms-related/` directory** (or
investigate/serialize access to `varinfo.yaml` if parallel throughput is
needed).

## 3. Results

### 3.1 Phase diagrams

![GEN, k-omega, k-epsilon c2xc4 phase diagrams](../figures/c2c4_phase_comparison_3closure_pycnocline.png)

*(Individual non-normalized versions: `figures/gen_c2c4_phase_pycnocline.png`,
`figures/kom_c2c4_phase_pycnocline.png`, `figures/ke_c2c4_phase_pycnocline.png`;
normalized versions with `_norm` suffix.)*

Numeric `t*_mix` grids, pycnocline length scale (rows = `C2`, columns = `c4/C1`):

**GEN** (`C1=1.0`, own default `C2=1.22`):
```
C2\c4/C1   0.1     0.3     0.55    0.8     1.0
0.52       2.28    2.38    NaN*    NaN*    NaN*
0.833      2.24    2.28    2.41    14.86   NaN*
1.22       2.21    2.23    2.28    2.39    3.41
1.6        2.19    2.21    2.23    2.28    2.35
1.92       2.18    2.19    2.21    2.24    2.27
```

**k-omega** (`C1=0.555`, own default `C2=0.833`):
```
C2\c4/C1   0.1     0.3     0.55    0.8     1.0
0.52       2.54    2.70    4.14    11.89   NaN*
0.833      2.42    2.46    2.55    2.81    3.97
1.22       2.36    2.38    2.41    2.46    2.52
1.6        2.33    2.34    2.36    2.38    2.41
1.92       2.32    2.32    2.34    2.35    2.37
```

**k-epsilon** (`C1=1.44`, own default `C2=1.92`, `CMU0=0.5477` -- see the
CMU0 fix note in section 3.5):
```
C2\c4/C1   0.1     0.3     0.55    0.8     1.0
0.52       2.50    8.63    NaN*    NaN*    NaN*
0.833      2.40    2.59    NaN*    NaN*    NaN*
1.22       2.35    2.40    2.92    NaN*    NaN*
1.6        2.33    2.36    2.43    3.30    8.91
1.92       2.31    2.33    2.38    2.50    3.48
```
(`*` = `phi_star` never dropped below the 5% threshold within the run --
see section 3.3.)

### 3.2 c2 and c4 are not separable, in any of the three closures

Previous work implicitly assumed `t*_mix ~= A(c4) * B(H0, zt)`, with `A(c4)`
a fixed function found once at the closure's default `C2`. **That assumption
breaks down once `C2` is allowed to vary, in all three closures**:

- **At high `C2` (near or above the closure's own default), `t*_mix` is
  essentially flat/independent of `c4`**, changing by only a few percent
  across the entire `c4 in (0, C1]` range. This is why the earlier `c4`-only
  sweeps (each run at its own closure's default `C2`) found relatively mild
  `c4`-sensitivity except very close to `c4=C1`.
- **At low `C2` (at or below the closure's own default), `t*_mix` grows
  sharply, non-linearly with `c4/C1`** -- an order of magnitude increase, or
  outright non-convergence, as `c4/C1` goes from 0.1 to 1.0, in all three
  closures.
- The interaction is qualitatively the **same shape** across all three
  closures once `c4` is normalized by each closure's own `C1` and the same
  absolute `C2` grid is used -- i.e. this is not a GEN-specific artifact of
  `P/M/N/SIGK/SIGP`, but a generic feature of the GLS psi-equation's
  production/destruction balance (`c4*Pd` vs `c2*eps`).

### 3.3 New finding: the *degree* of fragility is closure-dependent

With k-epsilon added, a genuine closure-dependence emerges in *how much*
danger zone each closure has, even though the qualitative shape is shared:

| Closure  | C1   | own default C2 | NaN cells (of 25) | t\*_mix at own C2, c4/C1=1.0 vs 0.1 |
|----------|------|-----------------|--------------------|--------------------------------------|
| k-omega  | 0.555| 0.833           | 1                  | 3.97 vs 2.42 (+64%)                  |
| GEN      | 1.0  | 1.22            | 3                  | 3.41 vs 2.21 (+54%)                  |
| k-epsilon| 1.44 | 1.92            | 8                  | 3.48 vs 2.31 (+51%)                  |

Two consistent patterns:

1. **Robustness at each closure's *own* literature `C2` is similar across
   all three closures**: none collapse even at `c4 = C1`, and the "spike" at
   `c4/C1=1.0` is a comparable +51-64% increase in `t*_mix` in all three.
   So if each closure is run at its textbook default `C2`, the structure-mixing
   parametrization behaves reasonably (though not negligibly-sensitive) all
   the way up to `c4 = C1`, regardless of which closure is chosen.
2. **k-epsilon is markedly more fragile than GEN or k-omega whenever `C2`
   is pushed below its own default** -- 8/25 grid cells collapse to the
   numerical floor (vs 3/25 for GEN, 1/25 for k-omega), and collapse sets in
   at moderate `c4/C1` (0.55) even at `C2=1.22`, a value that is GEN's own
   *safe* default. k-omega is the most forgiving: only the single most
   extreme corner (`C2=0.52`, `c4/C1=1.0`) fails to converge. This ranks the
   three closures by robustness to the c2/c4 interaction as
   **k-omega > GEN > k-epsilon**, likely reflecting k-epsilon's higher `C1`
   and `C2` and different `P/M/N` exponents (P=3.0, M=1.5 vs GEN's P=M=0,
   k-omega's P=-1, M=0.5), which shift where the `c4*Pd` production term is
   able to overwhelm the `c2*eps` destruction term at a given absolute `C2`.
3. **Practical corollary:** because the GLS literature default `C2` differs
   by closure (k-kl lowest, k-epsilon highest), *the same absolute `C2` value
   represents a very different risk level depending on which closure it's
   paired with* -- e.g. `C2=1.22` (GEN's own safe default) is already deep in
   k-epsilon's danger zone (3 collapsed cells at that row).

### 3.4 Mechanism: closure collapse to the numerical floor, not just slower mixing

Every `t*_mix = NaN` cell was checked directly against the raw model fields
(`tke`, `AKt`) at the end of the run, in all three closures:

| Case | max(tke) | max(AKt) |
|---|---|---|
| GEN, C2=0.52, c4/C1=0.55 (NaN) | 1.0e-8 (= `GLS_Kmin` floor) | ~1e-6 m2/s (background) |
| GEN, C2=0.833, c4/C1=1.0 (NaN) | 1.0e-8 (floor) | ~1e-6 m2/s |
| k-omega, C2=0.52, c4/C1=1.0 (NaN) | 1.0e-8 (floor) | ~1e-6 m2/s |
| k-epsilon, C2=0.52, c4/C1=0.55 (NaN) | 1.0e-8 (floor) | ~1e-6 m2/s |

In every non-converged case, `tke` is pinned at its numerical floor
(`GLS_Kmin`) and `AKt` never rises above its molecular-background value --
i.e. the flow **never mixes at all**; `phi_star` stays flat at 1.0 for the
entire run (see `figures/gen_c2low_c4series_pycnocline.png`,
`figures/gen_c2c4_anomaly_cases_pycnocline.png`). This is the same
"turbulence-collapses-to-floor" failure mode that
`mixing_timescale_analysis.md` section 4.3 documented for `c4 > C1`
(assumed to be a hard, `C1`-only boundary) -- **this sweep shows the same
collapse can happen well inside `c4 <= C1` if `C2` is low enough, and the
threshold for "low enough" is itself closure-dependent (section 3.3).**

The "spike" cells that *do* converge but slowly (e.g. GEN `C2=0.833,
c4/C1=0.8`: `t*_mix=14.86`; k-epsilon `C2=1.6, c4/C1=1.0`: `t*_mix=8.91`)
are genuinely still mixing -- `tke` and `AKt` reach normal-looking peak
values, just much more slowly, sitting right at the edge of the collapse
boundary rather than past it.

### 3.5 Root cause: a dissipation runaway in the psi-equation, confirmed via a derived epsilon field

Section 3.4 shows *that* the closure floors, but not *why*. The working
hypothesis was that the TKE dissipation rate `epsilon` becomes anomalously
large and dissipates the TKE faster than it can be produced. ROMS does not
output `epsilon` directly -- it outputs the two GLS prognostic fields `tke`
(=`k`) and `gls` (=the generic variable `psi`) -- but `epsilon` can be
recovered from them for *any* GLS closure via the standard relation (Umlauf
& Burchard 2003):

```
psi = CMU0^P * k^M * l^N            (defines the turbulence length scale l)
epsilon = CMU0^3 * k^(3/2) / l       (definition of epsilon)
=> epsilon = CMU0^(3+P/N) * k^(3/2+M/N) * psi^(-1/N)
```

This is now implemented as `utils.utils.compute_epsilon(ds, params)` and
reduces, per closure, exactly as anticipated:
- **k-epsilon** (P=3, M=1.5, N=-1): `epsilon = psi` exactly -- ROMS' `gls`
  field *is* epsilon for this closure.
- **GEN** (P=2, M=1, N=-2/3): `epsilon = psi^(3/2)`.
- **k-omega** (P=-1, M=0.5, N=-1): `epsilon = CMU0^4 * k * psi` (consistent
  with the standard `epsilon = Cmu*k*omega` relation, since
  `CMU0^4 ~= 0.090 ~= Cmu = 0.09` once `psi` is identified with `omega`).

Applying this to the very first hours of a collapsing run (`C2=0.52`, one
representative "danger-zone" `c4` per closure) versus a normal run at the
same `C2` confirms the hypothesis, **with an important refinement**:
`epsilon` does not need to become *astronomically* large in an absolute
sense -- it only needs to exceed the (roughly fixed, structure-drag-set)
production rate `Pd` by a small, sustained margin for `dk/dt = P + Pd -
epsilon` to be persistently negative, which then **snowballs** because of
the `psi/k` factor multiplying the net production/destruction imbalance in
the psi-equation itself:

```
d(psi)/dt = ... + (psi/k) * (c1*P + c4*Pd + c3*G - c2*eps) + diffusion
```

As `k` falls, `psi/k` grows, amplifying whatever imbalance exists and
driving `psi` (hence `epsilon`) up even faster -- a positive-feedback
runaway, not a one-off large dissipation event. Concretely (hourly output):

![Early-transient tke/epsilon collapse mechanism, normal vs collapsed, all three closures](../figures/epsilon_collapse_mechanism.png)

*(top row: `max(tke)` vs hours since start; bottom row: `epsilon` at the
depth of peak `tke`, with `max(Pd)` -- the roughly-fixed structure-drag
production rate -- shown as a dotted reference line; columns are GEN,
k-omega, k-epsilon, each at `C2=0.52`.)*

- **GEN** (`C2=0.52`, `c4/C1=0.55`, collapsed): `epsilon` at the depth of
  peak `tke` is `2.6-2.8e-6 m2/s3` for hours 1-6, only ~6-9% above
  `max(Pd)=2.65e-6 m2/s3` -- but sustained, this pushes `tke` down 5 orders
  of magnitude (`2.5e-3 -> 1.0e-8`) within 6 hours, at which point both
  `tke` and `psi` hit their floors (`GLS_Kmin`/`GLS_Pmin`) and the
  "epsilon" computed from the floored fields (`~1.8e-11`) becomes a
  numerical artifact, not a physical dissipation rate. The normal run
  (`c4/C1=0.1`) instead settles `epsilon ~= 1.55e-6 < Pd`, i.e.
  `dk/dt > 0`, and `tke` grows to and holds a steady, physically-sized
  value (`~3.4e-3`).
- **k-omega** (`C2=0.52`, `c4/C1=1.0`, collapsed): the runaway is far more
  dramatic and directly visible in `psi` (=`omega`) itself, which explodes
  from `0.14` to `964` over just 5 hours while `tke` crashes from `2.2e-4`
  to `3.2e-8` -- this is the classic "`k->0`, `omega->infinity`" numerical
  instability well documented for k-omega-type closures when destruction is
  too weak relative to production. The `psi/k` amplification factor itself
  was checked directly and grows from `O(1)` to `O(10^10)` over the same 4
  hours, confirming the feedback mechanism above.
- **k-epsilon** (`C2=0.52`, `c4/C1=0.55`, collapsed): the most extreme case
  -- `tke` and `epsilon` (`=psi` for this closure) are pinned at their
  floors from the very *first* output step onward; the closure never gets
  a chance to build up any real turbulence at all, consistent with
  k-epsilon's markedly larger "danger zone" found in section 3.3.

**Bottom line: yes, `phi_star` stays at 1.0 because `epsilon` (dissipation)
outruns production early in the run, in all three closures** -- but the
underlying trigger is not "epsilon becomes very large" in isolation; it's
that low `C2` weakens the psi-destruction term relative to the c4-driven
production term just enough that the resulting small, sustained
`epsilon > production` imbalance gets amplified through the `psi/k` factor
into an unstoppable runaway, which different closures resolve differently
(a comparatively gentle exponential TKE decay for GEN, an `omega`-blowup for
k-omega, or an instantaneous floor-lock for k-epsilon) but which all end at
the same numerical floor state.

**Bug found and fixed:** while implementing `compute_epsilon`,
`configs/variants/k-e.yaml` was found to be missing a `GLS.CMU0` override,
so all k-epsilon runs originally in this analysis (and presumably in
`mixing_timescale_analysis.md`'s earlier k-epsilon work) had silently used
`CMU0=0.527` (the GEN literature value, inherited from
`configs/baseline.yaml`) instead of k-epsilon's own literature value
`CMU0=0.5477` (`roms/ROMS/External/roms_basin.in:2027`). This was fixed in
`configs/variants/k-e.yaml` and `configs/variants/k-e_mixing_timescale.yaml`,
and all 25 k-epsilon `c2c4` runs were regenerated and rerun (all
`returncode=0`). As anticipated, the ~4% CMU0 difference did not change any
qualitative conclusion -- it does not appear in any of the `epsilon`
exponents (which reduce to 0 for k-epsilon, section 3.5), but does enter
`AKt`/`AKv` via the stability functions, producing small (a few percent)
shifts in `t*_mix` (see the updated numeric grid in section 3.1) and leaving
the NaN/danger-zone pattern (8/25 cells) unchanged. The two other
`configs/variants/k-e*.yaml` files used elsewhere in this repo
(`k-e_nostr.yaml`, `test_STRUCTURE_PRODUCTION.yaml`) have the same missing
`CMU0` override but were left untouched as out of scope for this note.

## 4. Practical implications

1. **`c4` calibration must be closure- *and* `C2`-aware.** A `c4` value
   found safe/efficient at one `C2` (e.g. a closure's literature default)
   is not guaranteed to be safe if `C2` is changed, or if a different
   closure with a different default `C2` is used with the same `c4`.
2. **Updated rule of thumb (pycnocline scale):** keep `c4 <~ 0.3*C1` if `C2`
   is at or below the closure's own literature default; above that default,
   `c4` can range up to `C1` with only a mild (<60%) effect on mixing
   time. Do not assume this holds "for free" when swapping closures --
   k-epsilon's safe region is noticeably smaller in absolute `C2` terms
   than GEN's or k-omega's.
3. **This c2/c4 interaction is a genuine, closure-generic GLS-physics
   effect** -- the same qualitative shape (flat at high `C2`, spike/collapse
   at low `C2`) appears in GEN, k-omega, and k-epsilon despite very
   different `P/M/N`/Schmidt-number settings, strongly suggesting the
   `c4*Pd` vs `c2*eps` balance in the psi-equation is the controlling
   mechanism, largely independent of the other GLS coefficients. A natural
   follow-up: re-parametrize the sweep by `c4/c2` directly (instead of `c4`
   and `C2` separately) to test whether the three closures' phase diagrams
   collapse onto a single 1D curve.
4. **But the *quantitative* danger threshold is closure-dependent**
   (k-epsilon most fragile, k-omega most robust; section 3.3) -- so a
   universal "safe `c4/C1`" number cannot be quoted independent of both the
   closure and its `C2` setting.

## 5. Artifacts

- Configs: `configs/variants/k-omega.yaml`,
  `configs/variants/k-omega_mixing_timescale.yaml`,
  `configs/variants/k-e_mixing_timescale.yaml` (reuses existing
  `configs/variants/k-e.yaml` coefficients)
- Sweep definitions: `templates/gen_c2c4_sweep.yaml`,
  `templates/k-omega_c2c4_sweep.yaml`, `templates/k-e_c2c4_sweep.yaml`
- Sweep outputs: `sweeps/gen_c2c4/`, `sweeps/k-omega_c2c4/`,
  `sweeps/k-e_c2c4/` (75 runs total, all `returncode=0`)
- Analysis tool: `analysis/plot_c2c4_phase.py` (supports
  `--length-scale {H0,pycnocline}`)
- New utility: `utils/utils.py:compute_epsilon(ds, params)` -- derives the
  TKE dissipation rate `epsilon(z,t)` from the `tke`/`gls` fields for any
  GLS closure (section 3.5)
- Figures (pycnocline length scale): `figures/gen_c2c4_phase_pycnocline.png`,
  `figures/gen_c2c4_phase_norm_pycnocline.png`,
  `figures/kom_c2c4_phase_pycnocline.png`,
  `figures/kom_c2c4_phase_norm_pycnocline.png`,
  `figures/ke_c2c4_phase_pycnocline.png`,
  `figures/ke_c2c4_phase_norm_pycnocline.png`,
  `figures/c2c4_phase_comparison_pycnocline.png` (GEN vs k-omega, 2-panel),
  `figures/c2c4_phase_comparison_3closure_pycnocline.png` (all three,
  3-panel), `figures/gen_c2c4_anomaly_cases_pycnocline.png`,
  `figures/gen_c2low_c4series_pycnocline.png`,
  `figures/epsilon_collapse_mechanism.png` (section 3.5: `tke`/`epsilon`
  early-transient time series, normal vs collapsed, all three closures)
- Superseded H0-based figures (kept for reference only, do not use for
  quantitative conclusions): `figures/gen_c2c4_phase.png`,
  `figures/gen_c2c4_phase_norm.png`, `figures/kom_c2c4_phase.png`,
  `figures/kom_c2c4_phase_norm.png`, `figures/c2c4_phase_comparison.png`,
  `figures/gen_c2c4_anomaly_cases.png`, `figures/gen_c2low_c4series.png`
