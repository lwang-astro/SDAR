# Reference note: sd-energy gauge (`etot_sd_ref_`) corruption — investigation record and implemented fix

Status: **reference note** — investigation record (2026-10-19→22), the
implemented fix mechanism, and its validation; open follow-up items live in
`docs/etot-sd-ref-tracking-corruption-plan.md`. Opened 2026-10-20 from the
PeTar restart-grind investigation (parent:
`PeTar/doc/restart-slowdown-cold-state-plan.md`, anomaly class B); root
causes fixed 2026-10-21; production-validated 2026-10-22 (N=1000 cluster
A/B, §Production validation). Headline: the original "tracking corruption"
hypothesis is disproven — the `etot_sd_ref_` bookkeeping is exact; the
residual is genuine committed integration error of one under-resolved
pericenter step, enabled by the `fix_step_option=later` policy freezing ds
(Findings 3–5).

## Symptom

A 4-body system (knife-edge hierarchical trio + one perturber) with the inner
binary pumped to e = 0.9955 (pericenter 1.9×10⁻⁹) by the encounter reaches a
pericenter plunge; during the plunge the AR integration aborts with
`Error! adaptive ds below time resolution` after the SD-channel Hamiltonian
turns NaN (`H = -nan` in the final record). The NaN is structural: halving ds
cannot cure it, so the NaN branch (unconditional halving) runs ds into the
time-resolution floor and aborts.

## Measured facts (gdb, abort frame + conditional watchpoint, 2026-10-20)

Abort-time members of `TimeTransformedSymplecticIntegrator` (symbolized replay
of the committed pre-collapse dump):

```
PHYS: ekin=3.9578e8   epot=-9.3595e8   etot_ref=-5.3920e8
      H_phys args: num=+9.3498e8  den=+9.3595e8        → finite ✓
SD:   ekin_sd=9433.55   epot_sd=-20379.40   etot_sd_ref=+9454.85   (!)
      true etot_sd = ekin_sd + epot_sd = -10945.8
      H_sd args: num = 9433.55 - 9454.85 = -21.30 < 0  → H_sd = ln(negative) = NaN
gt_kick_inv = gt_drift_inv = 20379.40 (= |epot_sd|)
tree: [0] inner 23.66+19.04, semi=4.18e-7, ecc=0.9955, κ=49575 (at cap)
      [1] root  42.71+0.084, semi=0.0182, ecc=0.9995, κ=1
```

- Both log denominators (`-epot`, `-epot_sd`) stay positive throughout —
  "hyperbolicity flips epot sign" is a wrong reading (epot stays negative; the
  log argument that flips sign is the *numerator*). The residual sits entirely
  in the reference vs. actual sd-energy difference.
- The residual R = (ekin_sd+epot_sd) − etot_sd_ref_ is violated by ~2×10⁴ ≈
  |epot_sd| — but see Findings: R is the integration-error gauge
  (`getDESlowDown`), not a bookkeeping invariant, and the violation was
  injected by a single committed step, not by tracking drift.
- Group-formation time R = 0 exactly (`Etot_SD_ref = -110.6`, consistent).
- The conditional watchpoint (`watch etot_sd_ref_ if etot_sd_ref_ > 0`)
  caught a write at **`src/AR/symplectic_integrator.h:1697`** flipping the
  reference −946.31 → +988.44 at t = 0.0012858. **2026-10-20 correction: that
  write is benign** — it is the legitimate κ-engagement jump (inner binding
  −1935 divided by κ as slowdown engages), correctly bridged by `de_sd`; at
  that moment R was already −20407. A sign-flip condition cannot identify the
  corrupting write when the residual predates the sign crossing (methodology
  lesson recorded in SDAR lessons #9 update).

## Findings (2026-10-20, measured — supersedes "Suspected root")

Method: full writer-attribution log (gdb breakpoint at the `de_sd` write,
`symplectic_integrator.h:1697`) — 701 writes from group formation to abort —
plus per-step controller tracing in the fatal window and a `later`-flip trace.
Tooling kept in `~/sdar_run/etot-research/` (`attrib.gdb`/`attrib.log`,
`stepchk.gdb`/`stepchk.log`, `setlater.gdb`, `analyze.py`, `analyze2.py`).

1. **The tracking is exact; the flip write is benign.** All 701 `de_sd` writes
   are force-sync (`force=1, rebuilt=0, intr=0`): both `ekin_sd_` (`calcEKin`)
   and `epot_sd_` (`calcAccPotAndGTKickInv`) are recomputed from the same live
   tree with identical per-node nested κ — the root-vs-node κ-mismatch
   hypothesis is disproven. Measured `de_sd ≡ S_post − S_pre` (S =
   ekin_sd_+epot_sd) to ≤0.01 absolute over ~10⁶–10⁷ value swings (≤1e-9
   relative). Consequently every write preserves the residual
   `R = S − etot_sd_ref_` exactly (R_pre = R_post at all 701 writes). The
   kick-writer (`etot_sd_ref_ += _dt*de`) is exonerated: Σ over the whole run
   = −5.7e-4.
2. **`etot_sd_ref_ ≈ ekin_sd_ + epot_sd_` is not a bookkeeping invariant** — it
   is the sd-energy integration-error gauge (`getDESlowDown`, printed as
   dE_cum). R ≈ 0 only while integration error is small. R-source decomposition
   (refined 2026-10-21): (i) κ CHANGES do not themselves inject R (de_sd
   bridges the rescale exactly) — they act as enabler (κ collapse shrinks
   P_eff ~8000× → frozen ds under-resolves the pericenter → error injection)
   and detonator (κ recovery shrinks |epot_sd| below the old residual);
   (ii) with κ FIXED and no perturbers (closed system) R is bounded symplectic
   oscillation — no secular accumulation; (iii) with perturbers, an
   OPEN-SYSTEM BOOKKEEPING MISMATCH accumulates R independently of κ changes:
   kicks add FULL-weight work to the reference
   (`kickEtotAndGTDriftTreeIter`: `etot_sd_ref_ += dt·Σ m(v·a_pert)`) while S
   absorbs only the per-node κ-weighted part → dR/dt ≈ −W·(1−1/κ), linear in
   the physical perturbation work W (measured g4e: −5.7e-4 total, weak
   perturber; cluster-core W≈0.2–14/time-unit per lessons 2026-10-19 #6) —
   precisely the "rate-type" drift the EMA subtracts by design, i.e. fed
   while deliberately ignored. Note this mismatch is specific to the n≥3 tree
   path; the 2-body path keeps `etot_sd_ref_ = etot_ref_·κinv` (consistent
   scaling, no mismatch). Candidate root fix (separate change, own
   validation): κ-weight the kick-time reference increment to match the tree
   decomposition. Justification (2026-10-21): the slowdown system is NOT a
   representation change but the modified dynamical system actually
   integrated — internal forces scaled 1/κ, internal drift at v/κ
   (driftPosTreeIter), perturbation kicks at FULL strength in real time —
   and its conserved energy S obeys dS/dt|pert = κinv·W_internal + W_cm
   (CM motion unscaled; ∂S/∂v_i is the κ-weighted gradient, exactly vel_sd
   in kickEtotAndGTDriftTreeIter). **IMPLEMENTED + validated 2026-10-21
   (TTL tree path only)**: `kickEtotAndGTDriftTreeIter` now accumulates
   `de_sd += m·(vel_sd·a_pert)` and credits `etot_sd_ref_ += dt·de_sd`
   (full-weight `de` still goes to `etot_ref_`); diagnostic accumulators
   `de_kick_full_cum_`/`de_kick_sd_cum_` (never serialized) record both
   credits. E1 identity audit (rp=0.15, one-off -g -O0 solver
   `~/sdar_run/etot-research/petar-merger-dbg`, 5881 sync samples):
   quiescent phase t=9.2e-6 → R = −8e-8 while the credit mismatch
   (de_full−de_sd) = 2.6e-5 (**~300× drift suppression**, R pinned at the
   truncation floor as predicted); t=1.8e-5 → R = −3.3e-7 vs mismatch 9.1e-5
   (~270×). Mid-dance phases remain dominated by committed integration
   error (controller domain; mismatch negligible there). E2 behavioral A/B
   (release binaries, rp=0.15/0.35/0.4 LogH + 0.4 BTLogH): 8/8 complete,
   rebases 0/0, 1/0, 0/0, 2/1 — never worse; rp=0.35's accumulated-flavor
   rebase gone, rp=0.4-BTLogH's degenerate rebase 2→1. Chaotic
   single-snapshot warnings cannot test the ÷κ rate-signal prediction
   statistically (that is E1's quiet window).
   **Residual localization + criterion root-fix (2026-10-21, completes the
   EMA-retirement path)**: per-step sampling on the fix build (`loc15.log`,
   6653 steps, rp=0.15) split the remaining persistent signal: κ-CHANGE
   steps (72%): |dH| med 2.9e-8, same-sign 0.97 — the H-gauge
   re-normalization artifact (H = ln(1+R/epot_sd); the sync bridge
   preserves R exactly but epot_sd is κ-rescaled → per-sync shift
   ~R·Δ(1/epot_sd)); κ-CONST steps: larger |dH| (1.4e-7) but ZERO-MEAN
   (signed median −9e-16) = ordinary truncation; correlation with kick-work
   increments ≈ 0 (C1-credit sampling excluded). Hence the root fix is a
   CRITERION change, not another bookkeeping term: the ds controller now
   uses the gauge-invariant per-step sd-energy error `control_error_abs =
   |ΔR|/max(|ref|,|epot_sd|)` (reduce/growth/recovery/ratio; robust
   denominator for dissolving groups) — the H signal and the EMA drift
   estimate remain PRINTED DIAGNOSTICS ONLY (EMA retired from control; its
   two jobs are root-fixed: C1 by the κ-weighted credit, κ re-normalization
   by construction of ΔR). Validation (`crit_*.log`, all four scenarios,
   EMA out of control): 4/4 complete, ZERO rebases (first time), warnings
   11/40/1/39 at the deliberately tiny ar-max-nstep=100 (vs 212+ with the
   old criterion and EMA off; vs 0–4 with EMA compensating) — the remaining
   warnings reflect real, ds-reducible truncation the controller now
   resolves honestly. **EMA fully deleted (2026-10-21, user decision)**: after
   confirming zero in-repo consumers of the printed fields
   (`|Int_err_eff/E|`, `|H_drift_rate|`), the members
   (`h_drift_rate_ema_`/`h_drift_nsample_`), the update block, the
   eff-subtraction computation, and the two printMessage columns were removed
   entirely (history notes left at the former declaration site and in the
   amended 2026-10-19 lessons entry); post-deletion validation identical
   (4/4 complete, 0 rebases, warnings 11/40/1/39 — bit-for-bit the crit-run
   outcome, as expected for a print-only removal). Full matrix (T1–T7 +
   functional + SDAR samples +
   2026-10-19 grind replay) still required before commit; LogH tree-path
   kickEtot credit remains a follow-up.
   Scope note: absolute W in
   these 4-body scenarios is small, so behavioral deltas are modest — the
   fix concentrates value in strong-perturbation long-lived environments
   (cluster cores, W·τ vs |epot_sd|); production-scale validation belongs to
   the full matrix. Follow-ups: LogH tree-path kickEtot needs the analogous
   tree walk for vel_sd (noted in code; PeTar hard uses TTL); decide
   keep-vs-strip of the audit accumulators at commit time; canonical-family
   T1–T7 + functional smoke.
   Empirical fingerprint of the mismatch:
   the 2026-10-19 measurements found per-step |ΔH| ≈ W·dt/|E| with the FULL
   W across both formation and equilibrium phases (consistent credit would
   read κinv·W). Reconciliation with the 2026-10-19/20 "changeover-only"
   finding (2026-10-21): the mismatch is generic but its per-step signal is
   set by the ORBIT-AVERAGED work W̄ — measured cluster rate ~0.1–1.2/unit →
   ~6.5e-9 per step at equilibrium dt, just below ar-max-error (invisible);
   the changeover k(r) band converts the smooth mean flow into an
   orbit-frequency MODULATED flow (members sweep the band at orbital speed
   each P_in), amplifying the instantaneous rate by orders of magnitude
   (standalone no-changeover reproduction reaches only 2e-7 vs PeTar-level
   3.6e-5, ~170×) — that modulation, not |W|, is what crossed the per-step
   threshold and ground the pre-EMA-fix controller (neutralization test: k=1
   everywhere = STRONGER force, grind GONE). Two stacked defects (refined
   2026-10-21, correcting the earlier "metric lacks the perturber potential
   term" phrasing): (a) the AR per-step H carries U_pert only IMPLICITLY
   through the work-accumulated reference (exact in continuous time; the H4
   audit layer by contrast has it explicitly: ekin+epot+epert−etot_ref) —
   the per-step |ΔH|≈W·dt/|E| signal is the DISCRETIZATION RESIDUAL of the
   kick's midpoint-velocity work sampling vs the exact state-recomputed
   energies; design-inherent, EMA-compensated (an explicit-U_pert metric was
   evaluated and rejected: changes etot_ref semantics and audit columns,
   and Upert jumps between calls have no work record); (b) the SD reference
   work recording uses the wrong weight — the genuine bug: the kick credits
   FULL-weight work to etot_sd_ref_ (same `de` as the physical channel)
   while S absorbs per-node κ-weighted work, even though inv_nest_sd/vel_sd
   are already computed in the same loop — fix shape: accumulate
   `de_sd += m·inv_nest_sd·(v·a_pert)` and credit that to etot_sd_ref_.
   Generic-existence of (b) ⟢ slow R fuel; changeover/κ-transients ⟢
   igniter — both needed for a controller-grade event.
3. **The residual is injected between writes by committed steps.** Through the
   plunge R oscillates ±(10²–3×10³) reversibly (pericenter passages, symplectic
   oscillation; e.g. dR −2967 then +3241, self-recovering). The fatal event:
   **one committed step** (step_count=524, t=0.00128094583533, dt=7.1e-12)
   injects dR = −2.75e4 (0.4% of |S|≈6.5e6); only +4.9e3 recovers; R freezes
   at ≈ −2.05e4 for the rest of the run.
4. **Why the controller did not catch it: `fix_step_option=later` froze ds for
   the group's entire life.** Set at group formation (424 re-assertions):
   `information.h:516` `else if (bin_root.stab<1.0) → later` (the n≥3 branch —
   the 2-body and post-merger grants never fired here: 424/424 hits on :516,
   0 on :514/:2936), with measured stab=1.35e-4 despite root ecc=0.9995.
   Formula reconstruction confirms the value:
   qst = 1.52·(sqrt(t_crit/P)/(1−e))^(1/6)·(fac·g)^(1/3) ≈ 13
   (t_crit=10000·P → the ^(1/6) factor is 7.68 at e=0.9995), q_apo =
   apo/a_in ≈ 9×10⁴ → 1.4e-4 ✓. Two deeper points: (i) the Kepler fit at
   formation is self-consistent — the third body sits at ~apo (r/apo=0.995,
   ecca≈π) moving at 1.6% of escape speed (deeply bound, NOT a hyperbolic
   misfit; vis-viva closes to <1%); (ii) even the peri-evaluated criterion
   would judge stable (q_peri ≈ 21 vs qst ≈ 13 → stab_peri ≈ 0.6 < 1; standard
   MA also says stable, a_crit ≈ 3.7·a_in vs 21·a_in). So the failure is not
   the Kepler estimation and not simply the apo/peri naming
   (`peri_out = semi*(ecc+1)` computes apo despite the name — real bug, but
   not decisive here): it is a **predicate mismatch** — `later` needs
   "resolution requirement adiabatic over the interval" (P_eff=P·κ constant;
   here pert_out swings (apo/peri)³≈6e10 per outer orbit → κ collapses 4
   decades → P_eff shrinks ~8000×), while `stab` answers "long-term survival
   against exchange/ionization" (a legitimately deep hierarchy at 22·a_in).
   MA-stable ≠ resolution-adiabatic. (Also suspicious, independently: the
   criterion's net e_out dependence ∝ (1−e)^(−1/6)/(1+e) makes e=0.9995
   ~4× "more stable" than a circular orbit of equal a — physically inverted,
   and these criteria are only calibrated to e≲0.9.)
   `later` disables BOTH error-reduce paths (only `none` reduces; `always`
   skips by design; `later` falls through silently) — per-step |H−H_bk| was
   1e-6–1e-5, i.e. 100–1000× over `ar-max-error`=1e-8, at every plunge step.
   ds re-estimation never rescued either: the quiescence gate (×3
   reproducibility) defers on the churning candidates; the cost valve
   requires `g_func_on` (off: ar-g-func=0); the safety valve
   (`n_live < 32/ds_scale`) never fires because gt explodes at pericenter
   making n_live huge. ds stayed 0.0087518 from formation to abort (step-100
   Warning already shows `Fix: later`). Note also: stab was never re-evaluated
   after formation (identical 1.352344e-4 across all 424 grants —
   `stableCheckIter` only refreshes on tree rebuild, and the tree never
   rebuilt, all DESD writes rebuilt=0).
5. **NaN chain (late detonation).** After the plunge κ[0] ramps back toward
   the cap (46057→49575) over ~8e-6 of time; |epot_sd| = |E_in|/κ shrinks with
   it; at t=0.0012895 |epot_sd| = 20379 finally drops below |R| = 20400 →
   num = ekin_sd−etot_sd_ref_ = R−epot_sd < 0 →
   `H = ln(num) − ln(−epot_sd)` = NaN (`calcH`, `sample/Hermite/ar_interaction.h`).
   The NaN branch halves ds, but the error is history-committed (not
   step-local), so halving cannot fix it → ds-floor abort (ds=1.27356e-13).
   Note also the log form compresses the signal: a 0.4%-of-|S| injection reads
   as |dH| ~ 1e-5.
6. **Archaeology — no bisect needed.** Not a tracking regression: the enabling
   policy is `d45b383` 2021-05-07 ("fix step for stable systems"); the freeze
   was hardened by the quiescence gate `598abd9` 2026-09-03 (and the
   hyperbolic-2body exception for `later` came later still, c3ce7b0
   2026-09-28 — n≥3 never got the analogous guard). `etot_sd_ref_` write paths
   last touched by `a197875` 2026-09-28 (snapshot ordering — working correctly
   here); the `881abf3` refactor is not implicated. Bit-identical aborts
   pre/post 2026-10-19 rate-fix: neither build ever reduces under `later`.
7. **Fix interplay warning.** Enabling `fix=none` alone may be insufficient in
   the rate-fixed build — the 2026-10-19 EMA subtraction can classify the
   smooth plunge dH (systematic, same-sign, dt-collapsing) as "rate" and erase
   it; its premise "genuine spikes ≫ rate·dt" fails for smooth
   under-resolution (though the drift rate itself explodes faster than the
   α=0.1 EMA tracks — undecided analytically). Fix validation must confirm the
   reduce path actually fires near the fatal step.
8. **BTLogH control experiment (2026-10-20, same binary, ar-g-func 0→1,
   tooling in `~/sdar_run/etot-research/btlogh/`)**: the reproducer binary is
   itself an `AR_G_FUNC_BTLOGH` build (`hard.hpp`: 0 = LogH, 1 = BTLogH
   runtime switch), so the gauge is a clean single-variable flip. Result:
   **BTLogH completes the window with no abort** (integrates to
   time_end=0.0078125; the trio breaks up physically, inner pair survives as
   a 2-body group; dE_SD/Etot_SD ≈ −1.9e-4, bounded) — at the cost of
   848,819 AR substeps. BUT the controller gap is NOT fixed, only masked:
   `Fix: later` is still granted and ds still frozen (ds=1.72, product-form
   scale); per-step |Int_err/E| ≈ 1.1e-7 sits ~10× over ar-max-error with
   the reduce path still disarmed — the product gauge simply retains enough
   margin that no O(epot)-scale one-sided injection occurs and the H numerator
   stays positive. Interpretation: (i) BTLogH is the right method for
   knife-edge hierarchical multiples (consistent with Wang 2025, ApJ 978, 65
   — LogH accuracy limits on hierarchical triples), and `ar-g-func 1` is the
   practical mitigation for such systems; (ii) the A2 fix remains necessary
   and gauge-independent — it is the safety net for the default (LogH)
   configuration and the controller hardening for BTLogH itself.

## Deterministic reproducer (fast, ~seconds to the abort)

State: the committed pre-collapse dump of the trio (n=3) plus its parameter
files, in `~/sdar_run/clean4/`:

```
dump:  g4e.dump  (group state at the step-count warning, before the collapse)
pars:  g4e.par g4e.par.hard g4e.par.bse g4e.par.rand
replay: OMP_NUM_THREADS=1 setarch $(uname -m) -R \
  ~/bin/petar.mpi.omp.avx512.hard.debug -p g4e.par g4e.dump
# → aborts with the same ds=1.27356e-13 / H=-nan, bit-identically with the
#   pre-fix binary ~/sdar_run/triple-stability-fixcheck/bin/old/petar...hard.debug
```

Solver-level reproduction (writes the dump itself): the 4-body plain IC
(`~/sdar_run/clean4/ic2.raw`, perturber at 0.4 on +x) with

```
OMP_NUM_THREADS=1 petar(plain,std) -u 0 -t 0.0312 -o 0.0312 -b 4 -s 0.0156 \
  --r 0.6 --r-ratio 0.15 --r-search-group 0.05 --G 0.00449830997959438 \
  --ar-max-nstep 100 input      # small cap forces an early warning → dumps
```

Placement scan (perturber x-position → outcome, both controllers identical):
0.03/0.05 clean · 0.15/0.3 rate-type grind (cured by the 2026-10-19 fix) ·
**0.2/0.25 divergent (chaos; one side hits this bug)** · **0.35/0.4/0.45/0.5
this bug (both abort)** · 1.5 (member-band corner) · 5.0 both abort.

gdb tooling already written (in `~/sdar_run/clean4/`):
- `nanprobe2.gdb` — break abort; frame 1; print phys/SD members + H args + tree.
- `watchref.gdb` — rbreak AR integrateToTime; `watch *(double*)&etot_sd_ref_ if
  *(double*)$ref > 0`; bt at the flip. (Template for watching any member.)

Research tooling (2026-10-20, `~/sdar_run/etot-research/`):
- `attrib.gdb` + `analyze.py`/`analyze2.py` — every `de_sd` write: flags,
  pre/post S, de_sd, ref, R, per-node κ; writer attribution by consecutive-
  write differencing (kick-writer sum = Δref − Σde_sd per interval).
- `stepchk.gdb` — per-step controller view (H, H_bk, R, ds, dt, fix option)
  inside a time window — use to watch reduce-path behavior when testing fixes.
- `setlater.gdb` — trace `fix_step_option=later` flips (stab/semi/ecc + call
  sites).

## Fix (final design — implemented 2026-10-21, after four review rounds +
Reviewer pass)

Four coordinated changes in `src/AR/symplectic_integrator.h` +
`src/AR/information.h` (+ `sample/AR/ar.cxx` help text):

a. *Grant-layer fix (user-directed revision, 2026-10-21 — supersedes the
   2026-10-20 symmetric-decision revision)*: `information.h`
   `calcDsAndStepOption` no longer auto-grants `later` for n≥3 stable
   roots (the `bin_root.stab<1.0` branch removed; history note in code);
   `later` keeps its ORIGINAL semantics (freeze ds after the initial
   adjustments) and is granted only to genuine elliptic two-body groups
   (auto) and the count_mass==2 post-merger grant — exactly the
   populations where locking ds is intended. All integrator feedback
   gates are RESTORED to HEAD (`==FixStepOption::none`); the 2026-10-20
   `!=FixStepOption::always` arming was reverted with it. `ar.cxx` help
   text updated (auto = later for elliptic two-body, none otherwise).
   n≥3 groups are now adaptive from formation — the plunge-freeze
   channel (Finding 4) is closed at the grant layer instead of
   compensated in the integrator. Four-scenario A/B re-verified at each
   revert stage (11/40/1/39 warnings, fin=1, 0 rebases — identical to
   the symmetric revision, as expected: n=3 groups receive `none`
   directly instead of `later`-treated-as-`none`).
   b. *Degenerate-gauge repair (final design 2026-10-21)*: single trigger —
      ds-floor escalation: at the loop-top ds-floor check, if the log
      argument is degenerate (`ekin_sd_−etot_sd_ref_ <
      1e-6·(−epot_sd_)`, covering both hard violation num≤0 and H<−13.8) →
      rebase + resume at the epoch ds; a genuine pure-resolution wall
      (healthy argument) still aborts. The floor IS the exact "halving is
      futile" criterion. **Physical-error gate (user-directed,
      2026-10-21)**: every repair is preceded by
      `|R|/max(|Ekin|,|Epot|) > energy_error_rel_max → abortSdGaugeRepair`
      (message + full column print + DATADUMP + abort): discarding the
      bookkeeping residual R is acceptable only while the PHYSICAL energy
      error is within tolerance — the pre-fix natural cases (|R|/|E|~1e-4,
      4 orders over) now fail loudly instead of continuing. Rebases are
      accounted into `de_sd_change_cum_` and warned.
      *The ISNAN(H) NaN branch was removed entirely (2026-10-21)* after a
      four-stage review: (1) transient/structural classification
      (2026-10-20); (2) trigger tightened to ISNAN(H) — the |H−H_bk| form
      also fires in a self-healing corner; (3) structural repair gated by
      the physical error; (4) removal for all builds. Grounds: TTL — the
      trajectory never consults num; the controller (control_error_abs,
      finite) handles real error, and a real-error transient carries
      |ΔR|=O(|epot_sd|) so the controller's restore+reduce IS the same
      cure (injection-verified: control_error_abs=1.68 fired, the restore
      recovered the reference bit-exactly from the intact backup); LogH —
      num is the drift gauge (dt=ds/num), so a domain exit makes dt≤0 and
      the PRE-EXISTING negative-dt tripwire catches it at the same
      post-step point (transients: same cure; committed violations:
      bounded ~1e3-step halving grind down to the ds-floor repair above).
      H is print-only everywhere now. Empirical basis (gauge-margin scans,
      `~/sdar_run/etot-research/bookfix/gauge_{pair,trio}.log`, 28,217
      steps): healthy runs sit ≥4.6e6× from a flip (log argument ≈1
      throughout; the κ cap bounds the denominator collapse — measured
      κ=1.011e8 ≈ cap), so the deleted detector guarded an event with no
      measured approach path in these configs. Residual accepted risks:
      an inflated reference inflates the control denominator (mild
      under-sensitivity); the suppressor corner lost its guaranteed
      backstop; LogH routing validated at compile level only (runtime
      negative-dt→floor validation on a LogH sample remains open, item 7).
   c. *EMA NaN guard*: moot — EMA deleted 2026-10-21 (see the fix-1 record
      above).
   d. *Not done (separate decision, unchanged)*: `stableCheckIter`
      apo-vs-peri / e-out-dependence fixes — global group-management
      semantics.
   **Validation (2026-10-20, solver-level rp scan, merger-variant binaries,
   strict A/B by temporary revert + rebuild; artifacts in
   `~/sdar_run/etot-research/solver-fix/`, `run_sym2_*.log`)**:
   - pre-fix: rp=0.4 aborts **bit-identically** (ds=1.27356e-13, gt=20379.4,
     t=0.0012895, step 735);
   - intermediate reduce-only revision: 4/4 passes but rejected as asymmetric;
     its BTLogH branch later showed the NaN-cascade hole (fixed by 3b);
   - final symmetric + exact-trigger revision (run_sym5_*): **4/4 complete**
     (rp=0.35/0.4/0.5 LogH, rp=0.4 BTLogH; exit=0 each). Rebases: 1/0/0/2 —
     the two LogH formation-plunge passages are now resolved properly by
     halving (0 rebases, unlike the streak=3 revision's premature 9.4e5
     rebase); rp=0.35's single late rebase (6330→−11005, t=0.00139) is the
     accumulated-residual flavor (within-tolerance accumulation + κ-recovery
     scale shrinkage); rp=0.4 BTLogH's two small rebases are the
     degenerate-gauge flavor. ds visibly adapts both ways under `later`
     (0.00875 → 2.2e-3 → 1.4e-4 → 4.3e-6 at plunges). The 2026-10-19
     rate-fix blind spot (Finding 7) was confirmed live on one branch
     (|Int_err|=1.4e-3, |Int_err_eff|=0) — the degenerate-gauge floor repair
     is the backstop for exactly that; revisit the EMA subtraction cap if it
     shows up in production matrices.
   Residual watch items: repairs now fire only via the gated ds-floor
   trigger — watch their frequency in production matrices; the post-sync
   estimator-channel clamp `ds=min(ds,info.ds)` remains one-way for ALL
   options (growth flows through the feedback enlargement — direct
   estimator-channel growth is still absent by design). Open follow-ups
   (T7 wiring + full matrix, LogH tree-path kick credit, LogH runtime
   routing validation, hyperbolic-root ds estimator) live in
   `etot-sd-ref-tracking-corruption-plan.md`.

## Production validation — N=1000 cluster A/B (2026-10-22, user-directed pre-commit gate)

Strict A/B: merger family, SDAR-HEAD baseline vs fixed working tree (provenance
by fix-string content check, banners identical by configure-time baking); IC:
mcluster N=1000, 100% binaries (500), Plummer R_h=1 pc, Q=0.5, Kroupa default
range, seed 42, `petar.init -s merger`; run: `-u 1 -b 500 -t 30 -o 0.5 -s
0.0009765625` (find.dt 0.001953125 halved), OMP=8, setarch -R; artifacts
`~/sdar_run/bench-n1k-arfix/` (commands.log, ab_wallclock.log, full logs).

* **Performance**: current 30 Myr in **4m24s** vs baseline stalled at
  0.6 Myr in a 2h cap (exit 124; extrapolated ≥500 h → >750×).
* **Energy**: current final global |dE/E| = **3.4e-9**, hard-energy events
  1.0/Myr, 4 H4-level info dumps (one event); baseline grind regime
  per-step |dE/E| = 3.66, cumulative H gauge destroyed (20.9), 16
  pathological AR dumps in 0.6 Myr.
* **Stability**: current 61/61 outputs, 0 aborts; baseline death spiral
  11 min after start, reproduced 3× (ds collapse 2.1e-17; trigger group
  = n3 containing a κ=4.1e5 tight binary — Finding 4 morphology at
  cluster scale). The 4 current-arm dumps = one legal 30-body/10-group
  resonant cascade at t≈5.9 Myr, replay-verified (`replay_c11.log`):
  2.16M AR steps, dissolves into 13 ordinary pairs, 0 rebase / 0 gate —
  legal-expensive, not pathological.

## Cross-references

- Open follow-ups: `docs/etot-sd-ref-tracking-corruption-plan.md`.
- Parent investigation: `PeTar/doc/restart-slowdown-cold-state-plan.md`
  (anomaly classes a/b; the rate-type controller fix and its T7 validation).
- SDAR lessons-learned 2026-10-20 entries #8/#9 (minimal-system localization;
  NaN three-step localization method: abort-frame member dump → conditional
  watchpoint → writer backtrace) — #9(b) updated 2026-10-20 with the measured
  root cause (`later`-frozen ds + one-step error injection) and the
  writer-attribution methodology correction.
- Prior same-family fix: commit a197875 (2026-09-28), snapshot-ordering —
  confirmed working here (Finding 1); its scenario differs (leak per sync
  event vs one-step committed error injection).
- Enabling commits: d45b383 (2021-05-07, `later` for stable n≥3), 598abd9
  (2026-09-03, quiescence gate hardening the freeze), c3ce7b0 (2026-09-28,
  the 2-body hyperbolic exception that n≥3 never received).
