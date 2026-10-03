# Group transition-layer unification plan (form/break/merge as inverse pure-state maps)

Status: proposed plan (2026-10-03), follows the diagnosis in
`PeTar/.github/plans/hard-large-energy-se-mass-loss.md` (Follow-up section).
Owner surface: `src/Hermite/hermite_integrator.h` (`addGroups`, `breakGroups`,
`writeBackGroupMembers`, slowdown bookkeeping), with PeTar-side replay gates.

## 1. Problem statement (measured)

Every group form/break (and form→break→re-form chatter around the criterion
boundary) injects an algorithmic energy error of ~1e-3 of the interacting-pair
energy into the trajectory:

- Pal5-IMF production: three replayed cases inject 16.4 / 1.12 / 0.138
  (6e-4..3.4e-2 relative) — BH pairs at peri of unbound or marginally bound
  (e > 0.9998) encounters; ~40 warning dumps/hour across the 7 runs.
- In-memory palindromes: no-event floor 5.6e-6 (dt-max-power 8); each event
  adds ~1e-3-level irreversibility; the core integrator (Hermite+AR+slowdown
  without events) is machine-symmetric (5.7e-12).

None of this is physical energy change; `de_change_cum` must stay reserved for
stellar-evolution/interrupt physics.

## 2. Goal and non-goals

**Goal**: group membership transitions become per-event invertible maps so the
injection per event drops by >= 20x (to <= 5e-5 of pair energy), with:

- P1 (per-event invertibility): `form` is a canonical change of variables
  (members <-> c.m. + relative coordinates + tree topology); `break` is its
  exact inverse to round-off. No cross-event memory, no lookahead: whether a
  group ever re-forms is irrelevant.
- P2 (purity): every auxiliary quantity (slowdown kappa and its perturbation
  sampling, ds, block dt, acc/jerk, energy references, vcm record) is a pure
  function of the current state variables, computed by the same code on both
  sides of the transition.
- P3 (sync alignment): transitions execute only where both representations
  have completed whole steps (Hermite block boundary; AR finishes its current
  transformed step before readout — no mid-step truncation).

**Non-goals** (accepted residuals, to be quantified, not fixed here):

- hybrid-scheme switching truncation (singles-Hermite vs AR+slowdown are
  different discretizations; each switch keeps the local truncation of the
  entering scheme's first step);
- the Hermite block-step non-self-adjointness floor (measured dpos 1.6e-3
  coarse / 5.6e-6 at dt-max-power 8);
- the slowdown frozen-perturbation approximation (accuracy set by the kappa
  update cadence, a state-determined choice, not an asymmetry).

If after Phase 2 the remaining per-event injection is dominated by these
residuals, the plan stops and the residual is documented as intrinsic; any
further change to the warning criterion semantics is a separate user decision.

## 3. Transition write audit (field-by-field)

| # | Quantity | form side | break side | Class | Action |
|---|----------|-----------|------------|-------|--------|
| 1 | member pos/vel | read from singles | `shiftToOriginFrame` + `writeBackMemberAll` | pure inverse | keep |
| 2 | tree topology | `generateBinaryTree` | `getTwoBranchParticleIndexOriginFromBinaryTree` split | pure | keep |
| 3 | kappa_ref, period, timescale | derived from tree | discarded (not needed by inverse) | pure | keep |
| 4 | kappa update schedule (`time_update_`, pert_in/out sampled per period) | anchored at formation time | discarded | **history** | re-anchor to the global block-step grid so the sampling times are state-determined; recompute pert_in/out at each AR sync from current state instead of caching |
| 5 | `vcm_record` + one-sided epert correction | CM velocity stored at formation, correction applied forward | no inverse term | **history, one-sided** | drop the record: recompute the correction increment per update interval from current-state CM acceleration (pure), or apply the exact reverse term at break — choose per Phase 1 budget |
| 6 | ds + AR transformed-time phase | sequence restarts at formation | readout at AR sync | phase | Phase 1 gate: only redesign if budget shows it matters (earlier estimate: 8 orders below injection) |
| 7 | returning singles' acc/jerk/dt/predictor | discarded at formation | fresh force evaluation at break instant | evaluation-point mismatch | verify the re-entry path uses the SAME function as an ordinary single force-update landing on that boundary; unify into one code path |
| 8 | `Etot_SD_ref`, epert refs | re-derived from fresh tree | implicit | pure (bookkeeping) | keep pure; dE semantics restored (PeTar re-baseline already withdrawn) |
| 9 | CM block dt | from CM force at formation | — | pure | keep |
| 10 | r_break_crit | max member r_group | — | pure | keep |

## 4. Phases

### Phase 0 — revert + bench (done 2026-10-03)

- PeTar re-baseline withdrawn (semantics: `de_change_cum` is physical-only).
- Bench = three Pal5 replay cases (`/data/lwang/Pal5_IMF/*/debug/`) + in-memory
  palindrome harness (`--reverse-at`, dt-max-power 8) + Plummer N=20 battery.

### Phase 1 — injection decomposition (instrumentation only; done 2026-10-03)

Add compile-time-gated per-event logging that decomposes the bookkeeping delta
at each form/break into: (a) vcm/epert correction terms, (b) slowdown sampling
(kappa update interval), (c) force re-evaluation (singles re-entry), (d) AR
restart phase, (e) reference re-derivation. Run on the bench.

Instrumentation: per-event decomposition trace under the standing
`ADJUST_GROUP_DEBUG` gate, all hooks in `hermite_integrator.h` — line types
`Group transition:` (rewrite purity, dE_rewrite), `Group formation:`,
`Group break:`, `Group re-entry:`, fields in the repo print convention
(`Pert_In/Pert_Out`, `SD/SD(org)`, `step/step(tsyn)`, `Etot_SD_ref`,
`Member_index`, `key: value`) (2026-10-03: originally a separate
`GROUP_EVENT_INJECTION_DEBUG` flag with a `hermite.geid` sample target;
consolidated into `ADJUST_GROUP_DEBUG` so every standard debug build —
including the PeTar `*.hard.debug` family — carries the decomposition trace).
Trace: `GEID ETRUE|BREAK|FORM|REENTRY` on stderr, physical time, event id.

Deliverable: ranked injection budget deciding which audit items (#4, #5, #6,
#7) actually matter; expected dominant terms per the 2026-09-30 attribution are
#4/#5.

Bench result (Plummer N=20 Henon units, dt-max-power 8, r-group 0.05, t=1000:
8715 events; 4389 breaks — 17 with kappa 5.5e3..1.2e4 from one persistent
binary, rest kappa=1 chatter):

| term | audit | measured |
|------|-------|----------|
| (a) `d_etot_sd` (break, slowdown shutdown) | #4/#5 | 0 at kappa=1 chatter; 2.8e-2..4.8e-2 (~30-50% of \|Etot\|) at every kappa>1 break — **dominant** |
| (a) `d_epert` (perturbation ref) | #5 | 0 — isolated bench has no perturbers; needs Pal5 replay |
| (a) `de_kin` (vcm one-sided) | #5 | <=1.4e-18 (roundoff) |
| (b) kappa/period/anchor state | #4 | logged both sides (kappa to 1.2e4) |
| (c) re-entry force jump | #7 | dacc_rel mean 6.4e-3, max 0.18; pot jumps 0.06..18 — secondary, every break |
| (d) AR sync gap at transition | #6 | 0 everywhere (P3 already holds) |
| (e) form-side `de_sd` booked | #8 | mirrors (a): 0 (kappa=1) / fresh-tree value (kappa>1) — the one-sided form/break pair is the asymmetry |
| state-rewrite purity \|dE_true\| pre/post | P1 | <=4.4e-16 — coordinate rewrite is exactly invertible |

Palindrome (seed 6, window [0, 2.10], mirrored event pair 2 -> 4 confirmed):
quiet floor dvel_rt 5.1e-6 (matches 2026-09-30); +1 chatter pair 5.3e-5 (10x);
same flyby Hermite-only (r-group 1e-7) O(1) — hybrid scheme is 4 orders more
reversible than the singles-only path through the same encounter.

Gate decision: injection is NOT >=80% unfixable hybrid truncation — at
bound-binary events it is dominated by the one-sided slowdown bookkeeping
(#4/#5). Proceed to Phase 2 with #4 and #5 primary, #7 secondary, #6 no
action. Chatter (kappa=1) events carry zero bookkeeping terms; their residual
(hybrid truncation + #7 re-entry) falls under Section 2 accepted residuals
unless the #7 unification also trims it.

Pending: Pal5 replays with the same gate compiled PeTar-side (measures the
(a) epert term with real perturbers; the three replay cases remain the
production-facing budget check).

Pal5 replay results (2026-10-03; `petar.mpi.omp.avx512.bse.galpy.hard.debug`
rebuilt with the gate as `build/...hard.debug.geid` from the workspace SDAR;
replayed `Hard Energy` finals match the run-side dump records bit-for-bit:
dE -0.401199 / +16.4051 / +0.138423):

- All three cases are the kappa=1 chatter class (form -> break, case 3 with
  re-form 16 us later). EVERY bookkeeping term is zero production-side too:
  `d_etot_sd` 0, `d_epert` 0 (real perturbers present but the epert reference
  path never engages for these short-lived groups), `de_kin` <=1.6e-14,
  `de_sd_booked` 0 at formation, sync gap 0.
- The groups live in extreme perturbation fields: pert_out/|pert_in| = 72
  (case 1), 73 (case 2), 4380 (case 3); 49-142 AR substeps per episode.
- The only nonzero discontinuities are the #7 re-entry terms: djerk_rel
  0.17-1.89, per-particle pot jumps 0.1-21.9, block dt reset (~6e-8 stash ->
  fresh re-init).

Production-side Phase-1 verdict (completes the budget):

- **Chatter class (all Pal5 dump cases, kappa=1)**: bookkeeping injection is
  identically zero; the flagged dE originates in the hybrid-scheme switch —
  the AR episode integrating the pair under pert_out/|pert_in| ~ 10^2-10^4
  tidal fields plus the #7 re-entry force/dt discontinuity. Actionable item
  for this class: #7 (pure re-entry/block-dt initialization), not #4/#5.
- **Bound-binary class (kappa>1)**: #4/#5 bookkeeping dominates (SDAR bench:
  3e-2 per break, Section above).
- State-rewrite purity holds production-side (|d_e| <= 2.8e-14 incl. the
  16-us re-form chatter).

Phase 2 scope confirmed on both classes: #4/#5 for bound binaries, #7 for
chatter; #6 no action.

Fixed en route: `sample/Hermite` checkpoint writer indexed groups by slot
i < getNGroup() instead of the sorted active index — aborted on
`reserveMem(0)` and wrote wrong group configs after slot masking; now uses
`HermiteIntegrator::getGroupIndexSorted()`.

### Phase 2 — purity refactor (SDAR)

Fix designs and verification record (2026-10-03, execution session; measured
evidence first, decision after):

**Fix 1 — slowdown bookkeeping (#4/#5), targeted at the bound-binary class.**

Measured (Plummer t=1000 run, the 17 kappa~5e3..1.2e4 form/break cycles of the
persistent hard binary):

- At formation AR books `de_sd = etot_sd_ref - etot_ref = +E_bin*(1-1/kappa)`;
  at break the Hermite side books `etot - etot_sd = -E_bin*(1-1/kappa)`.
  Example cycle: FORM +3.5375e-2 / BREAK -3.5335e-2 — **the two one-shot
  bookings cancel pairwise**; the net per cycle is the true binary energy
  drift (~4e-5 there). The Phase-1 "3e-2 per break" is therefore the
  definitional 1/kappa-scaling term, not a net injection; per-cycle net
  injection is ~1e-3 of E_bin, same order as the chatter class.
- Audit #4's prescription is already satisfied by current HEAD: slowdown
  state (pert_in/out, timescale, kappa) is recomputed from live perturber
  state at EVERY `integrateToTime` entry (block-grid sync;
  `syncTreeSlowDownAndDs(true, true)` at symplectic_integrator.h:2618), plus
  on interrupts and orbit updates; formation happens block-aligned
  (d_sync_gap = 0 measured everywhere). `SlowDown::time_update_` /
  `setUpdateTime` / `increaseUpdateTimeOnePeriod` are dead code (no callers)
  — the "schedule anchored at formation" assumed by the audit table no
  longer exists.
- #5 components measured negligible: `de_kin` (vcm_record term) <= 1.4e-18
  everywhere; `d_epert` = 0 in both benches and all three Pal5 replays.
  `vcm_record` is load-bearing for the SE mass-loss correction path
  (hermite_integrator.h interrupt handling) — removing it would touch
  physical bookkeeping for zero measured benefit.

Decision: **no code change for #4/#5.** The only #4/#5-compliant change left
(more frequent kappa recompute, e.g. per AR substep) would alter trajectories
of event-free persistent-group runs and violate the zero-event bit-identity
gate. The one-shot +/-E_bin transient bookings stay (they cancel; physical
channel `energy_.de_cum` never sees them). Actionable trajectory-level work
concentrates on #7.

**Fix 2 — re-entry initialization (#7), targeted at the chatter class (all
Pal5 production dumps).**

Verification (2026-10-03, code-path audit + measurements; no code change
warranted):

- The re-entry path is already the unified one: `initialIntegration` evaluates
  the re-entering singles with the same `calcAccJerkNBList` used by ordinary
  updates, stores the fresh acc/jerk/pot into `ptcl`, assigns dt through
  `calcDt2ndList` -> `step.calcBlockDt2nd` (the same estimator an ordinary
  particle gets on its initialization step), clears `initial_step_flag`, and
  `updateTimeNextList` re-enters them on the block grid. `predictAll` predicts
  from `ptcl` (with the fresh forces), so the early `pred_[k] = ptcl[k]`
  template copy cannot inject a zero-force predictor.
- The audited asymmetries (#4 cadence, #7 split init paths) were already
  removed by the 2026-09 consolidation; the audit table described
  pre-consolidation behavior.
- Empirical state-purity on both benches and all three Pal5 replays:
  |dE_true| pre/post rewrite <= 2.8e-14 (P1 holds everywhere); transitions
  land block-aligned (d_sync_gap = 0).
- The measured REENTRY jumps (djerk_rel 0.17..1.89, dpot 0.1..22) are the
  physical field evolution between formation and break — the same values any
  evaluation of that state returns — not an initialization artifact. The
  dt2nd-vs-dt4th estimator-order difference at re-init is the standard init
  convention, state-determined and applied symmetrically under reversal.

**Phase-2 outcome (both fixes verified, zero code changes).** With #4, #5, #7
all confirmed already-satisfied or negligible under the zero-event
bit-identity gate, the remaining per-event irreversibility (palindrome
5.3e-5 vs floor 5.1e-6 per chatter pair; Pal5 flagged dE 0.138..16.4) is
attributable to:

1. audit #6 family — the AR adaptive ds controller is error-driven
   (history-dependent, not self-adjoint): the mirrored run takes a different
   ds sequence through the same 50-140-substep episode under
   pert_out/|pert_in| ~ 10^2-10^4. This is now the last standing audited
   suspect; `fix_step_option` has no runtime switch, so the decisive
   attribution experiment needs a diagnostic build forcing
   `FixStepOption::always`.
2. the Section-2 accepted residual: hybrid-scheme switching truncation.

Per the plan's stopping rule, Phase 2 stops here unless the user opts to (a)
run the fixed-ds attribution experiment, (b) pursue a self-adjoint ds
controller (major AR redesign), or (c) revisit the warning-criterion
semantics (`dE` vs `dE_change` trigger) as an independent decision.

**Attribution experiments (2026-10-03, option (a) executed).**

1. Fixed-ds palindrome (seed-6 chatter window, `-e 1e10` disables the
   energy-error ds adaptation): **bit-identical** to the adaptive run
   (dvel_rt 5.280e-05 in both) — the adaptive ds controller never engaged in
   the window. Audit #6 is ruled out for the chatter class.
2. Matched-control weak-encounter palindrome (dedicated 3-body IC:
   m=0.45+0.45 hyperbolic encounter, peri 0.02, e=1.067 — the Pal5 marginal
   class — plus a distant 0.1 spectator; single form/break pair at
   t=0.93/0.94; window [0,1] reversed):

   | path through the same encounter | dvel_rt |
   |---|---|
   | pure block-Hermite singles (r-group 1e-7) | 9.5e-3 |
   | grouped (form -> AR episode -> break) | **3.0e-4** (32x better) |

   On current HEAD the group transition layer is no longer an error source
   relative to the no-group alternative — it is part of the solution. The
   pre-consolidation measurement (no-event 5.6e-6 vs 2-event 1.5e-3, 270x)
   does not reproduce with a properly encounter-matched control.

**Final attribution**: the remaining per-event irreversibility and the Pal5
flagged offsets (dE 0.138..16.4, reproduced bit-for-bit in replay) are the
intrinsic local truncation of switching discretizations mid-encounter
(Section-2 accepted residual), amplified by the encounter. The stopping rule
applies; the remaining user decision is (c) the warning-criterion semantics.

**Post-attribution finding (2026-10-03): the switch POSITION is the dominant
lever, not the transition implementation.** Palindrome scan over
`--r-group` on the dedicated hyperbolic IC (peri 0.02, e=1.067):

| r-group | weak perturber (isolated) | strong perturber (production-like) |
|---|---|---|
| 0.05 (switch near peri) | 3.0e-4 | 6.9e-4 |
| 0.1 | 7.4e-5 | — |
| 0.4 (switch in smooth field) | **7.5e-7** (400x) | 2.4e-4 (3x) |

Early switching moves both seams into weak-field regions and lets AR carry
the whole encounter (460 substeps vs 55; the reversal point falls inside the
group lifetime, so no new transition fires at all). Under a strong tidal
field the longer frozen-perturbation episode claws most of it back. The
optimal radius is environment-dependent; in production it couples to r_in
(changeover) and group-count/performance costs. Candidate follow-up: r-group
scan on the three Pal5 replays (edit the data.par.hard copy) before any
production recommendation. Caveat: r-group also feeds r_neighbor into
`calcAcc0OffsetSq`, so "no-group" runs at different r-group values are not
the same trajectory.

**Regression bench (2026-10-03)**: `sample/test/group_transition_bench.sh`
replaces the Pal5-specific replays for criterion/transition work — eight
deterministic no-SE cases in four families (weak/strong hyperbolic x late/
early switch; eccentric bound binary "breathing" late/early incl. a two-pair
chatter window; compact kappa>1 binary + light intruder for the >=3-member
capture/release path, gated structurally since the exchange is chaotic).
Gates: run completion, event mirroring (pal == 2x ref, counted from the
standing ADJUST_GROUP_DEBUG lines of the plain sample binary), palindrome
round-trip errors dV_RT and dE_RT/|E0| (per-case, 2-10x headroom), rewrite
purity (|dE_rewrite| <= 1e-12, same flag), and — for the triple family —
structure gates (a >=3-member event occurred; formation SD > 1). `--scan`
adds an r-group sweep with a lever-direction gate (dV_RT(0.4) <=
dV_RT(0.05)). Output carries an inline column legend. Depends only on the
plain `hermite` build + bash/python3; production (non-debug) builds never
define the flag.

**Criterion-optimization validation experiments (2026-10-03).**

1. Perturbation series (single dominant spectator, m = 1 / 3 / 10 at
   shrinking distance; grouped rg=0.05 vs ungrouped rg=1e-7, palindrome
   dE_RT/|E0|): grouped wins everywhere — 68x, 1.8e6x, 1.1e8x. In this
   family the "perturber" is really an active third body, so the ungrouped
   Hermite run goes through a genuine three-body tangle and fails; **no
   "forming is harmful" threshold exists for tangles**.
1b. Many-perturber ring (12 spectators on a near-circular ring at three
   violence levels, kappa_ref=1e-4 as in production; encounter pair at the
   center): grouped wins 7000x / 400x / 6.4x — the advantage shrinks
   monotonically with environmental violence but never flips; beyond the
   strongest ring both paths are already chaotic-garbage (dE ~ 1-7). A
   smooth far-field tide strong enough to matter self-gravitates into its
   own mini-cluster (ring members group among themselves), and production's
   pert ratios 72-4380 arise from discrete NEARBY perturbers, i.e. the
   single-spectator family, not a smooth field. Side observation: with the
   sample-default kappa_ref=1e-6 the criterion REJECTS formation entirely in
   these tidal environments (t1/t3/t5 formed nothing until kappa_ref=1e-4) —
   the kappa gate threshold and the kappa_ref convention are strongly
   coupled; production uses 1e-4.
   **Verdict: the "worth-forming gate" optimization direction is dead** —
   forming is always beneficial in every regime that matters; the remaining
   criterion lever is the switch radius (sample-side 400x/3x/flat), whose
   production validation requires a targeted run, not replay.
2. Pal5 replay r-group scan: **blocked by design** — the effective r_crit
   (2.02e-4 / 2.87e-4 / 6.68e-4 across the three cases = 3.7-12.3x the par
   value via per-particle mass scaling) is identical across par-file r-group
   factors 0.1-10; the replay takes the criterion radii from the dump's
   embedded per-particle changeover state, ignoring the par edit (x30 trips
   an r_search bound assert). All three formations sit exactly at the radius
   boundary (dr == r_crit_eff) with kappa_org 1.9e5-4.5e11 far above the
   1e-2 gate, i.e. membership for the dump class is radius-crossing timed
   and kappa-confirmed. Production-side switch-position validation therefore
   requires a targeted production run with modified options (or binary dump
   patching), not replay par edits.

Net: the switch-position lever is confirmed sample-side (400x weak / ~3x
strong / flat bound) but untestable via replay; the kappa gate is the active
membership constraint for the Pal5 dump class; the harmful-threshold
question needs a many-perturber bench. Both production-side questions reduce
to decisions about running targeted Pal5 configurations.

- Zero-event runs must stay bit-identical (no events -> no code path change).
- Unit tests (`sample/test`) and hermite sample builds green.

### Phase 3 — validation gates

- G1 palindrome: with-event dvel_rt within 2x of the no-event floor (>= 20x
  improvement from today's ~1e-3 per event).
- G2 three Pal5 replays: event-step dE <= 1e-5 (physical integration error
  scale), zero algorithmic component; warning dumps no longer fire for these.
- G3 Plummer N=20 battery: completes, bounded energy, form/break balanced.
- G4 PeTar functional smoke (5 cases, bse-galpy x3) green.
- G5 criterion-symmetry regression: event-time mirroring unchanged.

### Phase 4 — semantics guard

Regression tripwire: assert/diagnostic that `de_change_cum` receives
contributions only from the dm / interrupt paths, never from transition
rewrites.

## 5. Risks

- #4 re-anchoring changes kappa(t) sampling cadence slightly -> small physics
  drift vs current runs; acceptable (slowdown is an approximation either way),
  but must be re-validated by G2/G3.
- Chatter (form/break/form within one step) exercises the inverse twice in
  quick succession — the pure design makes it idempotent, which is also the
  strongest test of invertibility (case 1 and case 3 replays).
- PeTar `HardIntegrator::initial()` group-init path (r_break_crit, c.m.
  changeover sync) wraps the SDAR transitions; keep both sides' derivations
  shared, not duplicated.
