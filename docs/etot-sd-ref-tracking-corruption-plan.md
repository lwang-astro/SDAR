# Plan: sd-energy gauge fix — open follow-ups

Status: plan — open items only. The investigation and the implemented fix
(record, mechanism, validation) are archived in the reference note
`docs/etot_sd_ref_corruption_fix_notes.md`; prevention rules live in
`lessons-learned` (entries 2026-10-19 #5/#6, 2026-10-20 #8/#9, 2026-10-21,
and 2026-10-22 for item 5). Opened 2026-10-20; fix executed 2026-10-21;
production-validated 2026-10-22 (N=1000 cluster A/B — see the reference
note, §Production validation). Item 5 (BTLogH ds-decision defect) was
fixed and validated 2026-10-22 — its record is distilled into
`docs/hierarchical_blogh_impl_notes.md` (2026-10-22 entry); this file
keeps only its release-gating remainder. The LogH hyperbolic S-criterion
defect (merger NaN / negative-dt grind) was fixed and validated 2026-10-23
(S-cap in `calcLogHSumGaugeIter`), item 4 closed entirely the same day
(its g-func remainder was rejected on design review, note §5), and the
item-6 limit-cycle investigation closed; records are distilled into
`docs/step_size_effect_logh_notes.md` (§5 item-4 close-out, §6
implementation + validation, §7 residual mechanism).

1. **Regression test wiring**: ~~T7 abort-band placements~~ done
   2026-10-22 — record: PeTar `test/validation` T7 scenario
   `t7_triple_perturber_ar_ds_adjustment` (22 checks, former abort band
   pinned as complete-guards; commit `2ae102c`; the originally sketched
   κ-ramp invariant test was dropped as moot — `de_sd ≡ ΔS` proven at
   runtime to 1e-9 relative). ~~Full matrix~~ done 2026-10-22 on the
   fixed trees: T1–T8 all pipelines pass (9 reports), functional smoke
   green except `merger` (see item 6), SDAR sample builds + test binaries
   all clean.
2. **LogH tree-path κ-weighted kick credit** (code note at `kickEtot` in
   `src/AR/symplectic_integrator.h`): the TTL tree walk credits κ-weighted
   work via `vel_sd`; the flat LogH loop still credits full weight — needs
   a tree walk for `vel_sd` (PeTar hard uses TTL, so production is
   unaffected).
3. **LogH runtime validation** of the negative-dt → ds-floor-repair
   routing (compile-verified only; needs a LogH SDAR-sample injection
   run).
4. **ds estimator hyperbolic-root nonsense**: CLOSED 2026-10-23 — LogH
   side resolved by the S-criterion fix (semi-gauge u_min + t_cross
   floor + S≤0.25 cap in `calcLogHSumGaugeIter`; merger reproducer
   0 NaN / 0 negative-dt; controlled g5-vs-g6: warnings 168→6,
   wallclock 43→9.5 min, error same order); the proposed g-func
   remainder (mirror the t_cross floor into `calcEffectivePeriod`)
   REJECTED on design review — full rationale in
   `docs/step_size_effect_logh_notes.md` §5 (bundled-fix origin never
   isolated; the g-func product formula cancels an isolated hyperbolic
   pair's period, so the far-leg cost pathology is sum-gauge-specific;
   event-driven ds re-estimation is already load-bearing there).
5. **BTLogH ds-decision defect**: ~~fix + validation + commit~~ done
   2026-10-22 (commits SDAR `96f0d39`, PeTar `d14ce63`) — root cause,
   fix design, and validation numbers are archived in
   `docs/hierarchical_blogh_impl_notes.md` (2026-10-22 header entry);
   prevention rules in lessons-learned (2026-10-22).
6. **κ-huge hierarchy controller limit cycle (grind class) — attribution
   CORRECTED 2026-10-22 (both earlier attributions withdrawn after
   controlled verification); mechanism ROOT-CAUSED 2026-10-23 by dump
   replay**: the functional `merger` case stall and the
   N=1000 grind/energy-error events are a SHARED, pre-existing defect
   class (baseline, pre-etot-fix, already shows 10 warnings): compact
   binaries with capped κ (SD_org up to ~1e8, κ~6e4) inside hierarchies
   drive the step controller into a limit cycle with mega-grinds
   (7.6M–57M steps) and large dE_SD events (10–13%). Mechanism
   (2026-10-23 replay, CORRECTING the earlier "ds-independent error
   floor" reading — the per-step error scales with ds²): backup
   recovery (`countAndRecover`/`Reuse_backup_ds`) re-injects a
   too-large ds; the healthy halving cascade converges to an asymptote
   pinned at 1.00–1.02× ar-max-error; a quiet phase follows; then
   recovery fires again. Withdrawn
   attributions: (a) F1/semi-gauge "high-e under-resolution" — the g5
   catastrophic group (`h4_n35_g13_t1.90`: 14.3+13.9 κ=5.9e4 compact
   binary, ALL-elliptic) is F1-untouched; g5-vs-g6 cum differences are
   chaotic single-realization variance of this shared class; (b) item-5
   estimation removal — snapshot-restart A/B from `g6.6` (t=3→6, same
   trajectory): BOTH binaries grind on the same 83.7+21.6 group and the
   PRE-item-5 arm is WORSE in-window (cum −12,502 vs −252). What remains
   valid: the item-5 fix itself (deterministic bit-identical
   reproducer) and F1's merger-case NaN elimination (deterministic,
   S-criterion). Closed 2026-10-23: (i) replay root-cause — done
   (mechanism above; threshold scan 1e-6/1e-8/1e-10 → 710/127/118
   events, asymptote ds 5e-6/7.7e-8. Grind-group anatomy (g6): 4-body
   hierarchical quadruple — 105 Msun compact binary (a≈6.5e-7,
   κ≈4e5) + hardening light pair + wide e≈0.94 root — slowdown inverts
   the timescale ordering (light-pair collapse faster than the slowed
   inner binary), which the estimator's un-slowed-minimum-period
   structure cannot express, ~120× estimator-vs-controller divergence;
   record `docs/step_size_effect_logh_notes.md` §7); (ii) F1 isolated
   production verification — done (controlled g5-vs-g6, only F1
   differs: warnings 168→6, wallclock 43→9.5 min, error same order;
   deterministic merger reproducer 0 NaN / 0 negative-dt; gauge scope:
   BTLogH N=1000 (g3) 0 warnings / 10 min — the grind class is
   LogH-path and equals T8's r1/r2 "cascade" class, where BTLogH is
   38–934× cheaper per event). Open: the fix itself — bounded recovery
   (cap restored ds at ×2 last-good) vs survival window vs estimator
   re-ranking with slowed effective periods; decision pending.
   Evidence: `bench-n1k-arfix/{g3,g5,g6,rA6,rB6}.log` + dumps.

## Cross-references

- Investigation record + implemented fix (sd-energy gauge):
  `docs/etot_sd_ref_corruption_fix_notes.md` (this plan is its open-items
  remainder).
- Item-5 fix record (BTLogH ds estimation architecture):
  `docs/hierarchical_blogh_impl_notes.md` (2026-10-22 header entry +
  调用点表).
- LogH S-criterion fix (S-cap) + residual recovery-injection limit
  cycle: `docs/step_size_effect_logh_notes.md` (§6
  implementation/validation, §7 mechanism).
- Parent investigation: archived 2026-10-22 (`PeTar/doc/restart-slowdown-cold-state-plan.md`,
  distilled and deleted; record: `docs/etot_sd_ref_corruption_fix_notes.md` +
  PeTar lessons-learned 2026-10-19/2026-10-22).
