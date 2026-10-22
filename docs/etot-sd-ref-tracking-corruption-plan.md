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
keeps only its release-gating remainder.

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
4. **ds estimator hyperbolic-root nonsense**: after a tree scramble the
   estimator can produce nonsense `ds_init` from hyperbolic-root "periods"
   (measured ds_init ≈ 35, PeTar plan anomaly (a)); the first estimate
   from hyperbolic fits is *accepted* because no previous candidate
   exists. (Still open after the item-5 fix: the creation-time estimate
   accepts its first candidate unconditionally by design.)
5. **BTLogH ds-decision defect**: ~~fix + validation + commit~~ done
   2026-10-22 (commits SDAR `96f0d39`, PeTar `d14ce63`) — root cause,
   fix design, and validation numbers are archived in
   `docs/hierarchical_blogh_impl_notes.md` (2026-10-22 header entry);
   prevention rules in lessons-learned (2026-10-22).
6. **Functional `merger` case stalls in a 2-body contact grind** (found
   2026-10-22 during the item-1 matrix run): the collision-interrupt IC
   pair (10+10 Msun, marginally hyperbolic fit semi≈−1.8e-8, ecc≈1.0006)
   collapses ds (ds_init≈0.098 → 5.5e-7) and grinds at T≈3.6e-11 with
   NaN H columns until the runner's 1800 s step timeout. A/B with a
   pre-item-5-fix binary stalls IDENTICALLY (same state, same dump path)
   — pre-existing, unrelated to the item-5 fix; both binaries include
   the 2026-10-21 sd-gauge fix, so the regression window (if any) is
   that fix or earlier. Needs its own investigation (suspects:
   post-interrupt ds handling for near-contact hyperbolic pairs;
   interplay with plan item 4's junk hyperbolic-root first estimate).
   Evidence: `PeTar/test/out/functional_merger/`, A/B dirs `/tmp/mAB/`.

## Cross-references

- Investigation record + implemented fix (sd-energy gauge):
  `docs/etot_sd_ref_corruption_fix_notes.md` (this plan is its open-items
  remainder).
- Item-5 fix record (BTLogH ds estimation architecture):
  `docs/hierarchical_blogh_impl_notes.md` (2026-10-22 header entry +
  调用点表).
- Parent investigation: `PeTar/doc/restart-slowdown-cold-state-plan.md`.
