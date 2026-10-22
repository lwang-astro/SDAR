# Plan: sd-energy gauge fix — open follow-ups

Status: plan — open items only. The investigation and the implemented fix
(record, mechanism, validation) are archived in the reference note
`docs/etot_sd_ref_corruption_fix_notes.md`; prevention rules live in
`lessons-learned` (entries 2026-10-19 #5/#6, 2026-10-20 #8/#9, 2026-10-21).
Opened 2026-10-20; fix executed 2026-10-21; production-validated 2026-10-22
(N=1000 cluster A/B — see the reference note, §Production validation).

1. **Regression test wiring**: fold the rp≥0.35 placements into the renamed
   T7 scenario (`t7_triple_perturber_ar_ds_adjustment`, currently band-only
   r1/r2/r3) as abort→complete guards with step/energy budget asserts; the
   reproducer IC + `petar.init -s merger` recipe is in
   `~/sdar_run/etot-research/solver-fix/`. Then the full matrix on the
   canonical binary family: T1–T7 + functional smoke + SDAR samples
   (the merger variant needs `petar.init -s merger` input conversion; the
   g4e.dump replay needs a binary with matching defines). The originally
   sketched κ-ramp invariant test is moot (`de_sd ≡ ΔS` proven at runtime
   to 1e-9 relative).
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
   exists.

## Cross-references

- Investigation record + implemented fix:
  `docs/etot_sd_ref_corruption_fix_notes.md` (this plan is its open-items
  remainder).
- Parent investigation: `PeTar/doc/restart-slowdown-cold-state-plan.md`.
