# Legacy group form/break criteria (removed 2026-10-01, reference note)

The Hermite-layer group adjust (`HermiteIntegrator::adjustGroups` / `checkNewGroup`)
used a set of velocity-direction-gated criteria before the unified state-function
criterion (`groupedCriterion`) replaced them. This note preserves the legacy scheme
for reference (e.g. when comparing old runs or revisiting design decisions).

## Formation (checkNewGroup)

1. Distance gate: pair separation `d < r_crit`, `r_crit = max(getRGroup())` of the
   two objects (for a single joining a group: the group's `r_break_crit`).
2. **Direction gate**: merge only when incoming — `calcDrDv(pi, pj) < 0`
   (`drdv` is velocity-odd), or on the first adjust (`_start_flag`).
3. Perturbation gate: slowdown estimate `kappa_org = kappa_ref * pert_in/pert_out >= 1e-2`
   (`calcPertFromMR` for pert_in, `calcPertFromForcePot` on the CM acc/pot for pert_out);
   strongly perturbed candidates were rejected outright (no apoapsis consideration).

## Break (adjustGroups group loop)

All four legacy break paths required **outgoing** motion (`ecca > 0` eccentric-anomaly
sign, or `drdv > 0` radial velocity), i.e. velocity-odd gates:

1. *Binary escape*: elliptic root (`semi>0`), `ecca>0`, `r > r_break_crit`.
2. *Predicted escape*: elliptic, `apo > r_crit`, next-step prediction
   `rp = (drdv/r)*dt_cm + r > r_crit` (one-sided linear prediction in the outgoing
   direction); sub-branch: if currently deep inside (`r < 0.2 r_crit`) and the
   half-step prediction stays inside, halve the CM step instead of breaking.
3. *Hyperbolic escape*: `semi<0`, `drdv>0`, `r > r_crit` (with the same prediction
   variant and step-halving sub-branch).
4. *Strong perturbation*: `outgoing_flag` (from 1/3), `n_member==2`,
   `kappa_org < 1e-2`, and `apo > r_crit || semi < 0`.
5. (Compiled out in slowdown builds: inner-subtree `kappa_in_max > 5` break for
   few-body groups without inner AR slowdown.)

## Why replaced

Palindrome (time-reversal) tests showed the legacy scheme is structurally
asymmetric: formation and break used different gate types (incoming-only merge vs
outgoing-only break) plus one-sided next-step predictions, so the event sequence of
a reversed integration did not mirror the forward one (measured 2-13 forward vs 0-3
backward events on identical ICs; see `assets/lessons-learned.md` 2026-09-29/30
entries for the full attribution chain).

The replacement is a single state function used by BOTH sides
(`form <=> P(state)`, `break <=> !P(state)`):

```
grouped(d, apo, r_crit, kappa_org) = (d < r_crit) && (apo <= r_crit || kappa_org >= 1e-2)
```

Behavioural differences vs legacy (by design):
- no direction/prediction gates: events fire on the state alone (more events per
  orbit near the threshold; event times mirror under time reversal);
- tight-bound configurations (`apo <= r_crit`) merge/stay grouped regardless of the
  perturbation estimate (legacy rejected strongly perturbed candidates even when
  tightly bound);
- the step-halving "wait and see" sub-branches are gone (they relied on the
  outgoing direction).

The `kappa_org` input is computed inside `groupedCriterion` (not at the call
sites) so both sides share one recipe in the group-c.m. convention: specific
external acc/pot fed to `calcPertFromForcePot`, `pert_in = calcPertFromMR` on
the live separation, `kappa_ref` scaled by the mass ratio exactly like a newly
formed group. Before this centralisation the two sides silently disagreed in
three ways: formation fed force-unit mass-weighted fields that still contained
the internal mutual term while break used the group-c.m. specific fields (a
factor ~m_cm, plus the mutual term, in `pert_out`); formation scaled `kappa_ref`
by the mass ratio while break reused the stored (unscaled) group value; and the
break-side guard (`r > r_break_crit`) left the kappa branch unreachable, so
break reduced to `d > r_crit` — not the inverse of formation in the wide-apo
regime. Formation now passes the two-particle mass-weighted average of the
external field (mutual terms removed; for a candidate pair containing a group
c.m. the same formula applies since c.m. particles carry external-only fields);
break passes the group c.m. fields directly. The residual, physical difference
is O(d/L): two-point average vs one-point c.m. evaluation. The kappa test is
skipped on the break side for groups with more than two members and during the
first adjust (`_start_flag`): the root pair of an N>2 group is not the pair
formation judged. Validation: in-memory palindrome unchanged or improved
(dpos 2.5e-1 -> 1.7e-2 at the longest reversal; event counts 3/1 -> 4/3);
a targeted flyby IC (wide-apo e=0.9 binary + m=20 perturber) exercises the
now-live binary kappa break (d << r_crit, kappa_org ~ 4e-3 < 1e-2) with bounded
energy drift (max |dE/E| ~ 1e-3, one form/break pair per inner orbit).

