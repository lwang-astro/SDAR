# SDAR Minimal Question Sets

This file defines the minimum required inputs for each SDAR scenario.
When the skill is invoked, ask only the missing items from the relevant checklist.
If the user already provided a parameter, do not ask again.

## Rule

If enough information is already present in the user's request, do not ask follow-up questions.
Only ask for parameters that are both missing AND required for the scenario.

---

## AR Integration (LogH / TTL)

| # | Parameter | Flag / Where | Why it matters | Default |
|---|-----------|-------------|----------------|---------|
| 1 | Working directory | `mkdir ~/sdar_run/<name>` | Where outputs go; avoid polluting source tree | Must ask |
| 2 | Input file | positional | Particle data and optional tree structure | Must ask |
| 3 | AR variant | binary name | logh (accurate) vs ttl (fast), ±sd (slowdown), ±mpfrc (precision) | Must ask — see binary-scenario-map.md |
| 4 | End time `-t` or step count `-n` | `-t <t>` or `-n <n>` | Duration of integration | Must set one |
| 5 | Step size `-s` | `-s <ds>` | ≤0 = auto; >0 = fixed | Auto (0.0) — confirm |
| 6 | Gravitational constant `-G` | `-G <G>` | 1.0 (Henon) or 0.004498... (Msun/pc/Myr) | 1.0 — confirm if units matter |
| 7 | Energy error limit `-e` | `-e <de>` | Controls integration accuracy | 1e-10 — confirm |
| 8 | Output interval `-o` | `-o <dt>` | How often to print output; 0 = every step | 0.0 — confirm |
| 9 | Symplectic order `-k` | `-k <order>` | Even number; negative = Yoshida 2nd | -6 — confirm |
| 10 | Interrupt detection `-i` | `-i <0|1|2>` | 0=off, 1=modify, 2=record | 0 — confirm |
| 11 | Slowdown ref | `--slowdown-ref` | Perturbation ratio for slowdown (sd variants only) | 1e-6 |
| 12 | Slowdown timescale max | `--slowdown-timescale-max` | Max timescale for slowdown factor | time-end |

**Enforcement rule**: Before composing an AR command, confirm that at minimum items #1-4 are resolved. Items #5-12 have safe defaults but should be confirmed for production runs.

---

## Hermite+AR Hybrid Integration

| # | Parameter | Flag / Where | Why it matters | Default |
|---|-----------|-------------|----------------|---------|
| 1 | Working directory | `mkdir ~/sdar_run/<name>` | Where outputs go | Must ask |
| 2 | Input file | positional | N-particle data | Must ask |
| 3 | End time `-t` | `-t <t>` | Duration of integration | 1.0 — must confirm |
| 4 | Group radius `-r-group` | `-r-group <r>` | Controls AR/Hermite switching threshold | 1e-3 — must confirm |
| 5 | Neighbor radius factor | `-r-neighbor-over-group <f>` | Neighbor search = f × r_group | 2.0 — confirm |
| 6 | Gravitational constant `-G` | `-G <G>` | 1.0 (Henon) or physical | 1.0 — confirm |
| 7 | Energy error limit `-e` | `-e <de>` | AR accuracy within groups | 1e-10 — confirm |
| 8 | 4th-order time step coeff `-eta-4th` | `-eta-4th <eta>` | Hermite step accuracy | 0.1 — confirm |
| 9 | 2nd-order time step coeff `-eta-2nd` | `-eta-2nd <eta>` | Hermite step accuracy | 0.001 — confirm |
| 10 | Softening `-eps` | `-eps <eps>` | Force softening | 0.0 — confirm |
| 11 | Output interval `-o` | `-o <power>` | Power index of 0.5 for output | 2 — confirm |
| 12 | Interrupt detection `-i` | `-i <0|1|2>` | 0=off, 1=modify, 2=record | 0 — confirm |
| 13 | Slowdown ref | `--slowdown-ref` | Perturbation ratio for slowdown | 1e-6 |
| 14 | Slowdown timescale max | `--slowdown-timescale-max` | Max timescale for slowdown factor | time-end |

**Enforcement rule**: Items #1-4 are mandatory. For systems with N > 10, `-r-group` is particularly critical — a wrong value causes either missed groups or false positives. Items #5-14 have reasonable defaults but confirm for production.

---

## Kepler Binary Tree (keplertree)

| # | Parameter | Flag / Where | Why it matters | Default |
|---|-----------|-------------|----------------|---------|
| 1 | Input file | positional | Particle data (N lines) | Must ask |
| 2 | Inverse mode `-i` | `-i` | Read tree structure → output Cartesian | Off — ask |
| 3 | Unit system `-u` | `-u <0-4>` | Affects semi-major axis and period interpretation | 0 — confirm if not Henon |

---

## Kepler Orbit Conversion (keplerorbit)

| # | Parameter | Flag / Where | Why it matters | Default |
|---|-----------|-------------|----------------|---------|
| 1 | Input file | positional | Two particle lines per pair (or one Kepler line if -i) | Must ask |
| 2 | Inverse mode `-i` | `-i` | Kepler → Cartesian | Off — ask |
| 3 | Number of pairs `-n` | `-n <n>` | How many pairs to read | 1 — confirm |
| 4 | Unit system `-u` | `-u <0-4>` | Affects unit interpretation | 0 — confirm |

---

## Python Post-Processing

| # | Parameter | Where | Why it matters | Default |
|---|-----------|-------|----------------|---------|
| 1 | Output file to read | file path | AR vs Hermite output → different reader class | Must ask |
| 2 | Producing binary | user states | Determines reader: SDARData vs HermiteData | Must ask |
| 3 | Slowdown enabled? | compile flag | Pass `slowdown=True` to SDARData/HermiteData | False — confirm |
| 4 | Timing profile enabled? | compile flag | Pass `time_measure=True` to readers | False — confirm |
| 5 | Gravitational constant | `-G` value | Must match C++ run: 1.0 or 0.004498... | 1.0 — confirm |

**Enforcement rule**: Items #1-2 determine the reader class — getting them wrong silently produces garbage. Always confirm.

---

## Multi-Method Comparison

When the user wants to compare different AR methods (e.g., logh vs ttl, with/without slowdown):

| # | Parameter | Notes |
|---|-----------|-------|
| 1 | Input file | Same input for all methods |
| 2 | Methods to compare | From: logh, logh.ttl, logh.sd, logh.ttl.sd, logh.sd.kdkpert (current Makefile targets) |
| 3 | Common parameters | `-t`, `-n`, `-G`, `-e` apply to all |
| 4 | Binary availability | Check which variants are compiled before composing commands |

---

## Cross-Scenario Notes

- **AR → Hermite transition**: If the user starts with AR and later wants Hermite, the input format is identical (same N-particle file), but the parameter set changes significantly (`-r-group` vs `-s`). Remind the user of the different parameter semantics.
- **SDAR → PeTar**: When using SDAR to debug PeTar behavior, ensure the G constant and unit system match. PeTar's `--r-group` maps to SDAR Hermite's `-r-group`.
- **MPFRC variants**: Require `-lmpfr -lgmp` at link time. The compile flag `-D USE_MPFRC` must be set. Check `libmpfr-dev` availability before compiling.
