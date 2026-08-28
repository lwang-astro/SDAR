# SDAR Binary Scenario Map

This note documents SDAR integrator selection by scenario.

## Naming Scheme (2026-08-27)

Binary names follow `ar.<method>[.ttl][.sd][.kdkpert][.cm][.mpfrc]`:
- `<method>` is the g-function form: `logh`, `blogh`, `normblogh`, `mulall`, `btlogh`, `maxpot`, `addpot`
- `.ttl` marks the Time-Transformed Leapfrog implementation of the method (absent = LogH implementation)
- `.sd` is the tree-based hierarchical slowdown (the old `.sd.t` tag dropped its `.t`; `.sd.a` is deprecated)
- `.kdkpert` KDK perturbation splitting (old `kdk.pert`), `.cm` CM-frame build, `.mpfrc` MPFR precision

## AR (Algorithmic Regularization) Variants

All AR binaries are built from `sample/AR/ar.cxx` with different compile flags.

### Primary Integrators

| Binary | Method | Time Function | Slowdown | Use Case |
|--------|--------|---------------|----------|----------|
| `ar.logh` | LogH | `log(T+U)` | None | Isolated binary, best Kepler accuracy |
| `ar.logh.ttl` | LogH+TTL | `1/|U|` | None | Faster, simpler, for moderate eccentricity |
| `ar.logh.sd` | LogH | `log(T+U)` | Tree-mode (hierarchical) | Hierarchical triple/quadruple |
| `ar.logh.ttl.sd` | LogH+TTL | `1/|U|` | Tree-mode (hierarchical) | Hierarchical triple/quadruple (faster) |

> **Deprecated**: `ar.*.sd.a` (averaged slowdown) and the old naming (`ar.ttl.*`, `.sd.t`, `kdk.pert`)
> are removed from the current Makefile. Existing `~/bin/` binaries under old names come from older
> versions; do not use them.

### Extended AR Variants

| Binary | Additional Feature | Use Case |
|--------|-------------------|----------|
| `ar.logh.sd.kdkpert` | KDK perturbation splitting | Weakly perturbed binaries |
| `ar.logh.ttl.sd.kdkpert` | KDK perturbation splitting | Same, with TTL |
| `ar.blogh.ttl.sd[.cm]` | BLogH g-func (product of innermost pair potentials) | Systems with multiple inner binaries |
| `ar.normblogh.ttl.sd[.cm]` | Normalized BLogH (geometric mean) | Same, g ~ energy dimension |
| `ar.mulall.ttl.sd[.cm]` | Product over ALL pairs (needs `--s`) | Small-N systems |
| `ar.btlogh.ttl.sd[.cm]` | Tree-level product (BTLogH, inner x outer nodes) | Hierarchical quadruples+ (B-B) |
| `ar.maxpot.ttl.sd.cm` | Max innermost pair potential + CM | Alternative multi-binary handling |
| `ar.addpot.ttl.sd.cm` | Sum of innermost pair potentials + CM | Alternative multi-binary handling |
| `ar.logh.mpfrc` | MPFRC arbitrary precision | High-precision requirement |
| `ar.logh.ttl.mpfrc` | MPFRC arbitrary precision | High-precision with TTL |

### Decision Tree for AR

```
Is the system hierarchical (nested binaries)?
├── Yes → Use slowdown variant (ar.*.sd, N_sd = number of binary pairs)
│   ├── Weak perturbation? → ar.*.sd.kdkpert
│   └── Multiple inner binaries? → ar.blogh.ttl.sd (or normblogh/mulall)
│       └── Hierarchical quadruple+ (B-B)? → ar.btlogh.ttl.sd
└── No → Use plain variant
    ├── High eccentricity? → ar.logh
    └── Moderate eccentricity? → ar.logh.ttl (faster)
```

## Hermite+AR Hybrid

| Binary | Features | Use Case |
|--------|----------|----------|
| `hermite` | Standard 4th-order Hermite + AR | General few-body systems with close encounters |
| `hermite.kdkpert` | With KDK perturbation | Weakly perturbed subsystems |
| `hermite.mpfrc` | MPFRC arbitrary precision | High-precision few-body |

**When to use Hermite over AR:**
- System has both isolated particles and bound subsystems (e.g., a star cluster core with binaries)
- N > 10 where direct AR would be inefficient
- Need automatic group detection and switching between Hermite and AR

## Kepler Tools

| Binary | Purpose |
|--------|---------|
| `keplerorbit` | Convert between particle pairs and Kepler orbital parameters |
| `keplertree` | Build/dissolve hierarchical binary trees from particle data |

## Build Prerequisites

| Variant | Required Compile Flag | Required Libraries |
|---------|----------------------|-------------------|
| `sd` | `-D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE` | None |
| `kdkpert` | `-D AR_KDK_PERT` | None |
| `ttl` | `-D AR_TTL` | None |
| `blogh` | `-D AR_G_FUNC_BLOGH` | None |
| `normblogh` | `-D AR_G_FUNC_NORM_BLOGH` | None |
| `mulall` | `-D AR_G_FUNC_MUL_ALL_POT` | None |
| `btlogh` | `-D AR_G_FUNC_BTLOGH` | None |
| `maxpot` | `-D AR_G_FUNC_MAX_POT` | None |
| `addpot` | `-D AR_G_FUNC_ADD_INNER_POT` | None |
| `mpfrc` | `-D USE_MPFRC` | `libmpfr`, `libgmp` |

One g-func method macro per build (mutual exclusion enforced in `src/AR/g_func.h`).
`--g-func`: 0 = standard LogH, 1 = the method of this build, 2 = auto switch
(rejected for btlogh and mulall). The old `--g-func-switch` option and the
`AR_G_FUNC_MUL_POT` macro no longer exist.

## Skill Behavior

- If the user describes a system type, recommend the appropriate AR/Hermite variant.
- If the requested variant is not built, guide compilation before composing commands.
- If the user asks for slowdown but the system is non-hierarchical, warn that slowdown provides no benefit.
- For Hermite runs, ensure `-r-group` and `-r-neighbor-over-group` are appropriate for the particle configuration.
