# SDAR Binary Scenario Map

This note documents SDAR integrator selection by scenario.

## AR (Algorithmic Regularization) Variants

All AR binaries are built from `sample/AR/ar.cxx` with different compile flags.

### Primary Integrators

| Binary | Method | Time Function | Slowdown | Use Case |
|--------|--------|---------------|----------|----------|
| `ar.logh` | LogH | `log(T+U)` | None | Isolated binary, best Kepler accuracy |
| `ar.ttl` | TTL | `1/|U|` | None | Faster, simpler, for moderate eccentricity |
| `ar.logh.sd.t` | LogH | `log(T+U)` | Tree-mode (hierarchical) | Hierarchical triple/quadruple |
| `ar.ttl.sd.t` | TTL | `1/|U|` | Tree-mode (hierarchical) | Hierarchical triple/quadruple (faster) |

> **Deprecated**: `ar.*.sd.a` (averaged slowdown) was removed in current Makefile (`TARGET`).
> Existing `~/bin/ar.*.sd.a` binaries are from older versions; do not use them.

### Extended AR Variants

| Binary | Additional Feature | Use Case |
|--------|-------------------|----------|
| `ar.logh.sd.t.kdk.pert` | KDK perturbation splitting | Weakly perturbed binaries |
| `ar.ttl.sd.t.kdk.pert` | KDK perturbation splitting | Same, with TTL |
| `ar.ttl.sd.t.mulpot` | Multi-potential time function (product) | Systems with multiple inner binaries |
| `ar.ttl.sd.t.mulpot.cm` | Multi-potential + CM tracking | Drifting multi-binary systems |
| `ar.ttl.sd.t.maxpot.cm` | Max-potential time function + CM | Alternative multi-binary handling |
| `ar.ttl.sd.t.addpot.cm` | Additive-potential time function + CM | Alternative multi-binary handling |
| `ar.logh.mpfrc` | MPFRC arbitrary precision | High-precision requirement |
| `ar.ttl.mpfrc` | MPFRC arbitrary precision | High-precision with TTL |

### Decision Tree for AR

```
Is the system hierarchical (nested binaries)?
├── Yes → Use slowdown variant (ar.*.sd.t, N_sd = number of binary pairs)
│   ├── Weak perturbation? → ar.*.sd.t.kdk.pert
│   └── Multiple binaries? → ar.ttl.sd.t.mulpot
└── No → Use plain variant
    ├── High eccentricity? → ar.logh
    └── Moderate eccentricity? → ar.ttl (faster)
```

## Hermite+AR Hybrid

| Binary | Features | Use Case |
|--------|----------|----------|
| `hermite` | Standard 4th-order Hermite + AR | General few-body systems with close encounters |
| `hermite.kdk.pert` | With KDK perturbation | Weakly perturbed subsystems |
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
| `sd.t` | `-D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE` | None |
| `kdk.pert` | `-D AR_KDK_PERT` | None |
| `ttl` | `-D AR_TTL` | None |
| `mulpot` | `-D AR_TIME_FUNCTION_MUL_POT` | None |
| `maxpot` | `-D AR_TIME_FUNCTION_MAX_POT` | None |
| `addpot` | `-D AR_TIME_FUNCTION_ADD_POT` | None |
| `mpfrc` | `-D USE_MPFRC` | `libmpfr`, `libgmp` |

## Skill Behavior

- If the user describes a system type, recommend the appropriate AR/Hermite variant.
- If the requested variant is not built, guide compilation before composing commands.
- If the user asks for slowdown but the system is non-hierarchical, warn that slowdown provides no benefit.
- For Hermite runs, ensure `-r-group` and `-r-neighbor-over-group` are appropriate for the particle configuration.
