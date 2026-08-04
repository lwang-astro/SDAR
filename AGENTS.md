# SDAR Project Guidelines

SDAR (Symplectic, Drift, Acceleration, Regularization) is a C++11 header-only library for solving gravitational few-body problems.
It's the close-encounter and multiple-system integrator inside PeTar.
[Full documentation](https://lwang-astro.github.io/SDAR/docs/html/index.html)

## Build and Test

SDAR is a **header-only** library — `src/` contains only `.h` files. Executables are built from `sample/` subdirectories:

```bash
# Build all samples
make -C sample/Kepler
make -C sample/AR
make -C sample/Hermite
make -C sample/test

# Install Python tools (to ~/include/sdar/)
make -C tools
```

**Compile flags**: `-std=c++11 -I../../src -I./ -O2 -Wall`. Add `-g -O0` for debug builds.

**Feature flags** (enabled via `-D`):
- `AR_TTL` — Time-Transformed Leapfrog (vs LogH default)
- `AR_SLOWDOWN_TREE` — hierarchical slowdown for binaries
- `AR_KDK_PERT` — KDK perturbation scheme
- `USE_MPFRC` — arbitrary precision via MPFR (`-lmpfr -lgmp`)

## Architecture

Three namespaced components in `src/`, all header-only:

| Namespace | Directory | Purpose |
|-----------|-----------|---------|
| `COMM::` | `src/Common/` | Float type, List container, ParticleGroup, BinaryTree (Kepler conversion), Vector3/Matrix3, KDTree, I/O |
| `AR::` | `src/AR/` | Time-transformed symplectic integrator (LogH/TTL) with slowdown method |
| `H4::` | `src/Hermite/` | Hybrid 4th-order Hermite + AR integrator for global few-body systems |

`H4::HermiteIntegrator` wraps multiple `AR::TimeTransformedSymplecticIntegrator` instances — one per close subsystem detected by group analysis.

Both integrators are heavily templated: users provide particle, interaction, perturber, and information types.
See `sample/AR/` and `sample/Hermite/` for concrete instantiation patterns.

## Key Conventions

- **`COMM::List<T>` not `std::vector`**: Three modes — `copy` (own + track origin), `local` (own only), `link` (external data). Used throughout for particle management.
- **Compile-time feature selection**: Features are `#define`-gated, not runtime. Each `-D` flag produces a different binary.
- **Binary + ASCII I/O**: All data classes implement both `writeBinary`/`readBinary` and column-formatted ASCII output.
- **`Float` type**: Aliased in `src/Common/Float.h` — defaults to `double`, optionally `mpfr::mpreal`.
- **Python tools**: Post-processing package in `tools/`. Install with `make -C tools`, then `import sdar`. See `sample/data_analysis.ipynb` for usage patterns.

## Running Integrations

See [.github/skills/sdar-fewbody-integration/SKILL.md](.github/skills/sdar-fewbody-integration/SKILL.md) for detailed workflow guidance including input format, binary selection, and safety rules.

## Documentation

- [Doxygen docs](https://lwang-astro.github.io/SDAR/docs/html/index.html)
- [Local docs](docs/html/index.html)
- `docs/hierarchical_blogh_impl_notes.md` — BlogH implementation notes
- `docs/hierarchical_blogh_plan.md` — BlogH development plan
