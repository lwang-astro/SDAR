# SDAR Project Guidelines

SDAR (Symplectic, Drift, Acceleration, Regularization) is a C++11 header-only library for solving gravitational few-body problems. It is the close-encounter and multiple-system integrator inside PeTar.
[Full documentation](https://lwang-astro.github.io/SDAR/docs/html/index.html)

## Build and Test

`src/` contains only headers; executables are built from `sample/` subdirectories:

```bash
make -C sample/Kepler    # keplerorbit, keplertree
make -C sample/AR        # AR integrators (default flags build the .sd.cm variants)
make -C sample/Hermite   # hermite plus its variants
make -C sample/test      # unit tests
make -C tools            # installs the Python tools to ~/include/sdar/
```

**Compile flags**: `-std=c++11 -I../../src -I./ -O2 -Wall`. Add `-g -O0` for debug builds.
Feature selection is compile-time (`-D`), so each flag set produces a different binary — variant names, build targets, and per-variant prerequisites are in the [skill's binary-scenario-map](.github/skills/sdar-fewbody-integration/assets/binary-scenario-map.md).

### Feature flags

- `AR_TTL` — Time-Transformed Leapfrog (vs LogH default)
- `AR_SLOWDOWN_TREE` — hierarchical slowdown for binaries
- `AR_KDK_PERT` — KDK perturbation scheme
- `USE_MPFRC` — arbitrary precision via MPFR (requires `-lmpfr -lgmp`)
- `AR_G_FUNC_BLOGH` / `AR_G_FUNC_BTLOGH` — g-function method; exactly one per build (mutual exclusion enforced in `src/AR/g_func.h`). See `docs/hierarchical_blogh_impl_notes.md`; the last version implementing the removed variants is at git tag `gfunc-archive`.

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

The [sdar-fewbody-integration skill](.github/skills/sdar-fewbody-integration/SKILL.md) is the authority for standalone runs: input format, AR/Hermite variant selection, build and rebuild rules, run-safety gating, and Python post-processing. It is self-sufficient — no other file needs to be read to complete an integration.

## Lessons-Learned Capture

This process is canonical here and shared with the PeTar workspace (PeTar's `AGENTS.md` defers to it).

After each non-trivial task (integrator change, bug fix, simulation debugging, workflow change):

1. **Reflect**: did anything go wrong? A mistake, a misleading assumption, a silent failure, a confusing error message?
2. **If yes**: append the finding directly to [.github/skills/sdar-fewbody-integration/assets/lessons-learned.md](.github/skills/sdar-fewbody-integration/assets/lessons-learned.md) under the appropriate category, with **Mistake** / **Root cause** / **Prevention rule**, matching the entry style already present. **Title each entry `### <date>: <one-line rule>`** so that `grep -n '^### ' lessons-learned.md` is a complete index — do not split the file into per-lesson files (split trigger is in the maintenance doc below). Do not delegate this.
3. **If no**: no action needed.
4. **Periodically**: review entries and promote well-validated patterns into `SKILL.md` as hard rules.

Lessons rooted in SDAR code (`src/`, `sample/AR`, `sample/Hermite`, `tools/`) belong here — even when discovered inside a PeTar session.

**Before changing any customization file** (`AGENTS.md`, agent definitions, `SKILL.md`, assets) — in either repository — read [PeTar `.github/skills/README.md`](../PeTar/.github/skills/README.md): it owns the design goals, layering rules, update rules, and health-check procedure for the whole workspace.

## Repository Rules

- Do not commit or push without user confirmation.
- Bump `VERSION` in every commit: `bash get_version.sh` from the repo root (the script has no executable bit), then append the `e` experiment suffix and stage `VERSION`.
- Prefer repository-relative tracked files as inputs; never rely on untracked local files for test logic.

## Documentation

- [Doxygen docs](https://lwang-astro.github.io/SDAR/docs/html/index.html)
- [Local docs](docs/html/index.html)
- `docs/hierarchical_blogh_impl_notes.md` — BlogH implementation notes
- `docs/hierarchical_blogh_plan.md` — BlogH development plan
