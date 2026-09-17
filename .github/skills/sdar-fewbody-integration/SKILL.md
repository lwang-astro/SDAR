---
name: sdar-fewbody-integration
description: "Use when: setting up or running SDAR few-body integrations, including AR (LogH/TTL) symplectic integrator, Hermite+AR hybrid integrator, slowdown method for hierarchical systems, Kepler binary tree construction and orbit conversion, Python post-processing with tools/ modules, and reading SDAR/Hermite output snapshots."
---

# SDAR Few-Body Integration Skill

## Purpose

Provide strict, command-level guidance for SDAR workflows in this repository.
SDAR (Slow-Down Algorithmic Regularization) is a library for solving few-body problems,
serving as the close-encounter and multiple-system integrator inside PeTar.
This skill prioritizes correct integrator choice, input format compliance,
build correctness, and Python data analysis patterns.

## Non-Negotiable Rules

### Gate 1 — Working directory and inputs

- **Always `cd` into a dedicated working directory before running any SDAR integrator.** Never run inside `sample/` or the repository root — those hold source code and build artifacts only. If the user does not specify a working directory, ask rather than choosing; if they have no preference, suggest `~/sdar_run/<descriptive-name>` and confirm before creating.
- **Confirm which binary the scenario needs before composing a command** — AR or Hermite, LogH or TTL, with or without slowdown, MPFRC or not. Do not run whatever happens to be in `~/bin/`. See `assets/binary-scenario-map.md`.
- **Validate binary availability** first: `command -v <binary>` or `ls <path-to-build>`. If it is not built, guide the user through `make` in the correct `sample/` subdirectory.
- **Never mutate input files.** SDAR inputs are small ASCII tables — keep the original and copy one into the working directory if parameters need changing.
- **For hierarchical systems, prefer slowdown (`sd`) binaries** (`ar.logh.sd`, `ar.logh.ttl.sd`): the slowdown method dramatically accelerates weakly perturbed inner binaries compared with plain `ar.logh` / `ar.logh.ttl`.

### Gate 2 — Confirmation and traceability

- **Before any solver execution, present a structured parameter summary and wait for explicit confirmation.** Contents, in order: (1) working directory and input file; (2) full command line (binary, flags, redirect); (3) expected runtime hint (from particle count, end time, step size); (4) expected output files; (5) a clear prompt — *"Proceed? (y/n)"* — then stop and wait.
- **Record every command to `commands.log`** in the working directory:
  ```bash
  echo "# $(date): ~/bin/ar.logh.sd -t 1.0 -s 0.01 input.dat" >> commands.log
  ```
- **Redirect stdout of every major command to a file** — SDAR emits many diagnostic lines that cannot be recovered afterwards. Naming: `ar.logh >ar_logh.log`, `ar.logh.ttl >ar_logh_ttl.log`, `hermite >hermite.log`, `keplertree >keplertree.log`, `keplerorbit >keplerorbit.log`.

### Gate 3 — Units, constants, and precision

- **`-G` must match the unit system, and the C++ run and the Python analysis must agree.** The only two valid values are `G_HENON = 1.0` (unscaled / Henon N-body units — the C++ default) and `G_MSUN_PC_MYR = 0.00449830997959438` (Msun/pc/Myr, unit `-u 4`). Never mix them: a mismatch silently produces wrong energies and orbits.
- **Use the constants from `tools/particle.py`** — `G_MSUN_PC_MYR` and `G_HENON` — rather than hardcoding numeric values.
- **Use the `tools/` package via `import sdar`** after installing it with `make -C tools` (installs to `~/include/sdar/`).
- **For Python post-processing, always verify the data file format matches the reader class.** AR and Hermite outputs have different column layouts; the wrong reader class produces garbled results *without* an error message — see "Python Data Analysis Tools" below.
- **After editing any Jupyter notebook (`*.ipynb`), verify every modified cell for integrity**: (a) code blocks are complete — no truncated functions, missing loops, or orphaned `try`/`except`; (b) the cell executes without syntax or runtime errors; (c) every import and variable is defined in the cell or an earlier one. Run each modified cell before considering the edit done.
- **SDAR is the reference integrator for PeTar.** When debugging PeTar close-encounter behaviour, replicate the subsystem in standalone SDAR first to isolate SDAR-level from P3T-level issues.

## Algorithm Overview

SDAR provides three core components:

| Component | Location | Purpose |
|-----------|----------|---------|
| **BinaryTree** | `src/Common/BinaryTree.h` | Kepler orbit ↔ Cartesian coordinate conversion; hierarchical binary tree construction |
| **AR** | `src/AR/symplectic_integrator.h` | Time-transformed explicit symplectic integrator with slowdown method |
| **Hermite** | `src/Hermite/hermite_integrator.h` | Hybrid 4th-order Hermite + AR method for global few-body systems |

### AR: Time-Transformed Symplectic Integrator

The AR method combines three techniques:
1. **Time transformation** (LogH or TTL): decouples the fixed symplectic step `ds` from the
   variable physical time step `dt`, enabling accurate integration of eccentric Kepler orbits.
2. **Explicit symplectic integrator**: conserves Hamiltonian and angular momentum for
   long-term secular evolution.
3. **Slowdown method**: for perturbed binaries, artificially slows the orbital motion
   so that one numerical orbit represents the secular effect of many physical orbits.

Key references:
- Wang, Nitadori & Makino (2020, MNRAS, 493, 3398) — SDAR algorithm
- Wang (2025, ApJ, 978, 65) — BlogH hybrid and LogH accuracy limits

### AR Method Variants

Binary names follow `ar.<method>[.ttl][.sd][.kdkpert][.cm][.mpfrc]`: `<method>` is the g-function form, `.ttl` marks the Time-Transformed Leapfrog implementation (absent = LogH), and the remaining suffixes are feature builds.

Key semantics:

- `logh` — g = `log(f(T) − f(−U)) / (T+U)`; best for isolated binaries, follows exact Kepler with phase error only.
- `ttl` — g = `1/|U|`; simpler and faster per step, larger energy error at high eccentricity.
- `blogh` / `btlogh` — g-function methods (innermost-pair product / tree-level product); exactly one per build.
- `sd` — tree-based hierarchical slowdown; **the default for hierarchical systems** (old `.sd.t`; `.sd.a` deprecated).
- `kdkpert`, `cm`, `mpfrc` — KDK perturbation splitting, CM-frame build, and MPFR high-precision mode.

Runtime `--g-func`: 0 = standard LogH, 1 = the method of this build, 2 = auto switch between 1 and 0 (rejected for `btlogh`).

Full variant → scenario mapping, the decision tree, and per-variant compile prerequisites: `assets/binary-scenario-map.md`.

### Hermite: Hybrid Hermite+AR

The Hermite integrator combines:
- **4th-order Hermite** for global integration of single particles
- **AR** for subsystems (close binaries, triples, etc.)

Groups are formed via `-r-group` (group detection radius) and `-r-neighbor-over-group`.
Particles within a group are integrated with AR; inter-group and single-particle forces
use the Hermite scheme.

## Build System

Directory layout, feature flags, and the standard `make -C sample/<dir>` commands are in [SDAR/AGENTS.md](../../../AGENTS.md) — do not duplicate them here.

Build details this skill relies on:

- **One executable per base method** (`logh`, `logh.ttl`, `blogh.ttl`, `btlogh.ttl`); the feature suffix (order `.sd.cm.mpfrc.kdkpert`) is composed from the `use_sd` / `use_cm` / `use_mpfrc` / `use_kdkpert` flags (`use_kdkpert` forces `use_sd` on). The default flags build `ar.logh.sd.cm`, `ar.logh.ttl.sd.cm`, `ar.blogh.ttl.sd.cm`, `ar.btlogh.ttl.sd.cm`.
- **When to rebuild**: request a different feature combination with `use_*` flags (e.g. `make use_sd=no use_cm=no use_mpfrc=yes`) rather than editing Makefile rules. Only a method combination beyond the four base methods needs a new rule.
- **Install path**: `~/bin` by default; change `INSTALL_PATH` in each `sample/*/Makefile`.
- **Debug and profile flags**: `AR_DEBUG`, `BINARY_DEBUG`, `AR_DEEP_DEBUG`, `USE_OMP`; `SDAR_TIME_MEASURE` (timing profile) is on by default and must match the Python reader's `time_measure` kwarg.

## Input File Format

SDAR uses space-separated ASCII input files. The format depends on the executable.

### AR / Hermite Particle Input

```
<N>                              # number of particles
# For each particle (one line per particle):
<mass> <x> <y> <z> <vx> <vy> <vz> <radius>
# Optional: binary tree structure line (for hierarchical initialization)
```

Example (`sample/input/triple.stable.lowm3`):
```
3
        0.01195722020130        -0.00222235770313         0.00754906025388         0.00178321684717        50.66596846912270       -94.28549675106700       -42.30077368644540      0.0
        0.01563023273488        -0.00222253668503         0.00754825038571         0.00178306650888       -37.78302017271370        72.72213076044780        32.98784924095770      0.0
        0.00011903271872        -0.00236411166715         0.00742145332081         0.00163449907602        -3.15413398262258         4.22288109760738        -0.19173108631062      0.0
1 0 3 0 1 2
```

The last line is the binary tree specification used by `keplertree` and AR/Hermite
for hierarchical initialization.

### Kepler Tool Input

**keplerorbit** reads particle pairs and converts to Kepler orbital parameters:
```
# Forward mode (default): two particle lines per pair → one Kepler orbit line
# Inverse mode (-i): one Kepler orbit line → two particle lines
```

**keplertree** reads a full particle set and constructs the binary tree:
```
# Without -i: N particle lines → hierarchical binary tree output
# With -i: N particle lines + tree structure line → Cartesian particle coordinates
```

### Units

SDAR executables accept particles in various unit systems via the `-u` flag:

| `-u` | Mass | Position | Velocity | Semi-major axis | Period |
|------|------|----------|----------|-----------------|--------|
| 0 | unscaled | unscaled | unscaled | unscaled | unscaled |
| 1 | Msun | AU | AU/yr | AU | yr |
| 2 | Msun | AU | km/s | AU | days |
| 3 | Msun | pc | km/s | pc | days |
| 4 | Msun | pc | pc/Myr | pc | Myr |

For N-body simulations, unit 0 (unscaled, G=1, total mass=1) or unit 4 (Msun/pc/Myr, G=0.0044983...) are most common.

## Workflow Patterns

### Pattern 1: AR Integration of a Few-Body System

```
1. Prepare input file (N particles + optional tree structure)
2. Select AR variant based on system properties
3. Run AR integrator
4. Post-process output with Python tools
```

**Selecting the right AR variant**: decision tree, system → binary mapping, and build prerequisites are in `assets/binary-scenario-map.md`. Suffix features beyond the default `.sd.cm` (e.g. `.mpfrc`, `.kdkpert`) are opt-in via `use_*` flags at build time.

**Key AR command-line options** (full list with defaults: `assets/minimal-question-sets.md`):

| Flag | Argument | Description | Default |
|------|----------|-------------|---------|
| `-t` | float | End physical time | 0.0 (must set) |
| `-n` | int | Number of integration steps (overrides `-t`) | 0 |
| `-s` | float | Step size ds (≤0 = auto) | 0.0 |
| `-e` | float | Relative energy error limit | 1e-10 |
| `-G` | float | Gravitational constant — see Gate 3; must match the unit system | 1.0 (Henon) |
| `-p` | string | Load parameters from a file | "" |
| `--g-func` | int | g-function mode: 0 = standard LogH, 1 = method of this build, 2 = auto switch 1↔0 (rejected for btlogh; g-func builds only) | 0 |
| `--break-check` | flag | Record hyperbolic-escape break events to stderr (record only — does not stop integration) | off |

Other tunables — `-o`, `-r`, `-k`, `-i`, `--ds-scale`, `--dt-min`, `--slowdown-ref`, `--slowdown-timescale-max` — keep their defaults unless there is a specific reason; confirm any change with the user.

### Pattern 2: Hermite+AR Hybrid Integration

```
1. Prepare input file (N particles)
2. Run Hermite integrator
3. Post-process output with Python tools
```

**Key Hermite command-line options** (full list with defaults: `assets/minimal-question-sets.md`):

| Flag | Argument | Description | Default |
|------|----------|-------------|---------|
| `-t` | float | End physical time | 1.0 |
| `-r-group` | float | Group detection radius | 1e-3 |
| `-r-neighbor-over-group` | float | Neighbor radius = factor × r_group | 2.0 |
| `-G` | float | Gravitational constant — see Gate 3; must match the unit system | 1.0 (Henon) |
| `-e` | float | Relative energy error limit for AR | 1e-10 |
| `--g-func` | int | g-function mode, `hermite.btlogh` build only: 0 = standard LogH (default, bit-identical to the plain build), 1 = BTLogH, 2 = auto (rejected for BTLogH) | 0 |

Other tunables — `-eta-4th`, `-eta-2nd`, `-eps`, `-o`, `-i`, `-k`, `--dt-min-power`, `--dt-max-power`, `--n-neighbor-max` — keep their defaults unless there is a specific reason.

### Pattern 3: Kepler Binary Tree Construction

```
1. Prepare particle data
2. Use keplertree to build binary tree
3. Use keplerorbit to convert between Kepler ↔ Cartesian
```

**keplertree options:**
- `-i`: Inverse mode — read particles + tree structure, output Cartesian positions+velocities
- `-n`: Number of binaries/pairs
- `-u`: Unit system (0-4, see units table above)

**keplerorbit options:**
- `-i`: Inverse mode — read Kepler orbit parameters, output Cartesian
- `-n`: Number of pairs
- `-u`: Unit system

### Pattern 4: Python Post-Processing

```
1. Ensure tools/ is installed: make -C tools
2. Import sdar module in Python
3. Read output files with appropriate reader classes
4. Analyze and visualize
```

## Python Data Analysis Tools (`tools/`)

### Installation

```bash
cd /home/lwang/code/SDAR/tools
make      # installs to /home/lwang/include/sdar/
```

To use in Python:
```python
import sys
sys.path.append('/home/lwang/include')
import sdar
# or: from sdar import SDARData, HermiteData, findPair, ...
```

### Module Map

| Module | Key Classes/Functions | Purpose |
|--------|-----------------------|---------|
| `sdar.base` | `DictNpArrayMix` | Base class for all data containers — dictionary+array hybrid |
| `sdar.particle` | `SimpleParticle`, `ParticleGroup` | Basic particle data (mass, pos, vel) and particle groups |
| `sdar.functions` | `vecDot`, `vecRot`, `cantorPairing`, `calcTrh`, `calcTcr`, `calcRocheLobeRadius`, `calcTGW` | Utility functions for vector ops, relaxation time, Roche lobe, GW merger time |
| `sdar.ar` | `SDARParticle`, `SDARData`, `SDARBinary`, `SDARProfile`, `SDARInfo`, `SDARInterruptBinary`, `SlowDownGroup` | AR output readers |
| `sdar.hermite` | `HermiteBaseParticle`, `HermiteParticle`, `HermiteEnergy`, `HermiteProfile`, `HermiteData` | Hermite output readers |
| `sdar.group` | `findPair`, `findMultiple` | Binary/multiple system detection from particle data |


### Critical: `N_particle`, `time_measure`, `slowdown` Parameters

**All SDAR reader classes require matching keyword arguments at construction time, before `loadtxt`.**
These determine the expected column layout and enable built-in column-count validation via `readArray`.

```python
# CORRECT — kwargs at construction, then loadtxt:
data = SDARData(N_particle=2, time_measure=True)
data.loadtxt("output.log", skiprows=1)

# WRONG — kwargs passed to loadtxt go to np.loadtxt and will error:
data = SDARData()
data.loadtxt("output.log", N_particle=2, skiprows=1)  # np.loadtxt rejects N_particle

# ALSO WRONG — construct from ndarray skips readArray's column check:
arr = np.loadtxt("output.log", skiprows=1)
data = SDARData(arr, N_particle=2)  # no ncols validation
```

**Required kwargs by output type:**

| Output | Kwargs | Data cols | Expected cols |
|--------|--------|-----------|---------------|
| plain AR (any N) | `N_particle=N, time_measure=True` | 47 (N=2), 56 (N=3) | ✅ match |
| AR .sd (tree slowdown) | `slowdown=True, time_measure=True, N_particle=N, N_sd=M` | 75 (N=3, M=2) | ✅ match |
| Hermite | `time_measure=True, N_particle=N` | 114 (N=3) | ⚠️ 108, close |

> **Note on `N_sd`**: `N_sd` = total number of binary (slowdown) pairs in the system.
>
> **AR (tree slowdown, `.sd`):** For a fully-connected hierarchical tree of N particles,
> the number of binary pairs is `N_sd = N_particle - 1`. For example:
> - Hierarchical triple (N=3): inner binary (p0-p1) + outer binary ((p0+p1)-p2) → `N_sd=2`.
> - Hierarchical quadruple (N=4): 3 binary pairs → `N_sd=3`.
>
> **Hermite:** The Hermite+AR hybrid integrator embeds SDAR groups internally.
> `N_sd` depends on how many SDAR groups exist in the initial conditions, which can
> vary by case. If a standalone `hermite` binary output contains slowdown columns,
> determine `N_sd` from the column count: total columns minus 63 (non-slowdown base)
> divided by 3 (columns per slowdown pair). When in doubt, start with `slowdown=True`
> and adjust `N_sd` to eliminate the column mismatch warning.
>
> All SDAR sample binaries include `-D SDAR_TIME_MEASURE` by default.

### Hermite Output: Pre-Filtering Required

The standalone `~/bin/hermite` writes diagnostic messages (`Large_energy_error:`, `Step hist:`,
parameter summaries) to the same stdout stream as the data table. Before reading with `HermiteData`:

1. **Filter**: keep only lines starting with a digit, `-`, or `.`
2. **Skip column title**: the first filtered line is the column title

```python
# Step 1: Filter diagnostic lines from raw log
lines = open("hermite.log").readlines()
data_lines = [l.rstrip() for l in lines if l.strip() and l.strip()[0] in '0123456789-.']
with open("hermite_clean.dat", "w") as f:
    f.write("\n".join(data_lines))

# Step 2: Read with HermiteData
arr = np.loadtxt("hermite_clean.dat", skiprows=1)  # skip column-title line
data = HermiteData(arr, N_particle=3)
```

### Reader Class Quick Reference

**AR Output → `SDARData`:**
```python
from sdar import SDARData

# Construction with kwargs, then loadtxt
data = SDARData(N_particle=2, time_measure=True)
data.loadtxt("output.log", skiprows=1)

# With slowdown
data = SDARData(slowdown=True, time_measure=True, N_particle=3)
data.loadtxt("triple.logh.sd.log", skiprows=1)

# Access data
data.time        # physical time at each output
data.de          # energy error
data.ekin, data.epot  # kinetic, potential energy
data.particles   # ParticleGroup containing member particles
data.profile     # SDARProfile: step counts, timing
data.info        # SDARInfo: ds, time_offset, r_break_crit
```

**Hermite Output → `HermiteData`** (requires pre-filtering, see above):
```python
from sdar import HermiteData

# After pre-filtering (see "Hermite Output: Pre-Filtering Required" above):
data = HermiteData(time_measure=True, N_particle=3)
data.loadtxt("hermite_clean.dat", skiprows=1)

# Access
data.time              # physical time at each output
data.time_offset       # time offset
data.energy_phy        # HermiteEnergy: .de, .ekin, .epot, .epert, .de_change
data.energy_sd         # HermiteEnergy for slowdown component
data.profile           # HermiteProfile: .h4_step_single, .ar_step, .break_group, ...
data.particles         # ParticleGroup
```

**Particle Groups:**
```python
# Access individual particles in a group
for i in range(data.particles.n[0]):
    p = data.particles['p%d' % i]
    print(p.mass, p.pos, p.vel)

# Center of mass
cm = data.particles.cm
```

### Finding Binaries from Particle Data

```python
from sdar import SimpleParticle, findPair, findMultiple

# Read particle data into SimpleParticle
p = SimpleParticle()
p.mass = ...  # assign arrays
p.pos = ...
p.vel = ...

# Find binaries using KDTree (returns kdt, singles, binary)
G = 0.00449830997959438  # Msun, pc, Myr
rmax = 0.01  # pc
kdt, single, binary = findPair(p, G, rmax, use_kdtree=True)

# Without KDTree (PeTar status-column method, returns singles, binary):
# single, binary = findPair(p, G, rmax, use_kdtree=False)

# Find triples/quadruples from single+binary
triple, quad = findMultiple(single, binary, G, rmax)
```

### Important Unit Constants

```python
G_MSUN_PC_MYR = 0.00449830997959438   # Msun, pc, Myr
G_HENON = 1.0                          # Henon/N-body units
```

Always use these constants from `sdar.particle` rather than hardcoding numeric values.

## Parameter File (`-p`) Workflow

Both AR and Hermite support loading parameters from a file via `-p <filename>`.

**To save parameters from a run:**
Add `--save-param <filename>` to the command line.

**To reuse parameters:**
```bash
~/bin/ar.logh.sd -p saved.par input.dat
```

Command-line arguments override file values, so `-p` + `-t 2.0` uses `t=2.0` regardless of what `saved.par` contains.

## Common Pitfalls

1. **Step size too large for eccentric orbits.** The default auto step size may be insufficient for high-eccentricity binaries — use `-s <smaller_value>` or tighten `-e`.
2. **Unit system confusion in the Kepler tools.** `-u` changes the interpretation of semi-major axis, period, and velocities; with unit 4 (Msun/pc/Myr), `-G` must be `0.00449830997959438`.
3. **`findPair` returns a different tuple shape depending on `use_kdtree`** — `(kdt, singles, binary)` for `True`, `(singles, binary)` for `False`.
4. **`HermiteData` exposes energies via `data.energy_phy`**, not `data.energy`; the slowdown component is `data.energy_sd`.

The constructor-kwarg, pre-filtering, and column-count pitfalls are covered above and in `assets/data-readback-patterns.md`; the G-mismatch and output-redirection rules are Gates 2–3.

## Preferred References in This Repository

- `README.md` — user guide and algorithm introduction
- `docs/doc.h`, `docs/html/index.html` — algorithm derivation and full Doxygen documentation
- `sample/input/*.sh` — runnable examples (`binary_logh.sh`, `triple_logh_sd.sh`, `fewbody_hermite.sh`, `triple_compare_methods.sh`, `build_kepler_tree.sh`, `triple.stable.lowm3.sh`)
- `sample/{AR,Hermite,Kepler}/Makefile` — build targets and compile flags
- `sample/data_analysis.ipynb` — full analysis workflow examples
- Skill assets — `assets/binary-scenario-map.md`, `assets/data-readback-patterns.md`, `assets/minimal-question-sets.md`, `assets/lessons-learned.md`, `assets/HANDOFF.md`

## Reference Documents (Must Read)

The following asset files contain critical information not inlined in this document.
When a task falls into the corresponding category, **read the file explicitly** with `read_file` before proceeding — do not guess or rely on memory.

| When to read | File | What it contains |
|-------------|------|------------------|
| **Any SDAR scenario** (after scenario is identified) | `assets/minimal-question-sets.md` | Per-scenario required-ask lists, including the "do not ask if already known" rule |
| **Selecting AR variant** | `assets/binary-scenario-map.md` | Decision tree, binary-to-scenario mapping, compile flag requirements |
| **Python data analysis** (before writing any analysis code) | `assets/data-readback-patterns.md` | **MUST READ before writing any analysis code.** Contains verified readback patterns for SDARData, HermiteData, SimpleParticle, SDARBinary, and binary/multiple detection. |
| **Modifying integrator/time-sync code, A/B step-count comparisons, or adding defensive guards** | `assets/lessons-learned.md` | **MUST READ before touching `integrateToTime`/ds control or comparing runs.** Verified pitfalls: paired-evaluator consistency, same-binary A/B requirement, sorted-cck semantics, ds-floor gauge (float-resolution, not time_error). |
| **Continuing the BTLogH → Hermite/PeTar integration (or any cross-project handoff task)** | `assets/HANDOFF.md` | Current status, verified facts (compile-clean, LogH-inert flag, kappa-safety of the product form), remaining plumbing, risk list, test data paths, and the validation ladder. |
| **Full analysis workflow** | `sample/data_analysis.ipynb` | Complete Jupyter notebook demonstrating all analysis patterns with runnable code

Technical background (key algorithms):

- SDAR integrator: slow-down + time-transformed symplectic method for few-body systems
  (Wang et al. 2020, MNRAS, 493, 3398) — [https://doi.org/10.1093/mnras/staa480](https://doi.org/10.1093/mnras/staa480)
- BlogH hybrid method and LogH accuracy limits for hierarchical triples
  (Wang 2025, ApJ, 978, 65) — [https://doi.org/10.3847/1538-4357/ad98f3](https://doi.org/10.3847/1538-4357/ad98f3)
- PeTar code description (SDAR used as close-encounter integrator):
  Wang et al. 2020, MNRAS, 497, 536 — [https://doi.org/10.1093/mnras/staa1915](https://doi.org/10.1093/mnras/staa1915)

## Relationship to PeTar

SDAR is the few-body integrator embedded inside PeTar: when PeTar detects a close encounter or bound subsystem, it hands the particles to the AR/Hermite library in `src/`.

- **SDAR source (`src/`) is shared** — changes affect both the standalone executables and PeTar's internal SDAR integration.
- **PeTar's SDAR behaviour is reproducible standalone.** If a PeTar simulation shows unexpected close-binary evolution, extract the subsystem and replicate it here first (Gate 3).
- **Group detection maps across the boundary.** PeTar derives `--r-group` / `--r-search-group` from `r_in`; standalone SDAR takes `-r-group` directly. Check the equivalent parameter on the other side when debugging group formation.
- **Version consistency:** PeTar's `VERSION` encodes `PeTar_VERSION_SDAR_VERSION`. Comparing PeTar and standalone SDAR results requires matching versions.

## Scope Notes

- This skill is the authority for **standalone** SDAR usage (few-body systems, typically N ≤ 100) and is self-sufficient. For full N-body cluster simulations where SDAR is the close-encounter solver, use the [PeTar skill](../../../../PeTar/.github/skills/petar-nbody-simulation/SKILL.md).
- `src/`, build flags, and the `sdar` Python package internals are owned by [SDAR/AGENTS.md](../../../AGENTS.md) and this skill's assets — the PeTar skill does not restate them.
- The `BinaryTree` component is used both standalone (Kepler tools) and internally by AR/Hermite; standalone use is for constructing and analysing hierarchical systems.
- SDAR has no stellar evolution, external potentials, or MPI parallelism — those exist only in PeTar.
- High-precision (MPFRC) runs require `libmpfr` and `libgmp` at link time and `-D USE_MPFRC` at compile time.
