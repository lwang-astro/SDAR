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

- **Always `cd` to a dedicated working directory before running any SDAR integrator.**
  Do not run inside `sample/` or the repository root — these locations are for source code and build artifacts only.
  If the user does not specify a working directory, ask for one — do not decide the path yourself.
  If the user has no preference, suggest `~/sdar_run/<descriptive-name>` and confirm before creating.
- **Do not run whatever binary happens to be in `~/bin/` without confirming it matches the scenario.**
  Check what the user needs first (AR or Hermite, LogH or TTL, with/without slowdown, etc.), then select the correct binary.
- **Validate binary availability** before composing commands: `command -v <binary>` or `ls <path-to-build>`.
  If the binary is not built, guide the user through `make` in the appropriate `sample/` subdirectory.
- **Before any solver execution, present a structured parameter summary and wait for explicit confirmation.**
  The summary must include:
  1. Working directory and input file
  2. Full command line (binary, flags, output redirect)
  3. Expected runtime hint (based on N-body count, end time, and step size)
  4. Expected output files
  5. A clear prompt: *"Proceed? (y/n)"* — stop and wait for user response.
- **Redirect stdout of every major command to a file.** Use consistent naming:
  - `ar.logh ... >ar_logh.log` for AR LogH runs
  - `ar.ttl ... >ar_ttl.log` for AR TTL runs
  - `hermite ... >hermite.log` for Hermite runs
  - `keplertree ... >keplertree.log` for tree construction
  - `keplerorbit ... >keplerorbit.log` for orbit conversion
- **Record every command to a `commands.log` file in the working directory.**
  ```bash
  echo "# $(date): ~/bin/ar.logh.sd.t -t 1.0 -s 0.01 input.dat" >> commands.log
  ```
- **Never mutate input files.** SDAR input files are small ASCII tables. Always keep the original;
  if parameter changes are needed, copy the input file to the working directory first.
- **For Python post-processing, always verify the data file format matches the reader class.**
  AR output and Hermite output have different column layouts — using the wrong reader class
  produces garbled results without an error message. See "Python Data Analysis Tools" below.
- **After editing any Jupyter notebook (`*.ipynb`), verify every modified cell for integrity.**
  Check: (a) code blocks are complete — no truncated functions, missing loops, or orphaned
  try/except blocks; (b) the cell executes without syntax or runtime errors; (c) all imports
  and variables used in the cell are defined in previous cells or within the cell.
  Run each modified cell and confirm it produces expected output before considering the edit done.
- **Use the `tools/` module via `import sdar` after installing the tools to the Python path.**
  The tools install to `/home/lwang/include/sdar/` by the `tools/Makefile` (user-specific path).
  If not installed, run `make -C tools` first.
- **For gravitational constant, use `G_MSUN_PC_MYR = 0.00449830997959438` for Msun/pc/Myr units
  and `G_HENON = 1.0` for Henon/N-body units.** These constants are defined in `tools/particle.py`
  and must be consistent between the C++ integrator and Python analysis.
- **SDAR is the reference integrator for PeTar.** When debugging PeTar close-encounter behavior,
  replicate the subsystem in standalone SDAR first to isolate SDAR-level vs P3T-level issues.
- **For hierarchical systems, slowdown (`sd`) binaries are preferred.** The slowdown method
  dramatically accelerates weakly perturbed inner binaries. Use `ar.logh.sd.t` or `ar.ttl.sd.t`
  over plain `ar.logh`/`ar.ttl` for hierarchical systems.

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

| Suffix | Method | Time Transformation | Notes |
|--------|--------|---------------------|-------|
| `logh` | LogH | `g = log(f(T) - f(-U)) / (T+U)` | Best for isolated binaries. Numerical trajectory follows exact Kepler with phase error only. |
| `ttl` | TTL (Time-Transformed Leapfrog) | `g = 1/|U|` | Simpler, faster per step but larger energy error for high eccentricity. |
| `sd.t` | Slowdown, tree | With slowdown for inner+outer binaries | Hierarchical slowdown across binary tree levels. **Current default.** |
| `kdk.pert` | KDK with perturbation | LogH with KDK splitting | Alternative for weakly perturbed systems. |
| `mulpot` | Multi-potential | Product-of-pair-potentials time function | For systems with multiple binaries. |
| `cm` | Center-of-mass | With CM motion tracking | Useful when system drifts. |
| `mpfrc` | MPFRC | High-precision (mpfr::mpreal) | Arbitrary-precision mode. Requires `-lmpfr -lgmp`. |

### Hermite: Hybrid Hermite+AR

The Hermite integrator combines:
- **4th-order Hermite** for global integration of single particles
- **AR** for subsystems (close binaries, triples, etc.)

Groups are formed via `-r-group` (group detection radius) and `-r-neighbor-over-group`.
Particles within a group are integrated with AR; inter-group and single-particle forces
use the Hermite scheme.

## Build System

### Directory Layout

```
SDAR/
├── src/
│   ├── Common/      # Float.h, List.h, ParticleGroup.h, BinaryTree.h, io.h
│   ├── AR/          # symplectic_integrator.h, information.h
│   └── Hermite/     # hermite_integrator.h
├── sample/
│   ├── AR/          # AR standalone executables + Makefile
│   ├── Hermite/     # Hermite standalone executable + Makefile
│   ├── Kepler/      # keplerorbit, keplertree + Makefile
│   ├── input/       # Sample input files
│   └── test/        # Unit tests
├── tools/           # Python post-processing modules
└── docs/            # Doxygen documentation
```

### Compilation

Source files in `sample/` are the standalone executables. Each subdirectory has its own Makefile.

**To compile all AR variants:**
```bash
cd sample/AR
make              # build all targets to ./build/
make install      # install to ~/bin/ (or modify INSTALL_PATH in Makefile)
```

**AR build targets** (from `sample/AR/Makefile`):
```
ar.logh ar.logh.sd.t ar.ttl ar.ttl.sd.t ar.ttl.sd.t.mulpot
ar.ttl.sd.t.mulpot.cm ar.ttl.sd.t.maxpot.cm ar.ttl.sd.t.addpot.cm
ar.logh.sd.t.kdk.pert ar.ttl.sd.t.kdk.pert
```
Plus MPFRC variants: `ar.ttl.sd.t.cm ar.logh.mpfrc ar.ttl.mpfrc` etc.

**To compile Hermite:**
```bash
cd sample/Hermite
make              # builds hermite, hermite.mpfrc, hermite.kdk.pert
make install
```

**To compile Kepler tools:**
```bash
cd sample/Kepler
make              # builds keplerorbit, keplertree
make install
```

**Compile flags to be aware of:**
- `-D AR_SLOWDOWN_TREE` — enables hierarchical slowdown (required for `sd.t` variants)
- `-D AR_SLOWDOWN_TIMESCALE` — enables timescale-based slowdown control
- `-D AR_TTL` — use TTL time transformation instead of LogH
- `-D AR_KDK_PERT` — use KDK perturbation splitting
- `-D USE_MPFRC` — arbitrary-precision mode (requires `-lmpfr -lgmp`)
- `-D AR_DEBUG`, `-D BINARY_DEBUG`, `-D AR_DEEP_DEBUG` — debug flags
- `-D USE_OMP` — OpenMP parallelization
- `-D SDAR_TIME_MEASURE` — enable timing profile output (enabled by default)

**When to rebuild:** If the user requests a variant not listed in `sample/AR/Makefile`'s `TARGET`,
check whether the required `-D` flags are already enabled in the Makefile. If a new combination
is needed, advise adding a new target rule following the existing pattern.

**Install path:** Default is `~/bin`. To change, modify `INSTALL_PATH` in each `sample/*/Makefile`.

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

**Selecting the right AR variant:**

| System | Recommended Binary |
|--------|-------------------|
| Isolated binary | `ar.logh` |
| Isolated binary (faster, less accurate) | `ar.ttl` |
| Hierarchical triple/quadruple | `ar.logh.sd.t` or `ar.ttl.sd.t` |
| Weakly perturbed binary | `ar.logh.sd.t.kdk.pert` |
| System with multiple inner binaries | `ar.ttl.sd.t.mulpot` |
| High-precision requirement | `ar.logh.mpfrc` |

**Key AR command-line options:**

| Flag | Argument | Description | Default |
|------|----------|-------------|---------|
| `-t` | float | End physical time | 0.0 (must set) |
| `-n` | int | Number of integration steps (overrides `-t`) | 0 |
| `-s` | float | Step size ds (≤0 = auto) | 0.0 |
| `-e` | float | Relative energy error limit | 1e-10 |
| `-o` | float | Output time interval | 0.0 |
| `-r` | float | Distance criterion for stability check | 1e-3 |
| `-G` | float | Gravitational constant | 1.0 |
| `-k` | int | Symplectic integrator order (even; negative = Yoshida 2nd) | -6 |
| `-i` | int | Interrupt detection: 0=off, 1=modify orbits, 2=record only | 0 |
| `-p` | string | Load parameters from file | "" |
| `--ds-scale` | float | Step size scaling factor | 1.0 |
| `--dt-min` | float | Minimum physical time step | 1e-13 |
| `--slowdown-ref` | float | Slowdown perturbation ratio reference | 1e-6 |
| `--slowdown-timescale-max` | float | Max timescale for slowdown factor | time-end |

### Pattern 2: Hermite+AR Hybrid Integration

```
1. Prepare input file (N particles)
2. Run Hermite integrator
3. Post-process output with Python tools
```

**Key Hermite command-line options:**

| Flag | Argument | Description | Default |
|------|----------|-------------|---------|
| `-t` | float | End physical time | 1.0 |
| `-r-group` | float | Group detection radius | 1e-3 |
| `-r-neighbor-over-group` | float | Neighbor radius = factor × r_group | 2.0 |
| `-eta-4th` | float | Time step coefficient for 4th order | 0.1 |
| `-eta-2nd` | float | Time step coefficient for 2nd order | 0.001 |
| `-eps` | float | Softening parameter | 0.0 |
| `-G` | float | Gravitational constant | 1.0 |
| `-o` | int | Output interval (power index of 0.5) | 2 |
| `-e` | float | Relative energy error limit for AR | 1e-10 |
| `-i` | int | Interrupt detection option | 0 |
| `-k` | int | AR symplectic order | -6 |
| `--dt-min-power` | int | Power index for minimum Hermite step | 40 |
| `--dt-max-power` | int | Power index for maximum Hermite step | 2 |
| `--n-neighbor-max` | int | Max neighbors for group (-1 = same as N) | -1 |

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
| AR .sd.t (tree slowdown) | `slowdown=True, time_measure=True, N_particle=N, N_sd=M` | 75 (N=3, M=2) | ✅ match |
| Hermite | `time_measure=True, N_particle=N` | 114 (N=3) | ⚠️ 108, close |

> **Note on `N_sd`**: `N_sd` = total number of binary (slowdown) pairs in the system.
>
> **AR (tree slowdown, `.sd.t`):** For a fully-connected hierarchical tree of N particles,
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
data.loadtxt("triple.logh.sd.t.log", skiprows=1)

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
~/bin/ar.logh.sd.t -p saved.par input.dat
```

Command-line arguments override file values, so `-p` + `-t 2.0` uses `t=2.0` regardless of what `saved.par` contains.

## Common Pitfalls

1. **Missing `N_particle` or `time_measure` at construction.** Always pass `N_particle=N, time_measure=True`
   to `SDARData()` and `HermiteData()` constructors. Without them, `particles.n` is 0 and
   column-count validation fails. `loadtxt` does NOT accept these kwargs — they go to `np.loadtxt`
   which rejects them.

2. **Hermite output needs pre-filtering before `HermiteData`.** The raw log contains
   diagnostic messages mixed with data. Filter to numeric-starting lines first,
   then skip the column-title line, then pass to `HermiteData(arr, N_particle=N)`.

3. **Column counts vary by AR variant.** With N=3 particles: plain AR = 56 columns,
   `.sd.t` = 75 columns. Pass `slowdown=True` only for `.sd.t`; `N_sd` must match
   the number of binary pairs. For a fully-connected hierarchical tree,
   `N_sd = N_particle - 1`. For Hermite, `N_sd` varies by initial SDAR group count.
   When the column mismatch warning appears, adjust `N_sd` to resolve it.

4. **`findPair` signature depends on `use_kdtree`.** With `use_kdtree=True`, the return is
   `(kdt, singles, binary)` (3 values). With `use_kdtree=False`, it's `(singles, binary)` (2 values).

5. **HermiteData uses `data.energy_phy` not `data.energy`.** Access Hermite energy via
   `data.energy_phy.de`, `data.energy_phy.ekin`, etc. — not `data.energy`.

6. **G mismatch between C++ and Python.** The C++ default `-G 1.0` (Henon units) differs from
   `G_MSUN_PC_MYR = 0.00449830997959438`. Confirm the `-G` value used in the C++ run and
   use the matching constant in Python.

7. **Step size too large for eccentric orbits.** For high-eccentricity binaries, the default
   auto step size may be insufficient. Use `-s <smaller_value>` or `-e <tighter_tolerance>`.

8. **Lost output from terminal.** Always redirect stdout. SDAR can produce many output lines;
   without redirection, important diagnostics are lost.

9. **Unit system confusion.** The `-u` flag in Kepler tools changes interpretation of
   semi-major axis, period, and velocities. When unit 4 (Msun/pc/Myr) is used, G must be
   `0.00449830997959438`.

## Preferred References in This Repository

Use these as primary references for SDAR workflows:

- `README.md` — user guide and algorithm introduction
- `docs/doc.h` — Doxygen mainpage with complete algorithmic derivation
- `docs/html/index.html` — full Doxygen HTML documentation
- `sample/input/` — example input files
- `sample/input/triple.stable.lowm3.sh` — original AR method comparison script
- `sample/input/binary_logh.sh` — isolated binary AR LogH integration
- `sample/input/triple_logh_sd.sh` — hierarchical triple AR LogH + slowdown
- `sample/input/fewbody_hermite.sh` — Hermite+AR hybrid integration
- `sample/input/triple_compare_methods.sh` — multi-method comparison
- `sample/input/build_kepler_tree.sh` — Kepler binary tree construction
- `sample/AR/Makefile` — AR build targets and compile flags
- `sample/Hermite/Makefile` — Hermite build targets
- `sample/Kepler/Makefile` — Kepler tool build targets
- `tools/__init__.py` — Python module entry point
- `.github/skills/sdar-fewbody-integration/assets/binary-scenario-map.md` — AR/Hermite variant selection guide
- `.github/skills/sdar-fewbody-integration/assets/data-readback-patterns.md` — Python data readback patterns
- `.github/skills/sdar-fewbody-integration/assets/minimal-question-sets.md` — per-scenario minimal required inputs
- `sample/data_analysis.ipynb` — Jupyter notebook with full analysis workflow examples

## Reference Documents (Must Read)

The following asset files contain critical information not inlined in this document.
When a task falls into the corresponding category, **read the file explicitly** with `read_file` before proceeding — do not guess or rely on memory.

| When to read | File | What it contains |
|-------------|------|------------------|
| **Any SDAR scenario** (after scenario is identified) | `assets/minimal-question-sets.md` | Per-scenario required-ask lists, including the "do not ask if already known" rule |
| **Selecting AR variant** | `assets/binary-scenario-map.md` | Decision tree, binary-to-scenario mapping, compile flag requirements |
| **Python data analysis** (before writing any analysis code) | `assets/data-readback-patterns.md` | **MUST READ before writing any analysis code.** Contains verified readback patterns for SDARData, HermiteData, SimpleParticle, SDARBinary, and binary/multiple detection. |
| **Full analysis workflow** | `sample/data_analysis.ipynb` | Complete Jupyter notebook demonstrating all analysis patterns with runnable code

Technical background (key algorithms):

- SDAR integrator: slow-down + time-transformed symplectic method for few-body systems
  (Wang et al. 2020, MNRAS, 493, 3398) — [https://doi.org/10.1093/mnras/staa480](https://doi.org/10.1093/mnras/staa480)
- BlogH hybrid method and LogH accuracy limits for hierarchical triples
  (Wang 2025, ApJ, 978, 65) — [https://doi.org/10.3847/1538-4357/ad98f3](https://doi.org/10.3847/1538-4357/ad98f3)
- PeTar code description (SDAR used as close-encounter integrator):
  Wang et al. 2020, MNRAS, 497, 536 — [https://doi.org/10.1093/mnras/staa1915](https://doi.org/10.1093/mnras/staa1915)

## Relationship to PeTar

SDAR is the few-body integrator embedded inside PeTar. When PeTar detects a close encounter
or bound subsystem, it hands off the particles to SDAR (via the AR/Hermite library in `src/`).
Key relationship rules:

- **PeTar's SDAR behavior can be tested standalone.** If a PeTar simulation shows unexpected
  close-binary evolution, extract the subsystem and replicate it with standalone SDAR.
- **SDAR source (`src/`) is shared.** Changes to `SDAR/src/` affect both standalone SDAR
  executables and PeTar's internal SDAR integration.
- **PeTar's `--r-group`, `--r-search-group` map to SDAR's group detection.** When debugging
  PeTar SDAR group formation issues, check the equivalent SDAR standalone parameters.
- **Version consistency:** PeTar's VERSION file encodes `PeTar_VERSION_SDAR_VERSION`.
  When comparing PeTar and standalone SDAR results, ensure versions match.

## Scope Notes

- This skill focuses on standalone SDAR usage (few-body systems, typically N ≤ 100).
  For full N-body cluster simulations where SDAR is the close-encounter solver,
  use the PeTar skill instead.
- The `BinaryTree` component is used both standalone (Kepler tools) and internally
  by AR/Hermite. Its standalone use is for constructing/analyzing hierarchical systems.
- SDAR does not have stellar evolution, external potentials, or MPI parallelism —
  these features exist only in PeTar.
- For high-precision (MPFRC) runs, ensure `libmpfr` and `libgmp` are installed and
  the `-D USE_MPFRC` flag is set during compilation.
