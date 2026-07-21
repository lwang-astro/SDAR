# Python Data Readback Patterns (SDAR)

This file documents the verified patterns for reading SDAR output files from Python (`import sdar`).
These patterns have been derived from the `tools/` module source code and should be used when
an agent needs to read or analyze SDAR integration outputs.

**Reference implementation**: `tools/ar.py`, `tools/hermite.py`, `tools/particle.py`

## Prerequisites

```python
import sys
sys.path.insert(0, "/home/lwang/include")  # or wherever tools/Makefile installs
import sdar
import numpy as np
```

## Pattern 1: Read AR Output (LogH/TTL)

**CRITICAL**: 
1. All kwargs (`N_particle`, `time_measure`, `slowdown`) must be passed at **construction** time, NOT in `loadtxt`. `loadtxt` passes its kwargs to `np.loadtxt` which rejects them.
2. All SDAR sample binaries include `-D SDAR_TIME_MEASURE` by default, so **always** pass `time_measure=True`.

```python
from sdar import SDARData

# Correct pattern: construct with kwargs, then loadtxt
data = SDARData(N_particle=2, time_measure=True)
data.loadtxt("output.log", skiprows=1)

# With slowdown data (.sd.t variants only)
data = SDARData(slowdown=True, time_measure=True, N_particle=3)
data.loadtxt("triple.logh.sd.t.log", skiprows=1)

# .sd.a variants have fewer slowdown columns — use slowdown=False
data = SDARData(slowdown=False, time_measure=True, N_particle=3)
data.loadtxt("triple.logh.sd.a.log", skiprows=1)
data = SDARData(slowdown=True, time_measure=True)
arr = np.loadtxt("output.dat")
data = SDARData(arr)

# Access data
print("Times:", data.time)
print("Energy error:", data.de)
print("Kinetic energy:", data.ekin)
print("Potential energy:", data.epot)
print("Extended Hamiltonian:", data.H)

# Access particle data
particles = data.particles
n_particles = particles.n
for i in range(particles.n[0]):
    p = particles['p%d' % i]
    print(f"  Particle {i}: mass={p.mass}, pos={p.pos}, vel={p.vel}")

# CM particle
cm = particles.cm
print(f"CM: mass={cm.mass}, pos={cm.pos}, vel={cm.vel}")

# Access profile data
print("AR steps:", data.profile.n_step_sum)
print("Time sync steps:", data.profile.n_step_tsyn_sum)

# Access integrator info
print("Step ds:", data.info.ds)
print("Time offset:", data.info.time_offset)

# Access slowdown data (if slowdown=True)
if hasattr(data, 'de_sd'):
    print("Slowdown energy error:", data.de_sd)
    print("Slowdown factor:", data.sd.slowdown_factor)
```

## Pattern 2: Read Hermite Output

```python
from sdar import HermiteData

# Correct pattern: construct with kwargs, then loadtxt
data = HermiteData(time_measure=True, N_particle=3)
data.loadtxt("hermite_clean.dat", skiprows=1)

# Access energy data
print("Energy error:", data.energy_phy.de)
print("Kinetic:", data.energy_phy.ekin)
print("Potential:", data.energy_phy.epot)

# Access Hermite profile
print("H4 single steps:", data.profile.h4_step_single)
print("AR steps:", data.profile.ar_step)
```

## Pattern 3: Read Particle Data with Binary Tree

```python
from sdar import SimpleParticle

# Read particle data from text
p = SimpleParticle()
arr = np.loadtxt("particles.dat", skiprows=1)  # skip N line
p = SimpleParticle(arr)  # reads mass, pos, vel

# With MPFRC high-precision position
p = SimpleParticle(use_mpfrc=True)
arr = np.loadtxt("particles_mpfrc.dat", skiprows=1)
p = SimpleParticle(arr)

# Calculate derived quantities
p.calcR2()    # adds r2 member (distance squared)
p.calcEkin()  # adds ekin member (kinetic energy)

# Correct center-of-mass
cm_pos = np.array([0.0, 0.0, 0.0])
cm_vel = np.array([0.0, 0.0, 0.0])
p.correctCenter(cm_pos, cm_vel)
```

## Pattern 4: Find Binaries and Multiples

```python
from sdar import SimpleParticle, findPair, findMultiple

# Read particle data
p = SimpleParticle()
arr = np.loadtxt("particles.dat", skiprows=1)
p = SimpleParticle(arr)

G = 0.00449830997959438  # Msun, pc, Myr
rmax = 0.01  # pc

# Method 1: Using KDTree (slower, finds all pairs)
kdt, single, binary = findPair(p, G, rmax, use_kdtree=True)
print(f"Found {binary.size} binaries, {single.size} singles")

# Method 2: Using status column from PeTar (faster)
# Requires particle data to have a 'status' column (PeTar output)
single, binary = findPair(p, G, rmax, use_kdtree=False)

# Find triples and quadruples from singles + binaries
triple, quad = findMultiple(single, binary, G, rmax)

# Access binary parameters
if binary.size > 0:
    print("Semi-major axes:", binary.semi)
    print("Eccentricities:", binary.ecc)
    print("Periods:", binary.period)
```

## Pattern 5: Read SDAR Binary Parameters (from Interrupt Output)

```python
from sdar import SDARBinary, SDARInterruptBinary

# Read binary parameters from interrupt output
bin_data = SDARBinary()
arr = np.loadtxt("binary_params.dat")
bin_data = SDARBinary(arr)

print("Semi:", bin_data.semi)
print("Ecc:", bin_data.ecc)
print("Inclination:", bin_data.incline)
print("Omega (ascending node):", bin_data.rot_horizon)
print("omega (periapsis):", bin_data.rot_self)
print("Ecc anomaly:", bin_data.ecca)
print("Component masses:", bin_data.m1, bin_data.m2)
print("Relative distance:", bin_data.rrel)
print("Angular momentum:", bin_data.am)
print("Stability factor:", bin_data.stab)
```

## Unit Conversion

SDAR's Python tools do not include built-in unit conversion (unlike PeTar's `petar` module).
Use `astropy.units` for all conversions:

```python
import astropy.units as u

# SDAR gravitational constants
G_MSUN_PC_MYR = 0.00449830997959438   # Msun, pc, Myr
G_HENON = 1.0                          # Henon/N-body units

# pc → AU
semi_au = (binary.semi * u.pc).to(u.AU).value

# pc/Myr → km/s
vel_kms = (velocity * u.pc / u.Myr).to(u.km / u.s).value

# Time conversion
time_myr = time_henon  # if using Henon units, check scaling
```

## Common Pitfalls

1. **Missing `N_particle` argument.** Always pass `N_particle=N` to `SDARData` and `HermiteData`.
   Without it, `particles.n` will be 0 and particle member data (`p0`, `p1`, etc.) will not exist.
   This is the #1 cause of "data read successfully but no particles visible."

2. **Hermite output needs pre-filtering.** The raw `hermite.log` contains diagnostic
   messages (`Large_energy_error:`, `Step hist:`) mixed with data rows. Filter to
   numeric-starting lines first, then skip the column-title line, then pass to
   `HermiteData(arr, N_particle=N)`. The data rows consistently have 114 columns for N=3.

3. **Column counts vary by AR variant.** Plain `ar.logh` produces 56 columns, while
   `ar.logh.sd.t` produces 63 columns (extra slowdown data). Always pass `slowdown=True/False`
   to match the producing binary. Mismatch causes column alignment errors.

4. **Using SDARData on Hermite output or vice versa.** The column layouts are different.
   Always check which binary produced the output. Standalone Hermite output is NOT compatible
   with any SDAR Python reader class.

5. **`skiprows` depends on output variant.** AR output has a single column-title header line.
   Use the dynamic detection pattern shown in Pattern 1 rather than hardcoding a number.

6. **Mismatched G between C++ and Python.** If the C++ run used `-G 1.0`, use `G=1.0` in Python.
   If it used `-G 0.00449830997959438`, use `G_MSUN_PC_MYR`.

7. **ParticleGroup access.** After reading SDARData, particles are in `data.particles`.
   Individual members are `data.particles['p0']`, `data.particles['p1']`, etc.
   The count is in `data.particles.n` (a 1D array, use `data.particles.n[0]`).
   If `n` is 0, you forgot `N_particle`.
