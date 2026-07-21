#!/bin/bash
set -e

# AR LogH integration of an isolated binary
# Usage: ./binary_logh.sh [working_dir]
#   Creates a 2-body input file and runs ar.logh on it.
#   Runtime: < 1 second.

WORKDIR="${1:-/tmp/sdar_binary_logh}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Check binary exists
if ! command -v ~/bin/ar.logh >/dev/null 2>&1; then
    echo "Error: ~/bin/ar.logh not found. Run 'make -C sample/AR && make -C sample/AR install' first." >&2
    exit 1
fi

# Create a 2-body input file (N + particle lines: mass x y z vx vy vz radius)
cat > binary.dat << 'EOF'
2
 0.6  0.0  0.0  0.0   0.0  2.0  0.0  0.0
 0.4  0.02 0.0  0.0   0.0 -3.0  0.0  0.0
EOF

echo "# $(date): ~/bin/ar.logh -t 5.0 -o 1.0 -G 1.0 binary.dat" >> commands.log

# Run AR LogH integrator
# -t 5.0:   integrate for 5.0 time units (Henon)
# -o 1.0:   output snapshots every 1.0 time units
# -G 1.0:   gravitational constant (Henon units; use 0.00449830997959438 for Msun/pc/Myr)
~/bin/ar.logh -t 5.0 -o 1.0 -G 1.0 binary.dat > binary_logh.log

echo "Done. Output: $WORKDIR/binary_logh.log"
