#!/bin/bash
set -e

# Hermite+AR hybrid integration of a hierarchical triple
# Usage: ./fewbody_hermite.sh [working_dir]
#   Uses sample/input/triple.stable.lowm3 as input.
#   Runtime: < 0.1 second.
#   Note: Hermite standalone output is a human-readable log, NOT a fixed-column
#   table. Use grep to check for final time rather than HermiteData.read().

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT="$SCRIPT_DIR/triple.stable.lowm3"
WORKDIR="${1:-/tmp/sdar_hermite}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Check binary exists
if ! command -v ~/bin/hermite >/dev/null 2>&1; then
    echo "Error: ~/bin/hermite not found. Run 'make -C sample/Hermite && make -C sample/Hermite install' first." >&2
    exit 1
fi

echo "# $(date): ~/bin/hermite -t 1.0 -o 2 -G 1.0 $INPUT" >> commands.log

# Run Hermite+AR hybrid integrator
# -t 1.0: integrate for 1.0 time units (Henon)
# -o 2:   output interval as power index of 0.5 (2 → dt_out = 0.5^2 = 0.25)
# -G 1.0: gravitational constant (Henon units)
~/bin/hermite -t 1.0 -o 2 -G 1.0 "$INPUT" > hermite.log

echo "Done. Output: $WORKDIR/hermite.log"

# Quick sanity check: verify the simulation reached the target time
echo "--- Final time check ---"
grep "Step hist: time" hermite.log | tail -1
