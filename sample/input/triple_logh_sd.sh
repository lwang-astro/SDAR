#!/bin/bash
set -e

# AR LogH + slowdown (tree) integration of a hierarchical triple
# Usage: ./triple_logh_sd.sh [working_dir]
#   Uses sample/input/triple.stable.lowm3 as input.
#   Runtime: < 0.1 second.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT="$SCRIPT_DIR/triple.stable.lowm3"
WORKDIR="${1:-/tmp/sdar_triple_logh_sd}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Check binary exists
if ! command -v ~/bin/ar.logh.sd.t >/dev/null 2>&1; then
    echo "Error: ~/bin/ar.logh.sd.t not found. Run 'make -C sample/AR && make -C sample/AR install' first." >&2
    exit 1
fi

echo "# $(date): ~/bin/ar.logh.sd.t -t 1.0e-3 -n 1000 -G 1.0 $INPUT" >> commands.log

# Run AR LogH with slowdown (hierarchical tree mode)
# -t 1.0e-3: integrate for 0.001 time units (short for testing)
# -n 1000:   1000 integration steps
# -G 1.0:    gravitational constant (Henon units)
~/bin/ar.logh.sd.t -t 1.0e-3 -n 1000 -G 1.0 "$INPUT" > triple.logh.sd.t.log

echo "Done. Output: $WORKDIR/triple.logh.sd.t.log"
