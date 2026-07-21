#!/bin/bash
set -e

# Compare different AR methods for a hierarchical triple system
# Usage: ./triple_compare_methods.sh [working_dir]
#   Uses sample/input/triple.stable.lowm3 as input.
#   Runtime: < 0.5 second (6 methods × ~0.03s).

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT="$SCRIPT_DIR/triple.stable.lowm3"
WORKDIR="${1:-/tmp/sdar_triple_compare}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "# $(date): Comparing AR methods for $INPUT" >> commands.log

method_list='ttl.sd.t ttl logh.sd.t logh'

for suffix in $method_list
do
    BIN=~/bin/ar.$suffix
    if ! command -v "$BIN" >/dev/null 2>&1; then
        echo "Warning: $BIN not found, skipping" >&2
        continue
    fi
    echo "# $(date): $BIN -t 1.0e-3 -n 1000 $INPUT" >> commands.log
    echo "Running ar.$suffix ..."
    "$BIN" -t 1.0e-3 -n 1000 "$INPUT" > "triple.${suffix}.log"
done

echo "Done. Results in $WORKDIR/"
echo ""
echo "--- Quick comparison (last line of each log) ---"
for suffix in $method_list; do
    f="triple.${suffix}.log"
    if [ -f "$f" ]; then
        # Extract final time and energy error from the last data row
        echo -n "ar.$suffix: "
        tail -1 "$f" | awk '{printf "time=%.6g  dE=%.2e\n", $1, $2}'
    fi
done
