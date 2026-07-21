#!/bin/bash
set -e

# Build a Kepler binary tree and convert between Kepler ↔ Cartesian coordinates
# Usage: ./build_kepler_tree.sh [working_dir]
#   Demonstrates both keplertree (tree params → particles) and
#   keplerorbit (particle pairs → Kepler params).
#   Runtime: < 0.1 second.

WORKDIR="${1:-/tmp/sdar_kepler_tree}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

if ! command -v ~/bin/keplertree >/dev/null 2>&1; then
    echo "Error: ~/bin/keplertree not found. Run 'make -C sample/Kepler && make -C sample/Kepler install' first." >&2
    exit 1
fi

if ! command -v ~/bin/keplerorbit >/dev/null 2>&1; then
    echo "Error: ~/bin/keplerorbit not found. Run 'make -C sample/Kepler && make -C sample/Kepler install' first." >&2
    exit 1
fi

# --- Step 1: keplertree (binary tree parameters → Cartesian particles) ---
# Format: level branch m1 m2 semi ecc incline rot_horiz rot_self ecca
cat > binary_tree.dat << 'EOF'
0 0 0.6 0.4 0.5 0.0 0.0 0.0 0.0 0.0
EOF

echo "# $(date): ~/bin/keplertree binary_tree.dat" >> commands.log
~/bin/keplertree binary_tree.dat > keplertree_out.log
echo "--- keplertree output (tree → Cartesian) ---"
cat keplertree_out.log

# --- Step 2: keplerorbit (Cartesian particles → Kepler parameters) ---
# Note: keplerorbit expects particle lines WITHOUT an N header.
# Provide extra data lines to avoid a known EOF detection issue.
cat > kepler_particles.dat << 'EOF'
 0.6  0.0  0.0  0.0   0.0  2.0  0.0
 0.4  0.02 0.0  0.0   0.0 -3.0  0.0
 0.5  0.0  0.0  0.0   0.0  1.0  0.0
 0.5  0.03 0.0  0.0   0.0 -1.0  0.0
EOF

echo ""
echo "# $(date): ~/bin/keplerorbit -n 1 kepler_particles.dat" >> commands.log
echo "--- keplerorbit output (Cartesian → Kepler) ---"
~/bin/keplerorbit -n 1 kepler_particles.dat | head -1

echo ""
echo "Done. Results in $WORKDIR/"
