#!/bin/bash
# Group-transition regression bench (no stellar evolution).
#
# Purpose: self-contained replacement for the Pal5-specific replay cases when
# validating changes to the group form/break criterion or the transition layer
# (docs/transition_unification_plan.md). Each case is one deterministic
# hyperbolic encounter (the "chatter" class behind Pal5 hard_large_energy
# warnings: form -> AR episode -> break around peri of a marginally unbound
# pair, e=1.067, peri 0.02) with a controlled perturber field.
#
# Checks per case (kind=std):
#   1. runs complete
#   2. transitions exercised and mirrored: pal-run event count == 2 x ref
#      (events counted from the standing ADJUST_GROUP_DEBUG stderr lines
#      "Add new group"/"Break Group" of the plain sample build)
#   3. palindrome round-trip dvel_rt below the case gate (regression bound;
#      pre-2026-09-consolidation behavior was 100-1000x above these gates)
#   4. state-rewrite purity: max |dE_rewrite| on the "Group transition:" lines
#      <= 1e-12 (the injection-decomposition trace shares the
#      ADJUST_GROUP_DEBUG flag)
# Checks per case (kind=struct, for chaotic exchanges where round-trip error
# is intrinsically amplified): structure gates instead of 2-3 — a >=3-member
# group event must occur and the slowdown must engage (formation SD > 1);
# purity (4) still applies.
#
# Case families:
#   weak/strong hyperbolic encounter (Pal5 chatter class), late/early switch
#   bnd/chatter eccentric bound binary "breathing" through peri each orbit
#     (bnd: one form/break pair; chatter: two pairs = form-break-reform)
#   triple: compact kappa>1 binary + light intruder (3-member capture/release)
#
# --scan: r-group sweep on the weak hyperbolic IC; prints the switch-position
# table and gates on the lever direction (early switch dvel <= late dvel).
#
# Knobs: keep --slowdown-timescale-max pinned (it defaults to time-end, so
# ref/-t T and pal/-t 2T pairs would otherwise use different dynamics), and do
# not compare runs across different --r-group (r-group feeds r_neighbor into
# the block-dt controller via calcAcc0OffsetSq, so even no-group trajectories
# differ).
#
# Requires: bash, python3, make, g++ and a sample build with ADJUST_GROUP_DEBUG
# (the sample Makefile default); builds ../Hermite/build/hermite if missing.
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
HERMITE_DIR="$SCRIPT_DIR/../Hermite"
BIN="$HERMITE_DIR/build/hermite"
WORK="$(mktemp -d /tmp/sdar_gtb.XXXXXX)"
trap 'echo "artifacts: $WORK" >&2' EXIT

if [ ! -x "$BIN" ]; then
    echo "building $BIN ..." >&2
    make -C "$HERMITE_DIR" build/hermite >&2 || { echo "FAIL: build" >&2; exit 1; }
fi

# ---------------------------------------------------------------- ICs
gen_ic() {  # $1: file, $2: spectator mass, $3: sx, $4: sy, $5: svx, $6: svy
python3 - "$1" "$2" "$3" "$4" "$5" "$6" <<'PYEOF'
import sys
path, ms, sx, sy, svx, svy = sys.argv[1], *[float(v) for v in sys.argv[2:]]
# hyperbolic encounter pair: m1=m2=0.45, G=1 (mu=0.9), v_inf=3, peri q=0.02,
# e = 1+v_inf*q/mu = 1.0667, b = q*sqrt((e+1)/(e-1)); asymptotic placement at
# r0=3 so peri falls near t~0.94 (verified: single form/break pair)
import math
v_inf, q, mu, m = 3.0, 0.02, 0.9, 0.45
e = 1 + v_inf * q / mu
b = q * math.sqrt((e + 1) / (e - 1))
rows = [
    (m,  1.5,  b / 2, 0.0, -v_inf / 2, 0.0, 0.0),
    (m, -1.5, -b / 2, 0.0,  v_inf / 2, 0.0, 0.0),
    (ms, sx,   sy,    0.0,  svx,       svy,  0.0),
]
with open(path, "w") as f:
    f.write("3\n")
    for r in rows:
        f.write(" ".join(f"{v:.16e}" for v in r) + " 0.0\n")
    f.write("0\n")
PYEOF
}

gen_bound_ic() {  # $1: file — eccentric bound binary breathing through peri
python3 - "$1" <<'PYEOF'
import sys, math
# m1=m2=0.45, a=0.5, e=0.9 -> peri 0.05 < r_crit(0.055), apo 0.95 >> r_crit;
# placed at apo with perpendicular velocity; one form/break pair per orbit
# (peri passages at t~1.17, 3.51, ...; P=2.342); weak spectator
a, e, m = 0.5, 0.9, 0.45
mu = 2 * m
v_apo = math.sqrt(mu / a * (1 - e) / (1 + e))
apo = a * (1 + e)
rows = [
    (m,   +apo / 2, 0.0, 0.0, 0.0, -v_apo / 2, 0.0),
    (m,   -apo / 2, 0.0, 0.0, 0.0, +v_apo / 2, 0.0),
    (0.1, 5.0, 0.5, 0.0, 0.03, 0.02, 0.0),
]
with open(sys.argv[1], "w") as f:
    f.write("3\n")
    for r in rows:
        f.write(" ".join(f"{v:.16e}" for v in r) + " 0.0\n")
    f.write("0\n")
PYEOF
}

gen_triple_ic() {  # $1: file — compact kappa>1 binary + light intruder
python3 - "$1" <<'PYEOF'
import sys, math
# compact circular binary a=0.03, m=0.45+0.45 (apo 0.03 < r_crit 0.0553 ->
# slowdown engaged, SD ~ timescale/period); light intruder m=0.05 with
# v_inf=0.6, peri 0.05 arrives t~3: 2->3 member capture then release
m, a, m3, v_inf, q = 0.45, 0.03, 0.05, 0.6, 0.05
mu3 = 2 * m + m3
e = 1 + v_inf * q / mu3
b = q * math.sqrt((e + 1) / (e - 1))
rows = [
    (m,  +a / 2, 0.0, 0.0, 0.0, +math.sqrt(2 * m / a) / 2, 0.0),
    (m,  -a / 2, 0.0, 0.0, 0.0, -math.sqrt(2 * m / a) / 2, 0.0),
    (m3, 3.0, b, 0.0, -v_inf, 0.0, 0.0),
]
with open(sys.argv[1], "w") as f:
    f.write("3\n")
    for r in rows:
        f.write(" ".join(f"{v:.16e}" for v in r) + " 0.0\n")
    f.write("0\n")
PYEOF
}

# round-trip metrics vs the c.m.-shifted initial state from the pal(2T)
# checkpoint: velocity error, total-energy error (abs and relative to |Etot(0)|)
pal_metrics() {  # $1: ic file, $2: checkpoint; prints "dvel dE dE_rel"
python3 - "$1" "$2" <<'PYEOF'
import struct, sys, math
PREC, OFF_M, OFF_POS, OFF_VEL, OFF_ID = 184, 0, 8, 32, 64
G = 1.0
def read_ic(path):
    with open(path) as f:
        n = int(f.readline()); mass, pos, vel = [], [], []
        for _ in range(n):
            v = f.readline().split()
            mass.append(float(v[0])); pos.append([float(x) for x in v[1:4]]); vel.append([float(x) for x in v[4:7]])
    mt = sum(mass)
    for d in range(3):
        cp = sum(mass[i]*pos[i][d] for i in range(n))/mt; cv = sum(mass[i]*vel[i][d] for i in range(n))/mt
        for i in range(n): pos[i][d] -= cp; vel[i][d] -= cv
    return mass, pos, vel
def etot(mass, pos, vel):
    n = len(mass)
    ek = sum(0.5*mass[i]*sum(v*v for v in vel[i]) for i in range(n))
    ep = sum(-G*mass[i]*mass[j]/math.dist(pos[i], pos[j]) for i in range(n) for j in range(i+1, n))
    return ek + ep
with open(sys.argv[2], "rb") as f:
    data = f.read()
n, = struct.unpack_from("<i", data, 4)
by_id = {}
for i in range(n):
    base = 8 + i * PREC
    pid, = struct.unpack_from("<q", data, base + OFF_ID)
    by_id[pid] = (struct.unpack_from("<d", data, base + OFF_M)[0],
                  struct.unpack_from("<3d", data, base + OFF_POS),
                  struct.unpack_from("<3d", data, base + OFF_VEL))
mass0, pos0, vel0 = read_ic(sys.argv[1])
E0 = etot(mass0, pos0, vel0)
mass2 = [by_id[i+1][0] for i in range(n)]
pos2 = [by_id[i+1][1] for i in range(n)]
vel2 = [by_id[i+1][2] for i in range(n)]
E2 = etot(mass2, pos2, vel2)
dv = 0.0
for i in range(n):
    for d in range(3): dv = max(dv, abs(vel2[i][d] + vel0[i][d]))
dE = abs(E2 - E0)
print(f"{dv:.6e} {dE:.6e} {dE/max(abs(E0), 1e-300):.6e}")
PYEOF
}

n_events() {  # $1: stderr file; standing ADJUST_GROUP_DEBUG lines
    echo $(( $(grep -c '^Add new group' "$1") + $(grep -c '^Break Group' "$1") ))
}

# ---------------------------------------------------------------- runner
# name kind ic spectator(ms,sx,sy,svx,svy) r_group T dvel_gate dErel_gate
#   ic: hyp|bnd|tri; kind: std (mirror+dvel+dE gates) | struct (structure gates)
CASES=(
"weak_late    std  hyp 0.1 5.0 0.5 0.03 0.02 0.05 1   1e-3 5e-3"
"weak_early   std  hyp 0.1 5.0 0.5 0.03 0.02 0.4  1   1e-5 1e-4"
"strong_late  std  hyp 1.0 0.8 0.3 0.02 0.05 0.05 1   3e-3 5e-2"
"strong_early std  hyp 1.0 0.8 0.3 0.02 0.05 0.4  1   1e-3 2e-2"
"bnd_late     std  bnd  -    -   -   -    -    0.05 1.5 1e-3 1e-3"
"bnd_early    std  bnd  -    -   -   -    -    0.4  1.5 1e-3 1e-3"
"chatter      std  bnd  -    -   -   -    -    0.05 4   1e-3 1e-3"
"triple       struct tri -    -   -   -    -    0.05 3.5 -    -"
)
COMMON="-o 8 -G 1.0 --dt-max-power 8 --slowdown-timescale-max 100"

# --scan: r-group sweep on the weak hyperbolic IC (switch-position lever)
if [ "${1:-}" = "--scan" ]; then
    ic="$WORK/scan.dat"
    gen_ic "$ic" 0.1 5.0 0.5 0.03 0.02
    cat <<'LEGEND'
Columns:
  r-group     group-formation distance criterion (switch position: larger = switch earlier,
              in a smoother field)
  dV_RT       round-trip velocity error max|v(2T)+v(0)| after velocity reversal at T=1
  dE_RT_rel   round-trip total-energy error |Etot(2T)-Etot(0)|/|Etot(0)|
  result      event counts; gate = lever direction (dV_RT(0.4) <= dV_RT(0.05))
LEGEND
    printf "%-8s %-11s %-11s %s\n" r-group dV_RT dE_RT_rel result
    declare -A scan_dv
    fail=0
    for rg in 0.05 0.1 0.2 0.4; do
        tag="s$(echo "$rg" | tr -d '.')"
        ( cd "$WORK" && "$BIN" -t 1 $COMMON --r-group "$rg" -c "$tag.r.last" "$ic" >/dev/null 2>"$tag.r.err" \
          && "$BIN" -t 2 $COMMON --r-group "$rg" --reverse-at 1 -c "$tag.p.last" "$ic" >/dev/null 2>"$tag.p.err" ) || { echo "$rg: run FAIL"; fail=1; continue; }
        metrics=$(pal_metrics "$ic" "$WORK/$tag.p.last")
        dv=$(echo "$metrics" | awk '{print $1}')
        derel=$(echo "$metrics" | awk '{print $3}')
        scan_dv[$rg]="$dv"
        printf "%-8s %-11s %-11s %s\n" "$rg" "$dv" "$derel" "events ref $(n_events "$WORK/$tag.r.err") pal $(n_events "$WORK/$tag.p.err")"
    done
    # lever direction gate: switching in the smooth field (0.4) must not be
    # worse than switching near peri (0.05)
    awk -v a="${scan_dv[0.4]}" -v b="${scan_dv[0.05]}" 'BEGIN{exit !(a+0>b+0)}' \
        && { echo "switch-position lever inverted: dV_RT(0.4)=${scan_dv[0.4]} > dV_RT(0.05)=${scan_dv[0.05]}"; fail=1; }
    [ $fail -eq 0 ] && echo "scan: lever direction OK"
    exit $fail
fi

fail=0
cat <<'LEGEND'
Columns:
  CASE        test case (see header comment of this script)
  TYPE        std = mirrored-palindrome gates; struct = structure gates (chaotic exchange)
  EV_REF/PAL  group form+break event count, forward [0,T] run / velocity-reversed [0,2T] run
  MIRROR      event sequence time-mirrored (PAL == 2 x REF); struct cases are exempt
  REWRITE_dE  max |Etot before-after| across a single form/break state rewrite (gate 1e-12)
  dE_RT_rel   round-trip total-energy error |Etot(2T)-Etot(0)|/|Etot(0)|
  dV_RT       round-trip velocity error max|v(2T)+v(0)| (absolute, code units)
  RESULT      PASS/FAIL against the per-case gates (REWRITE_dE, MIRROR, dE_RT_rel, dV_RT)
LEGEND
printf "%-14s %-7s %-6s %-6s %-7s %-11s %-11s %-11s %s\n" CASE TYPE EV_REF EV_PAL MIRROR REWRITE_dE dE_RT_rel dV_RT RESULT
for case in "${CASES[@]}"; do
    read -r name kind icfam ms sx sy svx svy rg T gate degate <<<"$case"
    ic="$WORK/$name.dat"
    case "$icfam" in
        hyp) gen_ic "$ic" "$ms" "$sx" "$sy" "$svx" "$svy" ;;
        bnd) gen_bound_ic "$ic" ;;
        tri) gen_triple_ic "$ic" ;;
    esac
    T2=$(awk -v t="$T" 'BEGIN{printf "%.6g", 2*t}')
    ( cd "$WORK" && "$BIN" -t "$T" $COMMON --r-group "$rg" -c "$name.ref.last" "$name.dat" >/dev/null 2>"$name.ref.err" \
      && "$BIN" -t "$T2" $COMMON --r-group "$rg" --reverse-at "$T" -c "$name.pal.last" "$name.dat" >/dev/null 2>"$name.pal.err" ) || { echo "$name: run FAIL"; fail=1; continue; }
    nref=$(n_events "$WORK/$name.ref.err")
    npal=$(n_events "$WORK/$name.pal.err")
    de=$(awk '/^Group transition:/{d=$NF; if(d<0) d=-d; if(d>m) m=d} END{printf "%.2e", m+0}' "$WORK/$name.pal.err")
    purity="$de"
    ok=PASS
    if awk -v d="$de" 'BEGIN{exit !(d+0>1e-12)}'; then
        ok=FAIL; echo "  $name: rewrite impure |dE_rewrite|=$de"
    fi
    if [ "$kind" = std ]; then
        metrics=$(pal_metrics "$ic" "$WORK/$name.pal.last")
        dv=$(echo "$metrics" | awk '{print $1}')
        derel=$(echo "$metrics" | awk '{print $3}')
        [ "$nref" -ge 1 ] || { ok=FAIL; echo "  $name: transitions not exercised (ref events $nref)"; }
        [ "$npal" -eq $((2*nref)) ] || { ok=FAIL; echo "  $name: events not mirrored (ref $nref, pal $npal)"; }
        awk -v a="$dv" -v g="$gate" 'BEGIN{exit !(a+0>g+0)}' && { ok=FAIL; echo "  $name: dvel_rt $dv > gate $gate"; }
        awk -v a="$derel" -v g="$degate" 'BEGIN{exit !(a+0>g+0)}' && { ok=FAIL; echo "  $name: dE_rt/|E0| $derel > gate $degate"; }
    else
        # structure gates: >=3-member group event occurred, slowdown engaged
        n3=$(grep -cE '^Add new group.*Member_index: [0-9]+ [0-9]+ [0-9]+ ' "$WORK/$name.ref.err" "$WORK/$name.pal.err" | awk -F: '{s+=$2} END{print s+0}')
        sdmax=$(awk '/^Group formation:/{for(i=1;i<NF;i++) if($i=="SD:") v=$(i+1)+0; if(v>m) m=v} END{print m+0}' "$WORK/$name.ref.err")
        metrics=$(pal_metrics "$ic" "$WORK/$name.pal.last")
        dv="n/a"
        derel=$(echo "$metrics" | awk '{print $3}')
        [ "$n3" -ge 1 ] || { ok=FAIL; echo "  $name: no >=3-member group event"; }
        awk -v s="$sdmax" 'BEGIN{exit !(s>1)}' || { ok=FAIL; echo "  $name: slowdown not engaged (max formation SD=$sdmax)"; }
    fi
    [ "$ok" = FAIL ] && fail=1
    printf "%-14s %-7s %-6s %-6s %-7s %-11s %-11s %-11s %s\n" "$name" "$kind" "$nref" "$npal" \
        "$([ "$npal" -eq $((2*nref)) ] && echo yes || echo NO)" "$purity" "$derel" "$dv" "$ok"
done
exit $fail
