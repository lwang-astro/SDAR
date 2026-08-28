#!/usr/bin/env python
"""Unstable-triple BTLogH baseline analyzer (ds hyperbolic gauge plan, Phase 0).

Reads an AR sample log (unstable triple, g_func=4 auto) written by
``ar.ttl.sd.t.mulpot.cm`` and reproduces the manual baseline analysis of
``SDAR/docs/hierarchical_blogh_ds_hyperbolic_gauge_plan.md`` section 2.1:

NOTE (2026-08-27 g-func refactor): this script targets the ARCHIVED 2026-08-26
log format (mulpot binary, g_func column 0-4). To regenerate the log with the
current code use ``ar.btlogh.ttl.sd.cm --g-func 1``; its g_func column is the
active state 0/1 — map non-zero to 1 before comparing with the archived codes.

  1. ds timeline (unique ds states and their first-appearance times)
  2. dE envelope vs escaper separation r_esc
  3. steps-per-inner-orbit estimate from per-interval Nstep increments
  4. inner-binary osculating orbit (a, e, E_rel) decomposition
  5. slowdown factor (kappa) sequence

Error-channel attribution (v2 plan section 2.2, verified against the
2026-08-26 baseline log):
  - The close encounter of exchange #2 itself IS resolved (ds -> 4.09e-5,
    dE ~ 1e-6): the P_eff_min anchor (2026-08-15c) works.
  - The error accumulates AFTER the exchange, during the hyperbolic escape
    tail: the rebuilt tree freezes ds at 0.0127135 (bitwise constant to the
    end) because the hyperbolic root (a,e) are Kepler constants, while the
    runtime g(t) outer factor G*M_L*M_R/r_esc decays as 1/r.  The gauge
    mismatch ds(constant)/g(1/r) makes dt grow without bound; the inner
    binary resolution drops to ~7.9 steps/orbit (design 32) and dE climbs
    1e-6 -> 3.8e-2.
  - In the tail the ACTIVE slowdown factor is identically 1 (kappa_org
    ~ 1e-12, clamped to 1); only kappa_max = timescale/period ratchets
    5.5 -> 53.  Dominant channel: DKD operator-splitting error from the
    lost inner-binary resolution (channel A of the plan's remark); Fix A
    (restore the 32-steps/orbit invariant via r_ref = max(q, r_inst))
    repairs both channels.

Usage:
    python analysis_ustabtri_btlogh.py [logfile] [--dt-out DT]

Defaults target the ustabtri baseline:
    logfile = /home/lwang/localdata/SDAR_BLogH/ustabtri.btlogh_auto.s256.log
    dt_out  = 6.103515625e-05 (ar.cxx -o value; needed for steps/orbit)
"""

import argparse
import numpy as np
import sdar

DEFAULT_LOG = '/home/lwang/localdata/SDAR_BLogH/ustabtri.btlogh_auto.s256.log'
G_CONST = 1.0


def two_body_orbit(m1, pos1, vel1, m2, pos2, vel2, G=G_CONST):
    """Osculating two-body orbit from instantaneous states.

    Returns (a, e, E_rel) with E_rel the pair binding energy
    (E_rel = mu v^2/2 - G m1 m2 / r, negative for bound).
    """
    dr = np.asarray(pos2) - np.asarray(pos1)
    dv = np.asarray(vel2) - np.asarray(vel1)
    r = np.linalg.norm(dr)
    v2 = float(np.dot(dv, dv))
    mu = m1 * m2 / (m1 + m2)
    E_rel = 0.5 * mu * v2 - G * m1 * m2 / r
    L = mu * np.cross(dr, dv)
    L2 = float(np.dot(L, L))
    # a from energy: E = -G m1 m2 / (2a)
    a = -G * m1 * m2 / (2.0 * E_rel) if E_rel != 0.0 else np.inf
    # e from L: L^2 = G m1 m2 a(1-e^2) * mu ... use L^2 = mu^2 G(m1+m2) a (1-e^2)
    if a > 0:
        e = np.sqrt(max(0.0, 1.0 - L2 / (mu * mu * G * (m1 + m2) * a)))
    else:
        # hyperbolic: e from L and |a|
        e = np.sqrt(1.0 + L2 / (mu * mu * G * (m1 + m2) * (-a)))
    return a, e, E_rel


def split_triple(parts):
    """Split a 3-particle snapshot into (binary pair, escaper) by binding.

    Returns ((i, j), k): the most-bound pair indices and the escaper index.
    """
    assert parts['mass'].size == 3
    best, bidx = 0.0, None
    pairs = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]
    for i, j, k in pairs:
        _, _, E = two_body_orbit(parts['mass'][i], parts['pos'][i], parts['vel'][i],
                                 parts['mass'][j], parts['pos'][j], parts['vel'][j])
        if bidx is None or E < best:
            best, bidx = E, (i, j, k)
    return bidx


def analyze(logfile, dt_out):
    d = sdar.SDARData(g_func=True, N_particle=3, slowdown=True, N_sd=2,
                      time_measure=True)
    d.loadtxt(logfile, skiprows=1)
    n = d.time.size
    time = np.atleast_1d(d.time) + np.atleast_1d(d.info.time_offset)
    ds = np.atleast_1d(d.info.ds)
    de = np.abs(np.atleast_1d(d.de))
    gfunc = np.atleast_1d(d.g_func)
    nstep = np.atleast_1d(d.profile.n_step).astype(float)

    # --- particle snapshots ---
    n = d.time.size
    mass = np.array([np.atleast_1d(d.particles['p%d' % i].mass)[0:n] for i in range(3)]).T
    pos = np.array([np.atleast_1d(d.particles['p%d' % i].pos)[0:n, :] for i in range(3)])
    pos = np.transpose(pos, (1, 0, 2))     # (n, 3 particle, 3 xyz)
    vel = np.array([np.atleast_1d(d.particles['p%d' % i].vel)[0:n, :] for i in range(3)])
    vel = np.transpose(vel, (1, 0, 2))
    ids = np.array([np.atleast_1d(d.particles['p%d' % i].id)[0:n] for i in range(3)]).T

    # per-row: most-bound pair, escaper separation, inner orbit
    r_esc = np.zeros(n)
    a_in = np.zeros(n)
    e_in = np.zeros(n)
    E_in = np.zeros(n)
    esc_id = np.zeros(n, dtype=int)
    for it in range(n):
        (i, j, k) = split_triple({'mass': mass[it], 'pos': pos[it], 'vel': vel[it]})
        esc_id[it] = ids[it][k]
        r_esc[it] = np.linalg.norm(pos[it][k] -
                                   (mass[it][i] * pos[it][i] + mass[it][j] * pos[it][j]) /
                                   (mass[it][i] + mass[it][j]))
        a_in[it], e_in[it], E_in[it] = two_body_orbit(
            mass[it][i], pos[it][i], vel[it][i],
            mass[it][j], pos[it][j], vel[it][j])

    # kappa: ACTIVE slowdown factor (SD_factor column, kappa_ in slow_down.h).
    # kappa_max (SD_factor_max = timescale/period) is reported as well: in the
    # baseline escape tail the active kappa is identically 1 (kappa_org ~ 1e-12,
    # clamped), while kappa_max ratchets 5.5 -> 53 -- a symptom of the same
    # gauge mismatch, not the active error channel.
    sd0 = np.atleast_1d(d.sd.sd0.sd)
    sd1 = np.atleast_1d(d.sd.sd1.sd)
    kappa = np.maximum(sd0, sd1)
    kappa_max = np.maximum(np.atleast_1d(d.sd.sd0.sd_max),
                           np.atleast_1d(d.sd.sd1.sd_max))

    # dt per step within each output interval; steps per inner orbit
    dt_step = dt_out / np.maximum(nstep, 1.0)

    def _P_in(row):
        # total mass of the identified pair
        (i, j, k) = split_triple({'mass': mass[row], 'pos': pos[row], 'vel': vel[row]})
        mtot = mass[row][i] + mass[row][j]
        return 2.0 * np.pi * np.sqrt(abs(a_in[row])**3 / (G_CONST * mtot))

    P_in = np.array([_P_in(r_) for r_ in range(n)])
    steps_per_orbit = P_in / dt_step

    # ---------------- report ----------------
    print('=' * 100)
    print('file        :', logfile)
    print('rows        :', n, '  time span: [%.6g, %.6g]' % (time[0], time[-1]))
    print('g_func set  :', sorted(set(gfunc.tolist())))
    print('final |dE|  : %.6g    max |dE|: %.6g' % (de[-1], de.max()))
    print('final ds    : %.10g' % ds[-1])
    print()

    print('--- 1. ds timeline (unique states, first appearance) ---')
    u, idx = np.unique(np.round(ds, 12), return_index=True)
    order = np.argsort(idx)
    for oi in order:
        i = idx[oi]
        # last row with this value
        m = np.round(ds, 12) == u[oi]
        print('  t=[%.6g .. %.6g]  ds=%.10g  (Nstep_sum=%d)'
              % (time[i], time[np.where(m)[0][-1]], ds[i],
                 np.atleast_1d(d.profile.n_step_sum)[i]))
    print()

    print('--- 2. escape-tail sampling (dE envelope vs r_esc, kappa, resolution) ---')
    print('     t          ds          |dE|        r_esc    a_in        e_in     kappa  kap_max  st/orb')
    marks = [0, 1]
    for tq in (0.0022, 0.00232, 0.0365, 0.0598, 0.0600, 0.0627, 0.065, 0.0706, 0.0888):
        i = int(np.argmin(np.abs(time - tq)))
        if i not in marks: marks.append(i)
    marks.append(n - 1)
    for i in sorted(marks):
        print('  %-10.6g  %-11.5g  %-10.4g  %-7.4g  %-10.5g  %-6.4g  %-6.3g  %-7.3g  %-8.2f'
              % (time[i], ds[i], de[i], r_esc[i], a_in[i], e_in[i], kappa[i],
                 kappa_max[i], steps_per_orbit[i]))
    print()

    print('--- 3. tail diagnostics (t >= t_rebuild2 = first row with ds > 1e-2) ---')
    itail = np.where(ds > 1e-2)[0]
    if itail.size:
        i0 = itail[0]
        env = np.maximum.accumulate(de[i0:])
        print('  rebuild row %d  t=%.6g  ds=%.10g' % (i0, time[i0], ds[i0]))
        print('  ds ratio  last/first : %.6g   r_esc ratio last/first: %.6g'
              % (ds[-1] / ds[i0], r_esc[-1] / r_esc[i0]))
        print('  ds*r_esc first/last  : %.6g .. %.6g  (constant => ds ~ 1/r_esc)'
              % (ds[i0] * r_esc[i0], ds[-1] * r_esc[-1]))
        print('  min steps/orbit in tail: %.2f (design 32; baseline measured ~7.9'
              ' -- the 2.6 in the v2 plan used P_in=1.6e-4, but the rebuilt pair'
              ' a~2.24e-3, M~1.9 gives P_in~4.8e-4)'
              % steps_per_orbit[i0:].min())
        print('  active kappa in tail  : %.3g -> %.3g (identically 1 if slowdown'
              ' inactive; kappa_max %.3g -> %.3g)'
              % (kappa[i0], kappa[-1], kappa_max[i0], kappa_max[-1]))
        print('  dE envelope in tail  : %.3g -> %.3g (max %.3g)'
              % (de[i0], de[-1], de[i0:].max()))
    print()

    print('--- 4. total steps ---')
    nss = np.atleast_1d(d.profile.n_step_sum)[-1]
    nts = np.atleast_1d(d.profile.n_step_tsyn_sum)[-1]
    print('  Nstep_sum=%d  Nstep_tsyn_sum=%d' % (nss, nts))
    print('=' * 100)
    return 0


def main():
    ap = argparse.ArgumentParser(description='ustabtri BTLogH baseline analyzer')
    ap.add_argument('logfile', nargs='?', default=DEFAULT_LOG)
    ap.add_argument('--dt-out', type=float, default=6.103515625e-05,
                    help='ar.cxx -o output interval (default 6.103515625e-05)')
    args = ap.parse_args()
    analyze(args.logfile, args.dt_out)


if __name__ == '__main__':
    main()
