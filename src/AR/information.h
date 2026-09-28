#pragma once
#include "AR/g_func.h"
#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "AR/slow_down.h"

namespace AR {

    //! binary parameter with slowdown
    class BinarySlowDown: public COMM::Binary {
    public:
        SlowDown slowdown;
        Float stab_check_time;

        //! write class data to file with binary format
        /*! @param[in] _fp: FILE type file for output
         */
        void writeBinary(FILE *_fp) const {
            Binary::writeBinary(_fp);
            slowdown.writeBinary(_fp);
        }

        void printColumnBinary(std::ostream& _fout) const {
            Binary::printColumnBinary(_fout);
            slowdown.printColumnBinary(_fout);
        }

        //! read class data to file with binary format
        /*! @param[in] _fin FILE type file for reading
         */
        void readBinary(FILE *_fin) {
            Binary::readBinary(_fin);
            slowdown.readBinary(_fin);
        }

        void readBinary(std::istream& _fin) {
            Binary::readBinary(_fin);
            slowdown.readBinary(_fin);
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        static void printColumnTitleAscii(std::ostream & _fout, const int _width=20) {
            Binary::printColumnTitleAscii(_fout, _width);
            SlowDown::printColumnTitleAscii(_fout,_width);
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnAscii(std::ostream & _fout, const int _width=20){
            Binary::printColumnAscii(_fout, _width);
            slowdown.printColumnAscii(_fout,_width);
        }

    };
    
    //! define ar binary tree
    template <class Tparticle>
    using BinaryTree=COMM::BinaryTree<Tparticle,BinarySlowDown>;

    //! Fix step options for integration with adjusted step (not for time sychronizatio phase)
    /*! always: use the given step without change \n
        later: fix step after a few adjustment of initial steps due to energy error
        none: don't fix step
     */
    enum class FixStepOption {always, later, none};

    //! A class contains information (e.g. parameters, binary tree, indices) about the particle group
    /*! The member of this class should not be the data that must be recored and should be possible calculated any time based on the main class (TimeTransformedSymplecticIntegrator) data. 
      This class must be inherited when a different Information class is applied in the main class, since the binary tree is used in the integration for slowdown factor
      The basic members are used in the integration are \\
      ds: integration step size \\
      binarytree: the Kepler orbital parameters of the hierarchical systems and slowdown factors \\
      fix_step_option: option to control whether the adjustment of step sizes are used \\
     */
    template <class Tparticle, class Tpcm>
    class Information{
#ifdef AR_G_FUNC
    private:
        //! effective period (timescale) of one hierarchy level — shared by leaves and internal nodes
        /*! elliptic (semi>0): slowdown-effective period P*kappa;
            hyperbolic/degenerate (semi<=0): encounter timescale 2π |a|^{3/2} / √(G(m1+m2));
            invalid level (zero mass): returns 0.
         */
        Float calcEffectivePeriod(const BinaryTree<Tparticle>& _bin, const Float& _G) const {
            if (_bin.m1 <= 0 || _bin.m2 <= 0) return 0.0;
            if (_bin.semi > 0) {
                return _bin.slowdown.getEffectivePeriod();
            }
            // hyperbolic / degenerate / semi==0: encounter timescale
            Float abs_semi = -_bin.semi;
            return 2.0 * COMM::PI
                 * sqrt(pow(abs_semi, Float(3)) / (_G * (_bin.m1 + _bin.m2)));
        }

        //! Iteration to accumulate product of ds_i and periods for BLogH ds formula
        /*! For each innermost binary, computes per-orbit ds_i (no substep coeff)
            and accumulates ds_prod and period_prod. P_eff_min is the effective
            period of the FASTEST level — leaves (elliptic: P*kappa, hyperbolic:
            encounter timescale) AND internal nodes (their own orbit) — i.e. the
            true resolution driver of the whole hierarchy. Both use the shared
            calcEffectivePeriod helper.
            BLogH-family ds = prod(ds_i) * P_eff_min / prod(P_eff)
        */
        void calcBLogHDsIter(Float& _ds_prod, Float& _period_prod, int& _nbin,
                             Float& _P_eff_min,
                             BinaryTree<Tparticle>& _bin,
                             const int _int_order, const Float& _G) {
            if (_bin.getMemberN() > 2) {
                // internal node: its own orbit is a candidate for the fastest
                // level. During a close encounter of an outer member the node
                // period (NOT the leaf virtual period) is the true resolution
                // driver — without it the ds estimate ignores the fast plunging
                // orbit and stays far too large. Same effective-period convention
                // as the leaves, so healthy hierarchies are unaffected (the inner
                // leaf P_eff remains the minimum).
                Float P_node = calcEffectivePeriod(_bin, _G);
                if (P_node > 0 && P_node < _P_eff_min) _P_eff_min = P_node;
                for (int k=0; k<2; k++) {
                    if (_bin.isMemberTree(k)) {
                        calcBLogHDsIter(_ds_prod, _period_prod, _nbin, _P_eff_min,
                                       *_bin.getMemberAsTree(k), _int_order, _G);
                    }
                }
            } else {
                if (_bin.m1 > 0 && _bin.m2 > 0) {
                    // perturbation damping (well-defined for elliptic AND hyperbolic)
                    Float scale_factor = calcPertScale(_bin, _int_order);

                    Float ds_i;
                    // per-orbit / per-encounter ds, matching the LogH resolution
                    // convention after the global /32 step division:
                    //   elliptic  : 2π   -> 1/32 orbit per step
                    //   hyperbolic: 2π/8 -> 1/256 encounter per step
                    // (a hyperbolic peri-center crossing needs 8x finer resolution
                    // than orbit sampling; lost before when 2π was used for both)
                    if (_bin.semi > 0) {
                        ds_i = calcDsElliptic(_bin, _G, 2.0 * COMM::PI);
                    } else {
                        ds_i = calcDsHyperbolic(_bin, _G, 2.0 * COMM::PI / 8.0);
                    }
                    ds_i *= scale_factor;
                    _ds_prod *= ds_i;
                    // equivalent timescale / effective period (shared helper)
                    Float P_equiv = calcEffectivePeriod(_bin, _G);
                    _period_prod *= P_equiv;
                    if (P_equiv > 0 && P_equiv < _P_eff_min) _P_eff_min = P_equiv;
                    _nbin++;
                }
            }
        }

        //! Multiply non-leaf tree node potentials into ds for hierarchical BLogH
        /*! Walks binary tree recursively. For each node with >2 members (non-leaf),
            multiplies U_node = G * m1 * m2 / a (semi-major axis based, orbit-averaged
            potential) into ds, scaled by the node-level perturbation ratio.
            The instantaneous r_sep overestimates U at outer peri-center by
            1/(1-e_out) (10x at e=0.9), and when the hierarchy is transiently
            restructured r_sep is meaningless — both inflate ds and degrade accuracy.
            ds is a per-orbit resolution quantity, so the orbital average is the
            correct source.
        */
        void multiplyDsByNodePotentials(BinaryTree<Tparticle>& _bin, const Float& _G,
                                        const int _int_order) {
            if (_bin.getMemberN() > 2 && _bin.m1 > 0 && _bin.m2 > 0) {
                // orbit-averaged potential (semi-based); conservative fallbacks below
                Float U_node;
                if (_bin.semi > 0) {
                    // orbit-averaged (semi-based) gauge (restored 2026-09-01; the
                    // apoapsis gauge of 2026-08-29 reverted): on a Kepler orbit
                    // <1/r>_t = 1/a exactly, so U(a) = <U(r)>_t and the ds estimate
                    // matches the orbit-averaged step count N ~ 32/ds_scale per
                    // node orbit. With the LogH-family time transformation the
                    // steps are uniform in eccentric anomaly (U dt ~ dE), so there
                    // is no under-resolved apoapsis phase to compensate; the
                    // apoapsis gauge merely inflated the count to ~32*(1+e) per
                    // eccentric node (triple/quadruple tests: only the semi gauge
                    // reproduces the expected 32). The g-side floor (processOuterNode
                    // caps r at the apoapsis for elliptic nodes) is unchanged and
                    // remains the safety bound on dt.
                    U_node = _G * _bin.m1 * _bin.m2 / _bin.semi;
                } 
                else if (_bin.semi < 0) {
                    // hyperbolic outer orbit: no orbital average, use the
                    // peri-center (maximum) potential as a conservative estimate.
                    // Keplerian consistency: semi<0 implies ecc>1 (asserted) so a
                    // non-positive r_ref cannot corrupt ds.
                    // NOTE (2026-08-26): ds keeps this frozen peri-center gauge by
                    // design — ds must not change while the tree is unchanged
                    // (preserves the extended-phase-space structure and time
                    // symmetry). The hyperbolic escape-tail gauge mismatch
                    // (ds frozen vs runtime g ~ 1/r decaying) is handled on the
                    // g-function side (r cap in processOuterNode), not here.
                    ASSERT(_bin.ecc > 1.0);
                    Float r_ref = (-_bin.semi) * (_bin.ecc - 1.0);
                    ASSERT(r_ref > 0);
                    U_node = _G * _bin.m1 * _bin.m2 / r_ref;
                }
                else {
                    // semi==0 (parabolic/degenerate): keep the instantaneous-
                    // separation gauge. NOTE (2026-08-29 experiment): using the
                    // angular-momentum pericenter r_p = h^2/(2*G*(m1+m2)) here
                    // blows up ds during transient (three-body dance) fits where
                    // semi==0 occurs with r_p << r (U_node ~ G*m1*m2/r_p inflates
                    // the ds product by orders of magnitude, e.g. 150x spikes in
                    // the quintuple test). The g-side cap in processOuterNode
                    // still uses r_p (floor >= U(r_sep) is satisfied), which is
                    // the safety-relevant side.
                    auto* m0 = _bin.getMember(0);
                    auto* m1 = _bin.getMember(1);
                    Float dx = m0->pos[0] - m1->pos[0];
                    Float dy = m0->pos[1] - m1->pos[1];
                    Float dz = m0->pos[2] - m1->pos[2];
                    Float r_sep = sqrt(dx*dx + dy*dy + dz*dz);
                    ASSERT(r_sep > 0);
                    U_node = _G * _bin.m1 * _bin.m2 / r_sep;
                }

                // node-level perturbation scaling: when this hierarchy level is
                // strongly perturbed (pert_out >= pert_in) it cannot be resolved by
                // the multi-level structure, so its potential contribution must be
                // suppressed (same formula as the LogH leaf scaling). This makes ds
                // degrade to the conservative inner-binary dominated value instead of
                // jumping up during transient tree restructures.
                Float node_scale = calcPertScale(_bin, _int_order);
                // defensive checks (Fix B): r_ref > 0, node_scale in (0,1],
                // U_node finite and positive
                ASSERT(node_scale > 0.0 && node_scale <= 1.0);
                ASSERT(U_node > 0.0 && U_node < NUMERIC_FLOAT_MAX);
                ds *= U_node * node_scale;

                for (int k=0; k<2; k++) {
                    if (_bin.isMemberTree(k)) {
                        multiplyDsByNodePotentials(*_bin.getMemberAsTree(k), _G, _int_order);
                    }
                }
            }
        }
#endif

    public:
        Float ds;  ///> initial step size for integration
        Float peff_min;  ///> effective period of the fastest level at the last ds calculation (used by the quiescence-gate safety valve)
        Float ds_est_prev; ///> previous ds re-estimate candidate (quiescence gate v2: an epoch invariant must be reproducible across consecutive regens)
        bool ds_est_prev_valid; ///> whether ds_est_prev holds a usable candidate
        Float ds_pert_ratio_coff; ///> perturbation ds damping coefficient: scale = min(1,(coff*pert_in/pert_out)^(1/order)); never applied when perturbation data is missing (pert_out<=0 or pert_in<=0)
        Float time_offset; ///> offset of time to obtain real physical time (real time = TimeTransformedSymplecticIntegrator:time_ + info.time_offset)
        Float r_break_crit;    // group break radius criterion
        FixStepOption fix_step_option; ///> fix step option for integration
        COMM::List<BinaryTree<Tparticle>> binarytree; ///> a list of binary tree that contain the hierarchical orbital parameters of the particle group.
#ifdef AR_DEBUG_DUMP
        bool dump_flag; ///> for debuging dump
#endif

        //! initializer, set ds to zero, fix_step_option to none
        Information(): ds(0.0), peff_min(0.0), ds_est_prev(0.0), ds_est_prev_valid(false), ds_pert_ratio_coff(0.1), time_offset(0.0), r_break_crit(-1.0), fix_step_option(AR::FixStepOption::none), binarytree() {
#ifdef AR_DEBUG_DUMP
            dump_flag = false;
#endif
        }

        //! check whether parameters values are correct initialized
        /*! \return true: all correct
         */
        bool checkParams() {
            // must be strictly positive: coff=0 gives scale=0, i.e. ds=0
            ASSERT(ds_pert_ratio_coff>0.0);
            ASSERT(r_break_crit>=0.0);
            ASSERT(binarytree.getSize()>0);
            return true;
        }

        //! reserve memory of binarytree list
        void reserveMem(const int _nmax) {
            binarytree.setMode(COMM::ListMode::local);
            binarytree.reserveMem(_nmax);
        }

        //! get the root of binary tree
        BinaryTree<Tparticle>& getBinaryTreeRoot() const {
            int n = binarytree.getSize();
            ASSERT(n>0);
            return binarytree[n-1];
        }

        //! calculate ds for an elliptic orbit
        /*!
          For an ellpitic orbit, step ds=dt*G*m1*m2/r, the estimation of ds for one orbit is: 2*pi*sqrt(G*semi/(m1+m2))*m1*m2 
          @param[in] _bin: binary tree to check
          @param[in] _G: gravitational constant
          @param[in] _coff: coefficient for ds, default is 2*pi/32 (1/32 orbit)
         */
        inline Float calcDsElliptic(BinaryTree<Tparticle>& _bin, const Float& _G, const Float _coff = 0.19634954084) {
            return _coff*sqrt(_G*_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
        }

        //! calculate ds for a hyperbolic orbit
        /*!
          For an hyperbolic orbit, step ds=dt*G*m1*m2/r, the estimation of ds for one orbit is: 2*pi*sqrt(-G*semi/(m1+m2))*m1*m2 
          @param[in] _bin: binary tree to check
          @param[in] _G: gravitational constant
          @param[in] _coff: coefficient for ds, default is 2*pi/256 (1/256 orbit)
         */
        inline Float calcDsHyperbolic(BinaryTree<Tparticle>& _bin, const Float& _G, const Float _coff = 0.0245436926) {
            return _coff*sqrt(-_G*_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
        }

        //! compute the per-level ds damping scale from the perturbation ratio
        /*! scale = min(1, (ds_pert_ratio_coff * pert_in/pert_out)^(1/_int_order)):
            a level is damped only when ds_pert_ratio_coff*pert_ratio < 1, i.e. the
            outer tidal perturbation approaches the level binding (hierarchy
            breaking down).
            Data-presence guard (2026-09-01): when the perturbation measures are
            missing (pert_out<=0, e.g. the root of an isolated group with no
            external perturbers, or pert_in<=0) the scale is EXACTLY 1.0,
            independent of ds_pert_ratio_coff. A plain sentinel ratio of 1.0 fed
            through the coff multiplier would wrongly damp unperturbed levels by
            coff^(1/order) (0.56 at coff=0.1, order 4).
            The apo-based pert_in from the slowdown data is used for elliptic orbits
            (smooth, orbit-averaged measure). For hyperbolic orbits (semi<=0) it is
            negative/undefined, so the instantaneous tidal metric
            (COMM::Binary::calcPertFromMR) on the live member separation is used
            instead, giving a well-defined perturbation ratio everywhere.
         */
        Float calcPertScale(BinaryTree<Tparticle>& _bin, const int _int_order) {
            Float pert_out = _bin.slowdown.pert_out;
            if (pert_out <= 0) return 1.0; // no perturbation measure: no damping
            Float pert_in;
            if (_bin.semi > 0) {
                pert_in = _bin.slowdown.pert_in;
            } else {
                auto* pm0 = _bin.getMember(0);
                auto* pm1 = _bin.getMember(1);
                Float dr[3] = {pm1->pos[0]-pm0->pos[0],
                               pm1->pos[1]-pm0->pos[1],
                               pm1->pos[2]-pm0->pos[2]};
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                Float r = (r2 > 0) ? sqrt(r2) : 0.0;
                pert_in = (r > 0) ? COMM::Binary::calcPertFromMR(r, _bin.m1, _bin.m2) : 0.0;
            }
            if (pert_in <= 0) return 1.0; // no binding measure: no damping
            return std::min(Float(1.0),
                            pow(ds_pert_ratio_coff * pert_in / pert_out, 1.0 / Float(_int_order)));
        }

        //! iteration for the LogH sum-gauge ds: accumulate the perturbation-damped
        //! orbit-averaged potentials of ALL tree levels and track the smallest
        //! effective period
        /*!
          @param[out] _u_sum: sum over levels of chi_L * <U_L> (orbit average)
          @param[out] _P_eff_min: smallest effective period among all levels
          @param[in] _bin: current tree node
          @param[in] _int_order: symplectic integrator accurate order
          @param[in] _G: gravitational constant
         */
        void calcLogHSumGaugeIter(Float& _u_sum, Float& _P_eff_min, BinaryTree<Tparticle>& _bin, const int _int_order, const Float& _G) {
            if (_bin.m1 > 0 && _bin.m2 > 0) {
                // per-level potential gauge (same convention as the BLogH family):
                // elliptic: semi-major axis (orbit-averaged, <1/r>_t = 1/a);
                // hyperbolic: peri-center; degenerate: instantaneous
                Float u_min;
                if (_bin.semi > 0) {
                    u_min = _G * _bin.m1 * _bin.m2 / _bin.semi;
                }
                else if (_bin.semi < 0 && _bin.ecc > 1.0) {
                    u_min = _G * _bin.m1 * _bin.m2 / ((-_bin.semi) * (_bin.ecc - 1.0));
                }
                else {
                    auto* m0 = _bin.getMember(0);
                    auto* m1 = _bin.getMember(1);
                    Float dx = m0->pos[0] - m1->pos[0];
                    Float dy = m0->pos[1] - m1->pos[1];
                    Float dz = m0->pos[2] - m1->pos[2];
                    u_min = _G * _bin.m1 * _bin.m2 / sqrt(dx*dx + dy*dy + dz*dz);
                }
                Float scale = calcPertScale(_bin, _int_order);
                _u_sum += u_min * scale;

                // effective period: hyperbolic encounter timescale; elliptic uses
                // the UN-SLOWED Kepler period (2026-09-14 fix): multiplying by the
                // slowdown factor kappa (2026-08-29) inflated ds by kappa (10-500x
                // in the quad_sd2 slowdown B--B Hermite test), pushing the
                // per-interval energy error onto the -e check ceiling and producing
                // a linear semi-major-axis drift (1e-3 by t=40) plus giant initial
                // ds (60-700x tolerance trips at group creation). The pre-08-29
                // LogH estimator (calcDsKeplerBinaryTree) never scaled by kappa;
                // this restores that resolution. kappa==1 (no slowdown) is
                // bit-identical. The BLogH/BTLogH product path keeps its own
                // P*kappa gauge (getEffectivePeriod), unchanged.
                Float p_eff;
                if (_bin.semi > 0) {
                    p_eff = 2.0 * COMM::PI
                          * sqrt(pow(_bin.semi, Float(3)) / (_G * (_bin.m1 + _bin.m2)));
                }
                else if (_bin.semi < 0) {
                    p_eff = 2.0 * COMM::PI
                          * sqrt(pow(-_bin.semi, Float(3)) / (_G * (_bin.m1 + _bin.m2)));
                    // the energy-scale "period" of a hyperbolic level is the
                    // pericenter-region timescale; for a high-energy flyby
                    // seen at r >> |semi| it underestimates the encounter
                    // duration by ~(r/|semi|)^{3/2} (measured 2.4e5x after an
                    // SN-kick plunge), which pins ds at pericenter resolution
                    // over the whole in/out legs. Floor by the live crossing
                    // timescale of the current member separation.
                    const auto& pm0 = (_bin.isMemberTree(0) ? *_bin.getMemberAsTree(0) : *_bin.getMember(0));
                    const auto& pm1 = (_bin.isMemberTree(1) ? *_bin.getMemberAsTree(1) : *_bin.getMember(1));
                    Float dr[3] = {pm1.pos[0]-pm0.pos[0], pm1.pos[1]-pm0.pos[1], pm1.pos[2]-pm0.pos[2]};
                    Float dv[3] = {pm1.vel[0]-pm0.vel[0], pm1.vel[1]-pm0.vel[1], pm1.vel[2]-pm0.vel[2]};
                    Float r_now = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
                    Float v_now = sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
                    if (r_now > 0.0 && v_now > 0.0) {
                        Float t_cross = 2.0*r_now/v_now;
                        if (p_eff < t_cross) p_eff = t_cross;
                    }
                }
                else {
                    p_eff = 0.0; // degenerate level: no period contribution
                }
                if (p_eff > 0 && p_eff < _P_eff_min) _P_eff_min = p_eff;
            }
            for (int k = 0; k < 2; k++) {
                if (_bin.isMemberTree(k)) calcLogHSumGaugeIter(_u_sum, _P_eff_min, *_bin.getMemberAsTree(k), _int_order, _G);
            }
        }

        //! LogH ds from the total-potential orbit-averaged gauge (2026-08-29; semi gauge restored 2026-09-01)
        /*! ds = _ds_scale/32 * P_eff,min * sum_L chi_L * <U_L>, with <U_L> the
          orbit-averaged level potential (elliptic: G*mi*mj/a, since <1/r>_t = 1/a).
          The runtime LogH time transformation is g = sum over ALL pairs of
          G*mi*mj/rij (with slowdown); averaged over a Kepler orbit <g> ~ sum_L <U_L>,
          so each level receives ~32/_ds_scale substeps per effective period,
          distributed uniformly in eccentric anomaly (U dt ~ dE). The apoapsis
          lower-bound variant tried on 2026-08-29 was reverted: it over-resolves
          eccentric levels by ~(1+e) each (triple/quadruple tests). This replaces
          the old min-over-innermost-pairs form, which followed the tightest pair's
          OWN potential while g sums all pairs, over-resolving hierarchical systems
          by (sum_L <U_L>)/<U_tightest> (3-5x in typical hierarchies). For an
          isolated binary the two forms coincide.
         */
        Float calcDsLogHSumGauge(BinaryTree<Tparticle>& _bin, const int _int_order, const Float& _G, const Float& _ds_scale) {
            Float u_sum = 0.0;
            Float P_eff_min = NUMERIC_FLOAT_MAX;
            calcLogHSumGaugeIter(u_sum, P_eff_min, _bin, _int_order, _G);
            // degenerate group (post-merger interrupt: one member has zero
            // mass, no level contributes): no Kepler gauge exists; keep the
            // current ds - the caller drifts the single massive remnant out
            // of the interval, ds is not used for it
            if (!(u_sum > 0)) return ds;
            ASSERT(P_eff_min < NUMERIC_FLOAT_MAX);
            return (_ds_scale / 32.0) * P_eff_min * u_sum;
        }

        //! calculate ds from the inner most binary with minimum period, determine the fix step option
        /*! Estimate ds first from the inner most binary orbit (eccentric anomaly), set fix_step_option to later
          @param[in] _int_order: accuracy order of the symplectic integrator.
          @param[in] _G: gravitational constant
          @param[in] _ds_scale: scaling factor to determine ds
          @param[in] _g_func_on: true = the g-function method of this build is active (BLogH-family ds formula), false = standard LogH (default)
         */
        void calcDsAndStepOption(const int _int_order, const Float& _G, const Float& _ds_scale, const bool _g_func_on = false) {
            auto& bin_root = getBinaryTreeRoot();

#ifdef AR_G_FUNC
            if (_g_func_on) {
                // BLogH family: accumulate product of per-orbit ds_i and periods
                Float ds_prod = 1.0;
                Float period_prod = 1.0;
                Float P_eff_min = NUMERIC_FLOAT_MAX;
                int nbin = 0;
                calcBLogHDsIter(ds_prod, period_prod, nbin, P_eff_min, bin_root, _int_order, _G);

                // degenerate group (post-merger): no level contributes, keep
                // the current ds (see calcDsLogHSumGauge)
                if (nbin > 0) {
                    // plain product formula, ds ~ [energy^nbin·time]
                    // ds = Π(ds_i) * P_eff_min / Π(P_eff)
                    ds = ds_prod * P_eff_min / period_prod;
                    peff_min = P_eff_min;  // stored for the Fix-2 step-count ceiling (symplectic_integrator.h)
#ifdef AR_G_FUNC_BTLOGH
                    // with outer potential, eccentricity may affect ds determination that ds is not exact reach P_eff_min.
                    // node potentials are orbit-averaged (semi-based) and
                    // pert-ratio scaled, see multiplyDsByNodePotentials
                    multiplyDsByNodePotentials(bin_root, _G, _int_order);
#endif
                    // DKD integrator divides each orbit into n_sub substeps, default is 32 substeps, use _ds_scale to change it.
                    ds *= _ds_scale / 32.0;
                }
            } else {
                ds = calcDsLogHSumGauge(bin_root, _int_order, _G, _ds_scale);
            }
#else
            (void)_g_func_on; // no g-func method in this build; use the LogH sum-gauge form
            ds = calcDsLogHSumGauge(bin_root, _int_order, _G, _ds_scale);
#endif
            ASSERT(ds>0);

            const int n_particle = bin_root.getMemberN();

            // determine the fix step option
            fix_step_option = FixStepOption::none;
            // for two-body case, determine the step at begining then fix;
            // hyperbolic two-body (one-shot encounter) stays adaptive: the
            // pericenter needs a much smaller ds than the in/out legs, and
            // 'later' freezes the pericenter ds for the whole interval
            // (measured 2.3M/11.3M steps per hard block after an SN-kick
            // near-radial plunge)
            if (n_particle==2) {
                if (bin_root.semi > 0.0) fix_step_option = AR::FixStepOption::later;
            }
            else if (bin_root.stab<1.0) fix_step_option = AR::FixStepOption::later;
        }

        //! generate binary tree for the particle group
        /*! 
          - Construct the binary tree based on particle positions and velocities\n
          - If particle mass is 0.0, set to unused particles and put at the outermost orbits\n
          - If USE_CM_FRAME is used, switch each coordinate frame of the binary members to their center-of-the-mass frame\n

          @param[in] _particles: particle group 
          @param[in] _G: gravitational constant
         */
        void generateBinaryTree(COMM::ParticleGroup<Tparticle, Tpcm>& _particles, const Float _G) {
            // If the existing binarytree is not in original frame, first shift to original frame
#ifdef USE_CM_FRAME
            if (binarytree.getSize()>0) {
                auto& bin_root = getBinaryTreeRoot();
                if (!bin_root.isOriginFrame()) bin_root.shiftToOriginFrame();
            }
#endif

            const int n_particle = _particles.getSize();
            const int binary_tree_index = static_cast<int>(COMM::BinaryTreeMemberIndexTag::binarytree);
            ASSERT(n_particle>1);
            binarytree.resizeNoInitialize(n_particle-1);
            int particle_index_local[n_particle];
            int particle_index_unused[n_particle];
            int n_particle_real = 0;
            int n_particle_unused=0;
            for (int i=0; i<n_particle; i++) {
                if (_particles[i].mass>0.0) particle_index_local[n_particle_real++] = i;
                else particle_index_unused[n_particle_unused++] = i;
            }
            if (n_particle_real>1) 
                BinaryTree<Tparticle>::generateBinaryTree(binarytree.getDataAddress(), particle_index_local, n_particle_real, _particles.getDataAddress(), _G);

            // Add unused particles to the outmost orbit
            if (n_particle_real==0) {
                binarytree[0].setMembers(&(_particles[0]), &(_particles[1]), 0, 1);
                binarytree[0].mass = 0.0;
                binarytree[0].m1 = 0.0;
                binarytree[0].m2 = 0.0;  
                for (int i=2; i<n_particle; i++) {
                    binarytree[i-1].setMembers((Tparticle*)&(binarytree[i-2]), &( _particles[i]), binary_tree_index, i);
                    binarytree[i-1].mass = 0.0;
                    binarytree[i-1].m1 = 0.0;
                    binarytree[i-1].m2 = 0.0;
                }
            }
            else if (n_particle_real==1) {
                int i1 = particle_index_local[0];
                int i2 = particle_index_unused[0];
                binarytree[0].setMembers(&(_particles[i1]), &(_particles[i2]), i1 ,i2);
                binarytree[0].m1 = _particles[i1].mass;
                binarytree[0].m2 = 0.0;
                binarytree[0].mass = binarytree[0].m1;
                binarytree[0].pos[0] = _particles[i1].pos[0];
                binarytree[0].pos[1] = _particles[i1].pos[1];
                binarytree[0].pos[2] = _particles[i1].pos[2];
                binarytree[0].vel[0] = _particles[i1].vel[0];
                binarytree[0].vel[1] = _particles[i1].vel[1];
                binarytree[0].vel[2] = _particles[i1].vel[2];
                for (int i=1; i<n_particle_unused; i++) {
                    int k = particle_index_unused[i];
                    binarytree[i].setMembers((Tparticle*)&(binarytree[i-1]), &(_particles[k]), binary_tree_index, k);
                    binarytree[i].m1 = binarytree[i-1].mass;
                    binarytree[i].m2 = 0.0;
                    binarytree[i].mass = binarytree[i].m1;
                    binarytree[i].pos[0] = binarytree[i-1].pos[0];
                    binarytree[i].pos[1] = binarytree[i-1].pos[1];
                    binarytree[i].pos[2] = binarytree[i-1].pos[2];
                    binarytree[i].vel[0] = binarytree[i-1].vel[0];
                    binarytree[i].vel[1] = binarytree[i-1].vel[1];
                    binarytree[i].vel[2] = binarytree[i-1].vel[2];
                }
            }
            else {
                for (int i=0; i<n_particle_unused; i++) {
                    int ilast = n_particle_real-1+i;
                    ASSERT(ilast<n_particle-1);
                    int k = particle_index_unused[i];
                    binarytree[ilast].setMembers((Tparticle*)&(binarytree[ilast-1]), &(_particles[k]), binary_tree_index, k);
                    binarytree[ilast].m1 = binarytree[ilast-1].mass;
                    binarytree[ilast].m2 = 0.0;
                    binarytree[ilast].mass = binarytree[ilast].m1;
                    binarytree[ilast].pos[0] = binarytree[ilast-1].pos[0];
                    binarytree[ilast].pos[1] = binarytree[ilast-1].pos[1];
                    binarytree[ilast].pos[2] = binarytree[ilast-1].pos[2];
                    binarytree[ilast].vel[0] = binarytree[ilast-1].vel[0];
                    binarytree[ilast].vel[1] = binarytree[ilast-1].vel[1];
                    binarytree[ilast].vel[2] = binarytree[ilast-1].vel[2];
                }
            }

#ifdef USE_CM_FRAME
            getBinaryTreeRoot().shiftToCenterOfMassFrame();
#endif

        }

        //! check binary tree member pair id, if consisent, return ture. otherwise set the member pair id
        /*! 
          Note that if it is a quadruple system (B-B), and both binaries pre-exist with correct pair ids, the
          return flag will be true even it is a newly formed quadruple system.
          @param[in] _bin: binary tree to check
          @param[in] _reset_flag: if true, reset pair id to zero
        */
        bool checkAndSetBinaryPairIDIter(BinaryTree<Tparticle>& _bin, const bool _reset_flag) {
            bool return_flag=true;
            Tparticle* p[2] = {_bin.getLeftMember(), _bin.getRightMember()};

            for (int i=0; i<2; i++) {
                if (_bin.isMemberTree(i)) {
                    return_flag = return_flag & checkAndSetBinaryPairIDIter(*_bin.getMemberAsTree(i),_reset_flag);
                }
            }
            for (int i=0; i<2; i++) {
                if (!_bin.isMemberTree(i)) {
                    auto pair_id = p[1-i]->id;
                    return_flag = return_flag & (p[i]->getBinaryPairID()==pair_id);
                    if (_reset_flag) p[i]->setBinaryPairID(0);
                    else p[i]->setBinaryPairID(pair_id);
                }
            }
            if (p[0]->id<p[1]->id) _bin.id = -abs(p[0]->id);
            else _bin.id = -abs(p[1]->id);
            return return_flag;
        }

        //! get binary id for a pair
        /*! 
          @param[in] _p: particle
          \return binary id (negative value), 0 if not a binary pair
        */
        static int getBinaryID(const Tparticle& _p) {
            auto pair_id = _p.getBinaryPairID();
            if (pair_id==0) return 0;
            else {
                if (_p.id<pair_id) return -abs(_p.id);
                else return -abs(pair_id);
            }
        }

        //! clear function
        void clear() {
            ds=0.0;
            time_offset = 0.0;
            r_break_crit=-1.0;
            fix_step_option = FixStepOption::none;
            binarytree.clear();
#ifdef AR_DEBUG_DUMP
            dump_flag=false;
#endif
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitleAscii(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"ds";
            _fout<<std::setw(_width)<<"Time_offset";
            _fout<<std::setw(_width)<<"r_break_crit";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnAscii(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<ds;
            _fout<<std::setw(_width)<<time_offset;
            _fout<<std::setw(_width)<<r_break_crit;
        }

        //! write class data to file with binary format
        /*! @param[in] _fout: FILE type file for output
         */
        void writeBinary(FILE *_fout) const {
            fwrite(&ds, sizeof(int),1,_fout);
            fwrite(&time_offset, sizeof(Float),1,_fout);
            fwrite(&r_break_crit, sizeof(Float),1,_fout);
            fwrite(&fix_step_option, sizeof(FixStepOption),1,_fout);
        }

        void printColumnBinary(std::ostream& _fout) const {
            _fout.write(reinterpret_cast<const char*>(&ds), sizeof(int));
            _fout.write(reinterpret_cast<const char*>(&time_offset), sizeof(Float));
            _fout.write(reinterpret_cast<const char*>(&r_break_crit), sizeof(Float));
            _fout.write(reinterpret_cast<const char*>(&fix_step_option), sizeof(FixStepOption));
        }

        //! read class data to file with binary format
        /*! @param[in] _fin: FILE type file for reading
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(&ds, sizeof(int),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
            rcount = fread(&time_offset, sizeof(Float),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
            rcount = fread(&r_break_crit, sizeof(Float),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
            rcount = fread(&fix_step_option, sizeof(FixStepOption),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }    

        void readBinary(std::istream& _fin) {
            _fin.read(reinterpret_cast<char*>(&ds), sizeof(int));
            if (!_fin) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1.\n";
                abort();
            }
            _fin.read(reinterpret_cast<char*>(&time_offset), sizeof(Float));
            if (!_fin) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1.\n";
                abort();
            }
            _fin.read(reinterpret_cast<char*>(&r_break_crit), sizeof(Float));
            if (!_fin) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1.\n";
                abort();
            }
            _fin.read(reinterpret_cast<char*>(&fix_step_option), sizeof(FixStepOption));
            if (!_fin) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1.\n";
                abort();
            }
        }
    };

}
