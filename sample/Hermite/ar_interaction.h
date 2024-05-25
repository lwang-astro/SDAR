#pragma once

#include <cmath>
#include "Common/Float.h"
#include "AR/force.h"
#include "Hermite/neighbor.h"
#include "particle.h"

//! a sample interaction class with newtonian acceleration
class ARInteraction{
public:
    typedef H4::ParticleH4<Particle> H4Ptcl;

    Float eps_sq; // softening 
    Float gravitational_constant; ///> gravitational constant
    int interrupt_detection_option;    /// 2: record binary status when the pair distance is less than the sum of two members' radii, 1: merge when the pair distance is less than the sum of two members' radii; 0: no interruption
    ARInteraction(): eps_sq(Float(-1.0)), gravitational_constant(Float(-1.0)), interrupt_detection_option(0) {}

    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        ASSERT(eps_sq>=0.0);
        ASSERT(gravitational_constant>0.0);
        return true;
    }        

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"eps_sq: "<<eps_sq<<std::endl
             <<"G     : "<<gravitational_constant<<std::endl
             <<"Interrupt_option: "<<interrupt_detection_option<<std::endl;
    }    

    //! (Necessary) calculate inner member acceleration, potential and time transformation function gradient and factor for kick (two-body case)
    /*!
      @param[out] _f1: force for particle 1 to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard are overwritten, not accummulating old values)
      @param[out] _f2: force for particle 2
      @param[out] _epot: total inner potential energy
      @param[in] _p1: particle 1
      @param[in] _p2: particle 2
      @param[in] _pos_offset: position offset need to be added to calculate dr
      \return the time transformation factor (gt_kick_inv) for kick step
    */
    inline Float calcInnerAccPotAndGTKickInvTwo(AR::Force& _f1, AR::Force& _f2, Float& _epot, const Particle& _p1, const Particle& _p2, const Float* _pos_offset) {
        // acceleration
        const Float mass1 = _p1.mass;
        const Float* pos1 = _p1.pos;

        const Float mass2 = _p2.mass;
        const Float* pos2 = _p2.pos;

#ifdef USE_CM_FRAME
        Float dr[3] = {pos2[0] -pos1[0] + _pos_offset[0],
                       pos2[1] -pos1[1] + _pos_offset[1],
                       pos2[2] -pos1[2] + _pos_offset[2]};
#else
        Float dr[3] = {pos2[0] -pos1[0],
                       pos2[1] -pos1[1],
                       pos2[2] -pos1[2]};
#endif

        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
        ASSERT(r2>0.0);
        Float inv_r = 1.0/sqrt(r2);
        Float inv_r3 = inv_r*inv_r*inv_r;

        Float* acc1 = _f1.acc_in;
        Float* acc2 = _f2.acc_in;

        Float gmor3_1 = gravitational_constant*mass2*inv_r3;
        acc1[0] = gmor3_1 * dr[0];
        acc1[1] = gmor3_1 * dr[1];
        acc1[2] = gmor3_1 * dr[2];

        _f1.pot_in = -gravitational_constant*mass2*inv_r;

        Float gmor3_2 = gravitational_constant*mass1*inv_r3;
        acc2[0] = - gmor3_2 * dr[0];
        acc2[1] = - gmor3_2 * dr[1];
        acc2[2] = - gmor3_2 * dr[2];

        _f2.pot_in = -gravitational_constant*mass1*inv_r;

        Float gm1m2 = gravitational_constant*mass1*mass2;

#ifdef AR_TTL 
        // trans formation function gradient
#ifdef AR_TIME_FUNCTION_MUL_POT
        Float gm1m2or3 = inv_r*inv_r; // gt_kick_inv will be multiplied latter, thus only need 1/r^2
#else
        Float gm1m2or3 = gm1m2*inv_r3;
#endif
        Float* gtgrad1 = _f1.gtgrad;
        Float* gtgrad2 = _f2.gtgrad;

        gtgrad1[0] = gm1m2or3 * dr[0];
        gtgrad1[1] = gm1m2or3 * dr[1];
        gtgrad1[2] = gm1m2or3 * dr[2];

        gtgrad2[0] = - gtgrad1[0];
        gtgrad2[1] = - gtgrad1[1];
        gtgrad2[2] = - gtgrad1[2];
#endif

        // potential energy
        Float gm1m2or = gm1m2*inv_r;
        _epot = - gm1m2or;

        // transformation factor for kick
        Float gt_kick_inv = gm1m2or;

        return gt_kick_inv;
    }

    //! calculate inner member acceleration, potential and time transformation function gradient and factor for kick
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the time transformation factor (gt_kick_inv) for kick step
    */
    inline Float calcInnerAccPotAndGTKickInv(AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
        _epot = Float(0.0);
        Float gt_kick_inv = Float(0.0);
        for (int i=0; i<_n_particle; i++) {
            const Float massi = _particles[i].mass;
            const Float* posi = _particles[i].pos;
            Float* acci = _force[i].acc_in;
            acci[0] = acci[1] = acci[2] = Float(0.0);

#ifdef AR_TTL 
            Float* gtgradi = _force[i].gtgrad;
            gtgradi[0] = gtgradi[1] = gtgradi[2] = Float(0.0);
#endif

            Float poti = Float(0.0);
            Float gtki = Float(0.0);

            for (int j=0; j<_n_particle; j++) {
                if (i==j) continue;
                const Float massj = _particles[j].mass;
                const Float* posj = _particles[j].pos; 
                Float dr[3] = {posj[0] -posi[0],
                               posj[1] -posi[1],
                               posj[2] -posi[2]};
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                Float inv_r = 1.0/sqrt(r2);
                Float inv_r3 = inv_r*inv_r*inv_r;
                Float gmor3 = gravitational_constant*massj*inv_r3;
                acci[0] += gmor3 * dr[0];
                acci[1] += gmor3 * dr[1];
                acci[2] += gmor3 * dr[2];

#ifdef AR_TTL                     
                Float gmimjor3 = massi*gmor3;
                gtgradi[0] += gmimjor3 * dr[0];
                gtgradi[1] += gmimjor3 * dr[1];
                gtgradi[2] += gmimjor3 * dr[2];
#endif

                Float gmor = gravitational_constant*massj*inv_r;
                poti -= gmor;
                gtki += gmor;
                    
            }
            _epot += poti * massi;
            gt_kick_inv += gtki * massi;
        }
        _epot   *= 0.5;
        gt_kick_inv *= 0.5;

        return gt_kick_inv;
    }

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*!@param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
    */
    void calcAccPert(AR::Force* _force, const Particle* _particles, const int _n_particle, const H4Ptcl& _particle_cm, const H4::Neighbor<Particle>& _perturber, const Float _time) {
        static const Float inv3 = 1.0 / 3.0;

        const int n_pert = _perturber.neighbor_address.getSize();
        if (n_pert>0) {

            auto* pert_adr = _perturber.neighbor_address.getDataAddress();

            Float xp[n_pert][3], xcm[3], m[n_pert];

            for (int j=0; j<n_pert; j++) {
                H4::NBAdr<Particle>::Single* pertj;
                if (pert_adr[j].type==H4::NBType::group) pertj = &(((H4::NBAdr<Particle>::Group*)pert_adr[j].adr)->cm);
                else pertj = (H4::NBAdr<Particle>::Single*)pert_adr[j].adr;
                Float dt = _time - pertj->time;
                ASSERT(dt>=0.0);
                xp[j][0] = pertj->pos[0] + dt*(pertj->vel[0] + 0.5*dt*(pertj->acc0[0] + inv3*dt*pertj->acc1[0]));
                xp[j][1] = pertj->pos[1] + dt*(pertj->vel[1] + 0.5*dt*(pertj->acc0[1] + inv3*dt*pertj->acc1[1]));
                xp[j][2] = pertj->pos[2] + dt*(pertj->vel[2] + 0.5*dt*(pertj->acc0[2] + inv3*dt*pertj->acc1[2]));

                m[j] = pertj->mass;
            }

            Float dt = _time - _particle_cm.time;
            ASSERT(dt>=0.0);
            xcm[0] = _particle_cm.pos[0] + dt*(_particle_cm.vel[0] + 0.5*dt*(_particle_cm.acc0[0] + inv3*dt*_particle_cm.acc1[0]));
            xcm[1] = _particle_cm.pos[1] + dt*(_particle_cm.vel[1] + 0.5*dt*(_particle_cm.acc0[1] + inv3*dt*_particle_cm.acc1[1]));
            xcm[2] = _particle_cm.pos[2] + dt*(_particle_cm.vel[2] + 0.5*dt*(_particle_cm.acc0[2] + inv3*dt*_particle_cm.acc1[2]));

            Float acc_pert_cm[3]={0.0, 0.0, 0.0};
            Float pot_pert_cm=0.0;
            Float mcm = 0.0;
            if (_perturber.need_resolve_flag) {
                // calculate component perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    Float& pot_pert = _force[i].pot_pert;
                    const auto& pi = _particles[i];
                    acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
                    pot_pert = 0.0;

                    Float xi[3];
                    xi[0] = pi.pos[0] + xcm[0];
                    xi[1] = pi.pos[1] + xcm[1];
                    xi[2] = pi.pos[2] + xcm[2];

                    for (int j=0; j<n_pert; j++) {
                        Float dr[3] = {xp[j][0] - xi[0],
                                       xp[j][1] - xi[1],
                                       xp[j][2] - xi[2]};
                        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                        Float r  = sqrt(r2);
                        Float r3 = r*r2;
                        Float gm = gravitational_constant*m[j];
                        Float gmor3 = gm/r3;

                        acc_pert[0] += gmor3 * dr[0];
                        acc_pert[1] += gmor3 * dr[1];
                        acc_pert[2] += gmor3 * dr[2];

                        pot_pert += - gm/r;
                    }


                    acc_pert_cm[0] += pi.mass *acc_pert[0];
                    acc_pert_cm[1] += pi.mass *acc_pert[1];
                    acc_pert_cm[2] += pi.mass *acc_pert[2];

                    pot_pert_cm += pi.mass*pot_pert;
                    mcm += pi.mass;
                }
//#ifdef ARC_DEBUG
//                ASSERT(abs(mcm-_particles_cm.mass)<1e-10);
//#endif
                // get cm perturbation
                acc_pert_cm[0] /= mcm;
                acc_pert_cm[1] /= mcm;
                acc_pert_cm[2] /= mcm;

                pot_pert_cm /= mcm;

                // remove cm. perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    Float& pot_pert = _force[i].pot_pert;
                    acc_pert[0] -= acc_pert_cm[0]; 
                    acc_pert[1] -= acc_pert_cm[1];        
                    acc_pert[2] -= acc_pert_cm[2]; 

                    pot_pert -= pot_pert_cm;
                }
            }
            else {
                // first calculate c.m. acceleration and tidal perturbation
                for (int j=0; j<n_pert; j++) {
                    Float dr[3] = {xp[j][0] - xcm[0],
                                   xp[j][1] - xcm[1],
                                   xp[j][2] - xcm[2]};
                    Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                    Float r  = sqrt(r2);
                    Float r3 = r*r2;
                    Float gm = gravitational_constant*m[j];
                    Float gmor3 = gm/r3;

                    acc_pert_cm[0] += gmor3 * dr[0];
                    acc_pert_cm[1] += gmor3 * dr[1];
                    acc_pert_cm[2] += gmor3 * dr[2];

                    pot_pert_cm += - gm/r;
                }

                // calculate component perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    Float& pot_pert = _force[i].pot_pert;
                    const auto& pi = _particles[i];
                    acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
                    pot_pert = 0.0;

                    Float xi[3];
                    xi[0] = pi.pos[0] + xcm[0];
                    xi[1] = pi.pos[1] + xcm[1];
                    xi[2] = pi.pos[2] + xcm[2];

                    // remove c.m. perturbation acc
                    acc_pert[0] = -acc_pert_cm[0]; 
                    acc_pert[1] = -acc_pert_cm[1];        
                    acc_pert[2] = -acc_pert_cm[2]; 

                    pot_pert = - pot_pert_cm;

                    for (int j=0; j<n_pert; j++) {
                        Float dr[3] = {xp[j][0] - xi[0],
                                       xp[j][1] - xi[1],
                                       xp[j][2] - xi[2]};
                        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                        Float r  = sqrt(r2);
                        Float r3 = r*r2;
                        Float gm = gravitational_constant*m[j];
                        Float gmor3 = gm/r3;

                        acc_pert[0] += gmor3 * dr[0];
                        acc_pert[1] += gmor3 * dr[1];
                        acc_pert[2] += gmor3 * dr[2];

                        pot_pert += - gm/r;
                    }
                }

            }
        }
    }

    //! calculate perturbation from c.m. acceleration
    Float calcPertFromForce(const Float* _force, const Float _mp, const Float _mpert) {
        Float force2 = _force[0]*_force[0]+_force[1]*_force[1]+_force[2]*_force[2];
#ifdef AR_SLOWDOWN_PERT_R4
        return force2/(gravitational_constant*_mp*_mpert);
#else
        Float force = sqrt(force2)/gravitational_constant;
        return sqrt(force/(_mp*_mpert))*force;
#endif
    }

    //! calculate perturbation from binary tree
    static Float calcPertFromBinary(const AR::BinaryTree<Particle>& _bin) {
        Float apo = _bin.semi*(1.0+_bin.ecc);
        Float apo2 = apo*apo;
#ifdef AR_SLOWDOWN_PERT_R4
        return (_bin.m1*_bin.m2)/(apo2*apo2);
#else
        return (_bin.m1*_bin.m2)/(apo2*apo);
#endif
    }

    //! calculate perturbation from distance to perturber and masses of particle and perturber 
    static Float calcPertFromMR(const Float _r, const Float _mp, const Float _mpert) {
        Float r2 = _r*_r;
#ifdef AR_SLOWDOWN_PERT_R4
        return _mp*_mpert/(r2*r2);
#else
        return (_mp*_mpert)/(r2*_r);
#endif
    }

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
    //! calculate slowdown perturbation and timescale from particle j to particle i
    /*! 
      @param[out] _pert_out: perturbation from particle j
      @param[out] _t_min_sq: timescale limit from particle j
      @param[in] _pi: particle i (cm of binary)
      @param[in] _pj: particle j 
     */
    void calcSlowDownPertOne(Float& _pert_out, Float& _t_min_sq, const Particle& pi, const Particle& pj) {
        Float dr[3] = {pj.pos[0] - pi.pos[0],
                       pj.pos[1] - pi.pos[1],
                       pj.pos[2] - pi.pos[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float r = sqrt(r2);
        _pert_out += calcPertFromMR(r, pi.mass, pj.mass);
            
#ifdef AR_SLOWDOWN_TIMESCALE
        Float dv[3] = {pj.vel[0] - pi.vel[0],
                       pj.vel[1] - pi.vel[1],
                       pj.vel[2] - pi.vel[2]};

        Float v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
        Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];

        // identify whether hyperbolic or closed orbit
        Float gm = gravitational_constant*(pi.mass+pj.mass);
        Float semi = 1.0/(2.0/r - v2/gm);

        //hyperbolic, directly use velocity v
        if (semi<0) 
            _t_min_sq = std::min(_t_min_sq, r2/v2);
        else {
            if (r<semi) {
                // avoid decrese of vr once the orbit pass, calculate vr max at E=pi/2 (r==semi)
                // vr_max = sqrt(er*(drdv^2*er + r*vcr2^2))/(G(m1+m2)r)
                Float rv2 = r*v2;
                Float er = 2*gm - rv2;
                Float vcr2 = gm - rv2;
                Float vrmax_sq = er*(drdv*drdv*er + r*vcr2*vcr2)/(gm*gm*r2);
                _t_min_sq = std::min(_t_min_sq, semi*semi/vrmax_sq);
            }
            else {
                // r/vr
                Float rovr = r2/abs(drdv);
                _t_min_sq = std::min(_t_min_sq, rovr*rovr);
            }
        }

        // force dependent method
        // min sqrt(r^3/(G m))
        //Float gmor3 = (mp+mcm)*r*r2/(sdt->G*mp*mcm);
        //sdt->trf2_min =  std::min(sdt->trf2_min, gmor3);
#endif
    }

    //! (Necessary) calculate slowdown perturbation and timescale
    /*!
      @param[out] _pert_out: perturbation 
      @param[out] _t_min_sq: timescale limit 
      @param[in] _time: physical time for prediction
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
    */
    void calcSlowDownPert(Float& _pert_out, Float& _t_min_sq, const Float& _time, const H4Ptcl& _particle_cm, const H4::Neighbor<Particle>& _perturber) {
        static const Float inv3 = 1.0 / 3.0;

        const int n_pert = _perturber.neighbor_address.getSize();

        if (n_pert>0) {

            auto* pert_adr = _perturber.neighbor_address.getDataAddress();

            Float xp[3], xcm[3];
            Float dt = _time - _particle_cm.time;
            //ASSERT(dt>=0.0);
            xcm[0] = _particle_cm.pos[0] + dt*(_particle_cm.vel[0] + 0.5*dt*(_particle_cm.acc0[0] + inv3*dt*_particle_cm.acc1[0]));
            xcm[1] = _particle_cm.pos[1] + dt*(_particle_cm.vel[1] + 0.5*dt*(_particle_cm.acc0[1] + inv3*dt*_particle_cm.acc1[1]));
            xcm[2] = _particle_cm.pos[2] + dt*(_particle_cm.vel[2] + 0.5*dt*(_particle_cm.acc0[2] + inv3*dt*_particle_cm.acc1[2]));

            Float mcm = _particle_cm.mass;

#ifdef AR_SLOWDOWN_TIMESCALE
            // velocity dependent method 
            Float vp[3], vcm[3];

            vcm[0] = _particle_cm.vel[0] + dt*(_particle_cm.acc0[0] + 0.5*dt*_particle_cm.acc1[0]);
            vcm[1] = _particle_cm.vel[1] + dt*(_particle_cm.acc0[1] + 0.5*dt*_particle_cm.acc1[1]);
            vcm[2] = _particle_cm.vel[2] + dt*(_particle_cm.acc0[2] + 0.5*dt*_particle_cm.acc1[2]);
#endif

            for (int j=0; j<n_pert; j++) {
                H4::NBAdr<Particle>::Single* pertj;
                if (pert_adr[j].type==H4::NBType::group) pertj = &(((H4::NBAdr<Particle>::Group*)pert_adr[j].adr)->cm);
                else pertj = (H4::NBAdr<Particle>::Single*)pert_adr[j].adr;
                Float dt = _time - pertj->time;
                //ASSERT(dt>=0.0);
                xp[0] = pertj->pos[0] + dt*(pertj->vel[0] + 0.5*dt*(pertj->acc0[0] + inv3*dt*pertj->acc1[0]));
                xp[1] = pertj->pos[1] + dt*(pertj->vel[1] + 0.5*dt*(pertj->acc0[1] + inv3*dt*pertj->acc1[1]));
                xp[2] = pertj->pos[2] + dt*(pertj->vel[2] + 0.5*dt*(pertj->acc0[2] + inv3*dt*pertj->acc1[2]));

                Float mj = pertj->mass;

                Float dr[3] = {xp[0] - xcm[0],
                               xp[1] - xcm[1],
                               xp[2] - xcm[2]};

                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                Float r = sqrt(r2);
                _pert_out += calcPertFromMR(r, mcm, mj);

#ifdef AR_SLOWDOWN_TIMESCALE
                // velocity dependent method 
                vp[0] = pertj->vel[0] + dt*(pertj->acc0[0] + 0.5*dt*pertj->acc1[0]);
                vp[1] = pertj->vel[1] + dt*(pertj->acc0[1] + 0.5*dt*pertj->acc1[1]);
                vp[2] = pertj->vel[2] + dt*(pertj->acc0[2] + 0.5*dt*pertj->acc1[2]);

                Float dv[3] = {vp[0] - vcm[0],
                               vp[1] - vcm[1],
                               vp[2] - vcm[2]};

                Float v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
                Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];

                // identify whether hyperbolic or closed orbit
                Float gm = gravitational_constant*(mcm+mj);
                Float semi = 1.0/(2.0/r - v2/gm);

                //hyperbolic, directly use velocity v
                if (semi<0) 
                    _t_min_sq = std::min(_t_min_sq, r2/v2);
                else {
                    if (r<semi) {
                        // avoid decrese of vr once the orbit pass, calculate vr max at E=pi/2 (r==semi)
                        // vr_max = sqrt(er*(drdv^2*er + r*vcr2^2))/(G(m1+m2)r)
                        Float rv2 = r*v2;
                        Float er = 2*gm - rv2;
                        Float vcr2 = gm - rv2;
                        Float vrmax_sq = er*(drdv*drdv*er + r*vcr2*vcr2)/(gm*gm*r2);
                        _t_min_sq = std::min(_t_min_sq, semi*semi/vrmax_sq);
                    }
                    else {
                        // r/vr
                        Float rovr = r2/abs(drdv);
                        _t_min_sq = std::min(_t_min_sq, rovr*rovr);
                    }
                }
#endif
            }
        }
    }
#endif
        

#ifndef AR_TTL
    //! (Necessary) calcualte the time transformation factor for drift
    /*! The time transformation factor for drift only depends on (kinetic energy - total energy)
      @param[in] _ekin_minus_etot: ekin - etot
    */
    Float calcGTDriftInv(Float _ekin_minus_etot) {
        return _ekin_minus_etot;
    }

    //! (Necessary) calculate the time transformed Hamiltonian
    /*! calculate the time transformed Hamiltonian
      @param[in] _ekin_minus_etot: ekin - etot
    */
    Float calcH(Float _ekin_minus_etot, Float _epot) {
        return log(_ekin_minus_etot) - log(-_epot);
    }
#endif   

    //! modify one particle function
    template <class Tparticle>
    int modifyOneParticle(Tparticle& _p, const Float& _time_now, const Float& _time_end) {
        // sample of mass loss
        /*
        if (_p.time_check<_time_end) {
            _p.dm = -_p.mass*0.1;
            _p.mass *=0.9;
            _p.time_check = _time_end+1e-5;
            return true;
        }
        */
        return 0;
    }

    //! (Necessary) modify the orbits and interrupt check 
    /*! check the inner left binary whether their separation is smaller than particle radius sum and become close, if true, set one component stauts to merger with cm mass and the other unused with zero mass. Return the binary tree address 
      @param[in] _bin_interrupt: interrupt binary information: adr: binary tree address; time_now: current physical time; time_end: integration finishing time; status: interrupt status: change, merge,none
      @param[in] _bin: binarytree to check iteratively
     */
    void modifyAndInterruptIter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin) {
        if (_bin_interrupt.status==AR::InterruptStatus::none && interrupt_detection_option>0) {
            auto* p1 = _bin.getLeftMember();
            auto* p2 = _bin.getRightMember();

            auto merge = [&]() {
                _bin_interrupt.setBinaryTreeAddress(&_bin);
                _bin_interrupt.status = AR::InterruptStatus::merge;
                if (interrupt_detection_option==1) {
                    Float mcm = p1->mass + p2->mass;
                    for (int k=0; k<3; k++) {
                        p1->pos[k] = (p1->mass*p1->pos[k] + p2->mass*p2->pos[k])/mcm;
                        p1->vel[k] = (p1->mass*p1->vel[k] + p2->mass*p2->vel[k])/mcm;
                    }
                    p1->setBinaryInterruptState(BinaryInterruptState::none);
                    p2->setBinaryInterruptState(BinaryInterruptState::none);
                    p1->dm = mcm - p1->mass;
                    p2->dm = -p2->mass;
                    p1->mass = mcm;
                    p2->mass = 0.0;
                }
                else if (interrupt_detection_option == 2) {
                    p1->setBinaryInterruptState(BinaryInterruptState::collision);
                    p2->setBinaryInterruptState(BinaryInterruptState::collision);
                }
            };

            if(_bin.getMemberN()==2) {
                if (p1->getBinaryInterruptState()== BinaryInterruptState::delaycollision && 
                    p2->getBinaryInterruptState()== BinaryInterruptState::delaycollision &&
                    (p1->time_check<_bin_interrupt.time_end || p2->time_check<_bin_interrupt.time_end) &&
                    (p1->getBinaryPairID()==p2->id||p2->getBinaryPairID()==p1->id)) merge();
                else if (p1->getBinaryInterruptState() != BinaryInterruptState::collision && p2->getBinaryInterruptState() != BinaryInterruptState::collision) {
                    Float radius = p1->radius + p2->radius;
                    // slowdown case
                    if (_bin.slowdown.getSlowDownFactor()>1.0) {
                        Float drdv;
                        _bin.particleToSemiEcc(_bin.semi, _bin.ecc, _bin.r, drdv, *_bin.getLeftMember(), *_bin.getRightMember(), gravitational_constant);
                        Float peri = _bin.semi*(1 - _bin.ecc);
                        if (peri<radius && p1->getBinaryPairID()!=p2->id&&p2->getBinaryPairID()!=p1->id) {
                            Float ecc_anomaly  = _bin.calcEccAnomaly(_bin.r);
                            Float mean_anomaly = _bin.calcMeanAnomaly(ecc_anomaly, _bin.ecc);
                            Float mean_motion  = sqrt(gravitational_constant*_bin.mass/(fabs(_bin.semi*_bin.semi*_bin.semi))); 
                            Float t_peri = mean_anomaly/mean_motion;
                            if (drdv<0 && t_peri<_bin_interrupt.time_end-_bin_interrupt.time_now) merge();
                            else if (_bin.semi>0||(_bin.semi<0&&drdv<0)) {
                                p1->setBinaryPairID(p2->id);
                                p2->setBinaryPairID(p1->id);
                                p1->setBinaryInterruptState(BinaryInterruptState::delaycollision);
                                p2->setBinaryInterruptState(BinaryInterruptState::delaycollision);
                                p1->time_check = std::min(p1->time_check, _bin_interrupt.time_now + drdv<0 ? t_peri : (_bin.period - t_peri));
                                p2->time_check = std::min(p1->time_check, p2->time_check);
                            }
                        }
                    }
                    else { // no slowdown case, check separation directly
                        Float dr[3] = {p1->pos[0] - p2->pos[0], 
                                       p1->pos[1] - p2->pos[1], 
                                       p1->pos[2] - p2->pos[2]};
                        Float dr2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                        if (dr2<radius*radius) merge();
                    }
                }
            }
            else {
                for (int k=0; k<2; k++) 
                    if (_bin.isMemberTree(k)) modifyAndInterruptIter(_bin_interrupt, *_bin.getMemberAsTree(k));
            }
        }
    }


    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
        fwrite(this, sizeof(*this),1,_fp);
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this), 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }    
};

