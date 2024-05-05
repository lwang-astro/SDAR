#pragma once

#include <cmath>
#include "Common/Float.h"
#include "Common/particle_group.h"
#include "AR/force.h"
#include "AR/interrupt.h"
#include "AR/information.h"
#include "particle.h"
#include "perturber.h"

    
//! a sample interaction class with newtonian acceleration
class Interaction{
public:
    Float gravitational_constant; ///> gravitational constant
    int interrupt_detection_option;    /// 2: record binary status when the pair distance is less than the sum of two members' radii, 1: merge when the pair distance is less than the sum of two members' radii; 0: no interruption
    Interaction(): gravitational_constant(Float(-1.0)), interrupt_detection_option(0) {}

    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        ASSERT(gravitational_constant>0.0);
        return true;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
    }    

    //! (Necessary) calculate inner member acceleration, potential and time transformation function gradient and factor for kick (two-body case)
    /*!
      @param[out] _f1: force for particle 1 to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard are overwritten, not accummulating old values)
      @param[out] _f2: force for particle 2
      @param[out] _epot: total inner potential energy
      @param[in] _p1: particle 1
      @param[in] _p2: particle 2
      \return the time transformation factor (gt_kick_inv) for kick step
    */
    inline Float calcInnerAccPotAndGTKickInvTwo(AR::Force& _f1, AR::Force& _f2, Float& _epot, const Particle& _p1, const Particle& _p2) {
        // acceleration
        const Float mass1 = _p1.mass;
        const Float* pos1 = _p1.pos;

        const Float mass2 = _p2.mass;
        const Float* pos2 = _p2.pos;

        Float dr[3] = {pos2[0] -pos1[0],
                       pos2[1] -pos1[1],
                       pos2[2] -pos1[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
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
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                Float inv_r = 1.0/sqrt(r2);
                Float inv_r3 = inv_r*inv_r*inv_r;
                Float gmor3 = gravitational_constant*massj*inv_r3;
                acci[0] += gmor3 * dr[0];
                acci[1] += gmor3 * dr[1];
                acci[2] += gmor3 * dr[2];

#ifdef AR_TTL                     
                Float mimjor3 = gravitational_constant*massi*gmor3;
                gtgradi[0] += mimjor3 * dr[0];
                gtgradi[1] += mimjor3 * dr[1];
                gtgradi[2] += mimjor3 * dr[2];
#endif

                Float gmor = gravitational_constant*massj*inv_r;
                poti -= gmor;
                gtki += gmor;
                    
            }
            _epot += poti * massi;
            gt_kick_inv += gtki*massi;
        }
        _epot   *= 0.5;
        return 0.5*gt_kick_inv;
    }

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*!@param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
    */
    void calcAccPert(AR::Force* _force, const Particle* _particles, const int _n_particle, const Particle& _particle_cm, const Perturber& _perturber, const Float _time) {
        for (int i=0; i<_n_particle; i++) {
            Float* acc_pert = _force[i].acc_pert;
            acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
            _force[i].pot_pert = 0.0;
        }            
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

    //! (Necessary) calculate slowdown factor based on perturbers
    /*!
      @param[out] _pert_out: perturbation 
      @param[out] _t_min_sq: timescale limit 
      @param[in] _time: physical time for prediction
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
    */
    void calcSlowDownPert(Float& _pert_out, Float& _t_min_sq, const Float& _time, const Particle& _particle_cm, const Perturber& _perturber) {
        _pert_out = 0.0;
        _t_min_sq = 0.0;
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
        if (_ekin_minus_etot==0.0&&_epot==0.0) return 0;
        else return log(_ekin_minus_etot) - log(-_epot);
    }
#endif   

    //! (Necessary) modify the orbits and interrupt check 
    /*! check the inner left binary whether their separation is smaller than particle radius sum and become close, if true, set one component stauts to merger with cm mass and the other unused with zero mass. Return the binary tree address 
      @param[in] _bin_interrupt: interrupt binary information: adr: binary tree address; time_now: current physical time; time_end: integration finishing time; status: interrupt status: change, merge,none
      @param[in] _bin: binarytree to check iteratively
     */
    void  modifyAndInterruptIter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin) {
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
            };

            if(_bin.getMemberN()==2) {
                if (p1->getBinaryInterruptState()== BinaryInterruptState::collision && 
                    p2->getBinaryInterruptState()== BinaryInterruptState::collision &&
                    (p1->time_check<_bin_interrupt.time_end || p2->time_check<_bin_interrupt.time_end) &&
                    (p1->getBinaryPairID()==p2->id||p2->getBinaryPairID()==p1->id)) merge();
                else {
                    Float radius = p1->radius + p2->radius;
                    // slowdown case
                    if (_bin.slowdown.getSlowDownFactor()>1.0) {
                        Float semi, ecc, dr, drdv;
                        _bin.particleToSemiEcc(semi, ecc, dr, drdv, *_bin.getLeftMember(), *_bin.getRightMember(), gravitational_constant);
                        Float peri = semi*(1 - ecc);
                        if (peri<radius && p1->getBinaryPairID()!=p2->id&&p2->getBinaryPairID()!=p1->id) {
                            Float ecc_anomaly  = _bin.calcEccAnomaly(dr);
                            Float mean_anomaly = _bin.calcMeanAnomaly(ecc_anomaly, ecc);
                            Float mean_motion  = sqrt(gravitational_constant*_bin.mass/(fabs(_bin.semi*_bin.semi*_bin.semi))); 
                            Float t_peri = mean_anomaly/mean_motion;
                            if (drdv<0 && t_peri<_bin_interrupt.time_end-_bin_interrupt.time_now) merge();
                            else if (semi>0||(semi<0&&drdv<0)) {
                                p1->setBinaryPairID(p2->id);
                                p2->setBinaryPairID(p1->id);
                                p1->setBinaryInterruptState(BinaryInterruptState::collision);
                                p2->setBinaryInterruptState(BinaryInterruptState::collision);
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
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
    }    
};

