#pragma once

#include <cmath>
#include "Common/Float.h"
#include "Common/particle_group.h"
#include "AR/force.h"
#include "particle.h"
#include "perturber.h"

    
//! a sample interaction class with newtonian acceleration
class Interaction{
public:
    Float gravitational_constant; ///> gravitational constant

    Interaction(): gravitational_constant(Float(-1.0)) {}

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
      \return the time transformation factor (gt_kick) for kick step
    */
    inline Float calcInnerAccPotAndGTKickTwo(AR::Force& _f1, AR::Force& _f2, Float& _epot, const Particle& _p1, const Particle& _p2) {
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

        Float gmor3_2 = gravitational_constant*mass1*inv_r3;
        acc2[0] = - gmor3_2 * dr[0];
        acc2[1] = - gmor3_2 * dr[1];
        acc2[2] = - gmor3_2 * dr[2];

        Float gm1m2 = gravitational_constant*mass1*mass2;

#ifdef AR_TTL 
        // trans formation function gradient
        Float gm1m2or3 = gm1m2*inv_r3;
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
        Float gt_kick = 1.0/gm1m2or;

        return gt_kick;
    }

    //! calculate inner member acceleration, potential and time transformation function gradient and factor for kick
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the time transformation factor (gt_kick) for kick step
    */
    inline Float calcInnerAccPotAndGTKick(AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
        _epot = Float(0.0);
#ifdef AR_TTL_GT_MULTI
        Float gt_kick_inv = Float(1.0);
#else
        Float gt_kick_inv = Float(0.0);
#endif
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
#ifdef AR_TTL_GT_MULTI
            Float gtki = Float(1.0);
#else
            Float gtki = Float(0.0);
#endif

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
#ifdef AR_TTL_GT_MULTI
                Float inv_r2 = inv_r*inv_r;
                gtgradi[0] += inv_r2 * dr[0];
                gtgradi[1] += inv_r2 * dr[1];
                gtgradi[2] += inv_r2 * dr[2];
#else
                Float mimjor3 = gravitational_constant*massi*gmor3;
                gtgradi[0] += mimjor3 * dr[0];
                gtgradi[1] += mimjor3 * dr[1];
                gtgradi[2] += mimjor3 * dr[2];
#endif
#endif

                Float gmor = gravitational_constant*massj*inv_r;
                poti -= gmor;
#ifdef AR_TTL_GT_MULTI
                gtki *= inv_r;
#else
                gtki += gmor;
#endif
                    
            }
            _epot += poti * massi;
#ifdef AR_TTL_GT_MULTI
            gt_kick_inv *= gtki*massi;
#else
            gt_kick_inv += gtki*massi;
#endif
        }
        _epot   *= 0.5;
#ifdef AR_TTL_GT_MULTI
        gt_kick_inv = sqrt(gt_kick_inv);
        for (int i=0; i<_n_particle; i++) {
            Float* gtgradi = _force[i].gtgrad;
            gtgradi[0] *= gt_kick_inv;
            gtgradi[1] *= gt_kick_inv;
            gtgradi[2] *= gt_kick_inv;
        }
        return 1.0/gt_kick_inv;
#else
        return 2.0/gt_kick_inv;
#endif
    }

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*! The Force class acc_pert should be updated
      @param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[out] _epot: potential 
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
      \return time transformation factor for kick
    */
    Float calcAccPotAndGTKick(AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle, const Particle& _particle_cm, const Perturber& _perturber, const Float _time) {
        Float gt_kick;
        if (_n_particle==2) gt_kick = calcInnerAccPotAndGTKickTwo(_force[0], _force[1], _epot, _particles[0], _particles[1]);
        else gt_kick = calcInnerAccPotAndGTKick(_force, _epot, _particles, _n_particle);

        for (int i=0; i<_n_particle; i++) {
            Float* acc_pert = _force[i].acc_pert;
            acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
        }            
        return gt_kick;
    }    

    //! calculate perturbation from binary tree
    Float calcPertFromBinary(const COMM::BinaryTree<Particle>& _bin) {
        Float apo = _bin.semi*(1.0+_bin.ecc);
        Float apo2 = apo*apo;
#ifdef AR_SLOWDOWN_PERT_R4
        return (_bin.m1*_bin.m2)/(apo2*apo2);
#else
        return (_bin.m1*_bin.m2)/(apo2*apo);
#endif
    }

    //! calculate perturbation from distance to perturber and masses of particle and perturber 
    inline Float calcPertFromMR(const Float _r, const Float _mp, const Float _mpert) {
        Float r2 = _r*_r;
#ifdef AR_SLOWDOWN_PERT_R4
        return _mp*_mpert/(r2*r2);
#else
        return (_mp*_mpert)/(r2*_r);
#endif
    }

#ifdef AR_SLOWDOWN_INNER
    //! (Necessary) calculate slowdown factor for inner binary based on other particles and slowdown of system c.m.
    /*!
      @param[in,out] _slowdown: slowdown paramters for inner binary.
      @param[in] _slowdown_cm: slowdown paramters of system c.m..
      @param[in] _bin_root: binary tree root
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
     */
    void calcSlowDownInnerBinary(AR::SlowDown& _slowdown, const AR::SlowDown& _slowdown_cm, const COMM::BinaryTree<Particle>& _bin_root, const Particle* _particles, const int _n_particle) {
        _slowdown.pert_in = calcPertFromBinary(_bin_root);
        _slowdown.period = _bin_root.period;
        int imask[2] = {_bin_root.getMemberIndex(0), _bin_root.getMemberIndex(1)};
        const Float* xcm = _bin_root.pos;
        const Float  mcm = _bin_root.mass;
        Float pert_pot = 0.0;
#ifdef AR_SLOWDOWN_TIMESCALE
        const Float* vcm = _bin_root.vel;
        Float trf2_min = NUMERIC_FLOAT_MAX;
        Float mvor[3] = {0.0,0.0,0.0};
        Float mtot=0.0;
#endif
        for (int i=0; i<_n_particle; i++) {
            if (i==imask[0]||i==imask[1]) continue;
            const Float* xp = _particles[i].pos;

            const Float mj = _particles[i].mass;

            Float dr[3] = {xp[0] - xcm[0],
                           xp[1] - xcm[1],
                           xp[2] - xcm[2]};
            Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            Float r = sqrt(r2);
            pert_pot += calcPertFromMR(r, mcm, mj);

#ifdef AR_SLOWDOWN_TIMESCALE
            const Float* vp = _particles[i].vel;
            Float dv[3] = {vp[0] - vcm[0],
                           vp[1] - vcm[1],
                           vp[2] - vcm[2]};

            // velocity dependent method 
            // m_tot / |\sum m_j /|r_j| * v_j|
            Float mor = mj/r;
            mvor[0] += mor*dv[0];
            mvor[1] += mor*dv[1];
            mvor[2] += mor*dv[2];
            mtot += mj;

            // force dependent method
            // min sqrt(r^3/(G m))
            Float gmor3 = (mj+mcm)*r*r2/(gravitational_constant*mj*mcm);
            trf2_min =  std::min(trf2_min, gmor3);

#endif
        }            
        _slowdown.pert_out = pert_pot + _slowdown_cm.pert_out;
#ifdef AR_SLOWDOWN_TIMESCALE
        // velocity dependent method
        Float trv_ave = mtot/sqrt(mvor[0]*mvor[0] + mvor[1]*mvor[1] + mvor[2]*mvor[2]);
        // get min of velocity and force dependent values
        Float t_min = std::min(trv_ave, sqrt(trf2_min));
        _slowdown.timescale = 0.1*std::min(_slowdown.getTimescaleMax(), t_min);
#else
        _slowdown.timescale = _slowdown.getTimescaleMax();
#endif

        _slowdown.calcSlowDownFactor();
        //_slowdown.setSlowDownFactor(10.0);
    }
#endif

    //! (Necessary) calculate slowdown factor based on perturbers
    /*!
      @param[in,out] _slowdown: slowdown paramters.
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _bin_root: binary tree root
      @param[in] _perturber: pertuber container
    */
    void calcSlowDownPert(AR::SlowDown& _slowdown, const Particle& _particle_cm, const COMM::BinaryTree<Particle>& _bin_root, const Perturber& _perturber) {
        
        // slowdown inner perturbation: m1*m2/apo_in^4
        //Float apo_in = _bin_root.semi*(1.0+_bin_root.ecc);
        //Float apo_in2 = apo_in*apo_in;
        //_slowdown.pert_in = _bin_root.m1*_bin_root.m2/(apo_in2*apo_in2);
        //_slowdown.period  = _bin_root.period;
        _slowdown.pert_in = 0.0; // suppress slowdown 

        _slowdown.pert_out = 0.0;
        _slowdown.timescale = _slowdown.getTimescaleMax();
        _slowdown.calcSlowDownFactor();
    }



#ifndef AR_TTL
    //! (Necessary) calcualte the time transformation factor for drift
    /*! The time transformation factor for drift only depends on (kinetic energy - total energy)
      @param[in] _ekin_minus_etot: ekin - etot
    */
    Float calcGTDrift(Float _ekin_minus_etot) {
        return 1.0/_ekin_minus_etot;
    }

    //! (Necessary) calculate the time transformed Hamiltonian
    /*! calculate the time transformed Hamiltonian
      @param[in] _ekin_minus_etot: ekin - etot
    */
    Float calcH(Float _ekin_minus_etot, Float _epot) {
        return log(_ekin_minus_etot) - log(-_epot);
    }
#endif   

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

