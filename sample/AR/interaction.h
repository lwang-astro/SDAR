#pragma once

#include <cmath>
#include "Common/Float.h"
#include "Common/particle_group.h"
#include "AR/force.h"
#include "particle.h"
#include "perturber.h"

    
//! a sample interaction class with newtonian acceleration
class Interaction{
private:

    //! (Necessary) calculate inner member acceleration, potential and time transformation function gradient and factor for kick (two-body case)
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the time transformation factor (gt_kick) for kick step
    */
    inline Float calcInnerAccPotAndGTKickTwo(AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
        ASSERT(_n_particle==2);

        // acceleration
        const Float mass1 = _particles[0].mass;
        const Float* pos1 = _particles[0].pos;

        const Float mass2 = _particles[1].mass;
        const Float* pos2 = _particles[1].pos;

        Float dr[3] = {pos2[0] -pos1[0],
                       pos2[1] -pos1[1],
                       pos2[2] -pos1[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float inv_r = 1.0/sqrt(r2);
        Float inv_r3 = inv_r*inv_r*inv_r;

        Float* acc1 = _force[0].acc_in;
        Float* acc2 = _force[1].acc_in;

        Float mor3_1 = mass2*inv_r3;
        acc1[0] = mor3_1 * dr[0];
        acc1[1] = mor3_1 * dr[1];
        acc1[2] = mor3_1 * dr[2];

        Float mor3_2 = mass1*inv_r3;
        acc2[0] = - mor3_2 * dr[0];
        acc2[1] = - mor3_2 * dr[1];
        acc2[2] = - mor3_2 * dr[2];

        Float m1m2 = mass1*mass2;

#ifdef AR_TTL 
        // trans formation function gradient
        Float m1m2or3 = m1m2*inv_r3;
        Float* gtgrad1 = _force[0].gtgrad;
        Float* gtgrad2 = _force[1].gtgrad;

        gtgrad1[0] = m1m2or3 * dr[0];
        gtgrad1[1] = m1m2or3 * dr[1];
        gtgrad1[2] = m1m2or3 * dr[2];

        gtgrad2[0] = - gtgrad1[0];
        gtgrad2[1] = - gtgrad1[1];
        gtgrad2[2] = - gtgrad1[2];
#endif

        // potential energy
        Float m1m2or = m1m2*inv_r;
        _epot = - m1m2or;

        // transformation factor for kick
        Float gt_kick = 1.0/m1m2or;

        return gt_kick;
    }

    //! (Necessary) calculate inner member acceleration, potential and time transformation function gradient and factor for kick
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the time transformation factor (gt_kick) for kick step
    */
    inline Float calcInnerAccPotAndGTKick(AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
        _epot = Float(0.0);
        Float gt_kick = Float(0.0);
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
                Float mor3 = massj*inv_r3;
                acci[0] += mor3 * dr[0];
                acci[1] += mor3 * dr[1];
                acci[2] += mor3 * dr[2];

#ifdef AR_TTL                     
                Float mimjor3 = massi*mor3;
                gtgradi[0] += mimjor3 * dr[0];
                gtgradi[1] += mimjor3 * dr[1];
                gtgradi[2] += mimjor3 * dr[2];
#endif

                Float mor = massj*inv_r;
                poti -= mor;
                gtki += mor;
                    
            }
            _epot += poti * massi;
            gt_kick += gtki * massi;
        }
        _epot   *= 0.5;
        gt_kick = 2.0/gt_kick;

        return gt_kick;
    }

public:
    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        return true;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
    }    

#ifdef SLOWDOWN_INTEGRATE

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*! The AR::Force class acc_pert should be updated
      @param[in,out] _slowdown: slowdown class to store perturbation (member: pert: perturbation ([m]/[r^3]); timescale: limited timescale to calculate maximum slowdown factor)
      @param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[out] _epot: potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _bin_root: binary tree root
      @param[in] _perturber: pertuber container
      \return perturbation energy to calculate slowdown factor
    */
    Float calcAccPotGTKickAndSlowDownPert(AR::SlowDown& _slowdown, AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle, const Particle& _particle_cm, const COMM::BinaryTree<Particle>& _bin_root, const Perturber& _perturber) {

        Float gt_kick;
        if (_n_particle==2) gt_kick = calcInnerAccPotAndGTKickTwo(_force, _epot, _particles, _n_particle);
        else gt_kick = calcInnerAccPotAndGTKick(_force, _epot, _particles, _n_particle);

        // slowdown inner perturbation: m1*m2/apo_in^4
        //Float apo_in = _bin_root.semi*(1.0+_bin_root.ecc);
        //Float apo_in2 = apo_in*apo_in;
        //_slowdown.pert_in = _bin_root.m1*_bin_root.m2/(apo_in2*apo_in2);
        //_slowdown.period  = _bin_root.period;
        _slowdown.pert_in = 0.0; // suppress slowdown 

        for (int i=0; i<_n_particle; i++) {
            Float* acc_pert = _force[i].acc_pert;
            acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
        }            
        _slowdown.pert_out = 0;
        _slowdown.timescale = _slowdown.getTimescaleMax();

        return gt_kick;
    }

#else

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*! The Force class acc_pert should be updated
      @param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[out] _epot: potential 
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
      \return perturbation energy to calculate slowdown factor
    */
    Float calcAccPotAndGTKick(AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle, const Particle& _particle_cm, const Perturber& _perturber, const Float _time) {
        
        Float gt_kick;
        if (_n_particle==2) gt_kick = calcInnerAccPotAndGTKickTwo(_force, _epot, _particles, _n_particle);
        else gt_kick = calcInnerAccPotAndGTKick(_force, _epot, _particles, _n_particle);

        for (int i=0; i<_n_particle; i++) {
            Float* acc_pert = _force[i].acc_pert;
            acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
        }            
        return gt_kick;
    }    

    //! (Necessary) calculate slowdown perturbation and timescale
    /*!
      @param[in,out] _slowdown: slowdown paramters, (pert_out and timescale should be updated)
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
    }
#endif


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

