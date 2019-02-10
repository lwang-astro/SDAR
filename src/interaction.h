#pragma once

#include <cmath>
#include "Float.h"
#include "particle.h"
#include "particle_group.h"
#include "force.h"
#include "perturber.h"

namespace AR{
    
    //! a sample interaction class with newtonian acceleration
    class Interaction{
    public:

        //! (Necessary) calculate acceleration from perturber 
        /*! The Force class acc_pert should be updated
          @param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
          @param[in] _particles: member particle array
          @param[in] _n_particle: number of member particles
          @param[in] _particle_cm: center-of-mass particle
          @param[in] _perturber: pertuber container
         */
        void calcAccPert(Force* _force, const Particle* _particles, const int _n_particle, const Particle& _particle_cm, const Perturber& _perturber) {
            for (int i=0; i<_n_particle; i++) {
                Float* acc_pert = _force[i].acc_pert;
                acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
            }            
        }

        //! (Necessary) calculate slowdown factor 
        /*!
          \return slowdown factor
         */
        Float calcSlowDownPert(const Float _etot, const Particle* _particle_cm, const Perturber& _perturber)  {
            return 1.0;
        }

#ifdef AR_TTL 
        //! calculate inner member acceleration, potential and time transformation function gradient and factor for kick
        /*!
          @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
          @param[out] _epot: total inner potential energy
          @param[in] _particles: member particle array
          @param[in] _n_particle: number of member particles
          \return the time transformation factor (gt_kick) for kick step
         */
        Float calcAccPotGTGradAndGTKick(Force* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
            _epot = Float(0.0);
            Float gt_kick = Float(0.0);
            for (int i=0; i<_n_particle; i++) {
                const Float massi = _particles[i].mass;
                const Float* posi = _particles[i].pos;
                Float* acci = _force[i].acc_in;
                Float* gtgradi = _force[i].gtgrad;

                acci[0] = acci[1] = acci[2] = Float(0.0);
                gtgradi[0] = gtgradi[1] = gtgradi[2] = Float(0.0);


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
                    
                    Float mimjor3 = massi*mor3;
                    gtgradi[0] += mimjor3 * dr[0];
                    gtgradi[1] += mimjor3 * dr[1];
                    gtgradi[2] += mimjor3 * dr[2];

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
#else
        //! calculate inner member acceleration, potential and time transformation factor for kick
        /*!
          @param[out] _force: force array to store the calculation results (in acc_in[3], notice acc may need to reset zero to avoid accummulating old values)
          @param[out] _epot: total inner potential energy
          @param[in] _particles: member particle array
          @param[in] _n_particle: number of member particles
          \return the time transformation factor (gt_kick) for kick step
         */
        Float calcAccPotAndGTKick(Force* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
            _epot = Float(0.0);
            for (int i=0; i<_n_particle; i++) {
                const Float massi = _particles[i].mass;
                const Float* posi = _particles[i].pos;
                Float* acci = _force[i].acc_in;
                acci[0] = acci[1] = acci[2] = Float(0.0);

                Float poti = Float(0.0);
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

                    poti -= massj*inv_r;
                }
                _epot += poti * massi;
            }
            _epot   *= 0.5;
            // log H case
            Float gt_kick = -1.0/_epot;
            
            return gt_kick;
        }

        //! calcualte the time transformation factor for drift
        /*!
          @param[in] _ekin_etot: ekin - etot
         */
        Float calcGTDrift(Float _ekin_etot) {
            return 1.0/_ekin_etot;
        }
#endif   
    };



}
