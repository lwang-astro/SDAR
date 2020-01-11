#pragma once

#include <cmath>
#include "Common/Float.h"
#include "AR/force.h"
#include "Hermite/neighbor.h"
#include "particle.h"

//! a sample interaction class with newtonian acceleration
class ARInteraction{
public:
    typedef H4::ParticleAR<Particle> ARPtcl;
    typedef H4::ParticleH4<Particle> H4Ptcl;

    Float eps_sq; // softening 
    Float gravitational_constant; ///> gravitational constant

    ARInteraction(): eps_sq(Float(-1.0)), gravitational_constant(Float(-1.0)) {}

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
             <<"G     : "<<gravitational_constant<<std::endl;
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
    inline Float calcInnerAccPotAndGTKickTwo(AR::Force& _f1, AR::Force& _f2, Float& _epot, const ARPtcl& _p1, const Particle& _p2) {
        // acceleration
        const Float mass1 = _p1.mass;
        const Float* pos1 = _p1.pos;

        const Float mass2 = _p2.mass;
        const Float* pos2 = _p2.pos;

        Float dr[3] = {pos2[0] -pos1[0],
                       pos2[1] -pos1[1],
                       pos2[2] -pos1[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
        ASSERT(r2>0.0);
        Float inv_r = 1.0/sqrt(r2);
        Float inv_r3 = inv_r*inv_r*inv_r;

        Float* acc1 = _f1.acc_in;
        Float* acc2 = _f2.acc_in;

        Float mor3_1 = gravitational_constant*mass2*inv_r3;
        acc1[0] = mor3_1 * dr[0];
        acc1[1] = mor3_1 * dr[1];
        acc1[2] = mor3_1 * dr[2];

        Float mor3_2 = gravitational_constant*mass1*inv_r3;
        acc2[0] = - mor3_2 * dr[0];
        acc2[1] = - mor3_2 * dr[1];
        acc2[2] = - mor3_2 * dr[2];

        Float m1m2 = gravitational_constant*mass1*mass2;

#ifdef AR_TTL 
        // trans formation function gradient
        Float m1m2or3 = m1m2*inv_r3;
        Float* gtgrad1 = _f1.gtgrad;
        Float* gtgrad2 = _f2.gtgrad;

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

    //! calculate inner member acceleration, potential and time transformation function gradient and factor for kick
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the time transformation factor (gt_kick) for kick step
    */
    inline Float calcInnerAccPotAndGTKick(AR::Force* _force, Float& _epot, const ARPtcl* _particles, const int _n_particle) {
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
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                Float inv_r = 1.0/sqrt(r2);
                Float inv_r3 = inv_r*inv_r*inv_r;
                Float mor3 = gravitational_constant*massj*inv_r3;
                acci[0] += mor3 * dr[0];
                acci[1] += mor3 * dr[1];
                acci[2] += mor3 * dr[2];

#ifdef AR_TTL                     
                Float mimjor3 = massi*mor3;
                gtgradi[0] += mimjor3 * dr[0];
                gtgradi[1] += mimjor3 * dr[1];
                gtgradi[2] += mimjor3 * dr[2];
#endif

                Float mor = gravitational_constant*massj*inv_r;
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
    Float calcAccPotAndGTKick(AR::Force* _force, Float& _epot, const ARPtcl* _particles, const int _n_particle, const H4Ptcl& _particle_cm, const H4::Neighbor<Particle>& _perturber, const Float _time) {
        static const Float inv3 = 1.0 / 3.0;

        Float gt_kick;
        if (_n_particle==2) gt_kick = calcInnerAccPotAndGTKickTwo(_force[0], _force[1], _epot, _particles[0], _particles[1]);
        else gt_kick = calcInnerAccPotAndGTKick(_force, _epot, _particles, _n_particle);

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
            Float mcm = _particle_cm.mass;
            if (_perturber.need_resolve_flag) {
                // calculate component perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    const auto& pi = _particles[i];
                    acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);

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
                        Float mor3 = gravitational_constant*m[j]/r3;

                        acc_pert[0] += mor3 * dr[0];
                        acc_pert[1] += mor3 * dr[1];
                        acc_pert[2] += mor3 * dr[2];
                    }

                    acc_pert_cm[0] += pi.mass *acc_pert[0];
                    acc_pert_cm[1] += pi.mass *acc_pert[1];
                    acc_pert_cm[2] += pi.mass *acc_pert[2];
#ifdef ARC_DEBUG
                    mcm += pi.mass;
#endif
                }
#ifdef ARC_DEBUG
                ASSERT(abs(mcm-_particles_cm.mass)<1e-10);
#endif
                // get cm perturbation
                acc_pert_cm[0] /= mcm;
                acc_pert_cm[1] /= mcm;
                acc_pert_cm[2] /= mcm;

                // remove cm. perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    acc_pert[0] -= acc_pert_cm[0]; 
                    acc_pert[1] -= acc_pert_cm[1];        
                    acc_pert[2] -= acc_pert_cm[2]; 
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
                    Float mor3 = gravitational_constant*m[j]/r3;

                    acc_pert_cm[0] += mor3 * dr[0];
                    acc_pert_cm[1] += mor3 * dr[1];
                    acc_pert_cm[2] += mor3 * dr[2];
                }

                // calculate component perturbation
                for (int i=0; i<_n_particle; i++) {
                    Float* acc_pert = _force[i].acc_pert;
                    const auto& pi = _particles[i];
                    acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);

                    Float xi[3];
                    xi[0] = pi.pos[0] + xcm[0];
                    xi[1] = pi.pos[1] + xcm[1];
                    xi[2] = pi.pos[2] + xcm[2];

                    // remove c.m. perturbation acc
                    acc_pert[0] = -acc_pert_cm[0]; 
                    acc_pert[1] = -acc_pert_cm[1];        
                    acc_pert[2] = -acc_pert_cm[2]; 

                    for (int j=0; j<n_pert; j++) {
                        Float dr[3] = {xp[j][0] - xi[0],
                                       xp[j][1] - xi[1],
                                       xp[j][2] - xi[2]};
                        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] + eps_sq;
                        Float r  = sqrt(r2);
                        Float r3 = r*r2;
                        Float mor3 = gravitational_constant*m[j]/r3;

                        acc_pert[0] += mor3 * dr[0];
                        acc_pert[1] += mor3 * dr[1];
                        acc_pert[2] += mor3 * dr[2];
                    }
                }

            }
        }
        return gt_kick;
        
    }


    //! calculate perturbation from c.m. acceleration
    Float calcPertFromForce(const Float* _force, const Float _mp, const Float _mpert) {
        Float force2 = _force[0]*_force[0]+_force[1]*_force[1]+_force[2]*_force[2];
#ifdef AR_SLOWDOWN_PERT_R4
        return force2/(_mp*_mpert);
#else
        Float force = sqrt(force2);
        return sqrt(force/(_mp*_mpert))*force;
#endif
    }

    //! calculate perturbation from binary tree
    Float calcPertFromBinary(const COMM::BinaryTree<ARPtcl>& _bin) {
        Float apo = _bin.semi*(1.0+_bin.ecc);
        Float apo2 = apo*apo;
#ifdef AR_SLOWDOWN_PERT_R4
        return (_bin.m1*_bin.m2)/(apo2*apo2);
#else
        return (_bin.m1*_bin.m2)/(apo2*apo);
#endif
    }

    //! calculate perturbation from distance to perturber and masses of particle and perturber 
    Float calcPertFromMR(const Float _r, const Float _mp, const Float _mpert) {
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
    void calcSlowDownInnerBinary(AR::SlowDown& _slowdown, const AR::SlowDown& _slowdown_cm, const COMM::BinaryTree<ARPtcl>& _bin_root, const ARPtcl* _particles, const int _n_particle) {
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
            Float mor3 = (mj+mcm)*r*r2/(mj*mcm);
            trf2_min =  std::min(trf2_min, mor3);

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

    //! (Necessary) calculate slowdown perturbation and timescale
    /*!
      @param[in,out] _slowdown: slowdown paramters, (pert_out and timescale should be updated)
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _bin_root: binary tree root
      @param[in] _perturber: pertuber container
    */
    void calcSlowDownPert(AR::SlowDown& _slowdown, const H4Ptcl& _particle_cm, const COMM::BinaryTree<ARPtcl>& _bin_root, const H4::Neighbor<Particle>& _perturber) {
        static const Float inv3 = 1.0 / 3.0;

        // slowdown inner perturbation: 
        _slowdown.pert_in = calcPertFromBinary(_bin_root);
        _slowdown.period  = _bin_root.period;

        const Float time = _slowdown.getRealTime();

        const int n_pert = _perturber.neighbor_address.getSize();

        if (n_pert>0) {

            auto* pert_adr = _perturber.neighbor_address.getDataAddress();

            Float xp[3], xcm[3];
            Float dt = time - _particle_cm.time;
            //ASSERT(dt>=0.0);
            xcm[0] = _particle_cm.pos[0] + dt*(_particle_cm.vel[0] + 0.5*dt*(_particle_cm.acc0[0] + inv3*dt*_particle_cm.acc1[0]));
            xcm[1] = _particle_cm.pos[1] + dt*(_particle_cm.vel[1] + 0.5*dt*(_particle_cm.acc0[1] + inv3*dt*_particle_cm.acc1[1]));
            xcm[2] = _particle_cm.pos[2] + dt*(_particle_cm.vel[2] + 0.5*dt*(_particle_cm.acc0[2] + inv3*dt*_particle_cm.acc1[2]));

            Float pert_pot = 0.0; 
            Float mcm = _particle_cm.mass;

#ifdef AR_SLOWDOWN_TIMESCALE
            // velocity dependent method 
            Float vp[3], vcm[3];

            vcm[0] = _particle_cm.vel[0] + dt*(_particle_cm.acc0[0] + 0.5*dt*_particle_cm.acc1[0]);
            vcm[1] = _particle_cm.vel[1] + dt*(_particle_cm.acc0[1] + 0.5*dt*_particle_cm.acc1[1]);
            vcm[2] = _particle_cm.vel[2] + dt*(_particle_cm.acc0[2] + 0.5*dt*_particle_cm.acc1[2]);

            Float trf2_min = NUMERIC_FLOAT_MAX;
            Float mvor[3] = {0.0,0.0,0.0};
            Float mtot=0.0;
#endif

            for (int j=0; j<n_pert; j++) {
                H4::NBAdr<Particle>::Single* pertj;
                if (pert_adr[j].type==H4::NBType::group) pertj = &(((H4::NBAdr<Particle>::Group*)pert_adr[j].adr)->cm);
                else pertj = (H4::NBAdr<Particle>::Single*)pert_adr[j].adr;
                Float dt = time - pertj->time;
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
                pert_pot += calcPertFromMR(r, mcm, mj);

#ifdef AR_SLOWDOWN_TIMESCALE
                // velocity dependent method 
                vp[0] = pertj->vel[0] + dt*(pertj->acc0[0] + 0.5*dt*pertj->acc1[0]);
                vp[1] = pertj->vel[1] + dt*(pertj->acc0[1] + 0.5*dt*pertj->acc1[1]);
                vp[2] = pertj->vel[2] + dt*(pertj->acc0[2] + 0.5*dt*pertj->acc1[2]);

                Float dv[3] = {vp[0] - vcm[0],
                               vp[1] - vcm[1],
                               vp[2] - vcm[2]};

                // m_tot / |\sum m_j /|r_j| * v_j|
                Float mor = mj/r;
                mvor[0] += mor*dv[0];
                mvor[1] += mor*dv[1];
                mvor[2] += mor*dv[2];
                mtot += mj;

                // force dependent method
                // min sqrt(r^3/(G m))
                Float mor3 = (mj+mcm)*r*r2/(gravitational_constant*mj*mcm);
                trf2_min =  std::min(trf2_min, mor3);

#endif

            }
            _slowdown.pert_out = pert_pot;
#ifdef AR_SLOWDOWN_TIMESCALE
            // velocity dependent method
            Float trv_ave = mtot/sqrt(mvor[0]*mvor[0] + mvor[1]*mvor[1] + mvor[2]*mvor[2]);
            // get min of velocity and force dependent values
            Float t_min = std::min(trv_ave, sqrt(trf2_min));

            _slowdown.timescale = 0.1*std::min(_slowdown.getTimescaleMax(), t_min);
#else
            _slowdown.timescale = _slowdown.getTimescaleMax();
#endif
        }
        else{
            _slowdown.pert_out = 0.0;
            _slowdown.timescale = _slowdown.getTimescaleMax();
        }
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

