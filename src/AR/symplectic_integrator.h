#pragma once

#include "Common/list.h"
#include "Common/particle_group.h"
#include "AR/symplectic_step.h"
#include "AR/force.h"
#include "AR/slow_down.h"
#include "AR/profile.h"
#include "AR/information.h"

//! Algorithmic regularization chain (ARC) namespace
/*!
  All major ARC classes and related acceleration functions (typedef) are defined
*/
namespace AR {
    //! Symplectic integrator manager
    /*! Tmethod is the class contain the interaction function, see sample of interaction.h:\n
     */
    template <class Tmethod>
    class SymplecticManager {
    public:
        Float time_error_max_real; ///> maximum time error (absolute), should be positive and larger than round-off error 
        Float energy_error_relative_max; ///> maximum energy error requirement 
        Float time_step_real_min;        ///> minimum real time step allown
        Float slowdown_pert_ratio_ref;   ///> slowdown perturbation /inner ratio reference factor
        Float slowdown_mass_ref;         ///> slowdown mass factor reference
        Float slowdown_timescale_max;       ///> slowdown maximum timescale to calculate maximum slowdown factor
        long long unsigned int step_count_max; ///> maximum step counts
        
        Tmethod interaction; ///> class contain interaction function
        SymplecticStep step;  ///> class to manager kick drift step

        //! constructor
        SymplecticManager(): time_error_max_real(Float(-1.0)), energy_error_relative_max(Float(-1.0)), time_step_real_min(Float(-1.0)), slowdown_pert_ratio_ref(Float(-1.0)), slowdown_mass_ref(Float(-1.0)), step_count_max(-1), interaction(), step() {}

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            //ASSERT(time_error_max_real>ROUND_OFF_ERROR_LIMIT);
            ASSERT(time_error_max_real>0.0);
            ASSERT(energy_error_relative_max>ROUND_OFF_ERROR_LIMIT);
            //ASSERT(time_step_real_min>ROUND_OFF_ERROR_LIMIT);
            ASSERT(time_step_real_min>0.0);
            ASSERT(slowdown_pert_ratio_ref>0.0);
            ASSERT(slowdown_mass_ref>0.0);
            ASSERT(slowdown_timescale_max>0);
            ASSERT(step_count_max>0);
            ASSERT(step.getOrder()>0);
            ASSERT(interaction.checkParams());
            return true;
        }

        //! write class data with BINARY format
        /*! @param[in] _fout: file IO for write
         */
        void writeBinary(FILE *_fout) {
            size_t size = sizeof(*this) - sizeof(interaction) - sizeof(step);
            fwrite(this, size, 1,_fout);
            interaction.writeBinary(_fout);
            step.writeBinary(_fout);
        }

        //! read class data with BINARY format and initial the array
        /*! @param[in] _fin: file IO for read
         */
        void readBinary(FILE *_fin) {
            size_t size = sizeof(*this) - sizeof(interaction) - sizeof(step);
            size_t rcount = fread(this, size, 1, _fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
            interaction.readBinary(_fin);
            step.readBinary(_fin);
        }

        //! print parameters
        void print(std::ostream & _fout) const{
            _fout<<"time_error_max_real       : "<<time_error_max_real<<std::endl
                 <<"energy_error_relative_max : "<<energy_error_relative_max<<std::endl 
                 <<"time_step_real_min        : "<<time_step_real_min<<std::endl
                 <<"slowdown_pert_ratio_ref   : "<<slowdown_pert_ratio_ref<<std::endl
                 <<"slowdown_mass_ref         : "<<slowdown_mass_ref<<std::endl
                 <<"slowdown_timescale_max    : "<<slowdown_timescale_max<<std::endl
                 <<"step_count_max            : "<<step_count_max<<std::endl;
            interaction.print(_fout);
            step.print(_fout);
        }
    };

    //! Symplectic integrator class for a group of particles
    /*! The basic steps to use the integrator \n
      1. Add particles (particles.addParticle/particles.linkParticleList)  \n
      2. Initial system (initial) \n
      3. Integration (integrateOneStep/integrateToTime) \n
      Requirement for Tparticle class, public memebers: pos[3], vel[3], mass\n
      Template dependence: Tparticle: particle type; Tpcm: particle cm type  Tpert: perturber class type, Tmethod: interaction class;
    */
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo>
    class SymplecticIntegrator {
    private:
        //! slowdown with pair particle index
        struct SlowDownPair{
            SlowDown slowdown;
            COMM::BinaryTree<Tparticle>* bin;
        
            SlowDownPair(): slowdown(), bin(NULL) {}
        };

        // intergrated variables
        Float time_;   ///< integrated time (not real physical time if slowdown is on)
        Float etot_; ///< integrated system energy

        // calculated varaiables
        Float ekin_;   ///< kinetic energy
        Float epot_;   ///< potential

        // cumulative slowdown (inner + outer) energy change
        Float de_sd_cum_;  

#ifdef AR_TTL
        // transformation factors
        Float gt_drift_inv_;  ///< integrated time transformation factor for drift: dt(drift) = ds/gt_drift_inv_
        Float gt_kick_;       ///< time transformation factor for kick 
#ifdef AR_TTL_SLOWDOWN_INNER
        Float ekin_sdi_;  ///< slowdown (inner) kinetic energy
        Float epot_sdi_;  ///< slowdown (inner) potential energy
        Float etot_sdi_;  ///< slowdown (inner) total energy
#endif
#endif

        // force array
        COMM::List<Force> force_; ///< acceleration array 

    public:
        SymplecticManager<Tmethod>* manager; ///< integration manager
        COMM::ParticleGroup<Tparticle,Tpcm> particles; ///< particle group manager
#ifdef AR_TTL_SLOWDOWN_INNER
        COMM::List<SlowDownPair> slowdown_inner; /// inner binary slowdown
#endif
        SlowDown slowdown; ///< slowdown of the system
        Tpert    perturber; ///< perturber class 
        Tinfo    info;   ///< information of the system
        Profile  profile;  ///< profile to measure the performance
        
        //! Constructor
        SymplecticIntegrator(): time_(0), etot_(0), ekin_(0), epot_(0), de_sd_cum_(0),
#ifdef AR_TTL
                                gt_drift_inv_(0), gt_kick_(0), 
#ifdef AR_TTL_SLOWDOWN_INNER
                                ekin_sdi_(0), epot_sdi_(0), etot_sdi_(0), 
#endif
#endif
                                force_(), manager(NULL), particles(), 
#ifdef AR_TTL_SLOWDOWN_INNER
                                slowdown_inner(), 
#endif
                                slowdown(), perturber(), info(), profile() {}

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            ASSERT(manager!=NULL);
            ASSERT(manager->checkParams());
            ASSERT(perturber.checkParams());
            ASSERT(info.checkParams());
            return true;
        }

        //! reserve memory for force
        /*! The size of force depends on the particle data size.Thus particles should be added first before call this function
        */
        void reserveIntegratorMem() {
            // force array always allocated local memory
            int nmax = particles.getSizeMax();
            ASSERT(nmax>0);
            force_.setMode(COMM::ListMode::local);
            force_.reserveMem(nmax);
#ifdef AR_TTL_SLOWDOWN_INNER
            slowdown_inner.setMode(COMM::ListMode::local);
            slowdown_inner.reserveMem(nmax/2);
#endif
        }

        //! Clear function
        /*! Free dynamical memory space allocated
         */
        void clear() {
            time_ = 0.0;
            de_sd_cum_ = 0.0;
            force_.clear();
            particles.clear();
#ifdef AR_TTL_SLOWDOWN_INNER
            slowdown_inner.clear();
#endif
            slowdown.clear();
            perturber.clear();
            info.clear();
            profile.clear();
        }

        //! destructor
        ~SymplecticIntegrator() {
            clear();
        }

        //! operator = 
        /*! Copy function will remove the local data and also copy the particle data or the link
         */
        SymplecticIntegrator& operator = (const SymplecticIntegrator& _sym) {
            clear();
            time_   = _sym.time_;
            etot_   = _sym.etot_;

            ekin_   = _sym.ekin_;
            epot_   = _sym.epot_;
            de_sd_cum_= _sym.de_sd_cum_;
#ifdef AR_TTL
            gt_drift_inv_ = _sym.gt_drift_inv_;
            gt_kick_ = _sym.gt_kick_;
#ifdef AR_TTL_SLOWDOWN_INNER
            ekin_sdi_= _sym.ekin_sdi_;
            epot_sdi_= _sym.epot_sdi_;
            etot_sdi_= _sym.etot_sdi_;
#endif
#endif
            force_  = _sym.force_;
            manager = _sym.manager;
#ifdef AR_TTL_SLOWDOWN_INNER
            slowdown_inner = _sym.slowdown_inner;
#endif
            slowdown = _sym.slowdown;
            particles = _sym.particles;
            info = _sym.binarytree;
            profile = _sym.profile;

            return *this;
        }

    private:
        //! Calculate kinetic energy
        inline void calcEKin(){
            ekin_ = Float(0.0);
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            for (int i=0; i<num; i++) {
                const Float *vi=pdat[i].getVel();
                ekin_ += 0.5 * pdat[i].mass * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
            }
#ifdef AR_TTL_SLOWDOWN_INNER
            calcEkinSlowDownInner(ekin_);
#endif
        }

#ifdef AR_TTL_SLOWDOWN_INNER
        //! correct force, potential energy and gt_kick based on slowdown for inner binaries
        /*! 
          @param[in,out] _gt_kick: the time transformation factor for kick step (input), be corrected with slowdown (output)
          @param[in] _epot: total potential energy wihtout slowdown
         */
        inline void correctAccPotGTKickSlowDownInner(Float& _gt_kick, const Float& _epot) {
            int n = slowdown_inner.getSize();
            Float gt_kick_inv_cor = 0.0;
            Float de = 0.0;
            for (int i=0; i<n; i++) {
                auto& sdi = slowdown_inner[i];
                ASSERT(sdi.bin!=NULL);
                int i1 = sdi.bin->getMemberIndex(0);
                int i2 = sdi.bin->getMemberIndex(1);
                Float kappa = sdi.slowdown.getSlowDownFactor();
                if (i1>=0) {
                    ASSERT(i2>=0);
                    ASSERT(i1!=i2);
                    ASSERT(i1<particles.getSize());
                    ASSERT(i2<particles.getSize());
                    
                    // calculate pair interaction
                    Force fi[2];
                    Float epoti;
                    Float gt_kick_i = manager->interaction.calcInnerAccPotAndGTKickTwo(fi[0], fi[1], epoti, particles[i1], particles[i2]);
                    Float kappa_inv = 1.0/kappa;
                    // scale binary pair force with slowdown 
                    force_[i1].acc_in[0] += fi[0].acc_in[0]*kappa_inv - fi[0].acc_in[0];
                    force_[i1].acc_in[1] += fi[0].acc_in[1]*kappa_inv - fi[0].acc_in[1];
                    force_[i1].acc_in[2] += fi[0].acc_in[2]*kappa_inv - fi[0].acc_in[2];
                    force_[i2].acc_in[0] += fi[1].acc_in[0]*kappa_inv - fi[1].acc_in[0];
                    force_[i2].acc_in[1] += fi[1].acc_in[1]*kappa_inv - fi[1].acc_in[1];
                    force_[i2].acc_in[2] += fi[1].acc_in[2]*kappa_inv - fi[1].acc_in[2];

                    de += epoti*kappa_inv - epoti;
                    // scale gtgrad with slowdown
#ifdef AR_TTL_GT_BINARY_INNER
                    force_[i1].gtgrad[0] = fi[0].gtgrad[0]*kappa_inv;
                    force_[i1].gtgrad[1] = fi[0].gtgrad[1]*kappa_inv;
                    force_[i1].gtgrad[2] = fi[0].gtgrad[2]*kappa_inv;
                    force_[i2].gtgrad[0] = fi[1].gtgrad[0]*kappa_inv;
                    force_[i2].gtgrad[1] = fi[1].gtgrad[1]*kappa_inv;
                    force_[i2].gtgrad[2] = fi[1].gtgrad[2]*kappa_inv;

                    gt_kick_inv_cor += kappa_inv/gt_kick_i;
#else
                    force_[i1].gtgrad[0] += fi[0].gtgrad[0]*kappa_inv - fi[0].gtgrad[0];
                    force_[i1].gtgrad[1] += fi[0].gtgrad[1]*kappa_inv - fi[0].gtgrad[1];
                    force_[i1].gtgrad[2] += fi[0].gtgrad[2]*kappa_inv - fi[0].gtgrad[2];
                    force_[i2].gtgrad[0] += fi[1].gtgrad[0]*kappa_inv - fi[1].gtgrad[0];
                    force_[i2].gtgrad[1] += fi[1].gtgrad[1]*kappa_inv - fi[1].gtgrad[1];
                    force_[i2].gtgrad[2] += fi[1].gtgrad[2]*kappa_inv - fi[1].gtgrad[2];

                    // gt kick
                    gt_kick_inv_cor += (kappa_inv - 1.0)/gt_kick_i;
#endif
                }
            }
            epot_sdi_ = epot_ + de;
#ifdef AR_TTL_GT_BINARY_INNER
            if (gt_kick_inv_cor!=0.0) _gt_kick = 1.0/gt_kick_inv_cor;
#else
            if (gt_kick_inv_cor!=0.0) {
                Float gt_kick_inv = 1.0/_gt_kick;
                gt_kick_inv += gt_kick_inv_cor;
                _gt_kick = 1.0/gt_kick_inv;
            }
#endif
        }

        //! correct postion drift due to inner binary slowdown
        /*! @param[in] _dt: time step
         */
        inline void correctPosSlowDownInner(const Float _dt) {
            int n = slowdown_inner.getSize();
            for (int i=0; i<n; i++) {
                auto& sdi = slowdown_inner[i];
                ASSERT(sdi.bin!=NULL);
                Float kappa = sdi.slowdown.getSlowDownFactor();
                Float kappa_inv_m_one = 1.0/kappa - 1.0;
                Float* velcm = sdi.bin->getVel();
                for (int k=0; k<2; k++) {
                    int j = sdi.bin->getMemberIndex(k);
                    ASSERT(j>=0&&j<particles.getSize());
                    Float* pos = particles[j].getPos();
                    Float* vel = particles[j].getVel();

                    // only scale velocity referring to binary c.m.
                    Float vrel[3] = { vel[0] - velcm[0], 
                                      vel[1] - velcm[1], 
                                      vel[2] - velcm[2]}; 
                    pos[0] += _dt * vrel[0] * kappa_inv_m_one;
                    pos[1] += _dt * vrel[1] * kappa_inv_m_one;
                    pos[2] += _dt * vrel[2] * kappa_inv_m_one;
                }
            }
        }

        //! calculate gt_drift_inv in the case of inner binary slowdown.
        /*!
          @param[in,out] _dgt: d gt_drift_inv_ for correction, modified the original value (from all members) with the slowdown correction
          @param[in] _dt: time step
         */
        inline void correctDGTInvSlowDownInner(Float& _dgt, const Float _dt) {
            int n = slowdown_inner.getSize();
            Float dg = Float(0.0);
            for (int i=0; i<n; i++) {
                auto& sdi = slowdown_inner[i];
                ASSERT(sdi.bin!=NULL);
                Float kappa = sdi.slowdown.getSlowDownFactor();
                Float kappa_inv = 1.0/kappa;
                for (int k=0; k<2; k++) {
                    int j = sdi.bin->getMemberIndex(k);
                    if (j>=0) {
                        ASSERT(j<particles.getSize());
                        Float* gtgrad=force_[j].gtgrad;
                        Float* vel = particles[j].getVel();
#ifdef AR_TTL_GT_BINARY_INNER
                        dg +=  (vel[0] * kappa_inv* gtgrad[0] +
                                vel[1] * kappa_inv* gtgrad[1] +
                                vel[2] * kappa_inv* gtgrad[2]);
#else
                        dg +=  (vel[0] * (kappa_inv-1)* gtgrad[0] +
                                vel[1] * (kappa_inv-1)* gtgrad[1] +
                                vel[2] * (kappa_inv-1)* gtgrad[2]);
#endif
                    }
                }
            }
#ifdef AR_TTL_GT_BINARY_INNER
            _dgt = _dt*dg;
#else
            _dgt += _dt*dg;
#endif
        }

        //! update c.m. for binaries with slowdown inner
        inline void updateCenterOfMassForBinaryWithSlowDownInner() {
            int n = slowdown_inner.getSize();
            for (int i=0; i<n; i++) {
                auto& sdi = slowdown_inner[i];
                int i1 = sdi.bin->getMemberIndex(0);
                int i2 = sdi.bin->getMemberIndex(1);
                ASSERT(i1>=0&&i1<particles.getSize());
                ASSERT(i2>=0&&i2<particles.getSize());

                Float    m1 = particles[i1].mass;
                Float* pos1 = particles[i1].getPos();
                Float* vel1 = particles[i1].getVel();
                Float    m2 = particles[i2].mass;
                Float* pos2 = particles[i2].getPos();
                Float* vel2 = particles[i2].getVel();
                Float   mcm = m1+m2;

                // first obtain the binary c.m. velocity
                Float mcminv = 1.0/mcm;

                sdi.bin->mass = mcm;
                Float* pos = sdi.bin->getPos();
                pos[0] = (m1*pos1[0] + m2*pos2[0])*mcminv;
                pos[1] = (m1*pos1[1] + m2*pos2[1])*mcminv;
                pos[2] = (m1*pos1[2] + m2*pos2[2])*mcminv;

                Float* vel = sdi.bin->getVel();
                vel[0] = (m1*vel1[0] + m2*vel2[0])*mcminv;
                vel[1] = (m1*vel1[1] + m2*vel2[1])*mcminv;
                vel[2] = (m1*vel1[2] + m2*vel2[2])*mcminv;
            }
        }

        //! calculate kinetic energy with slowdown factor
        /*! @param[in] _ekin: total kinetic energy without slowdown
         */
        inline void calcEkinSlowDownInner(const Float& _ekin) {
            int n = slowdown_inner.getSize();
            Float de = Float(0.0);
            for (int i=0; i<n; i++) {
                auto& sdi = slowdown_inner[i];
                ASSERT(sdi.bin!=NULL);
                Float kappa = sdi.slowdown.getSlowDownFactor();
                Float kappa_inv_m_one = 1.0/kappa - 1.0;
                Float* velcm = sdi.bin->getVel();
                for (int k=0; k<2; k++) {
                    int j = sdi.bin->getMemberIndex(k);
                    ASSERT(j>=0&&j<particles.getSize());
                    Float* vel = particles[j].getVel();

                    // only scale velocity referring to binary c.m.
                    Float vrel[3] = { vel[0] - velcm[0], 
                                      vel[1] - velcm[1], 
                                      vel[2] - velcm[2]}; 

                    de += kappa_inv_m_one * particles[j].mass * (vrel[0]*vrel[0] + vrel[1]*vrel[1] + vrel[2]*vrel[2]);
                }
            }
            ekin_sdi_ = _ekin + 0.5*de;
        }

        //! find inner binaries for slowdown treatment iteration function
        static int findSlowDownInnerBinaryIter(COMM::List<SlowDownPair>& _slowdown_inner, const int& _c1, const int& _c2, COMM::BinaryTree<Tparticle>& _bin) {
            // find leaf binary
            if (_bin.getMemberN()==2 && _bin.semi>0.0) {
                _slowdown_inner.increaseSizeNoInitialize(1);
                auto& sdi_new = _slowdown_inner.getLastMember();
                sdi_new.bin = &_bin;
                return _c1+_c2+1;
            }
            else return _c1+_c2;
        }

        //! find inner binaries for slowdown treatment
        /*! record binary tree address and set update time to _time
          @param[in] _time: next slowdown update time
         */
        void findSlowDownInner(const Float _time) {
            auto& bin_root = info.getBinaryTreeRoot();
            slowdown_inner.resizeNoInitialize(0);
            int ncount[2]={0,0};
            int nsd = bin_root.processTreeIter(slowdown_inner, ncount[0], ncount[1], findSlowDownInnerBinaryIter);
            ASSERT(nsd==slowdown_inner.getSize());

            for (int i=0; i<nsd; i++) {
#ifdef SLOWDOWN_MASSRATIO
                const Float mass_ratio = manager->slowdown_mass_ref/slowdown_inner[i].bin->mass;
#else 
                const Float mass_ratio = 1.0;
#endif
                slowdown_inner[i].slowdown.initialSlowDownReference(mass_ratio*manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);
                slowdown_inner[i].slowdown.setUpdateTime(time_);
            }
        }

        //! calculate slowdown inner and correct gt_drift_inv_ and slowdown energy
        /*! Check whether the slowdown inner should be updated, and also initialize the new gt_drift_inv_ with new slowdown factor if updated
          @param[in] _time: current intergration time
         */
        void calcSlowDownInnerAndCorrectGTInvEnergy(const Float _time) {
            int n = slowdown_inner.getSize();
            bool modified_flag=false;
            for (int i=0; i<n; i++) {
                auto& sdi = slowdown_inner[i];
                //if (_time>=sdi.slowdown.getUpdateTime()) {
                sdi.bin->calcCenterOfMass();
                manager->interaction.calcSlowDownInnerBinary(sdi.slowdown, slowdown, *sdi.bin, particles.getDataAddress(), particles.getSize());
                sdi.slowdown.increaseUpdateTimeOnePeriod();
                modified_flag=true;
                //}
            }    
            // correct the initial gt_drift_inv_ for next integration step is important to be consistent with TTL method.
            if (modified_flag) {
                // back up old ekin_sdi and epot_sdi
                Float ekin_sdi_bk = ekin_sdi_;
                Float epot_sdi_bk = epot_sdi_;
                // initialize the gt_drift_inv_ with new slowdown factor
#ifdef AR_TTL_GT_BINARY_INNER
                Float gt_kick_sdi = 0.0;
#else
                Float gt_kick_sdi = manager->interaction.calcAccPotAndGTKick(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, slowdown.getRealTime());
#endif
                correctAccPotGTKickSlowDownInner(gt_kick_sdi, epot_);
                //gt_drift_inv_ = 1.0/gt_kick_sdi;
                gt_drift_inv_ += 1.0/gt_kick_sdi - 1.0/gt_kick_;
                gt_kick_ = gt_kick_sdi;

                // correct etot_sdi_ with new slowdown
                calcEkinSlowDownInner(ekin_);
                Float de_sdi = (ekin_sdi_ - ekin_sdi_bk) + (epot_sdi_ - epot_sdi_bk);
                etot_sdi_ += de_sdi;

                // add slowdown change to the global slowdown energy
                de_sd_cum_ += de_sdi/slowdown.getSlowDownFactor();
            }
        }
#endif

        //! kick velocity
        /*! First time step will be calculated, the velocities are kicked
          @param[in] _dt: time size
        */
        inline void kickVel(const Float _dt) {
            const int num = particles.getSize();            
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            const Float kappa = slowdown.getSlowDownFactor();
            for (int i=0; i<num; i++) {
                // kick velocity
                Float* vel = pdat[i].getVel();
                Float* acc = force[i].acc_in;
                Float* pert= force[i].acc_pert;
                // half dv 
                vel[0] += _dt * (acc[0] + kappa*pert[0]);
                vel[1] += _dt * (acc[1] + kappa*pert[1]);
                vel[2] += _dt * (acc[2] + kappa*pert[2]);
            }
        }

        //! drift time and position
        /*! First (real) time is drifted, then positions are drifted
          @param[in] _dt: time step
        */
        inline void driftTimeAndPos(const Float _dt) {
            // drift time 
            time_ += _dt;

            // update real time
            slowdown.driftRealTime(_dt);

            // drift position
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            for (int i=0; i<num; i++) {
                Float* pos = pdat[i].getPos();
                Float* vel = pdat[i].getVel();
                pos[0] += _dt * vel[0];
                pos[1] += _dt * vel[1];
                pos[2] += _dt * vel[2];
            }

#ifdef AR_TTL_SLOWDOWN_INNER
            // correct postion drift due to inner binary slowdown
            correctPosSlowDownInner(_dt);
#endif
        }

        //! calc force, potential and time transformation factor for kick
        /*!
          \return gt_kick: time transformation factor for kick
         */
        inline Float calcAccPotAndGTKick() {
            Float gt_kick = manager->interaction.calcAccPotAndGTKick(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, slowdown.getRealTime());            

//#ifdef AR_DEBUG
//            // check c.m. force 
//            Force fcm;
//            Float mcm=0.0;
//            for (int i=0; i<particles.getSize(); i++) {
//                for (int k=0; k<3; k++) {
//                    fcm.acc_in[k] += particles[i].mass * force_[i].acc_in[k];
//                }
//                mcm += particles[i].mass;
//            }
//            for (int k=0; k<3; k++) {
//                fcm.acc_in[k] /= mcm;
//                ASSERT(abs(fcm.acc_in[k])<1e-10);
//            }
//#endif
#ifdef AR_TTL_SLOWDOWN_INNER
            // inner slowdown binary acceleration
            correctAccPotGTKickSlowDownInner(gt_kick, epot_);
//#ifdef AR_DEBUG
//            // check c.m. force 
//            fcm.acc_in[0] = fcm.acc_in[1] = fcm.acc_in[2] = 0.0;
//            for (int i=0; i<particles.getSize(); i++) {
//                for (int k=0; k<3; k++) {
//                    fcm.acc_in[k] += particles[i].mass * force_[i].acc_in[k];
//                }
//            }
//            for (int k=0; k<3; k++) {
//                fcm.acc_in[k] /= mcm;
//                ASSERT(abs(fcm.acc_in[k])<1e-10);
//            }
//#endif
#endif
            return gt_kick;
        }


#ifdef AR_TTL
        //! kick energy and time transformation function for drift
        /*!
          @param[in] _dt: time step
        */

        inline void kickEtotAndGTDrift(const Float _dt) {
            Float de = Float(0.0);
            Float dg = Float(0.0);
            const Float kappa = slowdown.getSlowDownFactor();
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].getVel();
                Float* pert= force[i].acc_pert;
                Float* gtgrad=force[i].gtgrad;
                de += mass * (vel[0] * kappa*pert[0] +
                              vel[1] * kappa*pert[1] +
                              vel[2] * kappa*pert[2]);
                dg +=  (vel[0] * gtgrad[0] +
                        vel[1] * gtgrad[1] +
                        vel[2] * gtgrad[2]);
            }
            etot_ += _dt * de;
            Float dgt = _dt * dg;
#ifdef AR_TTL_SLOWDOWN_INNER
            etot_sdi_ += _dt * de;
            correctDGTInvSlowDownInner(dgt, _dt);
#endif
            gt_drift_inv_ += dgt;
        }

#else
        //! kick energy
        /*!
          @param[in] _dt: time step
        */
        inline void kickEtot(const Float _dt) {
            Float de = 0.0;
            const Float kappa = slowdown.getSlowDownFactor();
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].getVel();
                Float* pert= force[i].acc_pert;
                de += mass * (vel[0] * kappa*pert[0] +
                              vel[1] * kappa*pert[1] +
                              vel[2] * kappa*pert[2]);
            }
            etot_ += _dt * de;
        }
#endif
            
    public:

        //! initialization for integration
        /*! initialize the system. Acceleration, energy and time transformation factors are updated. If the center-of-mass is not yet calculated, the system will be shifted to center-of-mass frame.
          @param[in] _time_real: real physical time to initialize
        */
        void initialIntegration(const Float _time_real) {
            ASSERT(checkParams());

            // particle number and data address
            const int n_particle = particles.getSize();
            Tparticle* particle_data = particles.getDataAddress();

            // resize force array
            force_.resizeNoInitialize(n_particle);
            Force* force_data = force_.getDataAddress();

            // Initial intgrt value t (avoid confusion of real time when slowdown is used)
            time_ = Float(0.0);

            // slowdown initial time
            slowdown.setRealTime(_time_real);

            // check particle number
            ASSERT(particles.getSize()>=2);

            // reset particle modification flag
            particles.setModifiedFalse();


            // check the center-of-mass initialization
            if(particles.isOriginFrame()) {
                particles.calcCenterOfMass();
                particles.shiftToCenterOfMassFrame();
            }

            // set slowdown reference
#ifdef SLOWDOWN_MASSRATIO
            const Float mass_ratio = manager->slowdown_mass_ref/particles.cm.mass;
#else 
            const Float mass_ratio = 1.0;
#endif
            slowdown.initialSlowDownReference(mass_ratio*manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);

            // slowdown for the system
            manager->interaction.calcSlowDownPert(slowdown, particles.cm, info.getBinaryTreeRoot(), perturber);
 
#ifdef AR_TTL
            gt_kick_ = manager->interaction.calcAccPotAndGTKick(force_data, epot_, particle_data, n_particle, particles.cm,  perturber, _time_real);

            // initially gt_drift 
            gt_drift_inv_ = 1.0/gt_kick_;

#else
            manager->interaction.calcAccPotAndGTKick(force_data, epot_, particle_data, n_particle, particles.cm,  perturber, _time_real);
#endif

            // calculate kinetic energy
            calcEKin();

            // initial total energy
            etot_ = ekin_ + epot_;

#ifdef AR_TTL_SLOWDOWN_INNER
            // initialize the energy with slowdown inner
            ekin_sdi_ = ekin_;
            epot_sdi_ = epot_;
            etot_sdi_ = etot_;

            if (particles.getSize()>2) {
                findSlowDownInner(time_);
                // update c.m. of binaries 
                updateCenterOfMassForBinaryWithSlowDownInner();
                // calculate slowdown inner and correct gt and energy
                calcSlowDownInnerAndCorrectGTInvEnergy(time_);
            }
#endif
            // reset de_sd_cum, energy change due to the initializaiton of slowdown (inner/outer) is excluded.
            de_sd_cum_ = 0.0;
        }

        //! integration for one step
        /*!
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().
        */
        void integrateOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // time transformation factor for drift
#ifdef AR_TTL
                Float gt_drift = 1.0/gt_drift_inv_;
#else 
                Float gt_drift = manager->interaction.calcGTDrift(ekin_-etot_); // pt = -etot
#endif

                // drift
                Float dt_drift = ds_drift*gt_drift;

                // drift time and postion
                driftTimeAndPos(dt_drift);
                _time_table[i] = slowdown.getRealTime();

                // step for kick
                Float ds_kick = manager->step.getDK(i)*_ds;

                //! calc force, potential and time transformation factor for kick
                Float gt_kick = calcAccPotAndGTKick();

                // time step for kick
                Float dt_kick = ds_kick* gt_kick;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

#ifdef AR_TTL   
                // back up gt_kick 
                gt_kick_ = gt_kick;

                // kick total energy and time transformation factor for drift
                kickEtotAndGTDrift(dt_kick);
#else
                // kick total energy 
                kickEtot(dt_kick);
#endif
                // kick half step for velocity
                kickVel(0.5*dt_kick);
                
#ifdef AR_TTL_SLOWDOWN_INNER
                // update c.m. of binaries 
                updateCenterOfMassForBinaryWithSlowDownInner();
#endif

                // calculate kinetic energy
                calcEKin();
            }
        }

        //! integration for two body one step
        /*! For two-body problem the calculation can be much symplified to improve performance. 
          Besides, the slow-down factor calculation is embedded in the Drift (for time) and Kick (for perturbation). 
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().         
        */
        void integrateTwoOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

            // symplectic step coefficent group number
            const int nloop = manager->step.getCDPairSize();
            
            const int n_particle = particles.getSize();
            ASSERT(n_particle==2);

            const Float kappa = slowdown.calcSlowDownFactor();

            Tparticle* particle_data = particles.getDataAddress();
            Float mass1 = particle_data[0].mass;
            Float* pos1 = particle_data[0].getPos();
            Float* vel1 = particle_data[0].getVel();

            Float mass2 = particle_data[1].mass;
            Float* pos2 = particle_data[1].getPos();
            Float* vel2 = particle_data[1].getVel();

            Force* force_data = force_.getDataAddress();
            Float* acc1 = force_data[0].acc_in;
            Float* pert1= force_data[0].acc_pert;

            Float* acc2 = force_data[1].acc_in;
            Float* pert2= force_data[1].acc_pert;
#ifdef AR_TTL
            Float* gtgrad1 = force_data[0].gtgrad;
            Float* gtgrad2 = force_data[1].gtgrad;
#endif
            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds = manager->step.getCK(i)*_ds;
                // time transformation factor for drift
#ifdef AR_TTL
                Float gt = 1.0/gt_drift_inv_;
#else 
                Float gt = manager->interaction.calcGTDrift(ekin_-etot_); // pt = -etot
#endif
                // drift
                Float dt = ds*gt;
                ASSERT(!isnan(dt));
                
                // drift time 
                time_ += dt;

                // update real time
                slowdown.driftRealTime(dt);
                _time_table[i] = slowdown.getRealTime();

                // drift position
                pos1[0] += dt * vel1[0];
                pos1[1] += dt * vel1[1];
                pos1[2] += dt * vel1[2];

                pos2[0] += dt * vel2[0];
                pos2[1] += dt * vel2[1];
                pos2[2] += dt * vel2[2];

                // step for kick
                ds = manager->step.getDK(i)*_ds;

                gt = manager->interaction.calcAccPotAndGTKick(force_data, epot_, particle_data, n_particle, particles.cm, perturber, _time_table[i]);

                ASSERT(!isnan(epot_));

                // time step for kick
                dt = 0.5*ds*gt;

                // kick half step for velocity
                Float dvel1[3], dvel2[3];
                Float kpert1[3], kpert2[3];
                kpert1[0] = kappa*pert1[0];
                kpert1[1] = kappa*pert1[1];
                kpert1[2] = kappa*pert1[2];

                dvel1[0] = dt * (acc1[0] + kpert1[0]);
                dvel1[1] = dt * (acc1[1] + kpert1[1]);
                dvel1[2] = dt * (acc1[2] + kpert1[2]);

                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];

                kpert2[0] = kappa*pert2[0];
                kpert2[1] = kappa*pert2[1];
                kpert2[2] = kappa*pert2[2];

                dvel2[0] = dt * (acc2[0] + kappa*pert2[0]);
                dvel2[1] = dt * (acc2[1] + kappa*pert2[1]);
                dvel2[2] = dt * (acc2[2] + kappa*pert2[2]);

                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];

                // kick total energy and time transformation factor for drift
                etot_ += 2.0*dt * (mass1* (vel1[0] * kpert1[0] + 
                                           vel1[1] * kpert1[1] + 
                                           vel1[2] * kpert1[2]) +
                                   mass2* (vel2[0] * kpert2[0] + 
                                           vel2[1] * kpert2[1] + 
                                           vel2[2] * kpert2[2]));

#ifdef AR_TTL   
                // back up gt_kick
                gt_kick_ = gt;

                // integrate gt_drift_inv
                gt_drift_inv_ +=  2.0*dt* (vel1[0] * gtgrad1[0] +
                                           vel1[1] * gtgrad1[1] +
                                           vel1[2] * gtgrad1[2] +
                                           vel2[0] * gtgrad2[0] +
                                           vel2[1] * gtgrad2[1] +
                                           vel2[2] * gtgrad2[2]);
#endif

                // kick half step for velocity
                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];
                
                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];
                
                // calculate kinetic energy
                ekin_ = 0.5 * (mass1 * (vel1[0]*vel1[0]+vel1[1]*vel1[1]+vel1[2]*vel1[2]) +
                               mass2 * (vel2[0]*vel2[0]+vel2[1]*vel2[1]+vel2[2]*vel2[2]));
            }

#ifdef AR_TTL_SLOWDOWN_INNER
            // make consistent slowdown inner energy 
            etot_sdi_ = etot_;
            ekin_sdi_ = ekin_;
            epot_sdi_ = epot_;
#endif
        }
        
        // Integrate the system to a given time
        /*!
          @param[in] _ds: the integration step size
          @param[in] _time_end_real: the expected finishing real time 
         */
        void integrateToTime(const Float _time_end_real) {
            ASSERT(checkParams());

            // real full time step
            const Float dt_real_full = _time_end_real - slowdown.getRealTime();

            // time error 
            const Float time_error_real = manager->time_error_max_real;

            // energy error limit
            const Float energy_error_rel_max = manager->energy_error_relative_max;
            // expect energy_error using half step if energy_error_rel_max reached
            const Float energy_error_rel_max_half_step = energy_error_rel_max * manager->step.calcErrorRatioFromStepModifyFactor(0.5);
            const Float dt_real_min = manager->time_step_real_min;

            // backup data size
            const int bk_data_size = getBackupDataSize();
            
            Float backup_data[bk_data_size]; // for backup chain data
#ifdef AR_DEBUG_DUMP
            Float backup_data_init[bk_data_size]; // for backup initial data
#endif
            bool backup_flag=true; // flag for backup or restore

            // time table
            const int cd_pair_size = manager->step.getCDPairSize();
            Float time_table[cd_pair_size]; // for storing sub-integrated time 

            // step control
            Float ds[2] = {info.ds,info.ds}; // step with a buffer
            Float ds_backup = info.ds;  //backup step size
            int ds_switch=0;   // 0 or 1
            int n_step_wait=-1; // number of waiting step to change ds
            int n_step_end=0;  // number of steps integrated to reach the time end for one during the time sychronization sub steps

            // time end flag
            bool time_end_flag=false; // indicate whether time reach the end

            // step count
            long long unsigned int step_count=0; // integration step 
            long long unsigned int step_count_tsyn=0; // time synchronization step
            
            // particle data
            const int n_particle = particles.getSize();

#ifdef AR_DEBUG
            Float ekin_check = ekin_;
            calcEKin();
            ASSERT(abs(ekin_check-ekin_)<1e-10);
            ekin_ = ekin_check;
#endif

#ifdef AR_DEBUG_DUMP
            // back up initial data
            backupIntData(backup_data_init);
#endif
      
            // integration loop
            while(true) {
                // backup data
                if(backup_flag) {
                    // update slowdown
                    if (!time_end_flag) {
                        // system slowdown
                        //if (time_>=slowdown.getUpdateTime()) {
                            // backup old slowdown factor
                            Float sd_old = slowdown.getSlowDownFactor();
                            manager->interaction.calcSlowDownPert(slowdown, particles.cm, info.getBinaryTreeRoot(), perturber);
                            // add energy change due to slowdown change
                            Float sd_new = slowdown.getSlowDownFactor();
#ifdef AR_TTL_SLOWDOWN_INNER
                            de_sd_cum_ += (1.0/sd_new-1.0/sd_old)*(ekin_sdi_ + epot_sdi_);
#else
                            de_sd_cum_ += (1.0/sd_new-1.0/sd_old)*(ekin_ + epot_);
#endif
                            slowdown.increaseUpdateTimeOnePeriod();
                            //}
#ifdef AR_TTL_SLOWDOWN_INNER
                        // inner binary slowdown
                        calcSlowDownInnerAndCorrectGTInvEnergy(time_);
#endif
                    }

                    int bk_return_size = backupIntData(backup_data);
                    ASSERT(bk_return_size == bk_data_size);
                    (void)bk_return_size;

                }
                else { //restore data
                    int bk_return_size = restoreIntData(backup_data);
                    ASSERT(bk_return_size == bk_data_size);
                    (void)bk_return_size;
#ifdef AR_TTL_SLOWDOWN_INNER
                    // update c.m. of binaries 
                    // binary c.m. is not backup, thus recalculate to get correct c.m. velocity for position drift correction due to slowdown inner (the first drift in integrateonestep assume c.m. vel is up to date)
                    updateCenterOfMassForBinaryWithSlowDownInner();
#endif
                }

                // get real time 
                Float dt_real = slowdown.getRealTime();

                // integrate one step
                ASSERT(!isinf(ds[ds_switch]));
                if(n_particle==2) integrateTwoOneStep(ds[ds_switch], time_table);
                else integrateOneStep(ds[ds_switch], time_table);

                // real step size
                Float time_real = slowdown.getRealTime();
                dt_real =  time_real - dt_real;
                ASSERT(dt_real>0.0);
                
                step_count++;

                // energy check
#ifdef AR_TTL_SLOWDOWN_INNER
                Float energy_error_bk = getEnergyErrorSlowDownInnerFromBackup(backup_data);
                Float etot_bk = getEtotSlowDownInnerFromBackup(backup_data);
                Float energy_error = getEnergyErrorSlowDownInner();
#else
                Float energy_error_bk = getEnergyErrorFromBackup(backup_data);
                Float etot_bk = getEtotFromBackup(backup_data);
                Float energy_error = getEnergyError();
#endif
                Float energy_error_diff = energy_error - energy_error_bk;
                Float energy_error_rel_abs = abs(energy_error_diff/etot_bk);

                // time error
                Float time_diff_real_rel = (_time_end_real - time_real)/dt_real_full;

#ifdef AR_WARN
                // warning for large number of steps
                if(step_count>=manager->step_count_max) {
                    if(step_count%manager->step_count_max==0) {
                        std::cerr<<"Warning: step count is signficiant large "<<step_count<<std::endl;
                        std::cerr<<"Time(int): "<<time_
                                 <<" Time(real): "<<time_real
                                 <<" Time_end(real): "<<_time_end_real
                                 <<" Time_diff_rel(real): "<<time_diff_real_rel
                                 <<" ds(used): "<<ds[ds_switch]
                                 <<" ds(next): "<<ds[1-ds_switch]
                                 <<" Energy_error_rel_abs: "<<energy_error_rel_abs
                                 <<std::endl;
                        printColumnTitle(std::cerr);
                        std::cerr<<std::endl;
                        printColumn(std::cerr);
                        std::cerr<<std::endl;
                        DATADUMP("dump_large_step");
                    }
                }
#endif
          
                // When time sychronization steps too large, abort
                if(step_count_tsyn>manager->step_count_max) {
                    std::cerr<<"Error! step count after time synchronization is too large "<<step_count_tsyn<<std::endl;
                    std::cerr<<" Time(int): "<<time_
                             <<" Time(real): "<<time_real
                             <<" Time_end(real): "<<_time_end_real
                             <<" Time_error_rel(real): "<<time_diff_real_rel
                             <<" ds(used): "<<ds[ds_switch]
                             <<" ds(next): "<<ds[1-ds_switch]
                             <<" ds_backup: "<<ds_backup
                             <<" Energy_error_rel_abs: "<<energy_error_rel_abs
                             <<" Step_count: "<<step_count
                             <<std::endl;
                    printColumnTitle(std::cerr);
                    std::cerr<<std::endl;
                    printColumn(std::cerr);
                    std::cerr<<std::endl;
#ifdef AR_DEBUG_DUMP
                    restoreIntData(backup_data_init);
                    DATADUMP("dump_large_step");
#endif
                    abort();
                }


#ifdef AR_DEEP_DEBUG
                std::cerr<<" Time(int): "<<time_
                         <<" Time(real): "<<time_real
                         <<" Time_end(real): "<<_time_end_real
                         <<" Time_error_rel(real): "<<time_diff_real_rel
                         <<" ds(used): "<<ds[ds_switch]
                         <<" ds(next): "<<ds[1-ds_switch]
                         <<" ds_backup: "<<ds_backup
                         <<" Energy_error_rel_abs: "<<energy_error_rel_abs
                         <<" Step_count: "<<step_count
                         <<std::endl;
                std::cerr<<"Timetable: ";
                for (int i=0; i<cd_pair_size; i++) std::cerr<<" "<<time_table[manager->step.getSortCumSumCKIndex(i)];
                std::cerr<<std::endl;
#endif

                // modify step if energy error is large
                if(energy_error_rel_abs>energy_error_rel_max) {
                    if(info.fix_step_option!=FixStepOption::always) {

                        // energy error zero case, continue to avoid problem
                        if (energy_error_rel_abs==0.0) continue;

                        // estimate the modification factor based on the symplectic order
                        Float energy_error_ratio = energy_error_rel_max/energy_error_rel_abs;
                        // limit modify_factor to 0.125
                        Float modify_factor = std::max(manager->step.calcStepModifyFactorFromErrorRatio(energy_error_ratio), Float(0.125));
                        ASSERT(modify_factor>0.0);

                        // for initial steps
                        if(step_count<3) {
                            n_step_wait=-1;
                            ds[ds_switch] *= modify_factor;
                            ds[1-ds_switch] = ds[ds_switch];
                            ASSERT(!isinf(ds[ds_switch]));
                            backup_flag = false;
#ifdef AR_DEEP_DEBUG
                            std::cerr<<"Detected energy error too large, energy_error/max ="<<1.0/energy_error_ratio<<" energy_error_rel_abs ="<<energy_error_rel_abs<<" modify_factor ="<<modify_factor<<std::endl;
#endif
                            continue;
                        }
                        // for big energy error
                        else if (info.fix_step_option==FixStepOption::none && energy_error_ratio<0.1) {
                            if(backup_flag) ds_backup = ds[ds_switch];
                            ds[ds_switch] *= modify_factor;
                            ds[1-ds_switch] = ds[ds_switch];
                            ASSERT(!isinf(ds[ds_switch]));
                            backup_flag = false;
                            n_step_wait = 2*to_int(1.0/modify_factor);
#ifdef AR_DEEP_DEBUG
                            std::cerr<<"Detected energy error too large, energy_error/max ="<<1.0/energy_error_ratio<<" energy_error_rel_abs ="<<energy_error_rel_abs<<" modify_factor ="<<modify_factor<<" n_step_wait ="<<n_step_wait<<std::endl;
#endif
                            continue;
                        }
                    }
                }
// too much output
//#ifdef AR_WARN
//                if(energy_error_rel_abs>100.0*energy_error_rel_max) {
//                    std::cerr<<"Warning: symplectic integrator error > 100*criterion:"<<energy_error_rel_abs<<std::endl;
//                }
//#endif

                // abort when too small step found
                if(!time_end_flag&&dt_real<dt_real_min) {
                    std::cerr<<"Error! symplectic integrated time step ("<<dt_real<<") < minimum step ("<<dt_real_min<<")!\n";
                    std::cerr<<" Time(int): "<<time_
                             <<" Time(real): "<<time_real
                             <<" Time_end(real): "<<_time_end_real
                             <<" Time_error_rel(real): "<<(_time_end_real - time_real)/dt_real_full
                             <<" ds(used): "<<ds[ds_switch]
                             <<" ds(next): "<<ds[1-ds_switch]
                             <<" ds_backup: "<<ds_backup
                             <<" Energy_error_rel_abs: "<<energy_error_rel_abs
                             <<" Step_count: "<<step_count
                             <<std::endl;

#ifdef AR_DEBUG_DUMP
                    restoreIntData(backup_data_init);
                    DATADUMP("dump_negative_time");
#endif
                    abort();
                }
          
                // check integration time
                if(time_real < _time_end_real - time_error_real){
                    // step increase depend on n_step_wait or energy_error
                    if(info.fix_step_option==FixStepOption::none && !time_end_flag) {
                        // waiting step count reach
                        if(n_step_wait==0) {
#ifdef AR_DEEP_DEBUG
                            std::cerr<<"Recover to backup step ds_current="<<ds[ds_switch]<<" ds_next="<<ds[1-ds_switch]<<" ds_backup="<<ds_backup<<std::endl;
#endif
                            ds[1-ds_switch] = ds_backup;
                            ASSERT(!isinf(ds[1-ds_switch]));
                        }
                        // increase step size if energy error is small
                        else if(energy_error_rel_abs<energy_error_rel_max_half_step&&energy_error_rel_abs>0.0) {
                            Float energy_error_ratio = energy_error_rel_max/energy_error_rel_abs;
                            Float modify_factor = manager->step.calcStepModifyFactorFromErrorRatio(energy_error_ratio);
                            ASSERT(modify_factor>0.0);
                            ds[1-ds_switch] *= modify_factor;
                            ASSERT(!isinf(ds[1-ds_switch]));
#ifdef AR_DEEP_DEBUG
                            std::cerr<<"Energy error is small enought for increase step, energy_error_rel_abs="<<energy_error_rel_abs
                                     <<" energy_error_rel_max="<<energy_error_rel_max<<" step_modify_factor="<<modify_factor<<" new ds="<<ds[1-ds_switch]<<std::endl;
#endif
                        }
                        n_step_wait--;
                    }

                    // time sychronization on case, when step size too small to reach time end, increase step size
                    if(time_end_flag && ds[ds_switch]==ds[1-ds_switch]) {
                        step_count_tsyn++;

                        Float dt_real_end = _time_end_real - time_real;
                        if(n_step_end>1 && dt_real<0.3*dt_real_end) {
                            // dt_real should be >0.0
                            // ASSERT(dt_real>0.0);
                            ds[1-ds_switch] = ds[ds_switch] * dt_real_end/abs(dt_real);
                            ASSERT(!isinf(ds[1-ds_switch]));
#ifdef AR_DEEP_DEBUG
                            std::cerr<<"Time step dt(real) "<<dt_real<<" <0.3*(time_end-time)(real) "<<dt_real_end<<" enlarge step factor: "<<dt_real_end/dt_real<<" new ds: "<<ds[1-ds_switch]<<std::endl;
#endif
                        }
                        else n_step_end++;
                    }

                    // when used once, update to the new step
                    ds[ds_switch] = ds[1-ds_switch]; 
                    ASSERT(!isinf(ds[ds_switch]));
                    ds_switch = 1-ds_switch;

                    backup_flag = true;
                }
                else if(time_real > _time_end_real + time_error_real) {
                    time_end_flag = true;
                    backup_flag = false;

                    step_count_tsyn++;
                    n_step_end=0;

                    // check timetable
                    int i=-1,k=0; // i indicate the increasing time index, k is the corresponding index in time_table
                    for(i=0; i<cd_pair_size; i++) {
                        k = manager->step.getSortCumSumCKIndex(i);
                        if(_time_end_real<time_table[k]) break;
                    }
                    if (i==0) { // first step case
                        ASSERT(time_table[k]>0.0);
                        ds[ds_switch] *= manager->step.getSortCumSumCK(i)*_time_end_real/time_table[k];
                        ds[1-ds_switch] = ds[ds_switch];
                        ASSERT(!isinf(ds[ds_switch]));
#ifdef AR_DEEP_DEBUG
                        std::cerr<<"Time_end_real reach, time[k](real)= "<<time_table[k]<<" time(real)= "<<time_real<<" time_end/time[k](real)="<<_time_end_real/time_table[k]<<" CumSum_CK="<<manager->step.getSortCumSumCK(i)<<" ds(next) = "<<ds[ds_switch]<<" ds(next_next) = "<<ds[1-ds_switch]<<"\n";
#endif
                    }
                    else { // not first step case, get the interval time 
                        // previous integrated sub time in time table
                        Float time_real_prev = time_table[manager->step.getSortCumSumCKIndex(i-1)];
                        Float dt_real_k = time_table[k] - time_real_prev;
                        Float ds_tmp = ds[ds_switch];
                        // get cumsum CK factor for two steps near the time_end_real
                        Float cck_prev = manager->step.getSortCumSumCK(i-1);
                        Float cck = manager->step.getSortCumSumCK(i);
                        // in case the time is between two sub step, first scale the next step with the previous step CumSum CK cck(i-1)
                        ASSERT(!isinf(cck_prev));
                        ds[ds_switch] *= cck_prev;  
                        ASSERT(!isinf(ds[ds_switch]));
                        // then next next step, scale with the CumSum CK between two step: cck(i) - cck(i-1) 
                        ASSERT(dt_real_k>0.0);
                        ds[1-ds_switch] = ds_tmp*(cck-cck_prev)*std::min(Float(1.0),(_time_end_real-time_real_prev+time_error_real)/dt_real_k); 
                        ASSERT(!isinf(ds[1-ds_switch]));

#ifdef AR_DEEP_DEBUG
                        std::cerr<<"Time_end_real reach, time_prev(real)= "<<time_real_prev<<" time[k](real)= "<<time_table[k]<<" time(real)= "<<time_real<<" (time_end-time_prev)/dt(real)="<<(_time_end_real-time_real_prev)/dt_real<<" CumSum_CK="<<cck<<" CumSum_CK(prev)="<<cck_prev<<" ds(next) = "<<ds[ds_switch]<<" ds(next_next) = "<<ds[1-ds_switch]<<" \n";
#endif
                    }
                }
                else {
#ifdef AR_DEEP_DEBUG
                    std::cerr<<"Finish, time_diff_real_rel = "<<time_diff_real_rel<<" energy_error_rel_abs = "<<energy_error_rel_abs<<std::endl;
#endif
#ifdef AR_WARN
                    if (energy_error_rel_abs>energy_error_rel_max) {
                        std::cerr<<"AR large energy error, error_rel: "<<energy_error_rel_abs
                                 <<" N_members: "<<n_particle<<" Fix step option: "<<static_cast<typename std::underlying_type<FixStepOption>::type>(info.fix_step_option);
                        Float etot_init = getEtotFromBackup(backup_data_init);
                        std::cerr<<" etot_init: "<<etot_init
                                 <<" etot_end: "<<etot_
                                 <<" etot_diff: "<<etot_- etot_init;
#ifdef AR_TTL_SLOWDOWN_INNER
                        Float etot_sdi_init = getEtotSlowDownInnerFromBackup(backup_data_init);
                        std::cerr<<" etot_sdi_init: "<<etot_sdi_init
                                 <<" etot_sdi_end: "<<etot_sdi_
                                 <<" etot_sdi_diff: "<<etot_sdi_- etot_sdi_init;
#else
                        Float error_init = getEnergyErrorFromBackup(backup_data_init);
                        std::cerr<<" ekin+epot(init): "<<error_init+etot_init
                                 <<" ekin+epot(end): "<<etot_ + energy_error
                                 <<" ekin+epot(diff): "<<etot_ + energy_error - (error_init+etot_init);

#endif                      
                        std::cerr<<std::endl;
//#ifdef AR_DEBUG_DUMP
//                        restoreIntData(backup_data_init);
//                        DATADUMP("dump_large_error");
//#endif
//                        abort
                    }
#endif
                    break;
                }
            }

#ifdef AR_TTL_SLOWDOWN_INNER
            // find new inner slowdown binaries
            if (n_particle >2) {
                int nold = slowdown_inner.getSize();
                findSlowDownInner(time_);
                int nnew = slowdown_inner.getSize();
                // in case slowdown is disabled in the next step, gt_drift_inv_ should be re-initialized
                if (nold>0&&nnew==0) {
                    gt_kick_ = manager->interaction.calcAccPotAndGTKick(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, slowdown.getRealTime());
                    gt_drift_inv_ = 1.0/gt_kick_;
                }
            }
#endif
            // cumulative step count 
            profile.step_count = step_count;
            profile.step_count_tsyn = step_count_tsyn;
            profile.step_count_sum += step_count;
            profile.step_count_tsyn_sum += step_count_tsyn;
        }

        //! correct CM drift
        /*! calculate c.m. and correct the member data to the c.m. frame.
          This is used after the perturbation, in case the c.m. drift when members are in c.m. frame
         */
        void correctCenterOfMassDrift() {
            ASSERT(!particles.isOriginFrame());
            Float mcm=0.0, pos_cm[3]={0.0,0.0,0.0}, vel_cm[3]={0.0,0.0,0.0};
            auto* particle_data= particles.getDataAddress();
            for (int i=0; i<particles.getSize(); i++) {
                const Float *ri = particle_data[i].pos;
                const Float *vi = particle_data[i].getVel();
                const Float mi  = particle_data[i].mass;

                pos_cm[0] += ri[0] * mi;
                pos_cm[1] += ri[1] * mi;
                pos_cm[2] += ri[2] * mi;

                vel_cm[0] += vi[0] * mi;
                vel_cm[1] += vi[1] * mi;
                vel_cm[2] += vi[2] * mi;
                mcm += mi;
            }
            pos_cm[0] /= mcm; 
            pos_cm[1] /= mcm; 
            pos_cm[2] /= mcm; 
            vel_cm[0] /= mcm; 
            vel_cm[1] /= mcm; 
            vel_cm[2] /= mcm;

            for (int i=0; i<particles.getSize(); i++) {
                Float *ri = particle_data[i].pos;
                Float *vi = particle_data[i].getVel();

                ri[0] -= pos_cm[0]; 
                ri[1] -= pos_cm[1]; 
                ri[2] -= pos_cm[2]; 
                vi[0] -= vel_cm[0]; 
                vi[1] -= vel_cm[1]; 
                vi[2] -= vel_cm[2]; 
            }
        }

        //! write back particles with slowdown velocity
        /*! write back particles with slowdown velocity to original address
          @param[in] _particle_cm: center of mass particle to calculate the original frame, different from the particles.cm
         */
        template <class Tptcl>
        void writeBackSlowDownParticles(const Tptcl& _particle_cm) {
            ASSERT(particles.getMode()==COMM::ListMode::copy);
            ASSERT(!particles.isOriginFrame());
            auto* particle_adr = particles.getOriginAddressArray();
            auto* particle_data= particles.getDataAddress();
            const Float kappa_inv = 1.0/slowdown.getSlowDownFactor();
            for (int i=0; i<particles.getSize(); i++) {
                ASSERT(particle_adr[i]->mass == particle_data[i].mass);
                particle_adr[i]->mass = particle_data[i].mass;

                particle_adr[i]->pos[0] = particle_data[i].pos[0] + _particle_cm.pos[0];
                particle_adr[i]->pos[1] = particle_data[i].pos[1] + _particle_cm.pos[1];
                particle_adr[i]->pos[2] = particle_data[i].pos[2] + _particle_cm.pos[2];

                particle_adr[i]->vel[0] = particle_data[i].vel[0]*kappa_inv + _particle_cm.vel[0];
                particle_adr[i]->vel[1] = particle_data[i].vel[1]*kappa_inv + _particle_cm.vel[1];
                particle_adr[i]->vel[2] = particle_data[i].vel[2]*kappa_inv + _particle_cm.vel[2];
            }
        }

        //! write back particles to original address
        /*! If particles are in center-off-mass frame, write back the particle in original frame but not modify local copies to avoid roundoff error
         */
        template <class Tptcl>
        void writeBackParticlesOriginFrame() {
            ASSERT(particles.getMode()==COMM::ListMode::copy);
            auto* particle_adr = particles.getOriginAddressArray();
            auto* particle_data= particles.getDataAddress();
            if (particles.isOriginFrame()) {
                for (int i=0; i<particles.getSize(); i++) {
                    *(Tptcl*)particle_adr[i] = particle_data[i];
                }
            }
            else {
                for (int i=0; i<particles.getSize(); i++) {
                    Tptcl pc = particle_data[i];

                    pc.pos[0] = particle_data[i].pos[0] + particles.cm.pos[0];
                    pc.pos[1] = particle_data[i].pos[1] + particles.cm.pos[1];
                    pc.pos[2] = particle_data[i].pos[2] + particles.cm.pos[2];

                    pc.vel[0] = particle_data[i].vel[0] + particles.cm.vel[0];
                    pc.vel[1] = particle_data[i].vel[1] + particles.cm.vel[1];
                    pc.vel[2] = particle_data[i].vel[2] + particles.cm.vel[2];
                    
                    *(Tptcl*)particle_adr[i] = pc;
                }
            }
        }

        //! Get current physical time
        /*! \return current physical time
         */
        Float getTime() const {
            return slowdown.getRealTime();
        }

        //! Get current kinetic energy
        /*! \return current kinetic energy
         */
        Float getEkin() const {
            return ekin_;
        }

        //! Get current potential energy
        /*! \return current potetnial energy (negative value for bounded systems)
         */
        Float getEpot() const {
            return epot_;
        }

        //! Get current total integrated energy 
        /*! \return total integrated energy 
         */
        Float getEtot() const {
            return etot_;
        }

        //! get energy error 
        /*! \return energy error
         */
        Float getEnergyError() const {
            return ekin_ + epot_ - etot_;
        }

        //! get energy error from backup data
        Float getEnergyErrorFromBackup(Float* _bk) const {
            return -_bk[1] + _bk[2] + _bk[3];
        }

        //! get energy error from backup data
        Float getEtotFromBackup(Float* _bk) const {
            return _bk[1];
        }

        //! reset cumulative energy change due to slowdown change
        void resetDESlowDownCum() {
            de_sd_cum_ = 0.0;
        }

        //! get cumulative energy change due to slowdown change
        Float getDESlowDownCum() const {
            return de_sd_cum_;
        }


#ifdef AR_TTL_SLOWDOWN_INNER
        //! Get current kinetic energy with inner slowdown
        /*! \return current kinetic energy
         */
        Float getEkinSlowDownInner() const {
            return ekin_sdi_;
        }

        //! Get current potential energy with inner slowdown
        /*! \return current potetnial energy (negative value for bounded systems)
         */
        Float getEpotSlowDownInner() const {
            return epot_sdi_;
        }
        //! Get current total integrated energy with inner slowdown
        /*! \return total integrated energy 
         */
        Float getEtotSlowDownInner() const {
            return etot_sdi_;
        }
        //! get energy error with inner slowdown
        /*! \return energy error
         */
        Float getEnergyErrorSlowDownInner() const {
            return ekin_sdi_ + epot_sdi_ - etot_sdi_;
        }

        //! get energy error with inner slowdown from backup data 
        Float getEnergyErrorSlowDownInnerFromBackup(Float* _bk) const {
            return -_bk[7] + _bk[8] + _bk[9];
        }

        //! get energy error with inner slowdown from backup data
        Float getEtotSlowDownInnerFromBackup(Float* _bk) const {
            return _bk[7];
        }

#endif

        //! get backup data size
        int getBackupDataSize() const {
#ifdef AR_TTL
#ifdef AR_TTL_SLOWDOWN_INNER
            int bk_size = 10; 
#else
            int bk_size = 7; 
#endif
#else
            int bk_size = 5;
#endif
            bk_size += particles.getBackupDataSize();
            bk_size += slowdown.getBackupDataSize();
            return bk_size;
        }


        //! Backup integration data 
        /*! Backup #time_, #etot_, #ekin_, $epot_, #gt_drift_, $gt_kick_, #particles, $slowdown to one Float data array
          \return backup array size
        */
        int backupIntData(Float* _bk) {
            _bk[0] = time_;
            _bk[1] = etot_;
            _bk[2] = ekin_;
            _bk[3] = epot_;
            _bk[4] = de_sd_cum_;
#ifdef AR_TTL
            _bk[5] = gt_drift_inv_;
            _bk[6] = gt_kick_;
#ifdef AR_TTL_SLOWDOWN_INNER
            _bk[7] = etot_sdi_;
            _bk[8] = ekin_sdi_;
            _bk[9] = epot_sdi_;
            int bk_size = 10; 
#else
            int bk_size = 7; 
#endif
#else
            int bk_size = 5;
#endif
            bk_size += particles.backupParticlePosVel(&_bk[bk_size]); 
            bk_size +=  slowdown.backupSlowDownFactorAndTimeReal(&_bk[bk_size]); // time_real, slowdownfactor
            return bk_size;
        }

        //! Restore integration data
        /*! restore #time_, #etot_, #ekin_, $epot_, #gt_drift_, $gt_kick_, #particles, $slowdown from one Float data array
          \return backup array size
        */
        int restoreIntData(Float* _bk) {
            time_    = _bk[0];
            etot_    = _bk[1];
            ekin_    = _bk[2];
            epot_    = _bk[3];
            de_sd_cum_= _bk[4];
#ifdef AR_TTL
            gt_drift_inv_  = _bk[5];
            gt_kick_       = _bk[6];
#ifdef AR_TTL_SLOWDOWN_INNER
            etot_sdi_ = _bk[7];
            ekin_sdi_ = _bk[8];
            epot_sdi_ = _bk[9];
            int bk_size = 10; 
#else
            int bk_size = 7; 
#endif
#else
            int bk_size = 5; 
#endif
            bk_size += particles.restoreParticlePosVel(&_bk[bk_size]);
            bk_size +=  slowdown.restoreSlowDownFactorAndTimeReal(&_bk[bk_size]);
            return bk_size;
        }

#ifdef AR_TTL
        //! Get integrated inverse time transformation factor
        /*! In TTF case, it is calculated by integrating \f$ \frac{dg}{dt} = \sum_k \frac{\partial g}{\partial \vec{r_k}} \bullet \vec{v_k} \f$.
          Notice last step is the sub-step in one symplectic loop
          \return time transformation factor for drift
        */
        Float getGTDriftInv() const {
            return gt_drift_inv_;
        }
#endif

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"Time_int"
                 <<std::setw(_width)<<"dE"
                 <<std::setw(_width)<<"Etot"
                 <<std::setw(_width)<<"Ekin"
                 <<std::setw(_width)<<"Epot";
#ifdef AR_TTL
            _fout<<std::setw(_width)<<"GtD_inv";
#else
            _fout<<std::setw(_width)<<"H";
#endif
            _fout<<std::setw(_width)<<"dE_SD_cum";            

#ifdef AR_TTL_SLOWDOWN_INNER
            _fout<<std::setw(_width)<<"dE_SD_in" 
                 <<std::setw(_width)<<"Etot_SD_in" 
                 <<std::setw(_width)<<"Ekin_SD_in" 
                 <<std::setw(_width)<<"Epot_SD_in";
            _fout<<std::setw(_width)<<"N_SD_in";
            for (int i=0; i<slowdown_inner.getSize(); i++) {
                _fout<<std::setw(_width)<<"I1"
                     <<std::setw(_width)<<"I2";
                slowdown_inner[i].slowdown.printColumnTitle(_fout, _width);
            }
#endif
            slowdown.printColumnTitle(_fout, _width);
            perturber.printColumnTitle(_fout, _width);
            info.printColumnTitle(_fout, _width);
            profile.printColumnTitle(_fout, _width);
            particles.printColumnTitle(_fout, _width);
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<time_
                 <<std::setw(_width)<<getEnergyError()
                 <<std::setw(_width)<<etot_
                 <<std::setw(_width)<<ekin_
                 <<std::setw(_width)<<epot_;
#ifdef AR_TTL
            _fout<<std::setw(_width)<<gt_drift_inv_;
#else
            _fout<<std::setw(_width)<<manager->interaction.calcH(ekin_-etot_, epot_);
#endif
            _fout<<std::setw(_width)<<de_sd_cum_;
#ifdef AR_TTL_SLOWDOWN_INNER
            _fout<<std::setw(_width)<<getEnergyErrorSlowDownInner()
                 <<std::setw(_width)<<etot_sdi_ 
                 <<std::setw(_width)<<ekin_sdi_ 
                 <<std::setw(_width)<<epot_sdi_;
            _fout<<std::setw(_width)<<slowdown_inner.getSize();
            for (int i=0; i<slowdown_inner.getSize(); i++) {
                _fout<<std::setw(_width)<<slowdown_inner[i].bin->getMemberIndex(0)
                     <<std::setw(_width)<<slowdown_inner[i].bin->getMemberIndex(1);
                slowdown_inner[i].slowdown.printColumn(_fout, _width);
            }
#endif
            slowdown.printColumn(_fout, _width);
            perturber.printColumn(_fout, _width);
            info.printColumn(_fout, _width);
            profile.printColumn(_fout, _width);
            particles.printColumn(_fout, _width);
        }

        //! write class data with BINARY format
        /*! @param[in] _fout: file IO for write
         */
        void writeBinary(FILE *_fout) {
            fwrite(&time_, sizeof(Float), 1, _fout);
            fwrite(&etot_, sizeof(Float), 1, _fout);
            fwrite(&ekin_, sizeof(Float), 1, _fout);
            fwrite(&epot_, sizeof(Float), 1, _fout);
#ifdef AR_TTL
            fwrite(&gt_drift_inv_, sizeof(Float), 1, _fout);
#endif
            int size = force_.getSize();
            fwrite(&size, sizeof(int), 1, _fout);
            for (int i=0; i<size; i++) force_[i].writeBinary(_fout);
            
            particles.writeBinary(_fout);
            slowdown.writeBinary(_fout);
            perturber.writeBinary(_fout);
            info.writeBinary(_fout);
            profile.writeBinary(_fout);
        }

        //! read class data with BINARY format and initial the array
        /*! @param[in] _fin: file IO for read
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(&time_, sizeof(Float), 1, _fin);
            rcount += fread(&etot_, sizeof(Float), 1, _fin);
            rcount += fread(&ekin_, sizeof(Float), 1, _fin);
            rcount += fread(&epot_, sizeof(Float), 1, _fin);
            if (rcount<4) {
                std::cerr<<"Error: Data reading fails! requiring data number is 4, only obtain "<<rcount<<".\n";
                abort();
            }
#ifdef AR_TTL
            rcount = fread(&gt_drift_inv_, sizeof(Float), 1, _fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
#endif
            int size;
            rcount = fread(&size, sizeof(int),1, _fin);
            if(rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
            if(size<0) {
                std::cerr<<"Error: array size <0 "<<size<<"<=0!\n";
                abort();
            }
            if (size>0) {
                force_.setMode(COMM::ListMode::local);
                force_.reserveMem(size);
                force_.resizeNoInitialize(size);
                for (int i=0; i<size; i++) force_[i].readBinary(_fin);
            }
            
            particles.setMode(COMM::ListMode::local);
            particles.readBinary(_fin);
            slowdown.readBinary(_fin);
            perturber.readBinary(_fin);
            info.readBinary(_fin);
            profile.readBinary(_fin);
        }

    };
}
