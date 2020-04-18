#pragma once

#include <functional>
#include "Common/list.h"
#include "Common/particle_group.h"
#include "AR/symplectic_step.h"
#include "AR/force.h"
#include "AR/slow_down.h"
#include "AR/profile.h"
#include "AR/information.h"
#include "AR/interrupt.h"

//! Algorithmic regularization chain (ARC) namespace
/*!
  All major ARC classes and related acceleration functions (typedef) are defined
*/
namespace AR {
    //! Time Transformed Symplectic integrator manager
    /*! Tmethod is the class contain the interaction function, see sample of interaction.h:\n
     */
    template <class Tmethod>
    class TimeTransformedSymplecticManager {
    public:
        Float time_error_max; ///> maximum time error (absolute), should be positive and larger than round-off error 
        Float energy_error_relative_max; ///> maximum energy error requirement 
        Float time_step_min;        ///> minimum real time step allown
        Float slowdown_pert_ratio_ref;   ///> slowdown perturbation /inner ratio reference factor
#ifdef AR_SLOWDOWN_MASSRATIO
        Float slowdown_mass_ref;         ///> slowdown mass factor reference
#endif
        Float slowdown_timescale_max;       ///> slowdown maximum timescale to calculate maximum slowdown factor
        long long int step_count_max; ///> maximum step counts
        int interrupt_detection_option;    ///> 1: detect interruption; 0: no detection
        
        Tmethod interaction; ///> class contain interaction function
        SymplecticStep step;  ///> class to manager kick drift step

        //! constructor
        TimeTransformedSymplecticManager(): time_error_max(Float(-1.0)), energy_error_relative_max(Float(-1.0)), time_step_min(Float(-1.0)), slowdown_pert_ratio_ref(Float(-1.0)), 
#ifdef AR_SLOWDOWN_MASSRATIO
                                            slowdown_mass_ref(Float(-1.0)), 
#endif
                                            slowdown_timescale_max(0.0),
                                            step_count_max(-1), interrupt_detection_option(0), interaction(), step() {}

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            //ASSERT(time_error_max>ROUND_OFF_ERROR_LIMIT);
            ASSERT(time_error_max>0.0);
            ASSERT(energy_error_relative_max>ROUND_OFF_ERROR_LIMIT);
            //ASSERT(time_step_min>ROUND_OFF_ERROR_LIMIT);
            ASSERT(time_step_min>0.0);
            ASSERT(slowdown_pert_ratio_ref>0.0);
#ifdef AR_SLOWDOWN_MASSRATIO
            ASSERT(slowdown_mass_ref>0.0);
#endif
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
            _fout<<"time_error_max            : "<<time_error_max<<std::endl
                 <<"energy_error_relative_max : "<<energy_error_relative_max<<std::endl 
                 <<"time_step_min             : "<<time_step_min<<std::endl
                 <<"slowdown_pert_ratio_ref   : "<<slowdown_pert_ratio_ref<<std::endl
#ifdef AR_SLOWDOWN_MASSRATIO
                 <<"slowdown_mass_ref         : "<<slowdown_mass_ref<<std::endl
#endif
                 <<"slowdown_timescale_max    : "<<slowdown_timescale_max<<std::endl
                 <<"step_count_max            : "<<step_count_max<<std::endl;
            interaction.print(_fout);
            step.print(_fout);
        }
    };

    //! Time Transformed Symplectic integrator class for a group of particles
    /*! The basic steps to use the integrator \n
      1. Add particles (particles.addParticle/particles.linkParticleList)  \n
      2. Initial system (initial) \n
      3. Integration (integrateOneStep/integrateToTime) \n
      Requirement for Tparticle class, public memebers: pos[3], vel[3], mass\n
      Template dependence: Tparticle: particle type; Tpcm: particle cm type  Tpert: perturber class type, Tmethod: interaction class;
    */
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tinfo>
    class TimeTransformedSymplecticIntegrator {
    private:
        // intergrated variables
        Float time_;   ///< integrated time (not the physical time, it is always 0 initially)
        Float etot_ref_; ///< integrated system energy

        // calculated varaiables
        Float ekin_;   ///< kinetic energy
        Float epot_;   ///< potential

        // cumulative slowdown (inner + outer) energy change
        Float de_sd_change_cum_;  // slowdown energy change
        Float dH_sd_change_cum_;  // slowdown Hamiltonian change 

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
        Float ekin_sd_;  ///< slowdown (inner) kinetic energy
        Float epot_sd_;  ///< slowdown (inner) potential energy
        Float etot_sd_ref_;  ///< slowdown (inner) total energy
#endif

#ifdef AR_TTL
        // transformation factors
        Float gt_drift_inv_;  ///< integrated inverse time transformation factor for drift: dt(drift) = ds/gt_drift_inv_
        Float gt_kick_inv_;   ///< inverse time transformation factor for kick: dt(kick) = ds/gt_kick_inv_
#endif

        // force array
        COMM::List<Force> force_; ///< acceleration array 

    public:
        TimeTransformedSymplecticManager<Tmethod>* manager; ///< integration manager
        COMM::ParticleGroup<Tparticle,Tpcm> particles; ///< particle group manager
#ifdef AR_SLOWDOWN_ARRAY
        COMM::List<AR::BinaryTree<Tparticle>*> binary_slowdown; /// binary slowdown, first is root, then others are inner bianries
#endif
        Tpert    perturber; ///< perturber class 
        Tinfo    info;   ///< information of the system
        Profile  profile;  ///< profile to measure the performance
        
        //! Constructor
        TimeTransformedSymplecticIntegrator(): time_(0), etot_ref_(0), ekin_(0), epot_(0), de_sd_change_cum_(0), dH_sd_change_cum_(0),
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
                                ekin_sd_(0), epot_sd_(0), etot_sd_ref_(0), 
#endif
#ifdef AR_TTL
                                gt_drift_inv_(0), gt_kick_inv_(0), 
#endif
                                force_(), manager(NULL), particles(), 
#ifdef AR_SLOWDOWN_ARRAY
                                binary_slowdown(), 
#endif
                                perturber(), info(), profile() {}

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
#ifdef AR_SLOWDOWN_ARRAY
            binary_slowdown.setMode(COMM::ListMode::local);
            binary_slowdown.reserveMem(nmax/2+1);
#endif
        }

        //! Clear function
        /*! Free dynamical memory space allocated
         */
        void clear() {
            time_ = 0.0;
            etot_ref_ =0.0;
            ekin_ = 0.0;
            epot_ = 0.0;
            de_sd_change_cum_ = 0.0;
            dH_sd_change_cum_ = 0.0;
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            ekin_sd_ = 0.0;
            epot_sd_ = 0.0;
            etot_sd_ref_ = 0.0;
#endif
#ifdef AR_TTL
            gt_drift_inv_ = 0.0;
            gt_kick_inv_ = 0.0;
#endif
            force_.clear();
            particles.clear();
#ifdef AR_SLOWDOWN_ARRAY
            binary_slowdown.clear();
#endif
            perturber.clear();
            info.clear();
            profile.clear();
        }

        //! destructor
        ~TimeTransformedSymplecticIntegrator() {
            clear();
        }

        //! operator = 
        /*! Copy function will remove the local data and also copy the particle data or the link
         */
        TimeTransformedSymplecticIntegrator& operator = (const TimeTransformedSymplecticIntegrator& _sym) {
            clear();
            time_   = _sym.time_;
            etot_ref_   = _sym.etot_ref_;

            ekin_   = _sym.ekin_;
            epot_   = _sym.epot_;
            de_sd_change_cum_= _sym.de_sd_change_cum_;
            dH_sd_change_cum_= _sym.dH_sd_change_cum_;
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            ekin_sd_= _sym.ekin_sd_;
            epot_sd_= _sym.epot_sd_;
            etot_sd_ref_= _sym.etot_sd_ref_;
#endif
#ifdef AR_TTL
            gt_drift_inv_ = _sym.gt_drift_inv_;
            gt_kick_inv_ = _sym.gt_kick_inv_;
#endif
            force_  = _sym.force_;
            manager = _sym.manager;
#ifdef AR_SLOWDOWN_ARRAY
            binary_slowdown = _sym.binary_slowdown;
#endif
            particles = _sym.particles;
            info = _sym.binarytree;
            profile = _sym.profile;

            return *this;
        }

    private:

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
        //! iteration function to calculate perturbation and timescale information from binary tree j to binary i
        /*!
          @param[out] _pert_out: perturbation from particle j
          @param[out] _t_min_sq: timescale limit from particle j
          @param[in] _bini: binary i
          @param[in] _binj: binary tree j
         */
        void calcSlowDownPertInnerBinaryIter(Float& _pert_out, Float& _t_min_sq, AR::BinaryTree<Tparticle>& _bini, AR::BinaryTree<Tparticle>& _binj) {
            ASSERT(&_bini != &_binj);
            ASSERT(_bini.getMemberIndex(0)!=_binj.getMemberIndex(0));

            for (int k=0; k<2; k++) {
                if (_binj.isMemberTree(k)) {
                    auto* bink =  _binj.getMemberAsTree(k);
                    if (bink!=&_bini) {
                        if (bink->semi>0.0) manager->interaction.calcSlowDownPertOne(_pert_out, _t_min_sq, _bini, *bink);
                        else calcSlowDownPertInnerBinaryIter(_pert_out, _t_min_sq, _bini, *bink);
                    
                    }
                }
                else {
                    auto* pk =  _binj.getMember(k);
                    manager->interaction.calcSlowDownPertOne(_pert_out, _t_min_sq, _bini, *pk);
                }
            }
        }

        //! calculate slowdown factor for inner binary based on other particles and slowdown of system c.m.
        /*!
          @param[in] _bin: binary tree for calculating slowdown
        */
        void calcSlowDownInnerBinary(BinaryTree<Tparticle>& _bin) {
            _bin.slowdown.pert_in = manager->interaction.calcPertFromBinary(_bin);
            _bin.slowdown.period = _bin.period;

            Float pert_out = 0.0;
            Float t_min_sq = NUMERIC_FLOAT_MAX;
            auto& bin_root = info.getBinaryTreeRoot();

            calcSlowDownPertInnerBinaryIter(pert_out, t_min_sq, _bin, bin_root);

            _bin.slowdown.pert_out = pert_out + bin_root.slowdown.pert_out;

#ifdef AR_SLOWDOWN_TIMESCALE
            // velocity dependent method
            //Float trv_ave = sdtdat.mtot/sqrt(sdtdat.mvor[0]*sdtdat.mvor[0] + sdtdat.mvor[1]*sdtdat.mvor[1] + sdtdat.mvor[2]*sdtdat.mvor[2]);
            // get min of velocity and force dependent values
            //Float t_min = std::min(trv_ave, sqrt(sdtdat.trf2_min));
            _bin.slowdown.timescale = std::min(_bin.slowdown.getTimescaleMax(), sqrt(t_min_sq));
#else
            _bin.slowdown.timescale = _bin.slowdown.getTimescaleMax();
#endif

            _bin.slowdown.calcSlowDownFactor();
        }
#endif // AR_SLOWDOWN_TREE || AR_SLOWDOWN_ARRAY

#ifdef AR_SLOWDOWN_TREE
        //! update slowdown factor based on perturbation and record slowdown energy change
        /*! Update slowdown inner and global.
            @param [in] _update_energy_flag: Record cumulative slowdown energy change if true;
         */
        void updateSlowDownAndCorrectEnergy(const bool _update_energy_flag) {
            auto& bin_root = info.getBinaryTreeRoot();
            auto& sd_root = bin_root.slowdown;
            Float sd_backup = sd_root.getSlowDownFactor();
            //if (time_>=sd_root.getUpdateTime()) {
            sd_root.pert_in = manager->interaction.calcPertFromBinary(bin_root);
            sd_root.period = bin_root.period;

            Float t_min_sq= NUMERIC_FLOAT_MAX;
            sd_root.pert_out = 0.0;
            manager->interaction.calcSlowDownPert(sd_root.pert_out, t_min_sq, getTime(), particles.cm, perturber);
            
            sd_root.timescale = std::min(sd_root.getTimescaleMax(), sqrt(t_min_sq));
            sd_root.calcSlowDownFactor();

            //sd_root.increaseUpdateTimeOnePeriod();
            //}

            // inner binary slowdown
            bool inner_sd_change_flag=false;
            int n_bin = info.binarytree.getSize();
            for (int i=0; i<n_bin-1; i++) {
                auto& bini = info.binarytree[i];
                // only set slowdown if semi > 0 
                if (bini.semi>0.0) {
                    //if (time_>=bini.slowdown.getUpdateTime()) {
                    sdi->calcCenterOfMass();
                    calcSlowDownInnerBinary(bini);
                    //sdi->slowdown.increaseUpdateTimeOnePeriod();
                    inner_sd_change_flag=true;
                    //}
                }
            }

            if (_update_energy_flag) {
                Float ekin_sd_bk = ekin_sd_;
                Float epot_sd_bk = epot_sd_;
                Float etot_sd_ref_bk = etot_sd_ref_;
#ifdef AR_TTL
                Float gt_drift_inv_bk = gt_drift_inv_;
#endif
                if(inner_sd_change_flag) {
#ifdef AR_TTL
                    Float gt_kick_inv_new = calcAccPotAndGTKickInv();
                    gt_drift_inv_ += gt_kick_inv_new - gt_kick_inv_;
                    gt_kick_inv_ = gt_kick_inv_new;
#else                    
                    calcAccPotAndGTKickInv();
#endif
                    calcEkin();
                }
                else {
                    Float kappa_inv = 1.0/sd_root.getSlowDownFactor();
                    Float gt_kick_inv_new = gt_kick_inv_*sd_backup*kappa_inv;
#ifdef AR_TTL
                    gt_drift_inv_ += gt_kick_inv_sdi - gt_kick_inv_;
                    gt_kick_inv_ = gt_kick_inv_sdi;
#endif
                    ekin_sd_ = ekin_*kappa_inv;
                    epot_sd_ = epot_*kappa_inv;
                }
            }
            Float de_sd = (ekin_sd_ - ekin_sd_bk) + (epot_sd_ - epot_sd_bk);
            etot_sd_ref_ += de_sd;
#ifdef AR_TTL
            Float dH_sd = (ekin_sd_ + epot_sd_ - etot_sd_ref_)*gt_drift_inv_ - (ekin_sd_bk + epot_sd_bk - etot_sd_ref_bk)*gt_drift_inv_bk;
#else
            Float dH_sd = manager->interaction.calcH(ekin_sd_ - etot_sd_ref_, epot_sd_) 
                - manager->interaction.calcH(ekin_sd_bk - etot_sd_ref_bk, epot_sd_bk);
#endif

            // add slowdown change to the global slowdown energy
            de_sd_change_cum_ += de_sd;
            dH_sd_change_cum_ += dH_sd;
        }

        //! Calculate twice (slowdown) kinetic energy iteration function with binary tree
        /*! cumulative ekin_ and ekin_sd_. Notice these two values should be initialized to zero and reduce by two after iteration.
          @param[in] _vel_sd_up: upper cm sd vel
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary tree for kick etot and calc dgt_drift
        */
        void calcTwoEKinIter(const Float* _vel_sd_up, const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin){
            Float inv_nest_sd = _inv_nest_sd_up/_bin->slowdown.getSlowDownFactor();
            Float* vel_cm = _bin.getVel();
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    Float* vel = bink->getVel();
                    Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                       (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                       (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]}; 
                    calcEKinIter(vel_sd, inv_nest_sd, *bink);
                }
                else {
                    auto* pk = _bin.getMember(k);
                    Float* vel = pk->getVel();
                    Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                       (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                       (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]};
                    
                    ekin_ += pk->mass * (vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
                    ekin_sd_ += * pk->mass * (vel_sd[0]*vel_sd[0]+vel_sd[1]*vel_sd[1]+vel_sd[2]*vel_sd[2]);
                }
            }
        }

        //! Calculate (slowdown) kinetic energy
        void calcEkin() {
            ekin_ = ekin_sd_ = 0.0;
            auto& bin_root=info.getBinaryTreeRoot();
            Float vel_cm[3] = {0.0,0.0,0.0};
            Float sd_factor=1.0;
            
            calcTwoEKinIter(vel_cm, sd_factor, bin_root);
            ekin_ *= 0.5;
            ekin_sd_ *= 0.5;
        }


        //! update slowdown velocity iteration function with binary tree
        void updateBinaryVelIter(AR::BinaryTree<Tparticle>& _bin) {
            for (int k=0; k<2; k++) {
                _bin.vel[3]={0.0,0.0,0.0};
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    updateBinaryVelIter(*bink);
                    _bin.vel[0] += bink->mass*bink->vel[0];
                    _bin.vel[1] += bink->mass*bink->vel[1];
                    _bin.vel[2] += bink->mass*bink->vel[2];
                }
                else {
                    auto* pk = _bin.getMemberIndex(k);
                    _bin.vel[0] += pk->mass*pk->vel[0];
                    _bin.vel[1] += pk->mass*pk->vel[1];
                    _bin.vel[2] += pk->mass*pk->vel[2];
                }
            }
        }

        //! kick velocity
        /*! First time step will be calculated, the velocities are kicked
          @param[in] _dt: time size
        */
        void kickVel(const Float& _dt) {
            const int num = particles.getSize();            
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0; i<num; i++) {
                // kick velocity
                Float* vel = pdat[i].getVel();
                Float* acc = force[i].acc_in;
                Float* pert= force[i].acc_pert;
                // half dv 
                vel[0] += _dt * (acc[0] + pert[0]);
                vel[1] += _dt * (acc[1] + pert[1]);
                vel[2] += _dt * (acc[2] + pert[2]);
            }
            // update binary c.m. velocity interation
            updateBinaryVelIter(info.getBinaryTreeRoot());
        }

        //! drift position with slowdown tree
        /*!
          @param[in] _dt: drift time
          @param[in] _vel_sd_up: upper cm sd vel
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary to drift pos
        */
        void driftPosTreeIter(const Float& _dt, const Float* _vel_sd_up, const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin) { 
            // current nested sd factor
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();

            Float* pos_cm = _bin.getPos();
            Float* vel_cm = _bin.getVel();

            auto driftPos=[&](Float* pos, Float* vel, Float* vel_sd) {
                //scale velocity referring to binary c.m.
                vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                             (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                             (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]}; 
                pos[0] += _dt * vel_sd[0];
                pos[1] += _dt * vel_sd[1];
                pos[2] += _dt * vel_sd[2];
            };

            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* pj = _bin.getMemberAsTree(k);
                    Float* pos = pj->getPos();
                    Float* vel = pj->getVel();
                    Float vel_sd[3];
                    driftPos(pos, vel, vel_sd);
                    droftPosTreeIter(_dt, vel_sd, inv_nest_sd, *pj);
                }
                else {
                    auto* pj = _bin.getMember(k);
                    Float* pos = pj->getPos();
                    Float* vel = pj->getVel();
                    Float vel_sd[3];
                    driftPos(pos, vel, vel_sd);
                }
            }
        }

        //! drift time and position with slowdown tree
        /*! First (real) time is drifted, then positions are drifted
          @param[in] _dt: time step
        */
        void driftTimeAndPos(const Float& _dt) {
            // drift time 
            time_ += _dt;

            // the particle cm velocity is zero (assume in rest-frame)
            ASSERT(!particles.isOriginFrame());
            auto& bin_root=info.getBinaryTreeRoot();
            Float vel_cm[3] = {0.0,0.0,0.0};
            Float sd_factor=1.0;

            driftPosTreeIter(_dt, vel_cm, sd_factor, bin_root);
        }

        //! calc force, potential and inverse time transformation factor for one pair of particles
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _i: particle i index
          @param[in] _j: particle j index
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvTwo(const Float& _inv_nest_sd, const int _i, const int _j) {
            ASSERT(_i>=0&&_i<particles.getSize());
            ASSERT(_j>=0&&_j<particles.getSize());
            
            // calculate pair interaction
            Force fij[2];
            Float epotij;
            Float gt_kick_inv = manager->interaction.calcInnerAccPotAndGTKickInvTwo(fij[0], fij[1], epotij, particles[i], particles[j]);

            // scale binary pair force with slowdown 
            force_[_i].acc_in[0] += fij[0].acc_in[0]*_inv_nest_sd;
            force_[_i].acc_in[1] += fij[0].acc_in[1]*_inv_nest_sd;
            force_[_i].acc_in[2] += fij[0].acc_in[2]*_inv_nest_sd;
            force_[_j].acc_in[0] += fij[1].acc_in[0]*_inv_nest_sd;
            force_[_j].acc_in[1] += fij[1].acc_in[1]*_inv_nest_sd;
            force_[_j].acc_in[2] += fij[1].acc_in[2]*_inv_nest_sd;

            epot_    += epotij;
            epot_sd_ += epotij*_inv_nest_sd;

#ifdef AR_TTL
            // scale gtgrad with slowdown
            force_[_i].gtgrad[0] += fij[0].gtgrad[0]*_inv_nest_sd;
            force_[_i].gtgrad[1] += fij[0].gtgrad[1]*_inv_nest_sd;
            force_[_i].gtgrad[2] += fij[0].gtgrad[2]*_inv_nest_sd;
            force_[_j].gtgrad[0] += fij[1].gtgrad[0]*_inv_nest_sd;
            force_[_j].gtgrad[1] += fij[1].gtgrad[1]*_inv_nest_sd;
            force_[_j].gtgrad[2] += fij[1].gtgrad[2]*_inv_nest_sd;
#endif 
            
            return gt_kick_inv * _inv_nest_sd;
        }

        //! calc force, potential and inverse time transformation factor for one particle by walking binary tree
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _i: particle index
          @param[in] _bin: binary tree for walking
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvOneTreeIter(const Float& _inv_nest_sd, const int _i, AR::BinaryTree<Tparticle>& _bin) {
            Float gt_kick_inv=0.0;

            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) // particle - tree
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(_inv_nest_sd, _i, *(_bin.getMemberAsTree(k)));
                else  // particle - particle
                    gt_kick_inv += calcAccPotAndGTKickInvTwo(_inv_nest_sd, _i, _bin.getMemberIndex(k));
            }
            return gt_kick_inv;
        }

        //! calc crossing force, potential and inverse time transformation factor between binary tree i and binary tree j
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _bini: binary tree i for walking
          @param[in] _binj: binary tree j for walking
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvCrossTreeIter(const Float& _inv_nest_sd, AR::BinaryTree<Tparticle>& _bini, AR::BinaryTree<Tparticle>& _binj) {
            Float gt_kick_inv=0.0;
            for (int k=0; k<2; k++) { 
                if (_bini.isMemberTree(k))  // tree - tree
                    gt_kick_inv += calcAccPotAndGTKickInvCrossTreeIter(_inv_nest_sd, *(_bini.getMemberAsTree(k)), _binj);
                else  // particle - tree
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(_inv_nest_sd, _bini.getMemberIndex(k), _binj);
            }
            return gt_kick_inv;
        }

        //! calc force, potential and inverse time transformation factor for kick
        /*!
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary to drift pos
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        Float calcAccPotAndGTKickInvTreeIter(const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin) {
            // current nested sd factor
            Float inv_nest_sd = _inv_nest_sd_up/_bin->slowdown.getSlowDownFactor();
            Float gt_kick_inv = 0.0;

            // check left 
            if (_bin.isMemberTree(0)) { // left is tree
                auto* bin_left = _bin.getMemberAsTree(0);
                // inner interaction of left tree
                gt_kick_inv += calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_left);

                if (_bin.isMemberTree(1)) { // right is tree
                    auto* bin_right = _bin.getMemberAsTree(1);

                    // inner interaction of right tree
                    gt_kick_inv += calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_right);
                    
                    // cross interaction
                    gt_kick_inv += calcAccPotAndGTKickInvCrossTreeIter(inv_nest_sd, *bin_left, *bin_right);
                }
                else { // right is particle
                    // cross interaction from particle j to tree left
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(inv_nest_sd,  _bin.getMemberIndex(1), *_bin_left);
                }
            }
            else { // left is particle
                if (_bin.isMemberTree(1)) { // right is tree
                    auto* bin_right = _bin.getMemberAsTree(1);
                    // inner interaction of right tree
                    gt_kick_inv += calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_right);

                    // cross interaction from particle i to tree right
                    gt_kick_inv += calcAccPotAndGTKickInvOneTreeIter(inv_nest_sd, _bin.getMemberIndex(0), *_bin_right);
                }
                else { // right is particle
                    // particle - particle interaction
                    gt_kick_inv += calcAccPotAndGTKickInvTwo(inv_nest_sd, _bin.getMemberIndex(0), _bin.getMemberIndex(1));
                }
            }

            return gt_kick_inv;
        }
        

        //! calc force, potential and inverse time transformation factor for kick
        /*!
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        inline Float calcAccPotAndGTKickInv() {
            return calcAccPotAndGTKickInvTreeIter(1.0, info.getBinaryTreeRoot());
            // pertuber force
            manager->interaction.calcAccPert(force_.getDataAddress(), particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());
        }

#ifdef AR_TTL
        //! kick energy and time transformation function for drift of binary tree 
        /*!
          @param[in] _dt: time step
          @param[in] _vel_sd_up: upper cm sd vel
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary tree for kick etot and calc dgt_drift
          \return gt_drift_inv change 
        */
        Float kickEtotAndGTDriftTreeIter(const Float& _dt, const Float* _vel_sd_up, const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin) {
            // current nested sd factor
            Float inv_nest_sd = _inv_nest_sd_up/_bin->slowdown.getSlowDownFactor();
            Float dgt_drift_inv = 0.0;
            Float de = 0.0;

            Float* vel_cm = _bin.getVel();

            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    Float* vel = bink->getVel();
                    Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                       (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                       (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]}; 
                    dgt_drift_inv += kickEtotAndGTDriftTreeIter(_dt, vel_sd, inv_nest_sd, *bink);
                }
                else {
                    int i = _bin.getMemberIndex(k);
                    ASSERT(i>=0&&i<particles.getSize());
                    ASSERT(&particles[i]==_bin.getMember(k));
                    
                    Float* gtgrad = force_[i].gtgrad;
                    Float* vel = particles[i].getVel();
                    Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                       (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                       (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]};

                    de += mass * (vel[0] * pert[0] +
                                  vel[1] * pert[1] +
                                  vel[2] * pert[2]);

                    dgt_drift_inv +=  (vel_sd[0] * gtgrad[0] +
                                       vel_sd[1] * gtgrad[1] +
                                       vel_sd[2] * gtgrad[2]);
                }
            }
            etot_ref_     += _dt * de;
            etot_sd_ref_ += _dt * de;

            return dgt_drift_inv;
        }

        //! kick energy and time transformation function for drift
        /*!
          @param[in] _dt: time step
        */
        void kickEtotAndGTDrift(const Float _dt) {
            // the particle cm velocity is zero (assume in rest-frame)
            ASSERT(!particles.isOriginFrame());
            auto& bin_root=info.getBinaryTreeRoot();
            Float vel_cm[3] = {0.0,0.0,0.0};
            Float sd_factor=1.0;

            Float dgt_drift_inv = kickEtotAndGTDriftTreeIter(_dt, vel_cm, sd_factor, bin_root);
            gt_drift_inv_ += dgt_drift_inv*_dt;
        }

#else //! AR_TTL
        //! kick energy
        /*!
          @param[in] _dt: time step
        */
        inline void kickEtot(const Float _dt) {
            Float de = 0.0;
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].getVel();
                Float* pert= force[i].acc_pert;
                de += mass * (vel[0] * pert[0] +
                              vel[1] * pert[1] +
                              vel[2] * pert[2]);
            }
            etot_sd_ref_ += _dt * de;
            etot_ref_ += _dt * de;
        }
#endif // AR_TTL

#else //! AR_SLOWDOWN_TREE

#ifdef AR_SLOWDOWN_ARRAY
        //! correct force, potential energy and gt_kick_inv based on slowdown for inner binaries
        /*! 
          @param[in,out] _gt_kick_inv: the inverse time transformation factor for kick step (input), be corrected with slowdown (output)
         */
        inline void correctAccPotGTKickInvSlowDownInner(Float& _gt_kick_inv) {
            int n = binary_slowdown.getSize();
            Float gt_kick_inv_cor = 0.0;
            Float de = 0.0;
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                int i1 = sdi->getMemberIndex(0);
                int i2 = sdi->getMemberIndex(1);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv = 1.0/kappa;
                if (i1>=0) {
                    ASSERT(i2>=0);
                    ASSERT(i1!=i2);
                    ASSERT(i1<particles.getSize());
                    ASSERT(i2<particles.getSize());
                    
                    // calculate pair interaction
                    Force fi[2];
                    Float epoti;
                    Float gt_kick_inv_i = manager->interaction.calcInnerAccPotAndGTKickInvTwo(fi[0], fi[1], epoti, particles[i1], particles[i2]);
                                    // scale binary pair force with slowdown 
                    force_[i1].acc_in[0] += fi[0].acc_in[0]*kappa_inv - fi[0].acc_in[0];
                    force_[i1].acc_in[1] += fi[0].acc_in[1]*kappa_inv - fi[0].acc_in[1];
                    force_[i1].acc_in[2] += fi[0].acc_in[2]*kappa_inv - fi[0].acc_in[2];
                    force_[i2].acc_in[0] += fi[1].acc_in[0]*kappa_inv - fi[1].acc_in[0];
                    force_[i2].acc_in[1] += fi[1].acc_in[1]*kappa_inv - fi[1].acc_in[1];
                    force_[i2].acc_in[2] += fi[1].acc_in[2]*kappa_inv - fi[1].acc_in[2];

                    de += epoti*kappa_inv - epoti;
                    // scale gtgrad with slowdown
#ifdef AR_TTL
                    force_[i1].gtgrad[0] += fi[0].gtgrad[0]*kappa_inv - fi[0].gtgrad[0];
                    force_[i1].gtgrad[1] += fi[0].gtgrad[1]*kappa_inv - fi[0].gtgrad[1];
                    force_[i1].gtgrad[2] += fi[0].gtgrad[2]*kappa_inv - fi[0].gtgrad[2];
                    force_[i2].gtgrad[0] += fi[1].gtgrad[0]*kappa_inv - fi[1].gtgrad[0];
                    force_[i2].gtgrad[1] += fi[1].gtgrad[1]*kappa_inv - fi[1].gtgrad[1];
                    force_[i2].gtgrad[2] += fi[1].gtgrad[2]*kappa_inv - fi[1].gtgrad[2];
#endif

                    // gt kick
                    gt_kick_inv_cor += gt_kick_inv_i*(kappa_inv - 1.0);
                }
            }

            // global slowdown
            const Float kappa_inv_global = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            for (int i=0; i<force_.getSize(); i++) {
                force_[i].acc_in[0] *= kappa_inv_global;
                force_[i].acc_in[1] *= kappa_inv_global;
                force_[i].acc_in[2] *= kappa_inv_global;
#ifdef AR_TTL
                force_[i].gtgrad[0] *= kappa_inv_global;
                force_[i].gtgrad[1] *= kappa_inv_global;
                force_[i].gtgrad[2] *= kappa_inv_global;
#endif
            }

            epot_sd_ = (epot_ + de)*kappa_inv_global;
            _gt_kick_inv = (_gt_kick_inv + gt_kick_inv_cor)*kappa_inv_global;
        }

        //! correct postion drift due to inner binary slowdown
        /*! 
          @param[in] _dt: time step
          @param[in] _sd_global_inv: global slowdown factor inverse
         */
        inline void correctPosSlowDownInner(const Float _dt, const Float _sd_global_inv) {
            int n = binary_slowdown.getSize();
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv_m_one = (1.0/kappa - 1.0)*_sd_global_inv;
                Float* velcm = sdi->getVel();
                for (int k=0; k<2; k++) {
                    int j = sdi->getMemberIndex(k);
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

        //! update c.m. for binaries with slowdown inner
        inline void updateCenterOfMassForBinaryWithSlowDownInner() {
            int n = binary_slowdown.getSize();
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                int i1 = sdi->getMemberIndex(0);
                int i2 = sdi->getMemberIndex(1);
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

                sdi->mass = mcm;
                Float* pos = sdi->getPos();
                pos[0] = (m1*pos1[0] + m2*pos2[0])*mcminv;
                pos[1] = (m1*pos1[1] + m2*pos2[1])*mcminv;
                pos[2] = (m1*pos1[2] + m2*pos2[2])*mcminv;

                Float* vel = sdi->getVel();
                vel[0] = (m1*vel1[0] + m2*vel2[0])*mcminv;
                vel[1] = (m1*vel1[1] + m2*vel2[1])*mcminv;
                vel[2] = (m1*vel1[2] + m2*vel2[2])*mcminv;
            }
        }

        //! calculate kinetic energy with slowdown factor
        /*! @param[in] _ekin: total kinetic energy without slowdown
         */
        inline void calcEkinSlowDownInner(const Float& _ekin) {
            int n = binary_slowdown.getSize();
            if (n==0) return;
            const Float kappa_inv_global = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            Float de = Float(0.0);
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv_m_one = 1.0/kappa - 1.0;
                Float* velcm = sdi->getVel();
                for (int k=0; k<2; k++) {
                    int j = sdi->getMemberIndex(k);
                    ASSERT(j>=0&&j<particles.getSize());
                    Float* vel = particles[j].getVel();

                    // only scale velocity referring to binary c.m.
                    Float vrel[3] = { vel[0] - velcm[0], 
                                      vel[1] - velcm[1], 
                                      vel[2] - velcm[2]}; 

                    de += kappa_inv_m_one * particles[j].mass * (vrel[0]*vrel[0] + vrel[1]*vrel[1] + vrel[2]*vrel[2]);
                }
            }
            ekin_sd_ = (_ekin + 0.5*de)*kappa_inv_global;
        }

        //! find inner binaries for slowdown treatment iteration function
        static int findSlowDownInnerBinaryIter(COMM::List<AR::BinaryTree<Tparticle>*>& _binary_slowdown, const int& _c1, const int& _c2, AR::BinaryTree<Tparticle>& _bin) {
            // find leaf binary
            if (_bin.getMemberN()==2 && _bin.semi>0.0) {
                _binary_slowdown.increaseSizeNoInitialize(1);
                auto& sdi_new = _binary_slowdown.getLastMember();
                sdi_new = &_bin;
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
            binary_slowdown.resizeNoInitialize(1);
            int ncount[2]={0,0};
            int nsd = bin_root.processTreeIter(binary_slowdown, ncount[0], ncount[1], findSlowDownInnerBinaryIter);
            ASSERT(nsd==binary_slowdown.getSize()-1);

            for (int i=1; i<=nsd; i++) {
#ifdef AR_SLOWDOWN_MASSRATIO
                const Float mass_ratio = manager->slowdown_mass_ref/binary_slowdown[i].bin->mass;
#else 
                const Float mass_ratio = 1.0;
#endif
                binary_slowdown[i]->slowdown.initialSlowDownReference(mass_ratio*manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);
                binary_slowdown[i]->slowdown.setUpdateTime(time_);
            }
        }

#endif // AR_SLOWDOWN_ARRAY

        //! Calculate kinetic energy
        inline void calcEKin(){
            ekin_ = Float(0.0);
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            for (int i=0; i<num; i++) {
                const Float *vi=pdat[i].getVel();
                ekin_ += 0.5 * pdat[i].mass * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
            }
#ifdef AR_SLOWDOWN_ARRAY
            calcEkinSlowDownInner(ekin_);
#endif
        }

        //! kick velocity
        /*! First time step will be calculated, the velocities are kicked
          @param[in] _dt: time size
        */
        inline void kickVel(const Float _dt) {
            const int num = particles.getSize();            
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0; i<num; i++) {
                // kick velocity
                Float* vel = pdat[i].getVel();
                Float* acc = force[i].acc_in;
                Float* pert= force[i].acc_pert;
                // half dv 
                vel[0] += _dt * (acc[0] + pert[0]);
                vel[1] += _dt * (acc[1] + pert[1]);
                vel[2] += _dt * (acc[2] + pert[2]);
            }
#ifdef AR_SLOWDOWN_ARRAY
            // update c.m. of binaries 
            updateCenterOfMassForBinaryWithSlowDownInner();
#endif
        }

        //! drift time and position
        /*! First (real) time is drifted, then positions are drifted
          @param[in] _dt: time step
        */
        inline void driftTimeAndPos(const Float _dt) {
            // drift time 
            time_ += _dt;

            // drift position
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
#ifdef AR_SLOWDOWN_ARRAY
            const Float kappa_inv = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            const Float dt_sd = _dt * kappa_inv;

            for (int i=0; i<num; i++) {
                Float* pos = pdat[i].getPos();
                Float* vel = pdat[i].getVel();
                pos[0] += dt_sd * vel[0]; 
                pos[1] += dt_sd * vel[1];
                pos[2] += dt_sd * vel[2];
            }

            // correct postion drift due to inner binary slowdown
            correctPosSlowDownInner(_dt, kappa_inv);
#else
            for (int i=0; i<num; i++) {
                Float* pos = pdat[i].getPos();
                Float* vel = pdat[i].getVel();
                pos[0] += _dt * vel[0];
                pos[1] += _dt * vel[1];
                pos[2] += _dt * vel[2];
            }
#endif
        }

        //! calc force, potential and inverse time transformation factor for kick
        /*!
          \return gt_kick_inv: inverse time transformation factor for kick
         */
        inline Float calcAccPotAndGTKickInv() {
            Float gt_kick_inv = manager->interaction.calcAccPotAndGTKickInv(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());            

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
#ifdef AR_SLOWDOWN_ARRAY
            // slowdown binary acceleration
            correctAccPotGTKickInvSlowDownInner(gt_kick_inv);
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
            return gt_kick_inv;
        }


#ifdef AR_TTL
        //! kick energy and time transformation function for drift
        /*!
          @param[in] _dt: time step
        */

        inline void kickEtotAndGTDrift(const Float _dt) {
            Float de = Float(0.0);
            Float dg = Float(0.0);
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].getVel();
                Float* pert= force[i].acc_pert;
                Float* gtgrad=force[i].gtgrad;
                de += mass * (vel[0] * pert[0] +
                              vel[1] * pert[1] +
                              vel[2] * pert[2]);
                dg +=  (vel[0] * gtgrad[0] +
                        vel[1] * gtgrad[1] +
                        vel[2] * gtgrad[2]);
            }
            etot_ref_ += _dt * de;

#ifdef AR_SLOWDOWN_ARRAY
            etot_sd_ref_ += _dt * de;

            // correct gt_drift_inv
            const Float kappa_inv_global = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            int n = binary_slowdown.getSize();
            for (int i=1; i<n; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv = 1.0/kappa;
                Float* velcm = sdi->getVel();
                for (int k=0; k<2; k++) {
                    int j = sdi->getMemberIndex(k);
                    if (j>=0) {
                        ASSERT(j<particles.getSize());
                        Float* gtgrad=force_[j].gtgrad;
                        Float* vel = particles[j].getVel();
                        Float vrel[3] = { vel[0] - velcm[0], 
                                          vel[1] - velcm[1], 
                                          vel[2] - velcm[2]}; 
                        dg +=  (vrel[0] * (kappa_inv-1)* gtgrad[0] +
                                vrel[1] * (kappa_inv-1)* gtgrad[1] +
                                vrel[2] * (kappa_inv-1)* gtgrad[2]);
                    }
                }
            }
            gt_drift_inv_ += _dt * dg *kappa_inv_global;
#else
            gt_drift_inv_ += _dt * dg;
#endif
        }

#else //!AR_TTL
        //! kick energy
        /*!
          @param[in] _dt: time step
        */
        inline void kickEtot(const Float _dt) {
            Float de = 0.0;
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].getVel();
                Float* pert= force[i].acc_pert;
                de += mass * (vel[0] * pert[0] +
                              vel[1] * pert[1] +
                              vel[2] * pert[2]);
            }
#ifdef AR_SLOWDOWN_ARRAY
            etot_sd_ref_ += _dt * de;
#endif
            etot_ref_ += _dt * de;
        }
#endif //AR_TTL

#endif // AR_SLOWDOWN_TREE
            
    public:
#ifdef AR_SLOWDOWN_ARRAY
        //! update slowdown factor based on perturbation and record slowdown energy change
        /*! Update slowdown inner and global, update gt_inv
            @param [in] _update_energy_flag: Record cumulative slowdown energy change if true;
         */
        void updateSlowDownAndCorrectEnergy(const bool _update_energy_flag) {
            
            auto& sd_root = binary_slowdown[0]->slowdown;
            Float sd_backup = sd_root.getSlowDownFactor();
            //if (time_>=sd_root.getUpdateTime()) {
            sd_root.pert_in = manager->interaction.calcPertFromBinary(*binary_slowdown[0]);
            sd_root.period = binary_slowdown[0]->period;

            Float t_min_sq= NUMERIC_FLOAT_MAX;
            sd_root.pert_out = 0.0;
            manager->interaction.calcSlowDownPert(sd_root.pert_out, t_min_sq, getTime(), particles.cm, perturber);
            
            sd_root.timescale = std::min(sd_root.getTimescaleMax(), sqrt(t_min_sq));
            sd_root.calcSlowDownFactor();

            sd_root.increaseUpdateTimeOnePeriod();
            //}
            // inner binary slowdown
            int n = binary_slowdown.getSize();
            bool modified_flag=false;
            for (int i=1; i<n; i++) {
                auto* sdi = binary_slowdown[i];
                //if (time_>=sdi->slowdown.getUpdateTime()) {
                sdi->calcCenterOfMass();
                calcSlowDownInnerBinary(*sdi);
                sdi->slowdown.increaseUpdateTimeOnePeriod();
                modified_flag=true;
                //}
            }    

            if (_update_energy_flag) {
                Float ekin_sd_bk = ekin_sd_;
                Float epot_sd_bk = epot_sd_;
                Float etot_sd_ref_bk = etot_sd_ref_;
#ifdef AR_TTL
                Float gt_drift_inv_bk = gt_drift_inv_;
#endif

                if (modified_flag) {
                    // initialize the gt_drift_inv_ with new slowdown factor
                    Float gt_kick_inv_sdi = manager->interaction.calcAccPotAndGTKickInv(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());
                    correctAccPotGTKickInvSlowDownInner(gt_kick_inv_sdi);
#ifdef AR_TTL
                    gt_drift_inv_ += gt_kick_inv_sdi - gt_kick_inv_;
                    gt_kick_inv_ = gt_kick_inv_sdi;
#endif
                    // correct etot_sd_ref_ with new slowdown
                    calcEkinSlowDownInner(ekin_);
                }
                else {
                    Float kappa_inv = 1.0/sd_root.getSlowDownFactor();
                    Float gt_kick_inv_sdi = gt_kick_inv_*sd_backup*kappa_inv;
#ifdef AR_TTL
                    gt_drift_inv_ += gt_kick_inv_sdi - gt_kick_inv_;
                    gt_kick_inv_ = gt_kick_inv_sdi;
#endif
                    // only need to correct the total value
                    ekin_sd_ = ekin_*kappa_inv;
                    epot_sd_ = epot_*kappa_inv;
                }

                Float de_sd = (ekin_sd_ - ekin_sd_bk) + (epot_sd_ - epot_sd_bk);
                etot_sd_ref_ += de_sd;
#ifdef AR_TTL
                Float dH_sd = (ekin_sd_ + epot_sd_ - etot_sd_ref_)*gt_drift_inv_ - (ekin_sd_bk + epot_sd_bk - etot_sd_ref_bk)*gt_drift_inv_bk;
#else
                Float dH_sd = manager->interaction.calcH(ekin_sd_ - etot_sd_ref_, epot_sd_) 
                    - manager->interaction.calcH(ekin_sd_bk - etot_sd_ref_bk, epot_sd_bk);
#endif

                // add slowdown change to the global slowdown energy
                de_sd_change_cum_ += de_sd;
                dH_sd_change_cum_ += dH_sd;
            }
        }
#endif

        //! initialization for integration
        /*! initialize the system. Acceleration, energy and time transformation factors are updated. If the center-of-mass is not yet calculated, the system will be shifted to center-of-mass frame.
          @param[in] _time: real physical time to initialize
        */
        void initialIntegration(const Float _time) {
            ASSERT(checkParams());

            // particle number and data address
            const int n_particle = particles.getSize();

            // resize force array
            force_.resizeNoInitialize(n_particle);

            // Initial intgrt value t (avoid confusion of real time when slowdown is used)
            time_ = _time;

            // check particle number
            ASSERT(particles.getSize()>=2);

            // reset particle modification flag
            particles.setModifiedFalse();

            // check the center-of-mass initialization
            if(particles.isOriginFrame()) {
                particles.calcCenterOfMass();
                particles.shiftToCenterOfMassFrame();
                for (int i=0; i<info.binarytree.getSize(); i++) {
                    auto& bin = info.binarytree[i];
                    bin.pos[0] -= particles.cm.pos[0];
                    bin.pos[1] -= particles.cm.pos[1];
                    bin.pos[2] -= particles.cm.pos[2];
                    bin.vel[0] -= particles.cm.vel[0];
                    bin.vel[1] -= particles.cm.vel[1];
                    bin.vel[2] -= particles.cm.vel[2];
                }
            }
            ASSERT(info.getBinaryTreeRoot().pos[0]*info.getBinaryTreeRoot().pos[0]<1e-10);
            ASSERT(info.getBinaryTreeRoot().vel[0]*info.getBinaryTreeRoot().vel[0]<1e-10);

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
#ifdef AR_SLOWDOWN_TREE
            info.initialSlowDownReference(manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);
#else
            binary_slowdown.increaseSizeNoInitialize(1);
            binary_slowdown[0] = &info.getBinaryTreeRoot();

            // set slowdown reference
            SlowDown& slowdown_root = info.getBinaryTreeRoot().slowdown;

            // slowdown for the system
            slowdown_root.initialSlowDownReference(manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);

            if (particles.getSize()>2) {
                findSlowDownInner(time_);
                // update c.m. of binaries 
                //updateCenterOfMassForBinaryWithSlowDownInner();
            }
#endif

            updateSlowDownAndCorrectEnergy(false);

#ifdef AR_TTL
            gt_kick_inv_ = calcAccPotAndGTKickInv();

            // initially gt_drift 
            gt_drift_inv_ = gt_kick_inv_;

#else
            calcAccPotAndGTKickInv();
#endif

            calcEKin();

            etot_ref_ = ekin_ + epot_;
            etot_sd_ref_ = ekin_sd_ + epot_sd_;

            Float de_sd = etot_sd_ref_ - etot_ref_;
            etot_sd_ref_ += de_sd;

            // add slowdown change to the global slowdown energy
            de_sd_change_cum_ += de_sd;
            dH_sd_change_cum_ = 0.0;

#else // No SLOWDOWN
            Tparticle* particle_data = particles.getDataAddress();
            Force* force_data = force_.getDataAddress();

#ifdef AR_TTL
            gt_kick_inv_ = manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm,  perturber, _time);

            // initially gt_drift 
            gt_drift_inv_ = gt_kick_inv_;

#else
            manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm,  perturber, _time);
#endif

            // calculate kinetic energy
            calcEKin();

            // initial total energy
            etot_ref_ = ekin_ + epot_;

#endif
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

                // inverse time transformation factor for drift
#ifdef AR_TTL
                Float gt_drift_inv = gt_drift_inv_;
#else 
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
                Float gt_drift_inv = manager->interaction.calcGTDriftInv(ekin_sd_-etot_sd_ref_); // pt = -etot
#else
                Float gt_drift_inv = manager->interaction.calcGTDriftInv(ekin_-etot_ref_); // pt = -etot
#endif
#endif

                // drift
                Float dt_drift = ds_drift/gt_drift_inv;

                // drift time and postion
                driftTimeAndPos(dt_drift);
                _time_table[i] = time_;

                // step for kick
                Float ds_kick = manager->step.getDK(i)*_ds;

                //! calc force, potential and inverse time transformation factor for kick
                Float gt_kick_inv = calcAccPotAndGTKickInv();

                // time step for kick
                Float dt_kick = ds_kick/gt_kick_inv;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

#ifdef AR_TTL   
                // back up gt_kick 
                gt_kick_inv_ = gt_kick_inv;
                // kick total energy and inverse time transformation factor for drift
                kickEtotAndGTDrift(dt_kick);
#else
                // kick total energy 
                kickEtot(dt_kick);
#endif
                // kick half step for velocity
                kickVel(0.5*dt_kick);

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

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            const Float kappa_inv = 1.0/info.getBinaryTreeRoot().slowdown.calcSlowDownFactor();
#endif

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
                // inverse time transformation factor for drift
#ifdef AR_TTL
                Float gt_inv = gt_drift_inv_;
#else 
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
                Float gt_inv = manager->interaction.calcGTDriftInv(ekin_sd_-etot_sd_ref_); // pt = -etot_sd
#else
                Float gt_inv = manager->interaction.calcGTDriftInv(ekin_-etot_ref_); // pt = -etot
#endif
#endif
                // drift
                Float dt = ds/gt_inv;
                ASSERT(!isnan(dt));
                
                // drift time 
                time_ += dt;

                // update real time
                _time_table[i] = time_;

                Float dt_sd = dt*kappa_inv;

                // drift position
                pos1[0] += dt_sd * vel1[0];
                pos1[1] += dt_sd * vel1[1];
                pos1[2] += dt_sd * vel1[2];

                pos2[0] += dt_sd * vel2[0];
                pos2[1] += dt_sd * vel2[1];
                pos2[2] += dt_sd * vel2[2];

                // step for kick
                ds = manager->step.getDK(i)*_ds;

                gt_inv = manager->interaction.calcAccPotAndGTKickInv(force_data, epot_, particle_data, n_particle, particles.cm, perturber, _time_table[i]);

                ASSERT(!isnan(epot_));

                // kick half step for velocity
                Float dvel1[3], dvel2[3];

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
                // time step for kick
                gt_inv *= kappa_inv;

                dt = 0.5*ds/gt_inv;

                dvel1[0] = dt * (acc1[0]*kappa_inv + pert1[0]);
                dvel1[1] = dt * (acc1[1]*kappa_inv + pert1[1]);
                dvel1[2] = dt * (acc1[2]*kappa_inv + pert1[2]);

                dvel2[0] = dt * (acc2[0]*kappa_inv + pert2[0]);
                dvel2[1] = dt * (acc2[1]*kappa_inv + pert2[1]);
                dvel2[2] = dt * (acc2[2]*kappa_inv + pert2[2]);
#else
                dt = 0.5*ds/gt_inv;

                dvel1[0] = dt * (acc1[0] + pert1[0]);
                dvel1[1] = dt * (acc1[1] + pert1[1]);
                dvel1[2] = dt * (acc1[2] + pert1[2]);

                dvel2[0] = dt * (acc2[0] + pert2[0]);
                dvel2[1] = dt * (acc2[1] + pert2[1]);
                dvel2[2] = dt * (acc2[2] + pert2[2]);
#endif

                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];

                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];

                // kick total energy and time transformation factor for drift
                etot_ref_ += 2.0*dt * (mass1* (vel1[0] * pert1[0] + 
                                               vel1[1] * pert1[1] + 
                                               vel1[2] * pert1[2]) +
                                       mass2* (vel2[0] * pert2[0] + 
                                               vel2[1] * pert2[1] + 
                                               vel2[2] * pert2[2]));

#ifdef AR_TTL   
                // back up gt_kick_inv
                gt_kick_inv_ = gt_inv;

#if (defined AR_SLOWDOWN_ARRAY ) || (defined AR_SLOWDOWN_TREE)
                // integrate gt_drift_inv
                gt_drift_inv_ +=  2.0*dt*kappa_inv*kappa_inv* (vel1[0] * gtgrad1[0] +
                                                               vel1[1] * gtgrad1[1] +
                                                               vel1[2] * gtgrad1[2] +
                                                               vel2[0] * gtgrad2[0] +
                                                               vel2[1] * gtgrad2[1] +
                                                               vel2[2] * gtgrad2[2]);
#else
                // integrate gt_drift_inv
                gt_drift_inv_ +=  2.0*dt* (vel1[0] * gtgrad1[0] +
                                           vel1[1] * gtgrad1[1] +
                                           vel1[2] * gtgrad1[2] +
                                           vel2[0] * gtgrad2[0] +
                                           vel2[1] * gtgrad2[1] +
                                           vel2[2] * gtgrad2[2]);
#endif 

#endif // AR_TTL

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

#if (defined AR_SLOWDOWN_ARRAY ) || (defined AR_SLOWDOWN_TREE)
            // make consistent slowdown inner energy 
            etot_sd_ref_ = etot_ref_*kappa_inv;
            ekin_sd_ = ekin_*kappa_inv;
            epot_sd_ = epot_*kappa_inv;
#endif
        }
        
        // Integrate the system to a given time
        /*!
          @param[in] _ds: the integration step size
          @param[in] _time_end: the expected finishing real time 
          \return binary tree of the pair which triggers interruption condition
         */
        InterruptBinary<Tparticle> integrateToTime(const Float _time_end) {
            ASSERT(checkParams());

            // real full time step
            const Float dt_full = _time_end - time_;

            // time error 
            const Float time_error = manager->time_error_max;

            // energy error limit
            const Float energy_error_rel_max = manager->energy_error_relative_max;
            // expect energy_error using half step if energy_error_rel_max reached
            //const Float energy_error_rel_max_half_step = energy_error_rel_max * manager->step.calcErrorRatioFromStepModifyFactor(0.5);
            //const Float dt_min = manager->time_step_min;

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

            // two switch step control
            Float ds[2] = {info.ds,info.ds}; // step with a buffer
            Float ds_init   = info.ds;  //backup initial step
            int   ds_switch=0;   // 0 or 1

            // reduce ds control, three level
            const int n_reduce_level_max=2;

            struct DsBackupManager{
                Float ds_backup[n_reduce_level_max+1];
                int n_step_wait_recover_ds[n_reduce_level_max+1];
                int n_reduce_level; 

                DsBackupManager(const Float _ds) { initial(_ds), n_reduce_level = -1; }

                void initial(const Float _ds) {
                    for (int i=0; i<n_reduce_level_max; i++) {
                        ds_backup[i] = _ds;
                        n_step_wait_recover_ds[i] = 0;
                    }
                }

                // shift backup by one level, if successful, return true
                bool shiftReduceLevel() {
                    if (n_reduce_level>0) {
                        for (int i=0; i<n_reduce_level-1; i++) {
                            ds_backup[i] =ds_backup[i+1];
                            n_step_wait_recover_ds[i] = n_step_wait_recover_ds[i+1];
                        }
                        n_reduce_level--;
                        return true;
                    }
                    return false;
                }

                // record ds to backup
                void backup(const Float _ds, const Float _modify_factor) {
                    //if (n_reduce_level==n_reduce_level_max) shiftReduceLevel();  // not converse sometimes, becomes infinite small steps
                    if (n_reduce_level==n_reduce_level_max) 
                        n_step_wait_recover_ds[n_reduce_level] *= 2*to_int(1.0/_modify_factor);
                    else {
                        n_reduce_level++;
                        ds_backup[n_reduce_level] = _ds;
                        n_step_wait_recover_ds[n_reduce_level] = 2*to_int(1.0/_modify_factor);
                    }
                }

                // count step (return false) and recover ds if necessary (return true)
                bool countAndRecover(Float &_ds, Float &_modify_factor) {
                    if (n_reduce_level>=0) {
                        if (n_step_wait_recover_ds[n_reduce_level] ==0) {
                            _modify_factor = ds_backup[n_reduce_level]/_ds;
                            _ds = ds_backup[n_reduce_level];
                            n_step_wait_recover_ds[n_reduce_level] = -1;
                            n_reduce_level--;
                            return true;
                        }
                        else {
                            n_step_wait_recover_ds[n_reduce_level]--;
                            return false;
                        }
                    }
                    return false;
                }

            } ds_backup(info.ds);

            int reduce_ds_count=0; // number of reduce ds (ignore first few steps)
            Float step_modify_factor=1.0; // step modify factor 
            Float previous_step_modify_factor=1.0; // step modify factor 
            Float previous_error_ratio=-1; // previous error ratio when ds is reduced
            bool previous_is_restore=false; // previous step is reduced or not

            // time end control
            int n_step_end=0;  // number of steps integrated to reach the time end for one during the time sychronization sub steps
            bool time_end_flag=false; // indicate whether time reach the end

            // step count
            long long int step_count=0; // integration step 
            long long int step_count_tsyn=0; // time synchronization step
            InterruptBinary<Tparticle> bin_interrupt = {NULL, time_, _time_end, InterruptStatus::none};
            InterruptBinary<Tparticle> bin_interrupt_return = bin_interrupt;
            
            // particle data
            const int n_particle = particles.getSize();

/* This must suppress since after findslowdowninner, slowdown inner is reset to 1.0, recalculate ekin_sdi give completely wrong value for energy correction for slowdown change later
#ifdef AR_DEBUG
            Float ekin_check = ekin_;
            calcEKin();
            ASSERT(abs(ekin_check-ekin_)<1e-10);
            ekin_ = ekin_check;
#endif
*/


#ifdef AR_DEBUG_DUMP
            // back up initial data
            backupIntData(backup_data_init);
#endif
      
#ifdef AR_SLOWDOWN_ARRAY
            // find new inner slowdown binaries, the binary tree data may be modified, thus it is safer to recheck slowdown inner binary at beginning to avoid memory issue (bin is pointer).
#ifdef AR_TTL
            if (n_particle >2) {
                int nold = binary_slowdown.getSize();
                findSlowDownInner(time_);
                int nnew = binary_slowdown.getSize();
                // in case slowdown is disabled in the next step, gt_drift_inv_ should be re-initialized
                if (nold>0&&nnew==0) {
                    gt_kick_inv_ = manager->interaction.calcAccPotAndGTKickInv(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize(), particles.cm, perturber, time_);
                    gt_drift_inv_ = gt_kick_inv_;
                }
            }
#else
            if (n_particle >2) findSlowDownInner(time_);
#endif
#endif
            // integration loop
            while(true) {
                // backup data
                if(backup_flag) {
                    // update slowdown and correct slowdown energy and gt_inv
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
                    if (!time_end_flag) updateSlowDownAndCorrectEnergy(true);
#endif

                    int bk_return_size = backupIntData(backup_data);
                    ASSERT(bk_return_size == bk_data_size);
                    (void)bk_return_size;

                }
                else { //restore data
                    int bk_return_size = restoreIntData(backup_data);
                    ASSERT(bk_return_size == bk_data_size);
                    (void)bk_return_size;
#ifdef AR_SLOWDOWN_ARRAY
                    // update c.m. of binaries 
                    // binary c.m. is not backup, thus recalculate to get correct c.m. velocity for position drift correction due to slowdown inner (the first drift in integrateonestep assume c.m. vel is up to date)
                    updateCenterOfMassForBinaryWithSlowDownInner();
#elif AR_SLOWDOWN_TREE
                    updateBinaryVelIter(info.getBinaryTreeRoot());
#endif
                }

                // get real time 
                Float dt = time_;

                // integrate one step
                ASSERT(!isinf(ds[ds_switch]));
                if(n_particle==2) integrateTwoOneStep(ds[ds_switch], time_table);
                else integrateOneStep(ds[ds_switch], time_table);

                // real step size
                dt =  time_ - dt;
//                ASSERT(dt>0.0);
                
                step_count++;

                // energy check
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
                Float energy_error_bk = getEnergyErrorSlowDownInnerFromBackup(backup_data);
                Float etot_ref_bk = getEtotSlowDownInnerRefFromBackup(backup_data);
                Float energy_error = getEnergyErrorSlowDownInner();
#else
                Float energy_error_bk = getEnergyErrorFromBackup(backup_data);
                Float etot_ref_bk = getEtotRefFromBackup(backup_data);
                Float energy_error = getEnergyError();
#endif
                Float energy_error_diff = energy_error - energy_error_bk;

                Float energy_error_rel_abs = abs(energy_error_diff/etot_ref_bk);

                // get integration error for extended Hamiltonian
#ifdef AR_TTL
                Float integration_error_rel_abs = abs(energy_error_rel_abs*gt_drift_inv_);
                Float integration_error_rel_cum_abs = abs(energy_error*gt_drift_inv_/etot_ref_bk);
#else
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
                Float gt_drift_inv = manager->interaction.calcGTDriftInv(ekin_sd_-etot_sd_ref_);
#else
                Float gt_drift_inv = manager->interaction.calcGTDriftInv(ekin_-etot_ref_);
#endif
                Float integration_error_rel_abs = energy_error_rel_abs*gt_drift_inv;
                Float integration_error_rel_cum_abs = abs(energy_error*gt_drift_inv/etot_ref_bk);
#endif

                Float integration_error_ratio = energy_error_rel_max/integration_error_rel_abs;
      
                // time error
                Float time_diff_rel = (_time_end - time_)/dt_full;

                // error message print
                auto printMessage = [&](const char* message) {
                    std::cerr<<message<<std::endl;
                    std::cerr<<"  T: "<<time_
                             <<"  dT_err/T: "<<time_diff_rel
                             <<"  ds: "<<ds[ds_switch]
                             <<"  ds_init: "<<ds_init
                             <<"  |Int_err/E|: "<<integration_error_rel_abs
                             <<"  |Int_err_cum/E|: "<<integration_error_rel_cum_abs
                             <<"  |dE/E|: "<<energy_error_rel_abs
                             <<"  dE_cum: "<<energy_error
                             <<"  Etot_sd: "<<etot_ref_bk
                             <<"  T_end_flag: "<<time_end_flag
                             <<"  Step_count: "<<step_count;
                    switch (info.fix_step_option) {
                    case FixStepOption::always:
                        std::cerr<<"  Fix:  always"<<std::endl;
                        break;
                    case FixStepOption::later:
                        std::cerr<<"  Fix:  later"<<std::endl;
                        break;
                    case FixStepOption::none:
                        std::cerr<<"  Fix:  none"<<std::endl;
                        break;
                    default:
                        break;
                    }
                };

#ifdef AR_COLLECT_DS_MODIFY_INFO
                auto collectDsModifyInfo = [&](const char* error_message) {
                    auto& bin = info.getBinaryTreeRoot();
                    std::cerr<<error_message<<": "
                             <<"ds_new "<<ds[1-ds_switch]<<" "
                             <<"ds_init "<<ds_init<<" "
                             <<"modify "<<step_modify_factor<<" "
                             <<"steps "<<step_count<<" "
                             <<"n_mods "<<reduce_ds_count<<" "
                             <<"err "<<integration_error_rel_abs<<" "
                             <<"err/max "<<1.0/integration_error_ratio<<" "
                             <<"dt "<<dt<<" "
                             <<"n_ptcl "<<n_particle<<" "
                             <<"semi "<<bin.semi<<" "
                             <<"ecc "<<bin.ecc<<" "
                             <<"period "<<bin.period<<" "
                             <<"m1 "<<bin.m1<<" "
                             <<"m2 "<<bin.m2<<" "
                             <<"ecca "<<bin.ecca<<" "
                             <<"sd "<<bin.slowdown.getSlowDownFactor()<<" "
                             <<"sd_org "<<bin.slowdown.getSlowDownFactorOrigin()<<" "
                             <<"pert_in "<<bin.slowdown.getPertIn()<<" "
                             <<"pert_out "<<bin.slowdown.getPertOut()<<" "
                             <<std::endl;
                };
#endif 

#ifdef AR_WARN
                // warning for large number of steps
                if(step_count>=manager->step_count_max) {
                    if(step_count%manager->step_count_max==0) {
                        printMessage("Warning: step count is signficiant large");
                        printColumnTitle(std::cerr);
                        std::cerr<<std::endl;
                        printColumn(std::cerr);
                        std::cerr<<std::endl;
#ifdef AR_DEBUG_DUMP
                        DATADUMP("dump_large_step");
#endif
                    }
                }
#endif
          
                // When time sychronization steps too large, abort
                if(step_count_tsyn>manager->step_count_max) {
                    printMessage("Error! step count after time synchronization is too large");
                    printColumnTitle(std::cerr);
                    std::cerr<<std::endl;
                    printColumn(std::cerr);
                    std::cerr<<std::endl;
#ifdef AR_DEBUG_DUMP
//                    restoreIntData(backup_data_init);
                    DATADUMP("dump_large_step");
#endif
                    abort();
                }


#ifdef AR_DEEP_DEBUG
                printMessage("");
                std::cerr<<"Timetable: ";
                for (int i=0; i<cd_pair_size; i++) std::cerr<<" "<<time_table[manager->step.getSortCumSumCKIndex(i)];
                std::cerr<<std::endl;
#endif

                ASSERT(!isnan(integration_error_rel_abs));

                // modify step if energy error is large
                if(integration_error_rel_abs>energy_error_rel_max && info.fix_step_option!=FixStepOption::always) {

                    bool check_flag = true;

                    // check whether already modified
                    if (previous_step_modify_factor!=1.0) {
                        ASSERT(previous_error_ratio>0.0);

                        // if error does not reduce much, do not modify step anymore
                        if (integration_error_ratio>0.5*previous_error_ratio) check_flag=false;
                    }

                    if (check_flag) {
                        // for initial steps, reduce step permanently 
                        if(step_count<3) {

                            // estimate the modification factor based on the symplectic order
                            // limit step_modify_factor to 0.125
                            step_modify_factor = std::max(manager->step.calcStepModifyFactorFromErrorRatio(integration_error_ratio), Float(0.125));
                            ASSERT(step_modify_factor>0.0);

                            previous_step_modify_factor = step_modify_factor;
                            previous_error_ratio = integration_error_ratio;

                            ds[ds_switch] *= step_modify_factor;
                            ds[1-ds_switch] = ds[ds_switch];
                            // permanently reduce ds
                            info.ds = ds[ds_switch];
                            ASSERT(!isinf(info.ds));
                            ds_backup.initial(info.ds);

                            backup_flag = false;
#ifdef AR_COLLECT_DS_MODIFY_INFO
                            collectDsModifyInfo("Large_energy_error");
#endif
                            continue;
                        }
                        // for big energy error, reduce step temparely
                        else if (info.fix_step_option==FixStepOption::none && integration_error_ratio<0.1) {

                            // estimate the modification factor based on the symplectic order
                            // limit step_modify_factor to 0.125
                            step_modify_factor = std::max(manager->step.calcStepModifyFactorFromErrorRatio(integration_error_ratio), Float(0.125));
                            ASSERT(step_modify_factor>0.0);

                            previous_step_modify_factor = step_modify_factor;
                            previous_error_ratio = integration_error_ratio;

                            ds_backup.backup(ds[ds_switch], step_modify_factor);
                            if(previous_is_restore) reduce_ds_count++;

                            ds[ds_switch] *= step_modify_factor;
                            ds[1-ds_switch] = ds[ds_switch];
                            ASSERT(!isinf(ds[ds_switch]));

                            // if multiple times reduction happens, permanently reduce ds
                            //if (reduce_ds_count>3) {
                            //    bool shift_flag = ds_backup.shiftReduceLevel();
                            //    if (!shift_flag) ds_backup.initial(ds[ds_switch]);
                            //    info.ds = ds_backup.ds_backup[0];
                            //    reduce_ds_count=0;
                            //}

                            backup_flag = false;
#ifdef AR_COLLECT_DS_MODIFY_INFO
                            collectDsModifyInfo("Large_energy_error");
#endif
                            continue;
                        }
                    }
                }
// too much output
//#ifdef AR_WARN
//                if(integration_error_rel_abs>100.0*energy_error_rel_max) {
//                    std::cerr<<"Warning: symplectic integrator error > 100*criterion:"<<integration_error_rel_abs<<std::endl;
//                }
//#endif

                // if negative step, reduce step size
                if(!time_end_flag&&dt<0) {
                    // limit step_modify_factor to 0.125
                    step_modify_factor = std::min(std::max(manager->step.calcStepModifyFactorFromErrorRatio(abs(_time_end/dt)), Float(0.0625)),Float(0.5)); 
                    ASSERT(step_modify_factor>0.0);
                    previous_step_modify_factor = step_modify_factor;
                    previous_error_ratio = integration_error_ratio;

                    ds[ds_switch] *= step_modify_factor;
                    ds[1-ds_switch] = ds[ds_switch];
                    ASSERT(!isinf(ds[ds_switch]));

                    // for initial steps, reduce step permanently
                    if (step_count<3) {
                        info.ds = ds[ds_switch];
                        ds_backup.initial(info.ds);
                    }
                    else { // reduce step temparely
                        ds_backup.backup(ds[ds_switch], step_modify_factor);
                    }

                    backup_flag = false;

#ifdef AR_COLLECT_DS_MODIFY_INFO
                    collectDsModifyInfo("Negative_step");
#endif
                    continue;
//                    std::cerr<<"Error! symplectic integrated time step ("<<dt<<") < minimum step ("<<dt_min<<")!\n";
//                    printMessage();
//#ifdef AR_DEBUG_DUMP
//                    DATADUMP("dump_negative_time");
//#endif
//                    abort();
                }

                // check integration time
                if(time_ < _time_end - time_error){
                    // check interrupt condiction
                    if (manager->interrupt_detection_option>0) {
                        auto& bin_root = info.getBinaryTreeRoot();
                        bin_interrupt.time_now = time_;
                        manager->interaction.modifyAndInterruptIter(bin_interrupt, bin_root);
                        //InterruptBinary<Tparticle>* bin_intr_ptr = &bin_interrupt;
                        //bin_intr_ptr = bin_root.processRootIter(bin_intr_ptr, Tmethod::modifyAndInterruptIter);
                        ASSERT(bin_interrupt.checkParams());
                        if (bin_interrupt.status!=InterruptStatus::none) {
                            // the mode return back to the root scope
                            if (manager->interrupt_detection_option==2) {
                                // cumulative step count 
                                profile.step_count = step_count;
                                profile.step_count_tsyn = step_count_tsyn;
                                profile.step_count_sum += step_count;
                                profile.step_count_tsyn_sum += step_count_tsyn;

                                return bin_interrupt;
                            }
                            else {
                                ASSERT(!bin_interrupt.adr->isMemberTree(0));
                                ASSERT(!bin_interrupt.adr->isMemberTree(1));
                                Tparticle* p1 = bin_interrupt.adr->getLeftMember();
                                Tparticle* p2 = bin_interrupt.adr->getRightMember();


                                Float ekin_bk = ekin_;
                                Float epot_bk = epot_;
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
                                Float ekin_sd_bk = ekin_sd_;
                                Float epot_sd_bk = epot_sd_;
                                Float etot_sd_ref_bk = etot_sd_ref_;
#endif

#ifdef AR_TTL
                                Float gt_kick_inv_new = calcAccPotAndGTKickInv();
                                Float gt_drift_inv_bk = gt_drift_inv_;
                                gt_drift_inv_ += gt_kick_inv_new - gt_kick_inv_;
                                gt_kick_inv_ = gt_kick_inv_new;
#else
                                calcAccPotAndGTKickInv();
#endif
                                // calculate kinetic energy
                                calcEKin();

                                // get energy change
                                Float de_sd = (ekin_ - ekin_bk) + (epot_ - epot_bk);
                                etot_ref_ += de_sd;

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
                                de_sd = (ekin_sd_ - ekin_sd_bk) + (epot_sd_ - epot_sd_bk);
                                etot_sd_ref_ += de_sd;
                                
#ifdef AR_TTL
                                Float dH_sd = (ekin_sd_ + epot_sd_ - etot_sd_ref_)*gt_drift_inv_ - (ekin_sd_bk + epot_sd_bk - etot_sd_ref_bk)*gt_drift_inv_bk;
#else // !AR_TTL
                                Float dH_sd = manager->interaction.calcH(ekin_sd_ - etot_sd_ref_, epot_sd_) 
                                    - manager->interaction.calcH(ekin_sd_bk - etot_sd_ref_bk, epot_sd_bk);
#endif // AR_TTL

#else // NO SLOWDOWN
#ifdef AR_TTL
                                Float dH_sd = (ekin_ + epot_ - etot_ref_)*gt_drift_inv_ - (ekin_bk + epot_bk - etot_ref_bk)*gt_drift_inv_bk;
#else // !AR_TTL
                                Float etot_ref_bk = etot_ref_ - de_sd;
                                Float dH_sd = manager->interaction.calcH(ekin_ - etot_ref_, epot_) 
                                    - manager->interaction.calcH(ekin_bk - etot_ref_bk, epot_bk);
#endif // AR_TTL
#endif //SLOWDOWN
                                // add slowdown change to the global slowdown energy
                                de_sd_change_cum_ += de_sd;
                                dH_sd_change_cum_ += dH_sd;

#ifdef AR_DEBUG_PRINT
                                std::cerr<<"Interrupt condition triggered!";
                                std::cerr<<" Time: "<<time_;
                                std::cerr<<" Energy change: dE_SD: "<<de_sd<<" dH_SD: "<<dH_sd;
                                std::cerr<<" Slowdown: "<<bin_root.slowdown.getSlowDownFactor()<<std::endl;
                                bin_interrupt.adr->printColumnTitle(std::cerr);
                                std::cerr<<std::endl;
                                bin_interrupt.adr->printColumn(std::cerr);
                                std::cerr<<std::endl;
                                Tparticle::printColumnTitle(std::cerr);
                                std::cerr<<std::endl;
                                for (int j=0; j<2; j++) {
                                    bin_interrupt.adr->getMember(j)->printColumn(std::cerr);
                                    std::cerr<<std::endl;
                                }
#endif

                                // check merger case
                                if (bin_interrupt.status==InterruptStatus::merge) {
                                    if (n_particle==2) {
                                        Tparticle* p = NULL;
                                        if (p1->mass==0.0) p=p2;
                                        else if (p2->mass==0.0) p=p1;
                                        if (p!=NULL){
                                            Float dt = _time_end - time_;
                                            p->pos[0] += dt * p->vel[0];
                                            p->pos[1] += dt * p->vel[1];
                                            p->pos[2] += dt * p->vel[2];

                                            // cumulative step count 
                                            profile.step_count = step_count;
                                            profile.step_count_tsyn = step_count_tsyn;
                                            profile.step_count_sum += step_count;
                                            profile.step_count_tsyn_sum += step_count_tsyn;

                                            time_ += dt;

                                            return bin_interrupt;
                                        }
                                    }
                                    bin_interrupt_return = bin_interrupt;
                                }
                            }
                            bin_interrupt.clear();
                        }
                    }

                    // step increase depend on n_step_wait_recover_ds
                    if(info.fix_step_option==FixStepOption::none && !time_end_flag) {
                        // waiting step count reach
                        previous_is_restore=ds_backup.countAndRecover(ds[1-ds_switch], step_modify_factor);
                        if (previous_is_restore) {
                            previous_error_ratio = -1;
                            previous_step_modify_factor = 1.0;
#ifdef AR_COLLECT_DS_MODIFY_INFO
                            collectDsModifyInfo("Reuse_backup_ds");
#endif
                        }
                        // increase step size if energy error is small, not works correctly, suppress
                        /*
                        else if(integration_error_rel_abs<energy_error_rel_max_half_step&&integration_error_rel_abs>0.0) {
                            Float integration_error_ratio = energy_error_rel_max/integration_error_rel_abs;
                            Float step_modify_factor = manager->step.calcStepModifyFactorFromErrorRatio(integration_error_ratio);
                            ASSERT(step_modify_factor>0.0);
                            ds[1-ds_switch] *= step_modify_factor;
                            info.ds = ds[1-ds_switch];
                            ASSERT(!isinf(ds[1-ds_switch]));
#ifdef AR_DEEP_DEBUG
                            std::cerr<<"Energy error is small enought for increase step, integration_error_rel_abs="<<integration_error_rel_abs
                                     <<" energy_error_rel_max="<<energy_error_rel_max<<" step_modify_factor="<<step_modify_factor<<" new ds="<<ds[1-ds_switch]<<std::endl;
#endif
                        }
                        */
                    }

                    // time sychronization on case, when step size too small to reach time end, increase step size
                    if(time_end_flag && ds[ds_switch]==ds[1-ds_switch]) {
                        step_count_tsyn++;

                        Float dt_end = _time_end - time_;
                        if(n_step_end>1 && dt<0.3*dt_end) {
                            // dt should be >0.0
                            // ASSERT(dt>0.0);
                            ds[1-ds_switch] = ds[ds_switch] * dt_end/abs(dt);
                            ASSERT(!isinf(ds[1-ds_switch]));
#ifdef AR_DEEP_DEBUG
                            std::cerr<<"Time step dt(real) "<<dt<<" <0.3*(time_end-time)(real) "<<dt_end<<" enlarge step factor: "<<dt_end/dt<<" new ds: "<<ds[1-ds_switch]<<std::endl;
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
                else if(time_ > _time_end + time_error) {
                    time_end_flag = true;
                    backup_flag = false;

                    step_count_tsyn++;
                    n_step_end=0;

                    // check timetable
                    int i=-1,k=0; // i indicate the increasing time index, k is the corresponding index in time_table
                    for(i=0; i<cd_pair_size; i++) {
                        k = manager->step.getSortCumSumCKIndex(i);
                        if(_time_end<time_table[k]) break;
                    }
                    if (i==0) { // first step case
                        ASSERT(time_table[k]>0.0);
                        ds[ds_switch] *= manager->step.getSortCumSumCK(i)*_time_end/time_table[k];
                        ds[1-ds_switch] = ds[ds_switch];
                        ASSERT(!isinf(ds[ds_switch]));
#ifdef AR_DEEP_DEBUG
                        std::cerr<<"Time_end reach, time[k]= "<<time_table[k]<<" time= "<<time_<<" time_end/time[k]="<<_time_end/time_table[k]<<" CumSum_CK="<<manager->step.getSortCumSumCK(i)<<" ds(next) = "<<ds[ds_switch]<<" ds(next_next) = "<<ds[1-ds_switch]<<"\n";
#endif
                    }
                    else { // not first step case, get the interval time 
                        // previous integrated sub time in time table
                        Float time_prev = time_table[manager->step.getSortCumSumCKIndex(i-1)];
                        Float dt_k = time_table[k] - time_prev;
                        Float ds_tmp = ds[ds_switch];
                        // get cumsum CK factor for two steps near the time_end
                        Float cck_prev = manager->step.getSortCumSumCK(i-1);
                        Float cck = manager->step.getSortCumSumCK(i);
                        // in case the time is between two sub step, first scale the next step with the previous step CumSum CK cck(i-1)
                        ASSERT(!isinf(cck_prev));
                        ds[ds_switch] *= cck_prev;  
                        ASSERT(!isinf(ds[ds_switch]));
                        // then next next step, scale with the CumSum CK between two step: cck(i) - cck(i-1) 
                        ASSERT(dt_k>0.0);
                        ds[1-ds_switch] = ds_tmp*(cck-cck_prev)*std::min(Float(1.0),(_time_end-time_prev+time_error)/dt_k); 
                        ASSERT(!isinf(ds[1-ds_switch]));

#ifdef AR_DEEP_DEBUG
                        std::cerr<<"Time_end reach, time_prev= "<<time_prev<<" time[k]= "<<time_table[k]<<" time= "<<time_<<" (time_end-time_prev)/dt="<<(_time_end-time_prev)/dt<<" CumSum_CK="<<cck<<" CumSum_CK(prev)="<<cck_prev<<" ds(next) = "<<ds[ds_switch]<<" ds(next_next) = "<<ds[1-ds_switch]<<" \n";
#endif
                    }
                }
                else {
#ifdef AR_DEEP_DEBUG
                    std::cerr<<"Finish, time_diff_rel = "<<time_diff_rel<<" integration_error_rel_abs = "<<integration_error_rel_abs<<std::endl;
#endif
//#ifdef AR_WARN
//                    if (integration_error_rel_cum_abs>energy_error_rel_max) {
//                        std::cerr<<"AR large energy error at the end! ";
//                        printMessage();
//#ifdef AR_DEBUG_DUMP
////                        restoreIntData(backup_data_init);
//                        DATADUMP("dump_large_error");
//#endif
//                    }
//#endif
                    break;
                }
            }

            // cumulative step count 
            profile.step_count = step_count;
            profile.step_count_tsyn = step_count_tsyn;
            profile.step_count_sum += step_count;
            profile.step_count_tsyn_sum += step_count_tsyn;

            return bin_interrupt_return;
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

#ifdef AR_SLOWDOWN_ARRAY
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
            const Float kappa_inv = 1.0/binary_slowdown[0]->slowdown.getSlowDownFactor();
            for (int i=0; i<particles.getSize(); i++) {
                //ASSERT(particle_adr[i]->mass == particle_data[i].mass);
                particle_adr[i]->mass = particle_data[i].mass;

                particle_adr[i]->pos[0] = particle_data[i].pos[0] + _particle_cm.pos[0];
                particle_adr[i]->pos[1] = particle_data[i].pos[1] + _particle_cm.pos[1];
                particle_adr[i]->pos[2] = particle_data[i].pos[2] + _particle_cm.pos[2];

                particle_adr[i]->vel[0] = particle_data[i].vel[0]*kappa_inv + _particle_cm.vel[0];
                particle_adr[i]->vel[1] = particle_data[i].vel[1]*kappa_inv + _particle_cm.vel[1];
                particle_adr[i]->vel[2] = particle_data[i].vel[2]*kappa_inv + _particle_cm.vel[2];
            }

            // correct inner slowdown velocity
            int nsd= binary_slowdown.getSize();
            for (int i=1; i<nsd; i++) {
                auto& sdi = binary_slowdown[i];
                ASSERT(sdi!=NULL);
                Float kappa = sdi->slowdown.getSlowDownFactor();
                Float kappa_inv_m_one = (1.0/kappa - 1.0)*kappa_inv;
                Float* velcm = sdi->getVel();
                for (int k=0; k<2; k++) {
                    int j = sdi->getMemberIndex(k);
                    ASSERT(j>=0&&j<particles.getSize());
                    Float* vel = particle_data[j].getVel();

                    // only scale velocity referring to binary c.m.
                    Float vrel[3] = { vel[0] - velcm[0], 
                                      vel[1] - velcm[1], 
                                      vel[2] - velcm[2]}; 
                    particle_adr[j]->vel[0] += vrel[0] * kappa_inv_m_one;
                    particle_adr[j]->vel[1] += vrel[1] * kappa_inv_m_one;
                    particle_adr[j]->vel[2] += vrel[2] * kappa_inv_m_one;
                }
            }
        }
#endif

#ifdef AR_SLOWDOWN_TREE

        //! write back particles with slowdown velocity
        /*! write back particles with slowdown velocity to original address
          @param[in] _particle_cm: center of mass particle to calculate the original frame, different from the particles.cm
         */
        template <class Tptcl>
        void writeBackSlowDownParticles(const Tptcl& _particle_cm) {
            //! iteration function using binarytree
            auto writeBackIter =[](const Tptcl& _particle_cm, const Float* _vel_sd_up, const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin) {
                Float inv_nest_sd = _inv_nest_sd_up/_bin->slowdown.getSlowDownFactor();
                Float* vel_cm = _bin.getVel();
                for (int k=0; k<2; k++) {
                    if (_bin.isMemberTree(k)) {
                        auto* bink = _bin.getMemberAsTree(k);
                        Float* vel = bink->getVel();
                        Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                           (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                           (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]}; 
                        writeBackIter(_particle_cm, vel_sd, inv_nest_sd, _bin);
                    }
                    else {
                        int i = _bin.getMemberIndex(k);
                        auto& pk = particles[i];
                        auto* pk_adr = particles.getMemberOriginAddress(i);
                        Float* vel = pk.getVel();
                        Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                           (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                           (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]};
                        pk_adr->mass = pk.mass;

                        pk_adr->pos[0] = pk.pos[0] + _particle_cm.pos[0];
                        pk_adr->pos[1] = pk.pos[1] + _particle_cm.pos[1];
                        pk_adr->pos[2] = pk.pos[2] + _particle_cm.pos[2];

                        pk_adr->vel[0] = vel_sd[0] + _particle_cm.vel[0];
                        pk_adr->vel[1] = vel_sd[1] + _particle_cm.vel[1];
                        pk_adr->vel[2] = vel_sd[2] + _particle_cm.vel[2];
                    }
                }
            };
            auto& bin_root=info.getBinaryTreeRoot();
            Float vel_cm[3] = {0.0,0.0,0.0};
            Float sd_factor=1.0;
            writeBackIter(_particle_cm, vel_cm, sd_factor, bin_root);
        }
        
#endif
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
            return time_;
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
        Float getEtotRef() const {
            return etot_ref_;
        }

        //! Get current total energy from ekin and epot
        /*! \return total integrated energy 
         */
        Float getEtot() const {
            return ekin_ + epot_;
        }

        //! get energy error 
        /*! \return energy error
         */
        Float getEnergyError() const {
            return ekin_ + epot_ - etot_ref_;
        }

        //! get energy error from backup data
        Float getEnergyErrorFromBackup(Float* _bk) const {
            return -_bk[1] + _bk[2] + _bk[3];
        }

        //! get integrated energy from backup data
        Float getEtotRefFromBackup(Float* _bk) const {
            return _bk[1];
        }

        //! get total energy from backup data (ekin+epot)
        Float getEtotFromBackup(Float* _bk) const {
            return _bk[2] + _bk[3];
        }

        //! reset cumulative energy change due to slowdown change
        void resetDESlowDownChangeCum() {
            de_sd_change_cum_ = 0.0;
            dH_sd_change_cum_ = 0.0;
        }

        //! get cumulative energy change due to slowdown change
        Float getDESlowDownChangeCum() const {
            return de_sd_change_cum_;
        }

        //! get cumulative hamiltonian change due to slowdown change
        Float getDHSlowDownChangeCum() const {
            return dH_sd_change_cum_;
        }


#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
        //! Get current kinetic energy with inner slowdown
        /*! \return current kinetic energy with inner slowdown
         */
        Float getEkinSlowDownInner() const {
            return ekin_sd_;
        }

        //! Get current potential energy with inner slowdown
        /*! \return current potetnial energy with inner slowdown (negative value for bounded systems)
         */
        Float getEpotSlowDownInner() const {
            return epot_sd_;
        }

        //! Get current total integrated energy with inner slowdown
        /*! \return total integrated energy with inner slowdown
         */
        Float getEtotSlowDownInnerRef() const {
            return etot_sd_ref_;
        }

        //! Get current total energy with inner slowdown from ekin_sdi and epot_sdi
        /*! \return total energy with inner slowdown
         */
        Float getEtotSlowDownInner() const {
            return ekin_sd_ + epot_sd_;
        }

        //! get energy error with inner slowdown
        /*! \return energy error with inner slowdown
         */
        Float getEnergyErrorSlowDownInner() const {
            return ekin_sd_ + epot_sd_ - etot_sd_ref_;
        }

        //! get energy error with inner slowdown from backup data 
        Float getEnergyErrorSlowDownInnerFromBackup(Float* _bk) const {
            return -_bk[6] + _bk[7] + _bk[8];
        }

        //! get integrated energy with inner slowdown from backup data
        Float getEtotSlowDownInnerRefFromBackup(Float* _bk) const {
            return _bk[6];
        }

        //! get energy with inner slowdown from backup data (ekin_sdi + epot_sdi)
        Float getEtotSlowDownInnerFromBackup(Float* _bk) const {
            return _bk[7] + _bk[8];
        }

#endif

        //! get backup data size
        int getBackupDataSize() const {
            int bk_size = 6;
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            bk_size += 3; 
            bk_size += SlowDown::getBackupDataSize();
#endif
#ifdef AR_TTL
            bk_size += 2;
#endif
            bk_size += particles.getBackupDataSize();
            return bk_size;
        }


        //! Backup integration data 
        /*! Backup #time_, #etot_, #ekin_, $epot_, #gt_drift_, $gt_kick_inv_, #particles, $slowdown to one Float data array
          \return backup array size
        */
        int backupIntData(Float* _bk) {
            int bk_size=0;
            _bk[bk_size++] = time_;
            _bk[bk_size++] = etot_ref_;
            _bk[bk_size++] = ekin_;
            _bk[bk_size++] = epot_;
            _bk[bk_size++] = de_sd_change_cum_;
            _bk[bk_size++] = dH_sd_change_cum_;

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            _bk[bk_size++] = etot_sd_ref_;
            _bk[bk_size++] = ekin_sd_;
            _bk[bk_size++] = epot_sd_;
#endif

#ifdef AR_TTL
            _bk[bk_size++] = gt_drift_inv_;
            _bk[bk_size++] = gt_kick_inv_;
#endif

            bk_size += particles.backupParticlePosVel(&_bk[bk_size]); 
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            bk_size += info.getBinaryTreeRoot().slowdown.backup(&_bk[bk_size]); // slowdownfactor
#endif
            return bk_size;
        }

        //! Restore integration data
        /*! restore #time_, #etot_, #ekin_, $epot_, #gt_drift_, $gt_kick_inv_, #particles, $slowdown from one Float data array
          \return backup array size
        */
        int restoreIntData(Float* _bk) {
            int bk_size = 0;
            time_     = _bk[bk_size++];
            etot_ref_ = _bk[bk_size++];
            ekin_     = _bk[bk_size++];
            epot_     = _bk[bk_size++];
            de_sd_change_cum_= _bk[bk_size++];
            dH_sd_change_cum_= _bk[bk_size++];

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            etot_sd_ref_ = _bk[bk_size++];
            ekin_sd_     = _bk[bk_size++];
            epot_sd_     = _bk[bk_size++];
#endif
#ifdef AR_TTL
            gt_drift_inv_  = _bk[bk_size++];
            gt_kick_inv_   = _bk[bk_size++];
#endif
            bk_size += particles.restoreParticlePosVel(&_bk[bk_size]);
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            bk_size += info.getBinaryTreeRoot().slowdown.restore(&_bk[bk_size]);
#endif
            return bk_size;
        }

#ifdef AR_TTL
        //! Get integrated inverse time transformation factor
        /*! In TTF case, it is calculated by integrating \f$ \frac{dg}{dt} = \sum_k \frac{\partial g}{\partial \vec{r_k}} \bullet \vec{v_k} \f$.
          Notice last step is the sub-step in one symplectic loop
          \return inverse time transformation factor for drift
        */
        Float getGTDriftInv() const {
            return gt_drift_inv_;
        }
#endif

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width 
          @param[in] _n_sd: slowdown inner group
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20, const int _n_sd=0) {
            _fout<<std::setw(_width)<<"Time"
                 <<std::setw(_width)<<"dE"
                 <<std::setw(_width)<<"Etot"
                 <<std::setw(_width)<<"Ekin"
                 <<std::setw(_width)<<"Epot"
                 <<std::setw(_width)<<"Gt_drift"
                 <<std::setw(_width)<<"H"
                 <<std::setw(_width)<<"dE_SDC_cum" 
                 <<std::setw(_width)<<"dH_SDC_cum"; 
            perturber.printColumnTitle(_fout, _width);
            info.printColumnTitle(_fout, _width);
            profile.printColumnTitle(_fout, _width);
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            _fout<<std::setw(_width)<<"dE_SD" 
                 <<std::setw(_width)<<"Etot_SD" 
                 <<std::setw(_width)<<"Ekin_SD" 
                 <<std::setw(_width)<<"Epot_SD";
            _fout<<std::setw(_width)<<"N_SD";
            for (int i=0; i<_n_sd; i++) {
                _fout<<std::setw(_width)<<"I1"
                     <<std::setw(_width)<<"I2";
                SlowDown::printColumnTitle(_fout, _width);
            }
#endif
            particles.printColumnTitle(_fout, _width);
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width 
          @param[in] _n_sd: slowdown inner group
        */
        void printColumn(std::ostream & _fout, const int _width=20, const int _n_sd=0){
            _fout<<std::setw(_width)<<getTime()
                 <<std::setw(_width)<<getEnergyError()
                 <<std::setw(_width)<<etot_ref_
                 <<std::setw(_width)<<ekin_
                 <<std::setw(_width)<<epot_;
#ifdef AR_TTL
            _fout<<std::setw(_width)<<1.0/gt_drift_inv_
                 <<std::setw(_width)<<getEnergyError()/gt_drift_inv_;
#else
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            _fout<<std::setw(_width)<<1.0/manager->interaction.calcGTDriftInv(ekin_sd_-etot_sd_ref_)
                 <<std::setw(_width)<<manager->interaction.calcH(ekin_sd_-etot_sd_ref_, epot_sd_);
#else
            _fout<<std::setw(_width)<<1.0/manager->interaction.calcGTDriftInv(ekin_-etot_ref_)
                 <<std::setw(_width)<<manager->interaction.calcH(ekin_-etot_ref_, epot_);
#endif
#endif
            _fout<<std::setw(_width)<<de_sd_change_cum_
                 <<std::setw(_width)<<dH_sd_change_cum_;
            perturber.printColumn(_fout, _width);
            info.printColumn(_fout, _width);
            profile.printColumn(_fout, _width);
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            _fout<<std::setw(_width)<<getEnergyErrorSlowDownInner()
                 <<std::setw(_width)<<etot_sd_ref_ 
                 <<std::setw(_width)<<ekin_sd_ 
                 <<std::setw(_width)<<epot_sd_;
            SlowDown sd_empty;
#ifdef AR_SLOWDOWN_ARRAY
            int n_sd_now = binary_slowdown.getSize();
            _fout<<std::setw(_width)<<n_sd_now;
            for (int i=0; i<_n_sd; i++) {
                if (i<n_sd_now) {
                    _fout<<std::setw(_width)<<binary_slowdown[i]->getMemberIndex(0)
                         <<std::setw(_width)<<binary_slowdown[i]->getMemberIndex(1);
                    binary_slowdown[i]->slowdown.printColumn(_fout, _width);
                }
                else {
                    _fout<<std::setw(_width)<<-1
                         <<std::setw(_width)<<-1;
                    sd_empty.printColumn(_fout, _width);
                }
            }
#else
            int n_sd_now = info.binarytree.getSize();
            _fout<<std::setw(_width)<<n_sd_now;
            for (int i=0; i<_n_sd; i++) {
                if (i<n_sd_now) {
                    _fout<<std::setw(_width)<<info.binarytree[i].getMemberIndex(0)
                         <<std::setw(_width)<<info.binarytree[i].getMemberIndex(1);
                    info.binarytree[i].slowdown.printColumn(_fout, _width);
                }
                else {
                    _fout<<std::setw(_width)<<-1
                         <<std::setw(_width)<<-1;
                    sd_empty.printColumn(_fout, _width);
                }
            }
#endif
#endif
            particles.printColumn(_fout, _width);
        }

        //! write class data with BINARY format
        /*! @param[in] _fout: file IO for write
         */
        void writeBinary(FILE *_fout) {
            fwrite(&time_, sizeof(Float), 1, _fout);
            fwrite(&etot_ref_, sizeof(Float), 1, _fout);
            fwrite(&ekin_, sizeof(Float), 1, _fout);
            fwrite(&epot_, sizeof(Float), 1, _fout);
#ifdef AR_TTL
            fwrite(&gt_drift_inv_, sizeof(Float), 1, _fout);
#endif
            int size = force_.getSize();
            fwrite(&size, sizeof(int), 1, _fout);
            for (int i=0; i<size; i++) force_[i].writeBinary(_fout);
            
            particles.writeBinary(_fout);
            perturber.writeBinary(_fout);
            info.writeBinary(_fout);
            profile.writeBinary(_fout);
        }

        //! read class data with BINARY format and initial the array
        /*! @param[in] _fin: file IO for read
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(&time_, sizeof(Float), 1, _fin);
            rcount += fread(&etot_ref_, sizeof(Float), 1, _fin);
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
            perturber.readBinary(_fin);
            info.readBinary(_fin);
            profile.readBinary(_fin);
        }

    };
}
