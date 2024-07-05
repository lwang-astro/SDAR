
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

//! Algorithmic regularization (time transformed explicit symplectic integrator) namespace
/*!
  All major AR classes and related acceleration functions (typedef) are defined
*/
namespace AR {

    //! print features
    void printFeatures(std::ostream & fout) {
#ifdef AR_TTL
        fout<<"Use AR TTL method\n";
#else
        fout<<"Use AR LogH method\n";
#endif
#ifdef AR_SLOWDOWN_TREE
        fout<<"Use slowdown Tree method\n";
#endif
#ifdef AR_SLOWDOWN_TIMESCALE
        fout<<"Use slowdown timescale criterion\n";         
#endif
#ifdef AR_SLOWDOWN_MASSRATIO
        fout<<"Use slowdown mass ratio criterion\n";
#endif
    }

    //! print debug features
    void printDebugFeatures(std::ostream & fout) {
#ifdef AR_DEBUG
        fout<<"Debug mode: AR\n";
#endif        
    }

    //! print reference to cite
    void printReference(std::ostream & fout, const int offset=4) {
        for (int i=0; i<offset; i++) fout<<" ";
        fout<<"SDAR: Wang L., Nitadori K., Makino J., 2020, MNRAS, 493, 3398"
            <<std::endl;
    }

    //! Time Transformed Symplectic integrator manager
    /*! Tmethod is the class contain the interaction function, see sample of interaction.h:\n
     */
    template <class Tmethod>
    class TimeTransformedSymplecticManager {
    public:
        Float time_error_max; ///> maximum time error (absolute), should be positive and larger than round-off error 
        Float energy_error_relative_max; ///> maximum energy error requirement 
        Float time_step_min;        ///> minimum real time step allown
        Float ds_scale;            ///> scaling factor to determine ds
        Float slowdown_pert_ratio_ref;   ///> slowdown perturbation /inner ratio reference factor
#ifdef AR_SLOWDOWN_MASSRATIO
        Float slowdown_mass_ref;         ///> slowdown mass factor reference
#endif
        Float slowdown_timescale_max;       ///> slowdown maximum timescale to calculate maximum slowdown factor
        long long unsigned int step_count_max; ///> maximum step counts
        
        Tmethod interaction; ///> class contain interaction function
        SymplecticStep step;  ///> class to manager kick drift step

        //! constructor
        TimeTransformedSymplecticManager(): time_error_max(Float(-1.0)), energy_error_relative_max(Float(-1.0)), time_step_min(Float(-1.0)), ds_scale(1.0), slowdown_pert_ratio_ref(Float(-1.0)), 
#ifdef AR_SLOWDOWN_MASSRATIO
                                            slowdown_mass_ref(Float(-1.0)), 
#endif
                                            slowdown_timescale_max(0.0),
                                            step_count_max(0), interaction(), step() {}

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            //ASSERT(time_error_max>ROUND_OFF_ERROR_LIMIT);
            ASSERT(time_error_max>0.0);
            ASSERT(energy_error_relative_max>ROUND_OFF_ERROR_LIMIT);
            //ASSERT(time_step_min>ROUND_OFF_ERROR_LIMIT);
            ASSERT(time_step_min>0.0);
            ASSERT(ds_scale>0.0);
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
                std::cerr<<"Error: TimeTransformedSymplecticManager parameter reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
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
                 <<"step_count_max            : "<<step_count_max<<std::endl
                 <<"ds_scale                  : "<<ds_scale<<std::endl;
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
        Float de_change_interrupt_;      // energy change due to interruption
        Float dH_change_interrupt_;      // hamiltonian change due to interruption

#ifdef AR_SLOWDOWN_TREE
        Float ekin_sd_;  ///< slowdown (inner) kinetic energy
        Float epot_sd_;  ///< slowdown (inner) potential energy
        Float etot_sd_ref_;  ///< slowdown (inner) total energy

        Float de_sd_change_cum_;  // slowdown energy change
        Float dH_sd_change_cum_;  // slowdown Hamiltonian change 
        Float de_sd_change_interrupt_;   // slowdown energy change due to interruption
        Float dH_sd_change_interrupt_;   // slowdown energy change due to interruption

#endif

#ifdef AR_TTL
        // transformation factors
        Float gt_drift_inv_;  ///< integrated inverse time transformation factor for drift: dt(drift) = ds/gt_drift_inv_
#endif

        struct GtKickInv{
            Float value;  ///< value of minimum gt_kick_inv with slowdown
#ifdef AR_TIME_FUNCTION_MAX_POT
            int i; ///< index of binary member 1 having minimum gt_kick_inv
            int j; ///< index of binary member 1 having minimum gt_kick_inv
            Float gtgrad[2][3]; ///< time transformation function gradient with slowdown
            Float max; ///< max of gt_kick_inv
            Float scale; ///< scale factor to smooth gt_kick_inv change 
            bool initial;
            int inew;
            int jnew;
#elif AR_TIME_FUNCTION_MUL_POT
            int nbin; ///< number of binaries included in gt_kick_inv
            Float mul_pot_no_pow; ///< production of potential with no power
#endif

            // initialization
            GtKickInv(): 
#ifdef AR_TIME_FUNCTION_MAX_POT
                value(0.0), i(-1), j(-1), gtgrad{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, max(0.0), scale(1.0), initial(false), inew(-1), jnew(-1)
#elif AR_TIME_FUNCTION_MUL_POT
                value(1.0), nbin(0)
#else
                value(0.0)
#endif
            {}

            // reset 
            void reset() {
#ifdef AR_TIME_FUNCTION_MAX_POT
                value = 0.0;
                max = 0.0;
                initial = false;
#elif AR_TIME_FUNCTION_MUL_POT
                value = 1.0;
                nbin = 0;
#else
                value = 0.0;
#endif
            }

            // clear up
            void clear() {
#ifdef AR_TIME_FUNCTION_MAX_POT
                value = 0.0;
                i = j = -1;
                max = 0.0;
                scale = 1.0;
                initial = false;
                inew = jnew = -1;
#elif AR_TIME_FUNCTION_MUL_POT
                value = 1.0;
                nbin = 0;
#else
                value = 0.0;
#endif
            }

        } gt_kick_inv_;

        // force array
        COMM::List<Force> force_; ///< acceleration array 

    public:
        TimeTransformedSymplecticManager<Tmethod>* manager; ///< integration manager
        COMM::ParticleGroup<Tparticle,Tpcm> particles; ///< particle group manager
        Tpert    perturber; ///< perturber class 
        Tinfo    info;   ///< information of the system
        Profile  profile;  ///< profile to measure the performance
        
        //! Constructor
        TimeTransformedSymplecticIntegrator(): time_(0), etot_ref_(0), ekin_(0), epot_(0), de_change_interrupt_(0), dH_change_interrupt_(0),
#ifdef AR_SLOWDOWN_TREE
                                               ekin_sd_(0), epot_sd_(0), etot_sd_ref_(0), 
                                               de_sd_change_cum_(0), dH_sd_change_cum_(0), de_sd_change_interrupt_(0), dH_sd_change_interrupt_(0),
#endif
#ifdef AR_TTL
                                               gt_drift_inv_(0),
#endif
                                               gt_kick_inv_(), 
                                               force_(), manager(NULL), particles(), 
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
        }

        //! Clear function
        /*! Free dynamical memory space allocated
         */
        void clear() {
            time_ = 0.0;
            etot_ref_ =0.0;
            ekin_ = 0.0;
            epot_ = 0.0;
            de_change_interrupt_ = 0.0;
            dH_change_interrupt_ = 0.0;
#ifdef AR_SLOWDOWN_TREE
            ekin_sd_ = 0.0;
            epot_sd_ = 0.0;
            etot_sd_ref_ = 0.0;
            de_sd_change_cum_ = 0.0;
            dH_sd_change_cum_ = 0.0;
            de_sd_change_interrupt_ = 0.0;
            dH_sd_change_interrupt_ = 0.0;
#endif
#ifdef AR_TTL
            gt_drift_inv_ = 0.0;
#endif
            gt_kick_inv_.clear();
            force_.clear();
            particles.clear();
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
            de_change_interrupt_= _sym.de_change_interrupt_;
            dH_change_interrupt_= _sym.dH_change_interrupt_;
#ifdef AR_SLOWDOWN_TREE
            ekin_sd_= _sym.ekin_sd_;
            epot_sd_= _sym.epot_sd_;
            etot_sd_ref_= _sym.etot_sd_ref_;
            de_sd_change_cum_= _sym.de_sd_change_cum_;
            dH_sd_change_cum_= _sym.dH_sd_change_cum_;
            de_sd_change_interrupt_= _sym.de_sd_change_interrupt_;
            dH_sd_change_interrupt_= _sym.dH_sd_change_interrupt_;
#endif
#ifdef AR_TTL
            gt_drift_inv_ = _sym.gt_drift_inv_;
#endif
            gt_kick_inv_ = _sym.gt_kick_inv_;
            force_  = _sym.force_;
            manager = _sym.manager;
            particles = _sym.particles;
            info = _sym.binarytree;
            profile = _sym.profile;

            return *this;
        }

    private:

#ifdef AR_SLOWDOWN_TREE
        //! iteration function to calculate perturbation and timescale information from binary tree j to binary i
        /*!
          @param[out] _pert_out: perturbation from particle j
          @param[out] _t_min_sq: timescale limit from particle j
          @param[in] _bini: binary i
          @param[in] _binj: binary tree j
         */
        void calcSlowDownPertInnerBinaryIter(Float& _pert_out, Float& _t_min_sq, AR::BinaryTree<Tparticle>& _bini, AR::BinaryTree<Tparticle>& _binj) {
            ASSERT(&_bini != &_binj);
//            ASSERT(_bini.getMemberIndex(0)!=_binj.getMemberIndex(0));

            for (int k=0; k<2; k++) {
                if (_binj.isMemberTree(k)) {
                    auto* bink =  _binj.getMemberAsTree(k);
                    int check_flag = _bini.isSameBranch(*bink);
                    if (check_flag==0) {// no relation
                        if (bink->semi>0.0) manager->interaction.calcSlowDownPertOne(_pert_out, _t_min_sq, _bini, *bink);
                        else calcSlowDownPertInnerBinaryIter(_pert_out, _t_min_sq, _bini, *bink);
                    }
                    else if (check_flag==-2) { // _binj is the upper root tree
                        calcSlowDownPertInnerBinaryIter(_pert_out, _t_min_sq, _bini, *bink);
                    }
                    // other cases (same or sub branch), stop iteration.
                }
                else {
                    auto* pk =  _binj.getMember(k);
                    if (pk->mass>0.0) manager->interaction.calcSlowDownPertOne(_pert_out, _t_min_sq, _bini, *pk);
                }
            }
        }

        //! calculate slowdown factor for inner binary based on other particles and slowdown of system c.m.
        /*!
          @param[in] _bin: binary tree for calculating slowdown
        */
        void calcSlowDownInnerBinary(BinaryTree<Tparticle>& _bin) {
            _bin.slowdown.pert_in = manager->interaction.calcPertFromBinary(_bin);

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

            // stablility criterion
            // The slowdown factor should not make the system unstable, thus the Qst/Q set the limitation of the increasing of inner semi-major axis.
            if (_bin.stab>0 && _bin.stab != NUMERIC_FLOAT_MAX) {
                Float semi_amplify_max =  std::max(Float(1.0),1.0/_bin.stab);
                Float period_amplify_max = pow(semi_amplify_max,3.0/2.0);
                Float timescale_stab = period_amplify_max*_bin.period;
                _bin.slowdown.timescale = std::min(_bin.slowdown.timescale, timescale_stab);
            }

#else
            _bin.slowdown.timescale = _bin.slowdown.getTimescaleMax();
#endif

            // only set slowdown if semi > 0 and stable
            bool set_sd_flag = true;
            if (_bin.semi>0) {
                for (int k=0; k<2; k++) {
                    if (_bin.isMemberTree(k)) {
                        auto* bink = _bin.getMemberAsTree(k);
                        if (bink->stab>1.0) set_sd_flag = false;
                    }
                }
            }
            else set_sd_flag = false;
            
            if (set_sd_flag) {
                _bin.slowdown.period = _bin.period;
                _bin.slowdown.calcSlowDownFactor();
            }
            else {
                _bin.slowdown.setSlowDownFactor(1.0);
            }
        }


        //! Calculate twice (slowdown) kinetic energy iteration function with binary tree
        /*! cumulative ekin_ and ekin_sd_. Notice these two values should be initialized to zero and reduce by two after iteration.
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary tree for kick etot and calc dgt_drift
        */
        void calcTwoEKinIter(const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin){
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();
#ifndef USE_CM_FRAME
            Float* vel_cm = _bin.getVel();
#endif
            for (int k=0; k<2; k++) {
                Float* vk;
                Float  mk;
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    vk = bink->getVel();
                    mk = bink->mass;
                    calcTwoEKinIter(inv_nest_sd, *bink);
                }
                else {
                    auto* pk = _bin.getMember(k);
                    vk = pk->getVel();
                    mk = pk->mass;
                }
#ifdef USE_CM_FRAME
                Float* vrel = vk;
#else
                Float vrel[3] = {vk[0] - vel_cm[0],
                                 vk[1] - vel_cm[1],
                                 vk[2] - vel_cm[2]};
#endif
                ekin_ += mk * (vrel[0]*vrel[0]+vrel[1]*vrel[1]+vrel[2]*vrel[2]);
                ekin_sd_ += mk * inv_nest_sd * (vrel[0]*vrel[0] + vrel[1]*vrel[1] + vrel[2]*vrel[2]);
            }
        }

        //! Calculate (slowdown) kinetic energy
        void calcEKin() {
            ekin_ = ekin_sd_ = 0.0;
            auto& bin_root=info.getBinaryTreeRoot();
            Float sd_factor=1.0;
            
#ifdef USE_CM_FRAME
            ASSERT(!bin_root.isOriginFrame());
#endif
            calcTwoEKinIter(sd_factor, bin_root);
            Float* vcm = bin_root.getVel();
            // notice the cm velocity may not be zero after interruption, thus need to be added 
            ekin_sd_ += bin_root.mass*(vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2]);

            ekin_ *= 0.5;
            ekin_sd_ *= 0.5;
        }


        //! iteraction function to kick velocity with binary tree
        /*! First time step will be calculated, the velocities are kicked
          @param[in,out] _bin: binary tree to process
          @param[in] _dt: time step
        */
        void kickVelIter(AR::BinaryTree<Tparticle>& _bin, const Float& _dt) {
            Float dvb[3] = {0.0,0.0,0.0};
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    kickVelIter(*bink, _dt);
                    dvb[0] += bink->mass*bink->vel[0];
                    dvb[1] += bink->mass*bink->vel[1];
                    dvb[2] += bink->mass*bink->vel[2];
                }
                else {
                    auto* pk = _bin.getMember(k);
                    
                    // kick velocity
                    int pki = _bin.getMemberIndex(k);
                    Float* acc = force_[pki].acc_in;
                    Float* pert = force_[pki].acc_pert;
                    Float* vel = pk->getVel();
                    vel[0] += _dt * (acc[0] + pert[0]);
                    vel[1] += _dt * (acc[1] + pert[1]);
                    vel[2] += _dt * (acc[2] + pert[2]);

                    // calculate new c.m. velocity
                    dvb[0] += pk->mass*vel[0];
                    dvb[1] += pk->mass*vel[1];
                    dvb[2] += pk->mass*vel[2];
                }
            }
            Float mcm_inv = 1.0/_bin.mass;
#ifdef AR_DEBUG
            ASSERT(_bin.mass == _bin.getMember(0)->mass + _bin.getMember(1)->mass);
#endif
            dvb[0] *= mcm_inv;
            dvb[1] *= mcm_inv;
            dvb[2] *= mcm_inv;
#ifdef USE_CM_FRAME
            // correct c.m. and member velocity
            _bin.vel[0] += dvb[0];
            _bin.vel[1] += dvb[1];
            _bin.vel[2] += dvb[2];
            for (int k=0; k<2; k++) {
                auto* pk = _bin.getMember(k);
                pk->vel[0] -= dvb[0];
                pk->vel[1] -= dvb[1];
                pk->vel[2] -= dvb[2];
            }
#else
            // calculate new binary c.m. velocity
            _bin.vel[0] = dvb[0];
            _bin.vel[1] = dvb[1];
            _bin.vel[2] = dvb[2];
#endif
        }

        //! kick velocity
        /*! First time step will be calculated, the velocities are kicked
          @param[in] _dt: time size
        */
        void kickVel(const Float& _dt) {
            // update binary c.m. velocity interation
            auto& bin_root=info.getBinaryTreeRoot();
#ifdef USE_CM_FRAME
            ASSERT(!bin_root.isOriginFrame());
#endif
            kickVelIter(bin_root, _dt);
        }


        //! drift position with slowdown tree
        /*!
          @param[in] _dt: drift time
          @param[in] _vel_sd_up: upper cm sd vel
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary to drift pos
        */
        void driftPosTreeIter(const Float& _dt, 
#ifndef USE_CM_FRAME
                              const Float* _vel_sd_up, 
#endif
                              const Float& _inv_nest_sd_up, 
                              AR::BinaryTree<Tparticle>& _bin) { 
            // current nested sd factor
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();
#ifdef USE_CM_FRAME
            // when members are in the c.m. frame, no need to include c.m. velocity
            auto driftPos=[&](Float* pos, Float* vel) {
                pos[0] += _dt * vel[0] * inv_nest_sd;
                pos[1] += _dt * vel[1] * inv_nest_sd;
                pos[2] += _dt * vel[2] * inv_nest_sd;
            };
#else
            Float* vel_cm = _bin.getVel();
            auto driftPos=[&](Float* pos, Float* vel, Float* vel_sd) {
                //scale velocity referring to binary c.m.
                vel_sd[0] = (vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0];
                vel_sd[1] = (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1]; 
                vel_sd[2] = (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2];

                pos[0] += _dt * vel_sd[0];
                pos[1] += _dt * vel_sd[1];
                pos[2] += _dt * vel_sd[2];
            };
#endif

#if (defined USE_CM_FRAME) && (defined AR_DEBUG)
            Float dpb[3] = {0.0, 0.0, 0.0}; // c.m. position difference
            Float mcm_inv = 1.0/_bin.mass;
            ASSERT(_bin.mass == _bin.getMember(0)->mass + _bin.getMember(1)->mass);
#endif
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* pj = _bin.getMemberAsTree(k);
#ifdef USE_CM_FRAME
                    driftPos(pj->getPos(), pj->getVel());
                    driftPosTreeIter(_dt, inv_nest_sd, *pj);

#ifdef AR_DEBUG
                    dpb[0] += pj->mass*pj->pos[0];
                    dpb[1] += pj->mass*pj->pos[1];
                    dpb[2] += pj->mass*pj->pos[2];
#endif
#else
                    Float vel_sd[3];
                    driftPos(pj->getPos(), pj->getVel(), vel_sd);
                    driftPosTreeIter(_dt, vel_sd, inv_nest_sd, *pj);
#endif
                }
                else {
                    auto* pj = _bin.getMember(k);
#ifdef USE_CM_FRAME
                    driftPos(pj->getPos(), pj->getVel());

#ifdef AR_DEBUG
                    dpb[0] += pj->mass*pj->pos[0];
                    dpb[1] += pj->mass*pj->pos[1];
                    dpb[2] += pj->mass*pj->pos[2];
#endif
#else
                    Float vel_sd[3];
                    driftPos(pj->getPos(), pj->getVel(), vel_sd);
#endif
                }

            }

#if (defined USE_CM_FRAME) && (defined AR_DEBUG)
            dpb[0] *= mcm_inv;
            dpb[1] *= mcm_inv;
            dpb[2] *= mcm_inv;

            Float dpbm = sqrt(dpb[0]*dpb[0] + dpb[1]*dpb[1] + dpb[2]*dpb[2]);
            ASSERT(dpbm<ROUND_OFF_ERROR_LIMIT*10);
            /*// correct c.m. position
            _bin.pos[0] += dpb[0];
            _bin.pos[1] += dpb[1];
            _bin.pos[2] += dpb[2];
            for (int k=0; k<2; k++) {
                auto* pk = _bin.getMember(k);
                pk->pos[0] -= dpb[0];
                pk->pos[1] -= dpb[1];
                pk->pos[2] -= dpb[2];
            }
            */
#endif
        }

        //! drift time and position with slowdown tree
        /*! First (real) time is drifted, then positions are drifted
          @param[in] _dt: time step
        */
        void driftTimeAndPos(const Float& _dt) {
            // drift time 
            time_ += _dt;

            // the particle cm velocity is zero (assume in rest-frame)
            auto& bin_root=info.getBinaryTreeRoot();
            Float sd_factor=1.0;
#ifdef USE_CM_FRAME
            ASSERT(!bin_root.isOriginFrame());
            driftPosTreeIter(_dt, sd_factor, bin_root);
#else
            ASSERT(!particles.isOriginFrame());
            Float vel_cm[3] = {0.0,0.0,0.0};
            driftPosTreeIter(_dt, vel_cm, sd_factor, bin_root);
#endif
        }

        //! calc force, potential and inverse time transformation factor for one pair of particles
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _i: particle i index
          @param[in] _j: particle j index
          @param[in] _pos_offset: position offset need to be added to calculate dr
          @param[in] _calc_gt: if true, calculate gtgrad for AR_TTL mode and gt_kick_inv;
         */
        void calcAccPotAndGTKickInvTwo(const Float& _inv_nest_sd, 
                                       const int _i, const int _j, 
                                       const Float* _pos_offset, 
                                       const bool _calc_gt = true) {
            ASSERT(_i>=0&&_i<particles.getSize());
            ASSERT(_j>=0&&_j<particles.getSize());
            
            // calculate pair interaction
            Force fij[2];
            Float epotij;
            Float gt_kick_inv = manager->interaction.calcInnerAccPotAndGTKickInvTwo(
                fij[0], fij[1], epotij, particles[_i], particles[_j], _pos_offset);
            Float gt_kick_inv_sd = gt_kick_inv * _inv_nest_sd;

            // scale binary pair force with slowdown 
            force_[_i].acc_in[0] += fij[0].acc_in[0]*_inv_nest_sd;
            force_[_i].acc_in[1] += fij[0].acc_in[1]*_inv_nest_sd;
            force_[_i].acc_in[2] += fij[0].acc_in[2]*_inv_nest_sd;
            force_[_j].acc_in[0] += fij[1].acc_in[0]*_inv_nest_sd;
            force_[_j].acc_in[1] += fij[1].acc_in[1]*_inv_nest_sd;
            force_[_j].acc_in[2] += fij[1].acc_in[2]*_inv_nest_sd;

            epot_    += epotij;
            epot_sd_ += epotij*_inv_nest_sd;

            if (_calc_gt) { 
#ifdef AR_TTL
#ifdef AR_TIME_FUNCTION_MAX_POT
                // update maximum value of time transformation function gradient (gt_kick_inv) and save two paritcle indices and gtgrad values
                // if gt_kick_inv_sd > saved value, update the information
                //if (gt_kick_inv_sd > gt_kick_inv_.value) {
                if (gt_kick_inv_.i == _i && gt_kick_inv_.j == _j) {
                    Float factor = gt_kick_inv_.scale*_inv_nest_sd;
                    gt_kick_inv_.value = gt_kick_inv_.scale*gt_kick_inv_sd;
                    gt_kick_inv_.gtgrad[0][0] = fij[0].gtgrad[0]*factor;
                    gt_kick_inv_.gtgrad[0][1] = fij[0].gtgrad[1]*factor;
                    gt_kick_inv_.gtgrad[0][2] = fij[0].gtgrad[2]*factor;
                    gt_kick_inv_.gtgrad[1][0] = fij[1].gtgrad[0]*factor;
                    gt_kick_inv_.gtgrad[1][1] = fij[1].gtgrad[1]*factor;
                    gt_kick_inv_.gtgrad[1][2] = fij[1].gtgrad[2]*factor;
                }
                Float rinv = abs(epotij)/(particles[_i].mass*particles[_j].mass);
                if (rinv > gt_kick_inv_.max) {
                    gt_kick_inv_.inew = _i;
                    gt_kick_inv_.jnew = _j;
                    gt_kick_inv_.max = rinv;
                    if (gt_kick_inv_.i == -1 || gt_kick_inv_.initial) {
                        gt_kick_inv_.initial = true;
                        gt_kick_inv_.i = _i;
                        gt_kick_inv_.j = _j;
                        Float factor = gt_kick_inv_.scale*_inv_nest_sd;
                        gt_kick_inv_.value = gt_kick_inv_.scale*gt_kick_inv_sd;
                        gt_kick_inv_.gtgrad[0][0] = fij[0].gtgrad[0]*factor;
                        gt_kick_inv_.gtgrad[0][1] = fij[0].gtgrad[1]*factor;
                        gt_kick_inv_.gtgrad[0][2] = fij[0].gtgrad[2]*factor;
                        gt_kick_inv_.gtgrad[1][0] = fij[1].gtgrad[0]*factor;
                        gt_kick_inv_.gtgrad[1][1] = fij[1].gtgrad[1]*factor;
                        gt_kick_inv_.gtgrad[1][2] = fij[1].gtgrad[2]*factor;
                    }
                }

#elif AR_TIME_FUNCTION_MUL_POT
                // Here gtgrad excludes gt_kick_inv and slowdown factor, these factors will be multipled when gt_drift_inv is calculated.
                force_[_i].gtgrad[0] = fij[0].gtgrad[0];
                force_[_i].gtgrad[1] = fij[0].gtgrad[1];
                force_[_i].gtgrad[2] = fij[0].gtgrad[2];
                force_[_j].gtgrad[0] = fij[1].gtgrad[0];
                force_[_j].gtgrad[1] = fij[1].gtgrad[1];
                force_[_j].gtgrad[2] = fij[1].gtgrad[2];

                // add binary count and multiply gt_kick_inv_sd
                gt_kick_inv_.nbin++;
                gt_kick_inv_.value *= gt_kick_inv_sd;
#else
                // scale gtgrad with slowdown
                force_[_i].gtgrad[0] += fij[0].gtgrad[0]*_inv_nest_sd;
                force_[_i].gtgrad[1] += fij[0].gtgrad[1]*_inv_nest_sd;
                force_[_i].gtgrad[2] += fij[0].gtgrad[2]*_inv_nest_sd;
                force_[_j].gtgrad[0] += fij[1].gtgrad[0]*_inv_nest_sd;
                force_[_j].gtgrad[1] += fij[1].gtgrad[1]*_inv_nest_sd;
                force_[_j].gtgrad[2] += fij[1].gtgrad[2]*_inv_nest_sd;

                gt_kick_inv_.value += gt_kick_inv_sd;
#endif
#else // NO AR_TTL
                gt_kick_inv_.value += gt_kick_inv_sd;
#endif
            }
        }

        //! calc force, potential and inverse time transformation factor for one particle by walking binary tree
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _i: particle index
          @param[in] _bin: binary tree for walking
          @param[in] _pos_offset: position offset need to be added to calculate dr
          @param[in] _calc_gt: if true, calculate gtgrad and gt_kick_inv for AR_TTL mode, not used for LogH mode
         */
        void calcAccPotAndGTKickInvOneTreeIter(const Float& _inv_nest_sd, 
                                               const int _i, 
                                               AR::BinaryTree<Tparticle>& _bin, 
                                               const Float* _pos_offset,
                                               const bool _calc_gt) {
#ifdef USE_CM_FRAME
            Float pos_offset[3] = {
                _pos_offset[0] + _bin.pos[0],
                _pos_offset[1] + _bin.pos[1],
                _pos_offset[2] + _bin.pos[2]};
#else
            const Float* pos_offset = _pos_offset;
#endif
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) // particle - tree
                    calcAccPotAndGTKickInvOneTreeIter(_inv_nest_sd, _i, *(_bin.getMemberAsTree(k)), pos_offset, _calc_gt);
                else  // particle - particle
                    calcAccPotAndGTKickInvTwo(_inv_nest_sd, _i, _bin.getMemberIndex(k), pos_offset, _calc_gt);
            }
        }

        //! calc crossing force, potential and inverse time transformation factor between binary tree i and binary tree j
        /*!
          @param[in] _inv_nest_sd: inverse nested slowdown factor
          @param[in] _bini: binary tree i for walking
          @param[in] _binj: binary tree j for walking
          @param[in] _pos_offset: position offset need to be added to calculate dr
          @param[in] _calc_gt: if true, calculate gtgrad and gt_kick_inv for AR_TTL mode, not used for LogH mode
         */
        void calcAccPotAndGTKickInvCrossTreeIter(const Float& _inv_nest_sd, 
                                                  AR::BinaryTree<Tparticle>& _bini, 
                                                  AR::BinaryTree<Tparticle>& _binj, 
                                                  const Float* _pos_offset,
                                                  const bool _calc_gt) {
            ASSERT(&_bini!=&_binj);
#ifdef USE_CM_FRAME
            Float pos_offset[3] = {
                _pos_offset[0] - _bini.pos[0],
                _pos_offset[1] - _bini.pos[1],
                _pos_offset[2] - _bini.pos[2]};
#else
            const Float* pos_offset = _pos_offset;
#endif
            for (int k=0; k<2; k++) { 
                if (_bini.isMemberTree(k)) { // tree - tree
                    calcAccPotAndGTKickInvCrossTreeIter(_inv_nest_sd, *(_bini.getMemberAsTree(k)), _binj, pos_offset, _calc_gt);
                }
                else  // particle - tree
                    calcAccPotAndGTKickInvOneTreeIter(_inv_nest_sd, _bini.getMemberIndex(k), _binj, pos_offset, _calc_gt);
            }
        }

        //! calculate force, potential, gtgrad (for AR_TTL mode) and inverse time transformation factor for kick (gt_kick_inv_)
        /*!
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary to drift pos
         */
        void calcAccPotAndGTKickInvTreeIter(const Float& _inv_nest_sd_up, AR::BinaryTree<Tparticle>& _bin) {
            // current nested sd factor
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();
#ifdef USE_CM_FRAME
            Float pos_offset[3] = {0.0, 0.0, 0.0};
#else
            Float* pos_offset = NULL;
#endif

#if (defined AR_TIME_FUNCTION_MUL_POT) || (defined AR_TIME_FUNCTION_ADD_POT) || (defined AR_TIME_FUNCTION_MAX_POT)
            bool calc_gt_cross = false;
#else
            bool calc_gt_cross = true;
#endif

            // check left 
            if (_bin.isMemberTree(0)) { // left is tree
                auto* bin_left = _bin.getMemberAsTree(0);
                // inner interaction of left tree
                calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_left);

                if (_bin.isMemberTree(1)) { // right is tree
                    auto* bin_right = _bin.getMemberAsTree(1);

                    // inner interaction of right tree
                    calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_right);

                    // cross interaction
                    calcAccPotAndGTKickInvCrossTreeIter(inv_nest_sd, *bin_left, *bin_right, pos_offset, calc_gt_cross);
                }
                else { // right is particle
                    // cross interaction from particle j to tree left
                    calcAccPotAndGTKickInvOneTreeIter(inv_nest_sd, _bin.getMemberIndex(1), *bin_left, pos_offset, calc_gt_cross);
                }
            }
            else { // left is particle
                if (_bin.isMemberTree(1)) { // right is tree
                    auto* bin_right = _bin.getMemberAsTree(1);
                    // inner interaction of right tree
                    calcAccPotAndGTKickInvTreeIter(inv_nest_sd, *bin_right);

                    // cross interaction from particle i to tree right
                    calcAccPotAndGTKickInvOneTreeIter(inv_nest_sd, _bin.getMemberIndex(0), *bin_right, pos_offset, calc_gt_cross);
                }
                else { // right is particle
                    // particle - particle interaction
                    calcAccPotAndGTKickInvTwo(inv_nest_sd, _bin.getMemberIndex(0), _bin.getMemberIndex(1), pos_offset, true);
                }
            }
        }

        //! calc force, potential and inverse time transformation factor for kick
        inline void calcAccPotAndGTKickInv() {
            epot_ = 0.0;
            epot_sd_ = 0.0;
            for (int i=0; i<force_.getSize(); i++) force_[i].clear();

            gt_kick_inv_.reset();
#ifdef USE_CM_FRAME
            ASSERT(!info.getBinaryTreeRoot().isOriginFrame());
#endif
            calcAccPotAndGTKickInvTreeIter(1.0, info.getBinaryTreeRoot());
//#ifdef AR_TIME_FUNCTION_MUL_POT
//            // use power in gt_kick_inv_
//            gt_kick_inv_.mul_pot_no_pow = gt_kick_inv_.value;
//            gt_kick_inv_.value = pow(gt_kick_inv_.mul_pot_no_pow, 1.0/gt_kick_inv_.nbin);
//#endif

            // pertuber force
            manager->interaction.calcAccPert(force_.getDataAddress(), particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());
        }

#ifdef AR_TTL
        //! kick energy and time transformation function for drift of binary tree 
        /*!
          @param[in] _dt: time step
          @param[in] _vel_up: upper cm vel
          @param[in] _vel_sd_up: upper cm sd vel
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary tree for kick etot and calc dgt_drift
          \return gt_drift_inv change 
        */
        Float kickEtotAndGTDriftTreeIter(const Float& _dt, 
#ifdef USE_CM_FRAME
                                         const Float* _vel_up,
#endif
                                         const Float* _vel_sd_up, 
                                         const Float& _inv_nest_sd_up, 
                                         AR::BinaryTree<Tparticle>& _bin) {
            // current nested sd factor
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();
            Float dgt_drift_inv = 0.0;
            Float de = 0.0;

#ifndef USE_CM_FRAME
            Float* vel_cm = _bin.getVel();
#endif
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    // get no sd velocity in original frame
#ifdef USE_CM_FRAME
                    Float* vel_rel = bink->getVel();
                    Float vel[3] = {vel_rel[0] + _vel_up[0], 
                                    vel_rel[1] + _vel_up[1], 
                                    vel_rel[2] + _vel_up[2]};
#else
                    Float* vel = bink->getVel();
                    Float vel_rel[3] = {vel[0] - vel_cm[0], 
                                        vel[1] - vel_cm[1], 
                                        vel[2] - vel_cm[2]}; 
#endif
                    // sd velocity in original frame
                    Float vel_sd[3] = {vel_rel[0] * inv_nest_sd + _vel_sd_up[0], 
                                       vel_rel[1] * inv_nest_sd + _vel_sd_up[1], 
                                       vel_rel[2] * inv_nest_sd + _vel_sd_up[2]}; 

                    dgt_drift_inv += kickEtotAndGTDriftTreeIter(_dt, 
#ifdef USE_CM_FRAME
                                                                vel,
#endif
                                                                vel_sd, 
                                                                inv_nest_sd, 
                                                                *bink);
                }
                else {
                    int i = _bin.getMemberIndex(k);
                    ASSERT(i>=0&&i<particles.getSize());
                    ASSERT(&particles[i]==_bin.getMember(k));
                    
#ifdef USE_CM_FRAME
                    Float* vel_rel = particles[i].getVel();
                    Float vel[3] = {vel_rel[0] + _vel_up[0], 
                                    vel_rel[1] + _vel_up[1], 
                                    vel_rel[2] + _vel_up[2]};
#else
                    Float* vel = particles[i].getVel();
                    Float vel_rel[3] = {vel[0] - vel_cm[0], 
                                        vel[1] - vel_cm[1], 
                                        vel[2] - vel_cm[2]}; 
#endif
                    Float vel_sd[3] = {vel_rel[0] * inv_nest_sd + _vel_sd_up[0], 
                                       vel_rel[1] * inv_nest_sd + _vel_sd_up[1], 
                                       vel_rel[2] * inv_nest_sd + _vel_sd_up[2]}; 

#ifndef AR_TIME_FUNCTION_MAX_POT
                    Float* gtgrad = force_[i].gtgrad;
#else
                    // use recored gtgrad in gt_kick_inv_max_info instead of force_[i].gtgrad (not calculated)
                    Float* gtgrad = NULL;
                    if (i == gt_kick_inv_.i) 
                        gtgrad = gt_kick_inv_.gtgrad[0];
                    else if (i == gt_kick_inv_.j)
                        gtgrad = gt_kick_inv_.gtgrad[1];
                    if (gtgrad != NULL)
#endif
                        dgt_drift_inv += (vel_sd[0] * gtgrad[0] + 
                                          vel_sd[1] * gtgrad[1] +
                                          vel_sd[2] * gtgrad[2]);

                    Float* pert   = force_[i].acc_pert;

                    de += particles[i].mass * (vel[0] * pert[0] +
                                               vel[1] * pert[1] +
                                               vel[2] * pert[2]);
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
            auto& bin_root=info.getBinaryTreeRoot();
            Float vel_cm[3] = {0.0,0.0,0.0};
            Float sd_factor=1.0;
#ifdef USE_CM_FRAME
            ASSERT(!bin_root.isOriginFrame());
            Float dgt_drift_inv = kickEtotAndGTDriftTreeIter(_dt, vel_cm, vel_cm, sd_factor, bin_root); 
#else
            ASSERT(!particles.isOriginFrame());
            Float dgt_drift_inv = kickEtotAndGTDriftTreeIter(_dt, vel_cm, sd_factor, bin_root);
#endif
#ifdef AR_TIME_FUNCTION_MUL_POT
            //dgt_drift_inv *= gt_kick_inv_.value/gt_kick_inv_.nbin;
            dgt_drift_inv *= gt_kick_inv_.value;
#endif
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

        //! Calculate kinetic energy
        inline void calcEKin(){
            ekin_ = Float(0.0);
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            for (int i=0; i<num; i++) {
                const Float *vi=pdat[i].getVel();
                ekin_ += 0.5 * pdat[i].mass * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
            }
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
            for (int i=0; i<num; i++) {
                Float* pos = pdat[i].getPos();
                Float* vel = pdat[i].getVel();
                pos[0] += _dt * vel[0];
                pos[1] += _dt * vel[1];
                pos[2] += _dt * vel[2];
            }
        }

        //! calc force, potential and inverse time transformation factor for kick
        inline void calcAccPotAndGTKickInv() {

#ifdef USE_CM_FRAME
            Float pos_offset[3] = {0.0, 0.0 ,0.0};
#else
            Float* pos_offset = NULL;
#endif
            if (particles.getSize()==2) 
                gt_kick_inv_.value = manager->interaction.calcInnerAccPotAndGTKickInvTwo(force_[0], force_[1], epot_, particles[0], particles[1], pos_offset);
            else 
                gt_kick_inv_.value = manager->interaction.calcInnerAccPotAndGTKickInv(force_.getDataAddress(), epot_, particles.getDataAddress(), particles.getSize());

            // pertuber force
            manager->interaction.calcAccPert(force_.getDataAddress(), particles.getDataAddress(), particles.getSize(), particles.cm, perturber, getTime());

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
            gt_drift_inv_ += _dt * dg;
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
            etot_ref_ += _dt * de;
        }
#endif //AR_TTL

#endif // AR_SLOWDOWN_TREE

        //! set all binary c.m. mass to zero
        void setBinaryCMZeroIter(AR::BinaryTree<Tparticle>& _bin) {
            _bin.mass = 0.0;
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    setBinaryCMZeroIter(*bink);
                }
            }
        }

        //! update binary semi, ecc and period iteratively for unstable multiple systems
        bool updateBinarySemiEccPeriodIter(AR::BinaryTree<Tparticle>& _bin, const Float& _G, const Float _time, const bool _check=false) {
            bool check = _check;
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    check=updateBinarySemiEccPeriodIter(*bink, _G, _time, _check);
                }
            }
            if ((_time>_bin.stab_check_time&&_bin.stab>1.0&&_bin.m1>0.0&&_bin.m2>0.0)||check) {
                _bin.calcSemiEccPeriod(_G);
                _bin.stab_check_time = _time + _bin.period;
                return true;
            }
            return false;
        }
            
    public:
#ifdef AR_SLOWDOWN_TREE
        //! update slowdown factor based on perturbation and record slowdown energy change
        /*! Update slowdown inner and global.
            @param [in] _update_energy_flag: Record cumulative slowdown energy change if true;
            @param [in] _stable_check_flag: check whether the binary tree is stable if true;
         */
        void updateSlowDownAndCorrectEnergy(const bool _update_energy_flag, const bool _stable_check_flag) {
            auto& bin_root = info.getBinaryTreeRoot();
            auto& sd_root = bin_root.slowdown;

#ifdef AR_TTL
            Float sd_backup = sd_root.getSlowDownFactor();
#endif

            // when the maximum inner slowdown is large, the outer should not be slowed down since the system may not be stable.
            //if (inner_sd_change_flag&&sd_org_inner_max<1000.0*manager->slowdown_pert_ratio_ref) sd_root.setSlowDownFactor(1.0);
            //if (time_>=sd_root.getUpdateTime()) {
            sd_root.pert_in = manager->interaction.calcPertFromBinary(bin_root);
            sd_root.pert_out = 0.0;
            Float t_min_sq= NUMERIC_FLOAT_MAX;
            manager->interaction.calcSlowDownPert(sd_root.pert_out, t_min_sq, getTime(), particles.cm, perturber);
            sd_root.timescale = std::min(sd_root.getTimescaleMax(), sqrt(t_min_sq));

            //Float period_amplify_max = NUMERIC_FLOAT_MAX;
            if (_stable_check_flag) {
                // check whether the system is stable for 10000 out period and the apo-center is below break criterion
                Float stab = bin_root.stableCheckIter(bin_root,10000*bin_root.period);
                Float apo = bin_root.semi*(1+bin_root.ecc);
                if (stab<1.0 && apo<info.r_break_crit) {
                    sd_root.period = bin_root.period;
                    sd_root.calcSlowDownFactor();
                }
                else sd_root.setSlowDownFactor(1.0);

                // stablility criterion
                // The slowdown factor should not make the system unstable, thus the Qst/Q set the limitation of the increasing of inner semi-major axis.
                //if (stab>0 && stab != NUMERIC_FLOAT_MAX) {
                //    Float semi_amplify_max =  std::max(Float(1.0),1.0/stab);
                //    period_amplify_max = pow(semi_amplify_max,3.0/2.0);
                //}
            }
            else if (bin_root.semi>0) {
                sd_root.period = bin_root.period;
                sd_root.calcSlowDownFactor();
            }
            else sd_root.setSlowDownFactor(1.0);

            //sd_root.increaseUpdateTimeOnePeriod();
            //}

            // inner binary slowdown
            Float sd_org_inner_max = 0.0;
            bool inner_sd_change_flag=false;
            int n_bin = info.binarytree.getSize();
            for (int i=0; i<n_bin-1; i++) {
                auto& bini = info.binarytree[i];
                //if (time_>=bini.slowdown.getUpdateTime()) {
#ifndef USE_CM_FRAME
                // this is already updated when USE_CM_FRAME is used, should not calculate twice
                bini.calcCenterOfMass();
#endif
                calcSlowDownInnerBinary(bini);

                //sdi->slowdown.increaseUpdateTimeOnePeriod();
                sd_org_inner_max = std::max(bini.slowdown.getSlowDownFactorOrigin(),sd_org_inner_max);
                inner_sd_change_flag=true;
                //}
            }


            if (_update_energy_flag) {
                Float ekin_sd_bk = ekin_sd_;
                Float epot_sd_bk = epot_sd_;
                Float H_sd_bk = getHSlowDown();
                if(inner_sd_change_flag) {
#ifdef AR_TTL
                    Float gt_kick_inv_bk = gt_kick_inv_.value;
                    calcAccPotAndGTKickInv();
                    gt_drift_inv_ += gt_kick_inv_.value - gt_kick_inv_bk;
#else                    
                    calcAccPotAndGTKickInv();
#endif
                    calcEKin();
                }
                else {
                    Float kappa_inv = 1.0/sd_root.getSlowDownFactor();
#ifdef AR_TTL
                    Float gt_kick_inv_new = gt_kick_inv_.value*sd_backup*kappa_inv;
                    gt_drift_inv_ += gt_kick_inv_new - gt_kick_inv_.value;
                    gt_kick_inv_.value = gt_kick_inv_new;
#endif
                    ekin_sd_ = ekin_*kappa_inv;
                    epot_sd_ = epot_*kappa_inv;
                }
                Float de_sd = (ekin_sd_ - ekin_sd_bk) + (epot_sd_ - epot_sd_bk);
                etot_sd_ref_ += de_sd;

                Float dH_sd = getHSlowDown() - H_sd_bk;

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
#ifndef USE_CM_FRAME
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
#endif

#ifdef AR_SLOWDOWN_TREE
//#ifdef USE_CM_FRAME
//            auto& bin_root = info.getBinaryTreeRoot();
//            if (bin_root.isOriginFrame()) bin_root.shiftToCenterOfMassFrame();
//#endif
            for (int i=0; i<info.binarytree.getSize(); i++) 
                info.binarytree[i].slowdown.initialSlowDownReference(manager->slowdown_pert_ratio_ref, manager->slowdown_timescale_max);

            updateSlowDownAndCorrectEnergy(false,true);

#endif // END AR_SLOWDOWN_TREE

#ifdef AR_TTL
            calcAccPotAndGTKickInv();

            // initially gt_drift 
            gt_drift_inv_ = gt_kick_inv_.value;

#else
            calcAccPotAndGTKickInv();
#endif

            // calculate kinetic energy
            calcEKin();

            etot_ref_ = ekin_ + epot_;

#ifdef AR_SLOWDOWN_TREE
            etot_sd_ref_ = ekin_sd_ + epot_sd_;

            Float de_sd = etot_sd_ref_ - etot_ref_;

            // add slowdown change to the global slowdown energy
            de_sd_change_cum_ += de_sd;
            dH_sd_change_cum_ = 0.0;
#endif
        }

#ifdef AR_MULTI_STEP
        //! integration with KDKDK
        void integrateMultiStepChin4th(const Float _ds, 
#ifndef USE_CM_FRAME
                                       const Float* _vel_up,
                                       const Float* _vel_sd_up, 
#endif
                                       const Float _inv_nest_sd_up,
                                       AR::BinaryTree<Tparticle>& _bin) {
            ASSERT(checkParams());
            ASSERT(!particles.isModified());
            ASSERT(_ds>0);
#ifdef USE_CM_FRAME
            ASSERT(!_bin.isOriginFrame());
#else
            ASSERT(!particles.isOriginFrame());
#endif

            // Kick 1/6
            calcAccPot(_bin);

            Float dt_kick = _ds/6.0*gt_kick_inv_.value;
            
            kickVelIter(_bin, 0.5*dt_kick);
#ifdef AR_TTL   
            // kick total energy and inverse time transformation factor for drift
            kickEtotAndGTDriftTreeIter(dt_kick, _vel_up, _vel_sd_up, _inv_nest_sd_up, _bin);
#else
            // kick total energy
            kickEtotIter(dt_kick);
#endif
            kickVelIter(_bin, 0.5*dt_kick);
            
            // Drift 1/2
            Float dt_drift = 0.5*_ds*gt_drift_inv_;
#ifdef USE_CM_FRAME
            driftPosTreeIter(_dt, _inv_nest_sd_up, _bin);
#else
            driftPosTreeIter(_dt, _vel_sd_up, _inv_nest_sd_up, _bin);
#endif

            // Kick 4/6 with Grad
            calcAccPot();
            calcGrad();
            kickVelIter(_bin, 4.0*dt_kick);

            // Drift 1/2
            Float dt_drift = 0.5*_ds*gt_drift_inv_;
#ifdef USE_CM_FRAME
            driftPosTreeIter(_dt, _inv_nest_sd_up, _bin);
#else
            driftPosTreeIter(_dt, _vel_sd_up, _inv_nest_sd_up, _bin);
#endif
        }

        //! integration for one step
        /*!
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().
        */
        void integrateOneStepAR(const Float _ds, #ifndef USE_CM_FRAME
#ifndef USE_CM_FRAME
                                const Float* _vel_up,
                                const Float* _vel_sd_up, 
#endif
                                const Float _inv_nest_sd_up,
                                AR::BinaryTree<Tparticle>& _bin,
                                Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

#ifdef AR_TIME_FUNCTION_MAX_POT
            if (gt_kick_inv_.inew != gt_kick_inv_.i || gt_kick_inv_.jnew != gt_kick_inv_.jnew) {
                // update i and j for calculate gt_kick_inv
                gt_kick_inv_.i = gt_kick_inv_.inew;
                gt_kick_inv_.j = gt_kick_inv_.jnew;

                calcAccPotAndGTKickInv();

                // initially gt_drift
                gt_kick_inv_.scale = gt_drift_inv_/gt_kick_inv_.value;
                //gt_kick_inv_ = gt_drift_inv_;
                //gt_drift_inv_ = gt_kick_inv_.value;
            }
#endif

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // inverse time transformation factor for drift
#ifdef AR_TTL
                Float gt_drift_inv = gt_drift_inv_;
#else 
#ifdef AR_SLOWDOWN_TREE
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
                calcAccPotAndGTKickInv();
                Float gt_kick_inv = gt_kick_inv_.value;

                // time step for kick
                Float dt_kick = ds_kick/gt_kick_inv;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

#ifdef AR_TTL   
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

#endif // END AR_MULTI_STEP

        //! integration for one step
        /*!
          @param[in] _ds: step size
          @param[out] _time_table: for high order symplectic integration, store the substep integrated (real) time, used for estimate the step for time synchronization, size should be consistent with step.getCDPairSize().
        */
        void integrateOneStep(const Float _ds, Float _time_table[]) {
            ASSERT(checkParams());

            ASSERT(!particles.isModified());
            ASSERT(_ds>0);

#ifdef AR_TIME_FUNCTION_MAX_POT
            if (gt_kick_inv_.inew != gt_kick_inv_.i || gt_kick_inv_.jnew != gt_kick_inv_.jnew) {
                // update i and j for calculate gt_kick_inv
                gt_kick_inv_.i = gt_kick_inv_.inew;
                gt_kick_inv_.j = gt_kick_inv_.jnew;

                calcAccPotAndGTKickInv();

                // initially gt_drift
                gt_kick_inv_.scale = gt_drift_inv_/gt_kick_inv_.value;
                //gt_kick_inv_ = gt_drift_inv_;
                //gt_drift_inv_ = gt_kick_inv_.value;
            }
#endif

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // inverse time transformation factor for drift
#ifdef AR_TTL
                Float gt_drift_inv = gt_drift_inv_;
#else 
#ifdef AR_SLOWDOWN_TREE
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
                calcAccPotAndGTKickInv();
                Float gt_kick_inv = gt_kick_inv_.value;

                // time step for kick
                Float dt_kick = ds_kick/gt_kick_inv;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

#ifdef AR_TTL   
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

#ifdef AR_SLOWDOWN_TREE
            const Float kappa_inv = 1.0/info.getBinaryTreeRoot().slowdown.getSlowDownFactor();
#endif

#ifdef USE_CM_FRAME
            const Float pos_offset[3] = {0.0, 0.0, 0.0};
#else
            Float* pos_offset = NULL;
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

#ifdef AR_DEBUG_PRINT_DKD
            std::cout<<"K "<<time_<<" "
                     <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                     <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                     <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds = manager->step.getCK(i)*_ds;
                // inverse time transformation factor for drift
#ifdef AR_TTL
                Float gt_inv = gt_drift_inv_;
#else 
#ifdef AR_SLOWDOWN_TREE
                Float gt_inv = manager->interaction.calcGTDriftInv(ekin_sd_-etot_sd_ref_); // pt = -etot_sd
#else
                Float gt_inv = manager->interaction.calcGTDriftInv(ekin_-etot_ref_); // pt = -etot
#endif
#endif
                // drift
                Float dt = ds/gt_inv;
                ASSERT(!ISNAN(dt));
                
                // drift time 
                time_ += dt;

                // update real time
                _time_table[i] = time_;

#ifdef AR_SLOWDOWN_TREE
                Float dt_sd = dt*kappa_inv;

                // drift position
                pos1[0] += dt_sd * vel1[0];
                pos1[1] += dt_sd * vel1[1];
                pos1[2] += dt_sd * vel1[2];

                pos2[0] += dt_sd * vel2[0];
                pos2[1] += dt_sd * vel2[1];
                pos2[2] += dt_sd * vel2[2];
#else
                // drift position
                pos1[0] += dt * vel1[0];
                pos1[1] += dt * vel1[1];
                pos1[2] += dt * vel1[2];

                pos2[0] += dt * vel2[0];
                pos2[1] += dt * vel2[1];
                pos2[2] += dt * vel2[2];
#endif

                // step for kick
                ds = manager->step.getDK(i)*_ds;

                gt_inv = manager->interaction.calcInnerAccPotAndGTKickInvTwo(force_data[0], force_data[1], epot_, particle_data[0], particle_data[1], pos_offset);

                // pertuber force
                manager->interaction.calcAccPert(force_data, particle_data, n_particle, particles.cm, perturber, _time_table[i]);

                ASSERT(!ISNAN(epot_));

#ifdef AR_DEBUG_PRINT_DKD
                if (i>0)
                    std::cout<<"K "<<time_<<" "
                             <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                             <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                             <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick half step for velocity
                Float dvel1[3], dvel2[3];

#ifdef AR_SLOWDOWN_TREE
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

#ifdef AR_DEBUG_PRINT_DKD
                std::cout<<"D "<<time_<<" "
                         <<pos2[0]-pos1[0]<<" "<<pos2[1]-pos1[1]<<" "<<pos2[2]-pos1[2]<<" "
                         <<vel2[0]-vel1[0]<<" "<<vel2[1]-vel1[1]<<" "<<vel2[2]-vel1[2]<<" "
                         <<ekin_<<" "<<epot_<<" "<<etot_ref_<<std::endl;
#endif

                // kick total energy and time transformation factor for drift
                etot_ref_ += 2.0*dt * (mass1* (vel1[0] * pert1[0] + 
                                               vel1[1] * pert1[1] + 
                                               vel1[2] * pert1[2]) +
                                       mass2* (vel2[0] * pert2[0] + 
                                               vel2[1] * pert2[1] + 
                                               vel2[2] * pert2[2]));

#ifdef AR_TTL   
                // back up gt_kick_inv
                gt_kick_inv_.value = gt_inv;

#ifdef AR_SLOWDOWN_TREE
                // integrate gt_drift_inv
                Float dgt_drift_inv = 2.0*dt*kappa_inv*kappa_inv* (vel1[0] * gtgrad1[0] +
                                                                   vel1[1] * gtgrad1[1] +
                                                                   vel1[2] * gtgrad1[2] +
                                                                   vel2[0] * gtgrad2[0] +
                                                                   vel2[1] * gtgrad2[1] +
                                                                   vel2[2] * gtgrad2[2]);
#ifdef AR_TIME_FUNCTION_MUL_POT
                dgt_drift_inv *= gt_kick_inv_.value;
#endif
                gt_drift_inv_ += dgt_drift_inv;

#else // NO Slowdown
                // integrate gt_drift_inv
                gt_drift_inv_ +=  2.0*dt* (vel1[0] * gtgrad1[0] +
                                           vel1[1] * gtgrad1[1] +
                                           vel1[2] * gtgrad1[2] +
                                           vel2[0] * gtgrad2[0] +
                                           vel2[1] * gtgrad2[1] +
                                           vel2[2] * gtgrad2[2]);
#endif // END SLOWDOWN

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

#ifdef AR_SLOWDOWN_TREE
            // make consistent slowdown inner energy 
            etot_sd_ref_ = etot_ref_*kappa_inv;
            ekin_sd_ = ekin_*kappa_inv;
            epot_sd_ = epot_*kappa_inv;
#endif
        }
        
        // Integrate the system to a given time
        /*!
          @param[in] _time_end: the expected finishing time without offset
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
            
            Float backup_data[bk_data_size]; // for backup particle data
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

            // reduce ds control, the maximum level
            const int n_reduce_level_max=10;

            // ds manager
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
                        n_step_wait_recover_ds[n_reduce_level] += 2*to_int(1.0/_modify_factor);
                    else {
                        n_reduce_level++;
                        ds_backup[n_reduce_level] = _ds;
                        n_step_wait_recover_ds[n_reduce_level] = 2*to_int(1.0/_modify_factor);
                    }
                }

                // count step (return false) and recover ds if necessary (return true)
                bool countAndRecover(Float &_ds, Float &_modify_factor, const bool _recover_flag) {
                    if (n_reduce_level>=0) {
                        if (n_step_wait_recover_ds[n_reduce_level] ==0) {
                            if (_recover_flag) {
                                _modify_factor = ds_backup[n_reduce_level]/_ds;
                                _ds = ds_backup[n_reduce_level];
                                n_step_wait_recover_ds[n_reduce_level] = -1;
                                n_reduce_level--;
                                return true;
                            }
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
            long long unsigned int step_count=0; // integration step 
            long long unsigned int step_count_tsyn=0; // time synchronization step

            // interrupted binary recorder
            InterruptBinary<Tparticle> bin_interrupt;
            bin_interrupt.time_now=time_ + info.time_offset;
            bin_interrupt.time_end=_time_end + info.time_offset;
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

            // warning print flag
            bool warning_print_once=true;

#ifdef AR_DEBUG_DUMP
            // back up initial data
            backupIntData(backup_data_init);
#endif
      
#ifdef AR_SLOWDOWN_TREE
            // update slowdown and correct slowdown energy and gt_inv
            updateSlowDownAndCorrectEnergy(true, true);
#endif

            // reset binary stab_check_time
            for (int i=0; i<info.binarytree.getSize(); i++)
                info.binarytree[i].stab_check_time = time_;

            // integration loop
            while(true) {
                // backup data
                bool binary_update_flag=false;
                auto& bin_root = info.getBinaryTreeRoot();
                auto& G = manager->interaction.gravitational_constant;
                
                if(backup_flag) {
                    // check interrupt condiction, ensure that time end not reach
                    if (manager->interaction.interrupt_detection_option>0 && !time_end_flag) {
                        bin_interrupt.time_now = time_ + info.time_offset;
                        bin_interrupt.time_end = _time_end + info.time_offset;
                        // calc perturbation energy
                        //Float epert=0.0;
                        //for (int i=0; i<n_particle; i++) {
                        //    epert += force_[i].pot_pert*particles[i].mass;
                        //}
                        manager->interaction.modifyAndInterruptIter(bin_interrupt, bin_root);
                        //InterruptBinary<Tparticle>* bin_intr_ptr = &bin_interrupt;
                        //bin_intr_ptr = bin_root.processRootIter(bin_intr_ptr, Tmethod::modifyAndInterruptIter);
                        ASSERT(bin_interrupt.checkParams());
                        if (bin_interrupt.status!=InterruptStatus::none) {
                            if (manager->interaction.interrupt_detection_option==2) {

                                // return the first one of detection
                                if (bin_interrupt_return.status == InterruptStatus::none) 
                                    bin_interrupt_return = bin_interrupt;
                            }
                            else {
                                // check whether destroy appears (all masses becomes zero)
                                if (bin_interrupt.status==InterruptStatus::destroy) {
                                    // all particles become zero masses
#ifdef AR_DEBUG
                                    for (int j=0; j<n_particle; j++) {
                                        ASSERT(particles[j].mass==0.0);
                                    }
#endif
                                    de_change_interrupt_ -= etot_ref_;
                                    dH_change_interrupt_ -= getH();
                                    ekin_ = epot_ = etot_ref_ = 0.0;
#ifdef AR_SLOWDOWN_TREE
                                    de_sd_change_cum_ -= etot_sd_ref_;
                                    dH_sd_change_interrupt_ -= getHSlowDown();
                                    ekin_sd_ = epot_sd_ = etot_sd_ref_ = 0.0;
#endif                                    
#ifdef AR_DEBUG_PRINT
                                    std::cerr<<"Interrupt condition triggered! Destroy";
                                    std::cerr<<" Time: "<<time_;
                                    auto bin_adr = bin_interrupt.getBinaryTreeAddress();
                                    bin_adr->printColumnTitle(std::cerr);
                                    std::cerr<<std::endl;
                                    bin_adr->printColumn(std::cerr);
                                    std::cerr<<std::endl;
                                    Tparticle::printColumnTitle(std::cerr);
                                    std::cerr<<std::endl;
                                    for (int j=0; j<2; j++) {
                                        bin_adr->getMember(j)->printColumn(std::cerr);
                                        std::cerr<<std::endl;
                                    }
#endif

                                    // set binary tree mass to zero
                                    setBinaryCMZeroIter(bin_root);

                                    // cumulative step count 
                                    profile.step_count = step_count;
                                    profile.step_count_tsyn = step_count_tsyn;
                                    profile.step_count_sum += step_count;
                                    profile.step_count_tsyn_sum += step_count_tsyn;

                                    Float dt = _time_end - time_;
                                    time_ += dt;

                                    return bin_interrupt;
                                }

                                Float ekin_bk = ekin_;
                                Float epot_bk = epot_;
                                Float H_bk = getH();

#ifdef AR_SLOWDOWN_TREE
                                Float ekin_sd_bk = ekin_sd_;
                                Float epot_sd_bk = epot_sd_;
                                Float H_sd_bk = getHSlowDown();
#endif
                                
                                // update binary tree, put unused zero-mass particles out of the tree
                                info.generateBinaryTree(particles, G);
                                binary_update_flag = true;
                                
#ifdef AR_TTL
                                Float gt_kick_inv_bk = gt_kick_inv_.value;
                                calcAccPotAndGTKickInv();
                                Float d_gt_kick_inv = gt_kick_inv_.value - gt_kick_inv_bk;
                                // when the change is large, initialize gt_drift_inv_ to avoid large error
                                if (fabs(d_gt_kick_inv)/std::max(fabs(gt_kick_inv_bk),fabs(gt_kick_inv_.value)) >1e-3) 
                                    gt_drift_inv_ = gt_kick_inv_.value;
                                else 
                                    gt_drift_inv_ += d_gt_kick_inv;
#else
                                calcAccPotAndGTKickInv();
#endif
                                // calculate kinetic energy
                                calcEKin();

                                // Notice initially etot_ref_ does not include epert. The perturbation effect is accumulated in the integration. Here instance change of mass does not create any work. So no need to add de_pert
                                // get perturbation energy change due to mass change
                                //Float epert_new = 0.0;
                                //for (int i=0; i<n_particle; i++) {
                                //    epert_new += force_[i].pot_pert*particles[i].mass;
                                //}
                                //Float de_pert = epert_new - epert; // notice this is double perturbation potential

                                // get energy change
                                Float de = (ekin_ - ekin_bk) + (epot_ - epot_bk); //+ de_pert;
                                etot_ref_ += de;
                                de_change_interrupt_ += de;
                                dH_change_interrupt_ += getH() - H_bk;

#ifdef AR_SLOWDOWN_TREE
                                Float de_sd = (ekin_sd_ - ekin_sd_bk) + (epot_sd_ - epot_sd_bk);// + de_pert;
                                etot_sd_ref_ += de_sd;

                                Float dH_sd = getHSlowDown() - H_sd_bk;

                                // add slowdown change to the global slowdown energy
                                de_sd_change_interrupt_ += de_sd;
                                dH_sd_change_interrupt_ += dH_sd;
                                de_sd_change_cum_ += de_sd;
                                dH_sd_change_cum_ += dH_sd;
#endif //SLOWDOWN

#ifdef AR_DEBUG_PRINT
                                std::cerr<<"Interrupt condition triggered!";
                                std::cerr<<" Time: "<<time_;
#ifdef AR_SLOWDOWN_TREE
                                std::cerr<<" Energy change: dE_SD: "<<de_sd<<" dH_SD: "<<dH_sd;
                                std::cerr<<" Slowdown: "<<bin_root.slowdown.getSlowDownFactor()<<std::endl;
#endif
                                auto bin_adr = bin_interrupt.getBinaryTreeAddress();
                                bin_adr->printColumnTitle(std::cerr);
                                std::cerr<<std::endl;
                                bin_adr->printColumn(std::cerr);
                                std::cerr<<std::endl;
                                Tparticle::printColumnTitle(std::cerr);
                                std::cerr<<std::endl;
                                for (int j=0; j<2; j++) {
                                    bin_adr->getMember(j)->printColumn(std::cerr);
                                    std::cerr<<std::endl;
                                }
#endif

                                // change fix step option to make safety if energy change is large
                                //info.fix_step_option=FixStepOption::none;
                                
                                // if time_end flag set, reset it to be safety
                                //time_end_flag = false;

                                // check merger case
                                if (bin_interrupt.status==InterruptStatus::merge) {
                                    // count particle having mass
                                    int count_mass=0;
                                    int index_mass_last=-1;
                                    for (int j=0; j<n_particle; j++) {
                                        if (particles[j].mass>0.0) {
                                            count_mass++;
                                            index_mass_last=j;
                                        }
                                    }
                                    // only one particle has mass, drift directly
                                    if (count_mass==1) {
                                        ASSERT(index_mass_last<n_particle&&index_mass_last>=0);
                                        auto& p = particles[index_mass_last];
                                        Float dt = _time_end - time_;
                                        p.pos[0] += dt * p.vel[0];
                                        p.pos[1] += dt * p.vel[1];
                                        p.pos[2] += dt * p.vel[2];

                                        // cumulative step count 
                                        profile.step_count = step_count;
                                        profile.step_count_tsyn = step_count_tsyn;
                                        profile.step_count_sum += step_count;
                                        profile.step_count_tsyn_sum += step_count_tsyn;

                                        time_ += dt;

                                        return bin_interrupt;
                                    }
                                    // if only two particles have mass, switch off auto ds adjustment
                                    if (count_mass==2) {
                                        info.fix_step_option=FixStepOption::later;
                                    }
                                    //else {
                                    //    info.generateBinaryTree(particles, G);
                                    //}
                                }

#ifdef AR_SLOWDOWN_TREE
                                updateSlowDownAndCorrectEnergy(true, true);
#endif

                                info.ds = info.calcDsKeplerBinaryTree(*bin_interrupt.getBinaryTreeAddress(), manager->step.getOrder(), G, manager->ds_scale);
                                Float ds_max = manager->step.calcStepModifyFactorFromErrorRatio(2.0)*ds_init;
                                Float ds_min = manager->step.calcStepModifyFactorFromErrorRatio(0.5)*ds_init;
                                if (info.ds>ds_max || info.ds<ds_min) {
#ifdef AR_DEBUG_PRINT
                                    std::cerr<<"Change ds after interruption: ds(init): "<<ds_init<<" ds(new): "<<info.ds<<" ds(now): "<<ds[0]<<std::endl;
#endif
                                    ASSERT(info.ds>0);
                                    ds[0] = std::min(ds[0], info.ds);
                                    ds[1] = std::min(ds[1], info.ds);
                                    ds_backup.initial(info.ds);
                                    ds_init = info.ds;
                                }
                                else info.ds = ds_init;

                                // return one should be the top root
                                if (bin_interrupt_return.status!=InterruptStatus::none) {
                                    if (bin_interrupt_return.getBinaryTreeAddress()!= bin_interrupt.getBinaryTreeAddress()) {
                                        // give root address if interrupted binaries are different from previous one
                                        bin_interrupt_return.setBinaryTreeAddress(&(info.getBinaryTreeRoot()));
                                    }
                                    if (bin_interrupt.status==InterruptStatus::merge) 
                                        bin_interrupt_return.status = InterruptStatus::merge;
                                }
                                else bin_interrupt_return = bin_interrupt;
                            }
                            bin_interrupt.clear();
                        }
                    }


                    // update binary orbit and ds if unstable
                    if (!time_end_flag&&!binary_update_flag) {
                        bool update_flag=updateBinarySemiEccPeriodIter(bin_root, G, time_);

#ifdef AR_SLOWDOWN_TREE
                        updateSlowDownAndCorrectEnergy(true, true);
#endif

                        if (update_flag) {
                    // update slowdown and correct slowdown energy and gt_inv

#ifdef AR_DEBUG_PRINT
                            std::cerr<<"Update binary tree orbits, time= "<<time_<<"\n";
#endif
                            info.ds = info.calcDsKeplerBinaryTree(bin_root, manager->step.getOrder(), G, manager->ds_scale);
                            if (abs(ds_init-info.ds)/ds_init>0.1) {
#ifdef AR_DEBUG_PRINT
                                std::cerr<<"Change ds after update binary orbit: ds(init): "<<ds_init<<" ds(new): "<<info.ds<<" ds(now): "<<ds[0]<<std::endl;
#endif
                                ASSERT(info.ds>0);
                                ds[0] = std::min(ds[0], info.ds);
                                ds[1] = std::min(ds[1], info.ds);
                                ds_backup.initial(info.ds);
                                ds_init = info.ds;
                            }
                        }
                    }

                    int bk_return_size = backupIntData(backup_data);
                    ASSERT(bk_return_size == bk_data_size);
                    (void)bk_return_size;

                }
                else { //restore data
                    int bk_return_size = restoreIntData(backup_data);
                    ASSERT(bk_return_size == bk_data_size);
                    (void)bk_return_size;
//ifdef AR_SLOWDOWN_TREE
#ifndef USE_CM_FRAME
                    info.getBinaryTreeRoot().calcCenterOfMassIter();
#endif
//#endif
                }

                // get real time 
                Float dt = time_;

                // integrate one step
                ASSERT(!ISINF(ds[ds_switch]));
                if(n_particle==2) integrateTwoOneStep(ds[ds_switch], time_table);
                else integrateOneStep(ds[ds_switch], time_table);
                //info.generateBinaryTree(particles, G);

                // real step size
                dt =  time_ - dt;
//                ASSERT(dt>0.0);
                
                step_count++;

                // energy check
#ifdef AR_SLOWDOWN_TREE
                Float energy_error_bk = getEnergyErrorSlowDownFromBackup(backup_data);
                Float etot_ref_bk = getEtotSlowDownRefFromBackup(backup_data);
                Float energy_error = getEnergyErrorSlowDown();
                Float H_bk = getHSlowDownFromBackup(backup_data);
                Float H = getHSlowDown();
#else
                Float energy_error_bk = getEnergyErrorFromBackup(backup_data);
                Float etot_ref_bk = getEtotRefFromBackup(backup_data);
                Float energy_error = getEnergyError();
                Float H_bk = getHFromBackup(backup_data);
                Float H = getH();
#endif
                Float energy_error_diff = energy_error - energy_error_bk;

                Float energy_error_rel_abs = abs(energy_error_diff/etot_ref_bk);

                // get integration error for extended Hamiltonian
                Float integration_error_rel_abs = abs(H-H_bk);
                // H should be zero initially
                Float integration_error_rel_cum_abs = abs(H);

                Float integration_error_ratio = energy_error_rel_max/integration_error_rel_abs;
      
                // time error
                Float time_diff_rel = (_time_end - time_)/dt_full;

                //! regular block time step modification factor
                auto regularStepFactor = [](const Float _fac) {
                    Float fac = 1.0;
                    if (_fac<1) while (fac>_fac) fac *= 0.5;
                    else {
                        while (fac<=_fac) fac *= 2.0;
                        fac *= 0.5;
                    }
                    return fac;
                };

                Float error_increase_ratio_regular = manager->step.calcErrorRatioFromStepModifyFactor(2.0);

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
                    std::cerr<<error_message<<": "
                             <<"time "<<time_<<" " 
                             <<"ds_new "<<ds[1-ds_switch]<<" "
                             <<"ds_init "<<ds_init<<" "
                             <<"modify "<<step_modify_factor<<" "
                             <<"steps "<<step_count<<" "
                             <<"n_mods "<<reduce_ds_count<<" "
                             <<"err "<<integration_error_rel_abs<<" "
                             <<"err/max "<<1.0/integration_error_ratio<<" "
                             <<"errcum/E "<<integration_error_rel_cum_abs<<" "
                             <<"dt "<<dt<<" "
                             <<"n_ptcl "<<n_particle<<" ";
                    for (int i=0; i<info.binarytree.getSize(); i++) {
                        auto& bini = info.binarytree[i];
                        std::cerr<<"semi "<<bini.semi<<" "
                                 <<"ecc "<<bini.ecc<<" "
                                 <<"period "<<bini.period<<" "
                                 <<"m1 "<<bini.m1<<" "
                                 <<"m2 "<<bini.m2<<" "
                                 <<"stab "<<bini.stab<<" "
                                 <<"sd "<<bini.slowdown.getSlowDownFactor()<<" "
                                 <<"sd_org "<<bini.slowdown.getSlowDownFactorOrigin()<<" "
                                 <<"pert_in "<<bini.slowdown.getPertIn()<<" "
                                 <<"pert_out "<<bini.slowdown.getPertOut()<<" ";
                    }
                    std::cerr<<std::endl;
                };
#endif 

//#ifdef AR_WARN
                // warning for large number of steps
                if(warning_print_once&&step_count>=manager->step_count_max) {
                    if(step_count%manager->step_count_max==0) {
                        printMessage("Warning: step count is signficiant large");
                        for (int i=0; i<info.binarytree.getSize(); i++){
                            auto& bin = info.binarytree[i];
                            std::cerr<<"  Binary["<<i<<"]: "
                                     <<"  i1="<<bin.getMemberIndex(0)
                                     <<"  i2="<<bin.getMemberIndex(1)
                                     <<"  m1="<<bin.m1
                                     <<"  m2="<<bin.m2
                                     <<"  semi= "<<bin.semi
                                     <<"  ecc= "<<bin.ecc
                                     <<"  period= "<<bin.period
                                     <<"  stab= "<<bin.stab
                                     <<"  SD= "<<bin.slowdown.getSlowDownFactor()
                                     <<"  SD_org= "<<bin.slowdown.getSlowDownFactorOrigin()
                                     <<"  Tscale= "<<bin.slowdown.timescale
                                     <<"  pert_in= "<<bin.slowdown.pert_in
                                     <<"  pert_out= "<<bin.slowdown.pert_out;
                            std::cerr<<std::endl;
                            warning_print_once = false;
                        }
                        //printColumnTitle(std::cerr,20,info.binarytree.getSize());
                        //std::cerr<<std::endl;
                        //printColumn(std::cerr,20,info.binarytree.getSize());
                        //std::cerr<<std::endl;
#ifdef AR_DEBUG_DUMP
                        if (!info.dump_flag) {
                            DATADUMP("dump_large_step");
                            info.dump_flag=true;
                        }
#endif

//                        // increase step size if energy error is small, not works correctly, suppress
//                        if(integration_error_rel_abs<energy_error_rel_max) {
//                            Float integration_error_ratio = energy_error_rel_max/integration_error_rel_abs;
//                            Float step_modify_factor = manager->step.calcStepModifyFactorFromErrorRatio(integration_error_ratio);
//                            ASSERT(step_modify_factor>0.0);
//                            ds[ds_switch] *= step_modify_factor;
//                            info.ds = ds[ds_switch];
//                            ds[1-ds_switch] = ds[ds_switch];
//                            ds_backup.initial(info.ds);
//                            ds_init = info.ds;
//                            ASSERT(!ISINF(ds[ds_switch]));
//#ifdef AR_DEBUG_PRINT
//                            std::cerr<<"Energy error is small enough for increase step, integration_error_rel_abs="<<integration_error_rel_abs
//                                     <<" energy_error_rel_max="<<energy_error_rel_max<<" step_modify_factor="<<step_modify_factor<<" new ds="<<ds[1-ds_switch]<<std::endl;
//#endif
//                        }
                    }
                }
//#endif
          
                // When time sychronization steps too large, abort
                if(step_count_tsyn>manager->step_count_max) {
                    printMessage("Error! step count after time synchronization is too large");
                    printColumnTitle(std::cerr,20,info.binarytree.getSize());
                    std::cerr<<std::endl;
                    printColumn(std::cerr,20,info.binarytree.getSize());
                    std::cerr<<std::endl;
//                    restoreIntData(backup_data_init);
#ifdef AR_DEBUG_DUMP
                    if (!info.dump_flag) {
                        DATADUMP("dump_large_step");
                        info.dump_flag=true;
                    }
#endif
                    abort();
                }


#ifdef AR_DEEP_DEBUG
                printMessage("");
                std::cerr<<"Timetable: ";
                for (int i=0; i<cd_pair_size; i++) std::cerr<<" "<<time_table[manager->step.getSortCumSumCKIndex(i)];
                std::cerr<<std::endl;
#endif

                ASSERT(!ISNAN(integration_error_rel_abs));

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
                        if(step_count<5) {

                            // estimate the modification factor based on the symplectic order
                            // limit step_modify_factor to 0.125
                            step_modify_factor = std::max(regularStepFactor(manager->step.calcStepModifyFactorFromErrorRatio(integration_error_ratio)), Float(0.125));
                            ASSERT(step_modify_factor>0.0);

                            previous_step_modify_factor = step_modify_factor;
                            previous_error_ratio = integration_error_ratio;

                            ds[ds_switch] *= step_modify_factor;
                            ds[1-ds_switch] = ds[ds_switch];
                            // permanently reduce ds
                            // info.ds = ds[ds_switch];
                            // ASSERT(!ISINF(info.ds));
                            ds_backup.initial(info.ds);

                            backup_flag = false;
#ifdef AR_COLLECT_DS_MODIFY_INFO
                            collectDsModifyInfo("Large_energy_error");
#endif
                            continue;
                        }
                        // for big energy error, reduce step temparely
                        else if (info.fix_step_option==FixStepOption::none) {

                            // estimate the modification factor based on the symplectic order
                            // limit step_modify_factor to 0.125
                            step_modify_factor = std::max(regularStepFactor(manager->step.calcStepModifyFactorFromErrorRatio(integration_error_ratio)), Float(0.125));
                            ASSERT(step_modify_factor>0.0);

                            previous_step_modify_factor = step_modify_factor;
                            previous_error_ratio = integration_error_ratio;

                            ds_backup.backup(ds[ds_switch], step_modify_factor);
                            if(previous_is_restore) reduce_ds_count++;

                            ds[ds_switch] *= step_modify_factor;
                            ds[1-ds_switch] = ds[ds_switch];
                            ASSERT(!ISINF(ds[ds_switch]));

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
                    step_modify_factor = std::min(std::max(regularStepFactor(manager->step.calcStepModifyFactorFromErrorRatio(abs(_time_end/dt))), Float(0.0625)),Float(0.5)); 
                    ASSERT(step_modify_factor>0.0);
                    previous_step_modify_factor = step_modify_factor;
                    previous_error_ratio = integration_error_ratio;

                    ds[ds_switch] *= step_modify_factor;
                    ds[1-ds_switch] = ds[ds_switch];
                    ASSERT(!ISINF(ds[ds_switch]));

                    // for initial steps, reduce step permanently
                    if (step_count<5) {
                        //info.ds = ds[ds_switch];
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

                // if no modification, reset previous values
                previous_step_modify_factor = 1.0;
                previous_error_ratio = -1.0;

                // check integration time
                if(time_ < _time_end - time_error){
                    // step increase depend on n_step_wait_recover_ds
                    if(info.fix_step_option==FixStepOption::none && !time_end_flag) {
                        // waiting step count reach
                        previous_is_restore=ds_backup.countAndRecover(ds[1-ds_switch], step_modify_factor, integration_error_ratio>error_increase_ratio_regular);
                        if (previous_is_restore) {
                            //previous_error_ratio = -1;
                            //previous_step_modify_factor = 1.0;
#ifdef AR_COLLECT_DS_MODIFY_INFO
                            collectDsModifyInfo("Reuse_backup_ds");
#endif
                        }
                        // increase step size if energy error is small, not works correctly, integration error may not increase when ds becomes larger, then a very large ds may appear after several iterations. suppress
                        else if(integration_error_rel_abs<0.5*energy_error_rel_max && dt>0.0 && dt_full/dt>std::max(100.0,0.02*manager->step_count_max)) {
                            Float integration_error_ratio = energy_error_rel_max/integration_error_rel_abs;
                            Float step_modify_factor = std::min(Float(100.0),manager->step.calcStepModifyFactorFromErrorRatio(integration_error_ratio));
                            ASSERT(step_modify_factor>0.0);
                            ds[1-ds_switch] *= step_modify_factor;
                            info.ds = ds[1-ds_switch];
                            ASSERT(!ISINF(ds[1-ds_switch]));
#ifdef AR_DEBUG_PRINT
                            std::cerr<<"Energy error is small enough for increase step, integration_error_rel_abs="<<integration_error_rel_abs
                                     <<" energy_error_rel_max="<<energy_error_rel_max<<" step_modify_factor="<<step_modify_factor<<" new ds="<<ds[1-ds_switch]<<std::endl;
#endif
                        }
                    }

                    // time sychronization on case, when step size too small to reach time end, increase step size
                    if(time_end_flag && ds[ds_switch]==ds[1-ds_switch]) {
                        step_count_tsyn++;

                        Float dt_end = _time_end - time_;
                        if (dt<0) {
                            // limit step_modify_factor to 0.125
                            step_modify_factor = std::min(std::max(regularStepFactor(manager->step.calcStepModifyFactorFromErrorRatio(abs(_time_end/dt))), Float(0.0625)),Float(0.5)); 
                            ASSERT(step_modify_factor>0.0);

                            ds[ds_switch] *= step_modify_factor;
                            ds[1-ds_switch] = ds[ds_switch];
                            ASSERT(!ISINF(ds[ds_switch]));
                        }
                        else if (n_step_end>1 && dt<0.3*dt_end) {
                            // dt should be >0.0
                            // ASSERT(dt>0.0);
                            ds[1-ds_switch] = ds[ds_switch] * dt_end/dt;
                            ASSERT(!ISINF(ds[1-ds_switch]));
#ifdef AR_DEEP_DEBUG
                            std::cerr<<"Time step dt(real) "<<dt<<" <0.3*(time_end-time)(real) "<<dt_end<<" enlarge step factor: "<<dt_end/dt<<" new ds: "<<ds[1-ds_switch]<<std::endl;
#endif
                        }
                        else n_step_end++;
                    }

                    // when used once, update to the new step
                    ds[ds_switch] = ds[1-ds_switch]; 
                    ASSERT(!ISINF(ds[ds_switch]));
                    ds_switch = 1-ds_switch;

                    if (dt>0) backup_flag = true;
                    else backup_flag = false;
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
                        ASSERT(!ISINF(ds[ds_switch]));
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
                        ASSERT(!ISINF(cck_prev));
                        ds[ds_switch] *= cck_prev;  
                        ASSERT(!ISINF(ds[ds_switch]));
                        // then next next step, scale with the CumSum CK between two step: cck(i) - cck(i-1) 
                        ASSERT(dt_k>0.0);
                        ds[1-ds_switch] = ds_tmp*(cck-cck_prev)*std::min(Float(1.0),(_time_end-time_prev+time_error)/dt_k); 
                        ASSERT(!ISINF(ds[1-ds_switch]));

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

#ifdef AR_SLOWDOWN_TREE

        //! write back slowdown particle iteration function
        /*!
          @param[in] _pos_up: upper cm position
          @param[in] _vel_sd_up: upper slowdown cm velocity
          @param[in] _inv_nest_sd_up: upper inverse nested slowdown factor
          @param[in] _bin: current binary
         */
        void writeBackSlowDownParticlesIter(const Float* _pos_up,
                                            const Float* _vel_sd_up, 
                                            const Float& _inv_nest_sd_up, 
                                            AR::BinaryTree<Tparticle>& _bin) {
            Float inv_nest_sd = _inv_nest_sd_up/_bin.slowdown.getSlowDownFactor();
#ifndef USE_CM_FRAME
            Float* vel_cm = _bin.getVel();
#endif
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    Float* vel = bink->getVel();
#ifdef USE_CM_FRAME
                    Float* pos_rel = bink->getPos();
                    const Float pos[3] = {pos_rel[0] + _pos_up[0],
                                          pos_rel[1] + _pos_up[1],
                                          pos_rel[2] + _pos_up[2]};

                    Float vel_sd[3] = {vel[0] * inv_nest_sd + _vel_sd_up[0], 
                                       vel[1] * inv_nest_sd + _vel_sd_up[1], 
                                       vel[2] * inv_nest_sd + _vel_sd_up[2]}; 
#else
                    const Float* pos = _pos_up;
                    Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                       (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                       (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]}; 
#endif
                    writeBackSlowDownParticlesIter(pos, vel_sd, inv_nest_sd, *bink);
                }
                else {
                    int i = _bin.getMemberIndex(k);
                    auto& pk = particles[i];
                    auto* pk_adr = particles.getMemberOriginAddress(i);
                    Float* vel = pk.getVel();
                    pk_adr->mass = pk.mass;
#ifdef USE_CM_FRAME
                    Float vel_sd[3] = {vel[0] * inv_nest_sd + _vel_sd_up[0], 
                                       vel[1] * inv_nest_sd + _vel_sd_up[1], 
                                       vel[2] * inv_nest_sd + _vel_sd_up[2]}; 
#else
                    Float vel_sd[3] = {(vel[0] - vel_cm[0]) * inv_nest_sd + _vel_sd_up[0], 
                                       (vel[1] - vel_cm[1]) * inv_nest_sd + _vel_sd_up[1], 
                                       (vel[2] - vel_cm[2]) * inv_nest_sd + _vel_sd_up[2]};
#endif
                    pk_adr->pos[0] = pk.pos[0] + _pos_up[0];
                    pk_adr->pos[1] = pk.pos[1] + _pos_up[1];
                    pk_adr->pos[2] = pk.pos[2] + _pos_up[2];

                    pk_adr->vel[0] = vel_sd[0];
                    pk_adr->vel[1] = vel_sd[1];
                    pk_adr->vel[2] = vel_sd[2];
                }
            }
        }

        //! write back particles with slowdown velocity
        /*! write back particles with slowdown velocity to original address
          @param[in] _particle_cm: center of mass particle to calculate the original frame, different from the particles.cm
         */
        template <class Tptcl>
        void writeBackSlowDownParticles(const Tptcl& _particle_cm) {
            //! iteration function using binarytree
            auto& bin_root=info.getBinaryTreeRoot();
            const Float* vel_cm = &(_particle_cm.vel[0]);
            const Float* pos_cm = &(_particle_cm.pos[0]);
            Float sd_factor=1.0;
            writeBackSlowDownParticlesIter(pos_cm, vel_cm, sd_factor, bin_root);
        }
        
#endif

#ifdef USE_CM_FRAME

        //! write back slowdown particle iteration function
        /*!
          @param[in] _pos_up: upper cm position
          @param[in] _vel_up: upper cm velocity
          @param[in] _bin: current binary
         */
        template <class Tptcl>
        void writeBackParticlesOriginFrameIter(const Float* _pos_up,
                                               const Float* _vel_up, 
                                               AR::BinaryTree<Tparticle>& _bin) {
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    auto* bink = _bin.getMemberAsTree(k);
                    
                    Float* vel_rel = bink->getVel();
                    Float* pos_rel = bink->getPos();
                    Float pos[3] = {pos_rel[0] + _pos_up[0],
                                    pos_rel[1] + _pos_up[1],
                                    pos_rel[2] + _pos_up[2]};

                    Float vel[3] = {vel_rel[0] + _vel_up[0], 
                                    vel_rel[1] + _vel_up[1], 
                                    vel_rel[2] + _vel_up[2]}; 

                    writeBackParticlesOriginFrameIter(pos, vel, *bink);
                }
                else {
                    int i = _bin.getMemberIndex(k);
                    auto& pk = particles[i];
                    auto* pk_adr = particles.getMemberOriginAddress(i);

                    pk_adr->mass = pk.mass;
                    pk_adr->pos[0] = pk.pos[0] + _pos_up[0];
                    pk_adr->pos[1] = pk.pos[1] + _pos_up[1];
                    pk_adr->pos[2] = pk.pos[2] + _pos_up[2];

                    pk_adr->vel[0] = pk.vel[0] + _vel_up[0];
                    pk_adr->vel[1] = pk.vel[1] + _vel_up[1];
                    pk_adr->vel[2] = pk.vel[2] + _vel_up[2];
                }
            }
        }

        //! write back particles to original address
        /*! If particles are in center-of-the-mass frame, write back the particle in original frame but not modify local copies to avoid roundoff error
         */
        template <class Tptcl>
        void writeBackParticlesOriginFrame() {
            ASSERT(particles.getMode()==COMM::ListMode::copy);
            auto& bin_root=info.getBinaryTreeRoot();
            if (bin_root.isOriginFrame()) {
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
            else {
                //! iteration function using binarytree
                if (particles.isOriginFrame()) {
                    Float pos_cm = {0.0, 0.0, 0.0};
                    Float vel_cm = {0.0, 0.0, 0.0};
                    writeBackParticlesOriginFrameIter(pos_cm, vel_cm, bin_root);
                }
                else {
                    Float* vel_cm = &(particles.cm.vel[0]);
                    Float* pos_cm = &(particles.cm.pos[0]);
                    writeBackParticlesOriginFrameIter(pos_cm, vel_cm, bin_root);
                }
            }
        }

#else
        //! write back particles to original address
        /*! If particles are in center-of-the-mass frame, write back the particle in original frame but not modify local copies to avoid roundoff error
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
#endif

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

        //! get perturbation potential energy
        /*! \return perturbation potential energy
         */
        Float getEpert() const {
            Float epert=0.0;
            int n_particle = particles.getSize();
            for (int i=0; i<n_particle; i++) {
                epert += force_[i].pot_pert*particles[i].mass;
            }
            return epert;
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

        //! reset cumulative energy/hamiltonian change due to interruption
        void resetDEChangeBinaryInterrupt() {
            de_change_interrupt_ = 0.0;
            dH_change_interrupt_ = 0.0;
        }

        //! get cumulative energy change due to interruption
        Float getDEChangeBinaryInterrupt() const {
            return de_change_interrupt_;
        }

        //! get cumulative hamiltonian change due to interruption
        Float getDHChangeBinaryInterrupt() const {
            return dH_change_interrupt_;
        }

        //! get Hamiltonian
        Float getH() const {
#ifdef AR_TTL
            //return (ekin_ - etot_ref_)/gt_drift_inv_ + epot_/gt_kick_inv_;
            return (ekin_ + epot_ - etot_ref_)/gt_kick_inv_.value;
#else
            return manager->interaction.calcH(ekin_ - etot_ref_, epot_);
#endif
        }

        //! get Hamiltonian from backup data
        Float getHFromBackup(Float* _bk) const {
            Float& etot_ref =_bk[1];
            Float& ekin = _bk[2];
            Float& epot = _bk[3];
#ifdef AR_TTL
#ifdef AR_SLOWDOWN_TREE
            //Float& gt_drift_inv = _bk[13];
            Float& gt_kick_inv  = _bk[14];
#else
            //Float& gt_drift_inv = _bk[6];
            Float& gt_kick_inv  = _bk[7];
#endif
            return (ekin + epot - etot_ref)/gt_kick_inv;
            //return (ekin - etot_ref)/gt_drift_inv + epot/gt_kick_inv;
#else
            return manager->interaction.calcH(ekin - etot_ref, epot);
#endif
        }



#ifdef AR_SLOWDOWN_TREE
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

        //! reset cumulative energy change due to interruption
        void resetDESlowDownChangeBinaryInterrupt() {
            de_sd_change_interrupt_ = 0.0;
            dH_sd_change_interrupt_ = 0.0;
        }

        //! get cumulative energy change due to interruption
        Float getDESlowDownChangeBinaryInterrupt() const {
            return de_sd_change_interrupt_;
        }

        //! get cumulative hamiltonian change due to interruption
        Float getDHSlowDownChangeBinaryInterrupt() const {
            return dH_sd_change_interrupt_;
        }

        //! Get current kinetic energy with inner slowdown
        /*! \return current kinetic energy with inner slowdown
         */
        Float getEkinSlowDown() const {
            return ekin_sd_;
        }

        //! Get current potential energy with inner slowdown
        /*! \return current potetnial energy with inner slowdown (negative value for bounded systems)
         */
        Float getEpotSlowDown() const {
            return epot_sd_;
        }

        //! Get current total integrated energy with inner slowdown
        /*! \return total integrated energy with inner slowdown
         */
        Float getEtotSlowDownRef() const {
            return etot_sd_ref_;
        }

        //! Get current total energy with inner slowdown from ekin_sdi and epot_sdi
        /*! \return total energy with inner slowdown
         */
        Float getEtotSlowDown() const {
            return ekin_sd_ + epot_sd_;
        }

        //! get energy error with inner slowdown
        /*! \return energy error with inner slowdown
         */
        Float getEnergyErrorSlowDown() const {
            return ekin_sd_ + epot_sd_ - etot_sd_ref_;
        }

        //! get energy error with inner slowdown from backup data 
        Float getEnergyErrorSlowDownFromBackup(Float* _bk) const {
            return -_bk[6] + _bk[7] + _bk[8];
        }

        //! get slowdown Hamiltonian
        Float getHSlowDown() const {
#ifdef AR_TTL
            //return (ekin_sd_ - etot_sd_ref_)/gt_drift_inv_ + epot_sd_/gt_kick_inv_;
            return (ekin_sd_ + epot_sd_ - etot_sd_ref_)/gt_kick_inv_.value;
#else
            return manager->interaction.calcH(ekin_sd_ - etot_sd_ref_, epot_sd_);
#endif
        }

        //! get slowdown Hamiltonian from backup data
        Float getHSlowDownFromBackup(Float* _bk) const {
            Float& etot_sd_ref =_bk[6];
            Float& ekin_sd = _bk[7];
            Float& epot_sd = _bk[8];
#ifdef AR_TTL
            //Float& gt_drift_inv = _bk[13];
            Float& gt_kick_inv  = _bk[14];
            //return (ekin_sd - etot_sd_ref)/gt_drift_inv + epot_sd/gt_kick_inv;
            return (ekin_sd + epot_sd - etot_sd_ref)/gt_kick_inv;
#else
            return manager->interaction.calcH(ekin_sd - etot_sd_ref, epot_sd);            
#endif
        }

        //! get integrated energy with inner slowdown from backup data
        Float getEtotSlowDownRefFromBackup(Float* _bk) const {
            return _bk[6];
        }

        //! get energy with inner slowdown from backup data (ekin_sdi + epot_sdi)
        Float getEtotSlowDownFromBackup(Float* _bk) const {
            return _bk[7] + _bk[8];
        }

#endif

        //! get backup data size
        int getBackupDataSize() const {
            int bk_size = 6;
#ifdef AR_SLOWDOWN_TREE
            bk_size += 7; 
            //bk_size += SlowDown::getBackupDataSize();
#endif
#ifdef AR_TTL
            bk_size += 2;
#endif
            bk_size += particles.getSize()*6;
#ifdef USE_CM_FRAME
            bk_size += info.binarytree.getSize()*6;
#endif
            return bk_size;
        }

        //! Backup integration data 
        /*! Backup $time_, $etot_, $ekin_, $epot_, $gt_drift_, $gt_kick_inv_, #particles, $slowdown to one Float data array
          \return backup array size
        */
        int backupIntData(Float* _bk) {
            int bk_size=0;
            _bk[bk_size++] = time_;       //0
            _bk[bk_size++] = etot_ref_;   //1
            _bk[bk_size++] = ekin_;       //2
            _bk[bk_size++] = epot_;       //3 
            _bk[bk_size++] = de_change_interrupt_;     //4
            _bk[bk_size++] = dH_change_interrupt_;     //5
#ifdef AR_SLOWDOWN_TREE
            _bk[bk_size++] = etot_sd_ref_; //6
            _bk[bk_size++] = ekin_sd_;     //7
            _bk[bk_size++] = epot_sd_;     //8
            _bk[bk_size++] = de_sd_change_cum_; //9
            _bk[bk_size++] = dH_sd_change_cum_; //10
            _bk[bk_size++] = de_sd_change_interrupt_; //11
            _bk[bk_size++] = dH_sd_change_interrupt_; //12
#endif

#ifdef AR_TTL
            _bk[bk_size++] = gt_drift_inv_;  //13 / 6
            _bk[bk_size++] = gt_kick_inv_.value;   //14 / 7
#endif

            for (int i=0; i<particles.getSize(); i++) {
                Float* pos = particles[i].getPos();
                Float* vel = particles[i].getVel();
                _bk[bk_size++] = pos[0];
                _bk[bk_size++] = pos[1];
                _bk[bk_size++] = pos[2];
                _bk[bk_size++] = vel[0];
                _bk[bk_size++] = vel[1];
                _bk[bk_size++] = vel[2];
            }

#ifdef USE_CM_FRAME
            std::function<void(AR::BinaryTree<Tparticle>&)> backupCMPosVelIter = [&](AR::BinaryTree<Tparticle>& _bin) {
                _bk[bk_size++] = _bin.pos[0];
                _bk[bk_size++] = _bin.pos[1];
                _bk[bk_size++] = _bin.pos[2];
                _bk[bk_size++] = _bin.vel[0];
                _bk[bk_size++] = _bin.vel[1];
                _bk[bk_size++] = _bin.vel[2];
                for (int k=0; k<2; k++) {
                    if (_bin.isMemberTree(k)) {
                        auto* bink = _bin.getMemberAsTree(k);
                        backupCMPosVelIter(*bink);
                    } 
                }
            };
            backupCMPosVelIter(info.getBinaryTreeRoot());
#endif
//#ifdef AR_SLOWDOWN_TREE
//            bk_size += info.getBinaryTreeRoot().slowdown.backup(&_bk[bk_size]); // slowdownfactor
//#endif
            return bk_size;
        }

        //! Restore integration data
        /*! restore $time_, $etot_, $ekin_, $epot_, $gt_drift_, $gt_kick_inv_, #particles, $slowdown from one Float data array
          \return backup array size
        */
        int restoreIntData(Float* _bk) {
            int bk_size = 0;
            time_     = _bk[bk_size++];
            etot_ref_ = _bk[bk_size++];
            ekin_     = _bk[bk_size++];
            epot_     = _bk[bk_size++];
            de_change_interrupt_= _bk[bk_size++];
            dH_change_interrupt_= _bk[bk_size++];

#ifdef AR_SLOWDOWN_TREE
            etot_sd_ref_ = _bk[bk_size++];
            ekin_sd_     = _bk[bk_size++];
            epot_sd_     = _bk[bk_size++];

            de_sd_change_cum_= _bk[bk_size++];
            dH_sd_change_cum_= _bk[bk_size++];
            de_sd_change_interrupt_= _bk[bk_size++];
            dH_sd_change_interrupt_= _bk[bk_size++];
#endif
#ifdef AR_TTL
            gt_drift_inv_  = _bk[bk_size++];
            gt_kick_inv_.value   = _bk[bk_size++];
#endif

            //! restore member particle position and velocity
            for (int i=0; i<particles.getSize(); i++) {
                Float* pos = particles[i].getPos();
                Float* vel = particles[i].getVel();
                pos[0] = _bk[bk_size++];
                pos[1] = _bk[bk_size++];
                pos[2] = _bk[bk_size++];
                vel[0] = _bk[bk_size++];
                vel[1] = _bk[bk_size++];
                vel[2] = _bk[bk_size++];
            }

#ifdef USE_CM_FRAME
            std::function<void(AR::BinaryTree<Tparticle>&)> restoreCMPosVelIter = [&](AR::BinaryTree<Tparticle>& _bin) {
                _bin.pos[0] = _bk[bk_size++];
                _bin.pos[1] = _bk[bk_size++];
                _bin.pos[2] = _bk[bk_size++];
                _bin.vel[0] = _bk[bk_size++];
                _bin.vel[1] = _bk[bk_size++];
                _bin.vel[2] = _bk[bk_size++];
                for (int k=0; k<2; k++) {
                    if (_bin.isMemberTree(k)) {
                        auto* bink = _bin.getMemberAsTree(k);
                        restoreCMPosVelIter(*bink);
                    } 
                }
            };
            restoreCMPosVelIter(info.getBinaryTreeRoot());
#endif      

//#ifdef AR_SLOWDOWN_TREE
//            bk_size += info.getBinaryTreeRoot().slowdown.restore(&_bk[bk_size]);
//#endif
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

        //! print group information 
        /*! Message, Number of members, time, binary tree printing interation
          @param[in] _type: 0: new group (if pair id is same, no printing); 1: end group (always print and reset pair id)
          @param[in] _fout: FILE IO
          @param[in] _width: print width
          @param[in] _pcm: center of mass particle to calculate origin position and velocity, if NULL, assume cm pos and vel are zero
        */
        template<class Tptcl>
        void printGroupInfo(const int _type, std::ostream& _fout, const int _width, const Tptcl* _pcm=NULL) {
            auto& bin_root = info.getBinaryTreeRoot();
            //auto* p1 = bin_root.getLeftMember();
            //auto* p2 = bin_root.getRightMember();
            
            bool reset_flag = (_type==1 && bin_root.semi<0 && bin_root.ecca>0);

            if (info.checkAndSetBinaryPairIDIter(bin_root, reset_flag)) {
                if (_type==0) return; // if it is new but already existed binary, do not print
                else if (!reset_flag) return; // in the end case, if the system is still bound, do not print 
            }

            Float pos_cm[3], vel_cm[3];
            auto& pcm_loc = particles.cm;
            if (_pcm!=NULL) {
                pos_cm[0] = pcm_loc.pos[0] + _pcm->pos[0];
                pos_cm[1] = pcm_loc.pos[1] + _pcm->pos[1];
                pos_cm[2] = pcm_loc.pos[2] + _pcm->pos[2];
                vel_cm[0] = pcm_loc.vel[0] + _pcm->vel[0];
                vel_cm[1] = pcm_loc.vel[1] + _pcm->vel[1];
                vel_cm[2] = pcm_loc.vel[2] + _pcm->vel[2];
            }
            else {
                pos_cm[0] = pcm_loc.pos[0]; 
                pos_cm[1] = pcm_loc.pos[1]; 
                pos_cm[2] = pcm_loc.pos[2]; 
                vel_cm[0] = pcm_loc.vel[0]; 
                vel_cm[1] = pcm_loc.vel[1]; 
                vel_cm[2] = pcm_loc.vel[2]; 
            }
#pragma omp critical
            {
                _fout<<std::setw(_width)<<_type
                     <<std::setw(_width)<<bin_root.getMemberN()
                     <<std::setw(_width)<<time_ + info.time_offset;
                _fout<<std::setw(_width)<<pos_cm[0]
                     <<std::setw(_width)<<pos_cm[1]
                     <<std::setw(_width)<<pos_cm[2]
                     <<std::setw(_width)<<vel_cm[0]
                     <<std::setw(_width)<<vel_cm[1]
                     <<std::setw(_width)<<vel_cm[2];
                bin_root.printBinaryTreeIter(_fout, _width);
                _fout<<std::endl;
            }
            //if (_type==0) { // register pair id to avoid repeating printing
            //    p1->setBinaryPairID(p2->id);
            //    p2->setBinaryPairID(p1->id);
            //}
            //else { // break case reset pair id
            //    p1->setBinaryPairID(0);
            //    p2->setBinaryPairID(0);
            //}
        }

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
                 <<std::setw(_width)<<"dE_intr"
                 <<std::setw(_width)<<"dH_intr";
            perturber.printColumnTitle(_fout, _width);
            info.printColumnTitle(_fout, _width);
            profile.printColumnTitle(_fout, _width);
#ifdef AR_SLOWDOWN_TREE
            _fout<<std::setw(_width)<<"dE_SD" 
                 <<std::setw(_width)<<"Etot_SD" 
                 <<std::setw(_width)<<"Ekin_SD" 
                 <<std::setw(_width)<<"Epot_SD"
                 <<std::setw(_width)<<"dE_SDC_cum" 
                 <<std::setw(_width)<<"dH_SDC_cum"
                 <<std::setw(_width)<<"dE_SDC_intr" 
                 <<std::setw(_width)<<"dH_SDC_intr"; 
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
                 <<std::setw(_width)<<epot_
#ifdef AR_SLOWDOWN_TREE
#ifdef AR_TTL
                 <<std::setw(_width)<<1.0/gt_drift_inv_
#else
                 <<std::setw(_width)<<1.0/manager->interaction.calcGTDriftInv(ekin_sd_-etot_sd_ref_)
#endif
                 <<std::setw(_width)<<getHSlowDown()
#else
#ifdef AR_TTL
                 <<std::setw(_width)<<1.0/gt_drift_inv_
#else
                 <<std::setw(_width)<<1.0/manager->interaction.calcGTDriftInv(ekin_-etot_ref_)
#endif
                 <<std::setw(_width)<<getH()
#endif
                 <<std::setw(_width)<<de_change_interrupt_
                 <<std::setw(_width)<<dH_change_interrupt_;
            perturber.printColumn(_fout, _width);
            info.printColumn(_fout, _width);
            profile.printColumn(_fout, _width);
#ifdef AR_SLOWDOWN_TREE
            _fout<<std::setw(_width)<<getEnergyErrorSlowDown()
                 <<std::setw(_width)<<etot_sd_ref_ 
                 <<std::setw(_width)<<ekin_sd_ 
                 <<std::setw(_width)<<epot_sd_
                 <<std::setw(_width)<<de_sd_change_cum_
                 <<std::setw(_width)<<dH_sd_change_cum_
                 <<std::setw(_width)<<de_sd_change_interrupt_
                 <<std::setw(_width)<<dH_sd_change_interrupt_;
            SlowDown sd_empty;
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
#ifdef USE_CM_FRAME
            _fout<<std::setw(_width)<<particles.getSize();
            particles.cm.printColumn(_fout, _width);
            info.getBinaryTreeRoot().printMemberIter(_fout, _width);
#else                
            particles.printColumn(_fout, _width);
#endif            
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

#ifdef USE_CM_FRAME
            int n_particle = particles.getSize();
            fwrite(&n_particle, sizeof(int), 1, _fout);
            info.getBinaryTreeRoot().writeMemberBinaryIter(_fout);
#else
            particles.writeBinary(_fout);
#endif            
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
