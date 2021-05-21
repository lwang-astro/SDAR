#pragma once

#include "Common/Float.h"
#include "Common/list.h"
#include "AR/symplectic_integrator.h"
#include "Hermite/ar_information.h"
#include "Hermite/hermite_particle.h"
#include "Hermite/block_time_step.h"
#include "Hermite/neighbor.h"
#include "Hermite/profile.h"
#include <map>

namespace H4{

    //! print features
    void printFeatures(std::ostream & fout) {
#ifdef ADJUST_GROUP_PRINT
        fout<<"Print adjust group information\n";
#endif
    }

    //! print features
    void printDebugFeatures(std::ostream & fout) {
#ifdef HERMITE_DEBUG
        fout<<"Debug mode: HERMITE\n";
#endif
    }    

    //! Hermite manager class
    template <class Tmethod>
    class HermiteManager{
    public:
        //Float r_break_crit; ///> the distance criterion to break the groups
        //Float r_neighbor_crit; ///> the distance for neighbor search
        Tmethod interaction; ///> class contain interaction function
        BlockTimeStep4th step; ///> time step calculator
#ifdef ADJUST_GROUP_PRINT
        bool adjust_group_write_flag; ///> flag to indicate whether to output new/end group information
        std::ofstream fgroup; ///> pointer to a file IO to output new/end group information
#endif

#ifdef ADJUST_GROUP_PRINT
        HermiteManager(): interaction(), step(), adjust_group_write_flag(true), fgroup() {}
#endif


        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            //ASSERT(r_break_crit>=0.0);
            //ASSERT(r_neighbor_crit>=0.0);
            ASSERT(interaction.checkParams());
            ASSERT(step.checkParams());
#ifdef ADJUST_GROUP_PRINT
            ASSERT(!adjust_group_write_flag||(adjust_group_write_flag&&fgroup.is_open()));
#endif
            return true;
        }

        //! print parameters
        void print(std::ostream & _fout) const{
            //_fout<<"r_break_crit    : "<<r_break_crit<<std::endl
            //<<"r_neighbor_crit : "<<r_neighbor_crit<<std::endl;
            interaction.print(_fout);
            step.print(_fout);
        }


        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            interaction.printColumnTitle(_fout, _width);
            step.printColumnTitle(_fout, _width);
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            interaction.printColumn(_fout, _width);
            step.printColumn(_fout, _width);
        }

        //! write class data to file with binary format
        /*! @param[in] _fp: FILE type file for output
         */
        void writeBinary(FILE *_fp) const {
            //size_t size = sizeof(*this) - sizeof(interaction) - sizeof(step);
            //fwrite(this, size, 1, _fp);
            interaction.writeBinary(_fp);
            step.writeBinary(_fp);
#ifdef ADJUST_GROUP_PRINT
            fwrite(&adjust_group_write_flag, sizeof(bool),1,_fp);
#endif
        }

        //! read class data to file with binary format
        /*! @param[in] _fp: FILE type file for reading
         */
        void readBinary(FILE *_fin) {
            //size_t size = sizeof(*this) - sizeof(interaction) - sizeof(step);
            //size_t rcount = fread(this, size, 1, _fin);
            //if (rcount<1) {
            //    std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            //    abort();
            //}
            interaction.readBinary(_fin);
            step.readBinary(_fin);
#ifdef ADJUST_GROUP_PRINT
            size_t rcount = fread(&adjust_group_write_flag, sizeof(bool),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
#endif
        }

    };

    //! Hermite energy 
    struct HermiteEnergy{
    public:
        Float etot_ref;  ///< total energy for reference
        Float ekin;  ///< kinetic energy
        Float epot;  ///< potential energy
        Float epert; ///< perturbation extra energy
        Float de_cum; ///< cumulative energy change 
        Float de_binary_interrupt; ///< energy change due to interruption
        Float de_modify_single;    ///< energy change due to modify function

        HermiteEnergy(): etot_ref(0.0), ekin(0.0), epot(0.0), epert(0.0), de_cum(0.0), de_binary_interrupt(0.0), de_modify_single(0.0) {}

        //! clear function
        void clear() {
            etot_ref = ekin = epot = epert = 0.0;
            de_cum = de_binary_interrupt = de_modify_single = 0.0;
        }

        //! calc energy error 
        Float getEnergyError() const {
            return ekin + epot + epert - etot_ref;
        }

        //! calc total energy
        Float getEtot() const {
            return ekin + epot + epert;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"dE"
                 <<std::setw(_width)<<"Etot_ref"
                 <<std::setw(_width)<<"Ekin"
                 <<std::setw(_width)<<"Epot"
                 <<std::setw(_width)<<"Epert"
                 <<std::setw(_width)<<"dE_cum"
                 <<std::setw(_width)<<"dE_intr"
                 <<std::setw(_width)<<"dE_mod";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<getEnergyError()
                 <<std::setw(_width)<<etot_ref
                 <<std::setw(_width)<<ekin
                 <<std::setw(_width)<<epot
                 <<std::setw(_width)<<epert
                 <<std::setw(_width)<<de_cum
                 <<std::setw(_width)<<de_binary_interrupt
                 <<std::setw(_width)<<de_modify_single;
        }
    };

    //!Hermite integrator class
    template <class Tparticle, class Tpcm, class Tpert, class TARpert, class Tacc, class TARacc, class Tinfo>
    class HermiteIntegrator{
    private:
        typedef ParticleH4<Tparticle> H4Ptcl;
        typedef AR::TimeTransformedSymplecticIntegrator<Tparticle, H4Ptcl, TARpert, TARacc, ARInformation<Tparticle>> ARSym;

        // time 
        Float time_;   ///< integrated time 
        Float time_offset_; ///< offset to obtain the real time (real time = time_ + time_offset_)
        Float time_next_min_; ///< the next minimum time 
        Float dt_limit_;  // maximum step size allown for next integration step calculation

        // Energy
        Float energy_init_ref_; // initial energy reference
        HermiteEnergy energy_; // true energy
        HermiteEnergy energy_sd_;  // slowdown energy

        // active particle number
        int n_act_single_;    /// active particle number of singles
        int n_act_group_;    /// active particle number of groups

        // initial particle number
        int n_init_single_;  /// number of singles need to initial
        int n_init_group_;   /// number of groups need to initial

        // group offset
        int index_offset_group_; /// offset of group index in pred_, force_ and time_next_

        // interrupt group 
        int interrupt_group_dt_sorted_group_index_; /// interrupt group position in index_dt_sorted_group_
        AR::InterruptBinary<Tparticle> interrupt_binary_; /// interrupt binary tree address
        COMM::List<int> index_group_merger_; /// merger group index

        // flags
        bool initial_system_flag_; /// flag to indicate whether the system is initialized with all array size defined (reest in initialSystem)
        bool modify_system_flag_;  /// flag to indicate whether the system (group/single) is added/removed (reset in adjustSystemAfterModify)

        // arrays
        // sorted time index list for select active particles
        COMM::List<int> index_dt_sorted_single_; // index of particles with time next sorted order (small to large)
        COMM::List<int> index_dt_sorted_group_; // index list of active groups

        // index list to store the groups with members resolved/only cm in Hermite integration 
        COMM::List<int> index_group_resolve_;   // index of resolved group for integration
        COMM::List<int> index_group_cm_;        // index of group with cm for integration

        // particle prediction, force, time array
        COMM::List<Tparticle> pred_; // predictor
        COMM::List<ForceH4> force_;  // force
        COMM::List<Float> time_next_; // next integrated time of particles

        // table of mask to show which particles are removed from the integration
        COMM::List<int> index_group_mask_;   // index list to record the masked groups
        COMM::List<bool> table_group_mask_; // bool mask to indicate whether the particle of index (group) (same index of table) is masked (true) or used (false)
        COMM::List<bool> table_single_mask_; // bool mask to indicate whether the particle of index (single) (same index of table) is masked (true) or used (false)

    public:
        BlockTimeStep4th step; ///> time step calculator
        HermiteManager<Tacc>* manager; ///< integration manager
        AR::TimeTransformedSymplecticManager<TARacc>* ar_manager; ///< integration manager
        COMM::ParticleGroup<H4Ptcl, Tpcm> particles; // particles
        COMM::List<ARSym> groups; // integrator for sub-groups
        COMM::List<Neighbor<Tparticle>> neighbors; // neighbor information of particles
        Tpert perturber; // external perturber
        Tinfo info; ///< information of the system
        Profile profile; // profile to measure the status

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            ASSERT(manager!=NULL);
            ASSERT(manager->checkParams());
            ASSERT(ar_manager!=NULL);
            ASSERT(ar_manager->checkParams());
            ASSERT(perturber.checkParams());
            ASSERT(info.checkParams());
            return true;
        }

        //! constructor
        HermiteIntegrator(): time_(0.0), time_offset_(0.0), time_next_min_(0.0), dt_limit_(0.0),
                             energy_init_ref_(0.0),
                             energy_(), energy_sd_(),
                             n_act_single_(0), n_act_group_(0), 
                             n_init_single_(0), n_init_group_(0), 
                             index_offset_group_(0), 
                             interrupt_group_dt_sorted_group_index_(-1), interrupt_binary_(), index_group_merger_(),
                             initial_system_flag_(false), modify_system_flag_(false),
                             index_dt_sorted_single_(), index_dt_sorted_group_(), 
                             index_group_resolve_(), index_group_cm_(), 
                             pred_(), force_(), time_next_(), 
                             index_group_mask_(), table_group_mask_(), table_single_mask_(), step(),
                             manager(NULL), ar_manager(NULL), particles(), groups(), neighbors(), perturber(), info(), profile() {}

        //! clear function
        void clear() {
            time_ = 0.0;
            time_offset_ = 0.0;
            time_next_min_ = 0.0;
            dt_limit_ = 0.0;
            energy_init_ref_ = 0.0;
            energy_.clear();
            energy_sd_.clear();
            n_act_single_ = n_act_group_ = n_init_single_ = n_init_group_ = 0;
            index_offset_group_ = 0;
            interrupt_group_dt_sorted_group_index_ = -1;
            interrupt_binary_.clear();
            index_group_merger_.clear();
            initial_system_flag_ = false;
            modify_system_flag_ = false;

            particles.clear();
            groups.clear();
            index_dt_sorted_single_.clear();
            index_dt_sorted_group_.clear();
            index_group_resolve_.clear();
            index_group_cm_.clear();
            index_group_mask_.clear();
            table_group_mask_.clear();
            table_single_mask_.clear();
            pred_.clear();
            force_.clear();
            time_next_.clear();
            neighbors.clear();
            perturber.clear();
            info.clear();
            profile.clear();
        }

        //! reserve memory for system
        /*! The memory size depends in the particles and groups memory sizes, thus particles and groups are reserved first
         */
        void reserveIntegratorMem() {
            // reserve memory
            const int nmax = particles.getSizeMax();
            const int nmax_group = groups.getSizeMax();
            const int nmax_tot = nmax + nmax_group;
            ASSERT(nmax>0);

            index_dt_sorted_single_.setMode(COMM::ListMode::local);
            index_dt_sorted_group_.setMode(COMM::ListMode::local);
            index_group_merger_.setMode(COMM::ListMode::local);
            index_group_resolve_.setMode(COMM::ListMode::local);
            index_group_cm_.setMode(COMM::ListMode::local);
            index_group_mask_.setMode(COMM::ListMode::local);
            table_group_mask_.setMode(COMM::ListMode::local);
            table_single_mask_.setMode(COMM::ListMode::local);
            pred_.setMode(COMM::ListMode::local);
            force_.setMode(COMM::ListMode::local);
            time_next_.setMode(COMM::ListMode::local);
            neighbors.setMode(COMM::ListMode::local);

            index_group_merger_.reserveMem(nmax_group);
            
            index_dt_sorted_single_.reserveMem(nmax);
            index_dt_sorted_group_.reserveMem(nmax_group);

            index_group_resolve_.reserveMem(nmax_group);
            index_group_cm_.reserveMem(nmax_group);

            index_group_mask_.reserveMem(nmax_group);
            table_group_mask_.reserveMem(nmax_group);
            table_single_mask_.reserveMem(nmax);

            // total allocation size includes both single and group c.m. size
            index_offset_group_ = nmax;
            pred_.reserveMem(nmax_tot);
            force_.reserveMem(nmax_tot);
            time_next_.reserveMem(nmax_tot);

            // reserve for neighbor list
            neighbors.reserveMem(nmax);
            auto* nb_ptr = neighbors.getDataAddress();
            for (int i=0; i<nmax; i++) {
                nb_ptr[i].neighbor_address.setMode(COMM::ListMode::local);
                nb_ptr[i].neighbor_address.reserveMem(nmax_tot);
            }
        }

    private:
        //! Calculate 2nd order time step for lists particles 
        /*! Calculate 2nd order time step 
          @param[in] _index_single: active particle index for singles
          @param[in] _n_single: number of active singles
          @param[in] _index_group: active particle index for groups
          @param[in] _n_group: number of active groups
          @param[in] _dt_limit: time step limit
        */
        void calcDt2ndList(const int* _index_single,
                           const int _n_single,
                           const int* _index_group,
                           const int _n_group,
                           const Float _dt_limit) {
            ASSERT(manager!=NULL);
            // for single
            for (int i=0; i<_n_single; i++){
                const int k = _index_single[i];
                const Float* acc0 = force_[k].acc0;
                const Float* acc1 = force_[k].acc1;
                const Float dt =step.calcBlockDt2nd(acc0, acc1, _dt_limit);
                particles[k].dt = dt;
                neighbors[k].initial_step_flag = false;
            }
            // for group
            for (int i=0; i<_n_group; i++) {
                const int k = _index_group[i];
                const int kf = k + index_offset_group_;
                const Float* acc0 = force_[kf].acc0;
                const Float* acc1 = force_[kf].acc1;
                const Float dt =step.calcBlockDt2nd(acc0, acc1, std::min(_dt_limit, groups[k].info.dt_limit));
                groups[k].particles.cm.dt = dt;
                groups[k].perturber.initial_step_flag = false;
            }
        }



        //! predict particles to the time
        /*! @param[in] _time_pred: time for prediction
         */
        void predictAll(const Float _time_pred) {
            static thread_local const Float inv3 = 1.0 / 3.0;
            // single
            const int n_single = index_dt_sorted_single_.getSize();
            ASSERT(n_single <= particles.getSize());
            const auto* ptcl = particles.getDataAddress();
            auto* pred = pred_.getDataAddress();
            for (int k=0; k<n_single; k++){
                const int i = index_dt_sorted_single_[k];
                const Float dt = _time_pred - ptcl[i].time;
                ASSERT(i<particles.getSize());
                ASSERT(i>=0);
                pred[i].pos[0] = ptcl[i].pos[0] + dt*(ptcl[i].vel[0]  + 0.5*dt*(ptcl[i].acc0[0] + inv3*dt*ptcl[i].acc1[0]));
                pred[i].pos[1] = ptcl[i].pos[1] + dt*(ptcl[i].vel[1]  + 0.5*dt*(ptcl[i].acc0[1] + inv3*dt*ptcl[i].acc1[1]));
                pred[i].pos[2] = ptcl[i].pos[2] + dt*(ptcl[i].vel[2]  + 0.5*dt*(ptcl[i].acc0[2] + inv3*dt*ptcl[i].acc1[2]));

                pred[i].vel[0] = ptcl[i].vel[0] + dt*(ptcl[i].acc0[0] + 0.5*dt*ptcl[i].acc1[0]);
                pred[i].vel[1] = ptcl[i].vel[1] + dt*(ptcl[i].acc0[1] + 0.5*dt*ptcl[i].acc1[1]);
                pred[i].vel[2] = ptcl[i].vel[2] + dt*(ptcl[i].acc0[2] + 0.5*dt*ptcl[i].acc1[2]);

                pred[i].mass = ptcl[i].mass;
            }
            // group
            const int n_group = index_dt_sorted_group_.getSize();
            ASSERT(n_group <= groups.getSize());
            auto* group_ptr = groups.getDataAddress();
            for (int k=0; k<n_group; k++) {
                const int i = index_dt_sorted_group_[k];
                ASSERT(i<groups.getSize());
                ASSERT(i>=0);
                const auto& pcm = group_ptr[i].particles.cm;
                const Float dt = _time_pred - pcm.time;
                ASSERT(dt>=0.0);
                // group predictor is after single predictor with offset n_single
                auto& predcm = pred[i+index_offset_group_];
                predcm.pos[0] = pcm.pos[0] + dt*(pcm.vel[0]  + 0.5*dt*(pcm.acc0[0] + inv3*dt*pcm.acc1[0]));
                predcm.pos[1] = pcm.pos[1] + dt*(pcm.vel[1]  + 0.5*dt*(pcm.acc0[1] + inv3*dt*pcm.acc1[1]));
                predcm.pos[2] = pcm.pos[2] + dt*(pcm.vel[2]  + 0.5*dt*(pcm.acc0[2] + inv3*dt*pcm.acc1[2]));

                predcm.vel[0] = pcm.vel[0] + dt*(pcm.acc0[0] + 0.5*dt*pcm.acc1[0]);
                predcm.vel[1] = pcm.vel[1] + dt*(pcm.acc0[1] + 0.5*dt*pcm.acc1[1]);
                predcm.vel[2] = pcm.vel[2] + dt*(pcm.acc0[2] + 0.5*dt*pcm.acc1[2]);

                predcm.mass = pcm.mass;
            }
        }

        //! correct particle and calculate step 
        /*! Correct particle and calculate next time step
          @param[in] _pi: particle to corect
          @param[in] _fi: force for correction
          @param[in] _dt_limit: maximum step size allown
          @param[in] _inti_step_flag: true: only calculate dt 2nd instead of 4th
         */
        void correctAndCalcDt4thOne(H4Ptcl& _pi, ForceH4& _fi, const Float _dt_limit, const bool _init_step_flag) {
            static thread_local const Float inv3 = 1.0 / 3.0;
            const Float dt = _pi.dt;
            const Float h = 0.5 * dt;
            const Float hinv = 2.0 / dt;
            const Float A0p[3] = {_fi.acc0[0] + _pi.acc0[0], _fi.acc0[1] + _pi.acc0[1], _fi.acc0[2] + _pi.acc0[2]};
            const Float A0m[3] = {_fi.acc0[0] - _pi.acc0[0], _fi.acc0[1] - _pi.acc0[1], _fi.acc0[2] - _pi.acc0[2]};
            const Float A1p[3] = {(_fi.acc1[0] + _pi.acc1[0])*h, (_fi.acc1[1] + _pi.acc1[1])*h, (_fi.acc1[2] + _pi.acc1[2])*h};
            const Float A1m[3] = {(_fi.acc1[0] - _pi.acc1[0])*h, (_fi.acc1[1] - _pi.acc1[1])*h, (_fi.acc1[2] - _pi.acc1[2])*h};

            const Float vel_new[3] = {_pi.vel[0] + h*( A0p[0] - inv3*A1m[0] ), 
                                      _pi.vel[1] + h*( A0p[1] - inv3*A1m[1] ), 
                                      _pi.vel[2] + h*( A0p[2] - inv3*A1m[2] )};
            _pi.pos[0] += h*( (_pi.vel[0] + vel_new[0]) + h*(-inv3*A0m[0]));
            _pi.pos[1] += h*( (_pi.vel[1] + vel_new[1]) + h*(-inv3*A0m[1]));
            _pi.pos[2] += h*( (_pi.vel[2] + vel_new[2]) + h*(-inv3*A0m[2]));

            _pi.vel[0] = vel_new[0];
            _pi.vel[1] = vel_new[1];
            _pi.vel[2] = vel_new[2];

            _pi.acc0[0] = _fi.acc0[0];
            _pi.acc0[1] = _fi.acc0[1];
            _pi.acc0[2] = _fi.acc0[2];

            _pi.acc1[0] = _fi.acc1[0];
            _pi.acc1[1] = _fi.acc1[1];
            _pi.acc1[2] = _fi.acc1[2];

            _pi.time += dt;

            const Float acc3[3] = {(1.5*hinv*hinv*hinv) * (A1p[0] - A0m[0]),
                                   (1.5*hinv*hinv*hinv) * (A1p[1] - A0m[1]),
                                   (1.5*hinv*hinv*hinv) * (A1p[2] - A0m[2])};
            const Float acc2[3] = {(0.5*hinv*hinv) * A1m[0] + h*acc3[0], 
                                   (0.5*hinv*hinv) * A1m[1] + h*acc3[1], 
                                   (0.5*hinv*hinv) * A1m[2] + h*acc3[2]};

#ifdef HERMITE_DEBUG_ACC
            _pi.acc2[0] = acc2[0];
            _pi.acc2[1] = acc2[1];
            _pi.acc2[2] = acc2[2];

            _pi.acc3[0] = acc3[0];
            _pi.acc3[1] = acc3[1];
            _pi.acc3[2] = acc3[2];
#endif

            _pi.pot = _fi.pot;

            const Float dt_old = _pi.dt;
            if(_init_step_flag) {
                _pi.dt = step.calcBlockDt2nd(_pi.acc0, _pi.acc1, _dt_limit);
#ifdef HERMITE_DEBUG
                std::cerr<<"Initial step flag on: pi.id: "<<_pi.id<<" step size: "<<_pi.dt<<" time: "<<_pi.time<<std::endl;
#endif                
            }
            else _pi.dt = step.calcBlockDt4th(_pi.acc0, _pi.acc1, acc2, acc3, _dt_limit);

            ASSERT((dt_old > 0.0 && _pi.dt >0.0));
            (void)dt_old;
        }

        //! correct particle and calculate step 
        /*! Correct particle and calculate next time step
          @param[in] _index_single: active particle index for singles
          @param[in] _n_single: number of active singles
          @param[in] _index_group: active particle index for groups
          @param[in] _n_group: number of active groups
          @param[in] _dt_limit: time step maximum
        */
        void correctAndCalcDt4thList(const int* _index_single,
                                     const int _n_single,
                                     const int* _index_group,
                                     const int _n_group,
                                     const Float _dt_limit) {
            ASSERT(_n_single<=index_dt_sorted_single_.getSize());
            ASSERT(_n_group<=index_dt_sorted_group_.getSize());

            // single
            auto* ptcl = particles.getDataAddress();
            ForceH4* force = force_.getDataAddress();
            for(int i=0; i<_n_single; i++){
                const int k = _index_single[i];
                correctAndCalcDt4thOne(ptcl[k], force[k], _dt_limit, neighbors[k].initial_step_flag);
                neighbors[k].initial_step_flag = false;
            }
            // group
            auto* group_ptr = groups.getDataAddress();
            for(int i=0; i<_n_group; i++) {
                const int k = _index_group[i];
                auto& groupi = group_ptr[k];
                correctAndCalcDt4thOne(groupi.particles.cm, force[k+index_offset_group_], std::min(_dt_limit, groupi.info.dt_limit), groupi.perturber.initial_step_flag);
                groupi.perturber.initial_step_flag = false;
                //groupi.correctCenterOfMassDrift();
            }
        }

        //! check whether group need to resolve
        /*! @param[in] _n_group: number of groups in index_dt_sorted_group_ to check
         */
        void checkGroupResolve(const int _n_group) {
            if(_n_group==0) return;
            //AR::SlowDown sd;
            ASSERT(_n_group<=index_dt_sorted_group_.getSize());
            for (int i=0; i<_n_group; i++) {
                const int k = index_dt_sorted_group_[i];
                auto& groupk = groups[k];
                // remove soft perturbation 
                //sd.initialSlowDownReference(groupk.slowdown.getSlowDownFactorReference(), groupk.slowdown.getSlowDownFactorMax());
                //const Float pert_out = std::max(0.0, groupk.slowdown.getPertOut()-groupk.perturber.soft_pert_min);
                //const Float pert_in= groupk.slowdown.getPertIn();
                //const Float kappa = sd.calcSlowDownFactor(pert_in, pert_out);
                // calculate non soft perturbation
                groupk.perturber.checkGroupResolve();
            }
        }

        //! Generate j particle list from groups for force calculation
        /*! check group status, if the components are need to resolve, writeback the data with slowdown velocity, save the index to _index_group_resolve list, otherwisde save to _index_group_cm list
          @param[in]: _write_back_flag: if true, write back slowdown particles from groups, otherwise only record group index
         */
        void writeBackResolvedGroupAndCreateJParticleList(const bool _write_back_flag) {
            index_group_resolve_.resizeNoInitialize(0);
            index_group_cm_.resizeNoInitialize(0);
            const int n_group_org = index_dt_sorted_group_.getSize();
            auto* group_ptr = groups.getDataAddress();
            auto* ptcl = pred_.getDataAddress();
            for (int i=0; i<n_group_org; i++) {
                int k = index_dt_sorted_group_[i];
                ASSERT(k<groups.getSize());
                ASSERT(k>=0);
                auto& group_k = group_ptr[k];
                if (group_k.perturber.need_resolve_flag) {
                    if (_write_back_flag) group_k.writeBackSlowDownParticles(ptcl[k+index_offset_group_]);
                    index_group_resolve_.addMember(k);
                }
                else index_group_cm_.addMember(k);
            }
        }

        //! Insert particle index in an existed or new group
        /*! Used for checkNewGroup, if particle _i, _j are already in new groups found before, insert _i or _j into this group, otherwise make new group
          @param[in] _i: particle index _i
          @param[in] _j: particle index _j
          @param[in,out] _used_mask: mask table recored whether index is used in which group
          @param[in,out] _new_group_particle_index_origin: particle index (for hermite particle data) of new group members
          @param[in,out] _new_n_group_offset:  group member boundary, first value is defined already 
          @param[in,out] _new_n_particle: total number of new member particles 
          @param[in,out] _new_n_group:    number of new group
          @param[in,out] _break_group_index_with_offset:   break group index in groups
          @param[out] _n_break_no_add: number of break groups without adding new particles
          @param[in] _n_break:   number of break groups
         */
        void insertParticleIndexToGroup(const int _i, 
                                        const int _j, 
                                        int* _used_mask, 
                                        int* _new_group_particle_index_origin,
                                        int* _new_n_group_offset, 
                                        int& _new_n_particle,
                                        int& _new_n_group,
                                        int* _break_group_index_with_offset,
                                        int& _n_break_no_add,
                                        const int _n_break) {
            int insert_group=-1, insert_index=-1;
            if(_used_mask[_i]>=0) {
                insert_group = _used_mask[_i];
                insert_index = _j;
            }
            else if(_used_mask[_j]>=0) {
                insert_group = _used_mask[_j];
                insert_index = _i;
            }
            // the case of merging group
            if(insert_group>=0) {
                // check insert particle number
                int insert_n = 1;
                if (insert_index>=index_offset_group_) insert_n = groups[insert_index-index_offset_group_].particles.getSize();
                
                // shift the first insert_n index in last group to last 
                int last_group_offset = _new_n_group_offset[_new_n_group-1];
                for (int j=0; j<insert_n; j++) 
                    _new_group_particle_index_origin[_new_n_particle++] = _new_group_particle_index_origin[last_group_offset+j];

                // in the case insert_group is not the last group, shift the offset and also the first insert_n index of each group into right group.
                if(insert_group<_new_n_group-1) {
                    for (int k=_new_n_group-1; k>insert_group; k--) {
                        int k0_group_offset = _new_n_group_offset[k-1];
                        for (int j=0; j<insert_n; j++) 
                            _new_group_particle_index_origin[_new_n_group_offset[k]++] = _new_group_particle_index_origin[k0_group_offset+j];
                    }
                }
                // replace the first position of insert_group with insert ptcl
                int insert_group_offset=_new_n_group_offset[insert_group];
                if(insert_n>1) {
                    auto& group = groups[insert_index-index_offset_group_];
                    for (int j=0; j<insert_n; j++) 
                        _new_group_particle_index_origin[insert_group_offset++] = group.info.particle_index[j];
                    // add group index to break list
                    _break_group_index_with_offset[_n_break+_n_break_no_add] = insert_index;
                    _n_break_no_add++;
                }
                else _new_group_particle_index_origin[insert_group_offset] = insert_index;
                _used_mask[insert_index] = insert_group;
            }
            else {   // new group case
                _new_n_group_offset[_new_n_group] = _new_n_particle;
                if (_i>=index_offset_group_) {
                    auto& group = groups[_i-index_offset_group_];
                    int insert_n = group.particles.getSize();
                    for (int j=0; j<insert_n; j++)
                        _new_group_particle_index_origin[_new_n_particle++] = group.info.particle_index[j];
                    // add group index to break list
                    _break_group_index_with_offset[_n_break+_n_break_no_add] = _i;
                    _n_break_no_add++;
                }
                else _new_group_particle_index_origin[_new_n_particle++] = _i;
                // used mask regist the current staying group of particle _i
                _used_mask[_i] = _new_n_group;

                if (_j>=index_offset_group_) {
                    auto& group = groups[_j-index_offset_group_];
                    int insert_n = group.particles.getSize();
                    for (int j=0; j<insert_n; j++)
                        _new_group_particle_index_origin[_new_n_particle++] = group.info.particle_index[j];
                    // add group index to break list
                    _break_group_index_with_offset[_n_break+_n_break_no_add] = _j;
                    _n_break_no_add++;
                }
                else _new_group_particle_index_origin[_new_n_particle++] = _j;
                // used mask regist the current staying group of particle _j
                _used_mask[_j] = _new_n_group;
  
                _new_n_group++;
            }
        }

        //! calculate one particle interaction from all singles and groups
        /*! Neighbor information is also updated
          @param[out] _fi: acc and jerk of particle i
          @param[out] _nbi: neighbor information of i
          @param[in] _pi: i particle
          @param[in] _pid: i particle id
        */
        template <class Tpi>
        inline void calcOneSingleAccJerkNB(ForceH4 &_fi, 
                                           Neighbor<Tparticle> &_nbi,
                                           const Tpi &_pi,
                                           const int _pid) {
            // clear force
            _fi.clear();
            _nbi.resetNeighbor();

            // single list
            const int* single_list = index_dt_sorted_single_.getDataAddress();
            const int n_single = index_dt_sorted_single_.getSize();
            auto* ptcl = pred_.getDataAddress();
            for (int i=0; i<n_single; i++) {
                const int j = single_list[i];
                ASSERT(j<pred_.getSize());
                const auto& pj = ptcl[j];
                if (_pid==pj.id) continue;
                if (pj.mass==0) continue;
                Float r2 = manager->interaction.calcAccJerkPairSingleSingle(_fi, _pi, pj);
                ASSERT(r2>0.0);
                _nbi.checkAndAddNeighborSingle(r2, particles[j], neighbors[j], j);
            }

            auto* group_ptr = groups.getDataAddress();

            // resolved group list
            const int n_group_resolve = index_group_resolve_.getSize();
            for (int i=0; i<n_group_resolve; i++) {
                const int j =index_group_resolve_[i];
                auto& groupj = group_ptr[j];
                if (_pid==groupj.particles.cm.id) continue;
                if (groupj.particles.cm.mass==0) continue;
                Float r2 = manager->interaction.calcAccJerkPairSingleGroupMember(_fi, _pi, groupj);
                ASSERT(r2>0.0);
                _nbi.checkAndAddNeighborGroup(r2, groupj, j+index_offset_group_);
            }

            // cm group list
            const int n_group_cm = index_group_cm_.getSize();
            for (int i=0; i<n_group_cm; i++) {
                const int j = index_group_cm_[i];
                auto& groupj = group_ptr[j];
                if (_pid==groupj.particles.cm.id) continue;
                if (groupj.particles.cm.mass==0) continue;
                // used predicted particle instead of original cm
                const auto& pj = ptcl[j+index_offset_group_];
                ASSERT(j+index_offset_group_<pred_.getSize());
                ASSERT(_pi.id!=pj.id);
                Float r2 = manager->interaction.calcAccJerkPairSingleGroupCM(_fi, _pi, groupj, pj);
                ASSERT(r2>0.0);
                _nbi.checkAndAddNeighborGroup(r2, groupj, j+index_offset_group_);
            }
            ASSERT(_nbi.n_neighbor_group + _nbi.n_neighbor_single == _nbi.neighbor_address.getSize());
        }

        //! calculate one cm group interaction from all singles and groups
        /*! Neighbor information is also updated
          @param[out] _fi: acc and jerk of particle i
          @param[out] _groupi: group i 
          @param[in] _pi: predicted group i cm 
        */
        template <class Tpi, class Tgroupi>
        inline void calcOneGroupCMAccJerkNB(ForceH4 &_fi, 
                                             Tgroupi &_groupi, 
                                             const Tpi &_pi) {
            // clear force
            _fi.clear();
            auto& nbi = _groupi.perturber;
            nbi.resetNeighbor();

            // single list
            const int* single_list = index_dt_sorted_single_.getDataAddress();
            const int n_single = index_dt_sorted_single_.getSize();
            auto* ptcl = pred_.getDataAddress();
            for (int i=0; i<n_single; i++) {
                const int j = single_list[i];
                ASSERT(j<pred_.getSize());
                const auto& pj = ptcl[j];
                if (_pi.id==pj.id) continue;
                if (pj.mass==0.0) continue;
                Float r2 = manager->interaction.calcAccJerkPairGroupCMSingle(_fi, _groupi, _pi, pj);
                ASSERT(r2>0.0);
                nbi.checkAndAddNeighborSingle(r2, particles[j], neighbors[j], j);
            }

            auto* group_ptr = groups.getDataAddress();

            // resolved group list
            const int n_group_resolve = index_group_resolve_.getSize();
            for (int i=0; i<n_group_resolve; i++) {
                const int j =index_group_resolve_[i];
                auto& groupj = group_ptr[j];
                if(_pi.id==groupj.particles.cm.id) continue;
                if (groupj.particles.cm.mass==0.0) continue;
                Float r2 = manager->interaction.calcAccJerkPairGroupCMGroupMember(_fi, _groupi, _pi, groupj);
                ASSERT(r2>0.0);
                nbi.checkAndAddNeighborGroup(r2, groupj, j+index_offset_group_);
            }

            // cm group list
            const int n_group_cm = index_group_cm_.getSize();
            for (int i=0; i<n_group_cm; i++) {
                const int j = index_group_cm_[i];
                auto& groupj = group_ptr[j];
                if (_pi.id==groupj.particles.cm.id) continue;
                if (groupj.particles.cm.mass==0.0) continue;
                // used predicted particle instead of original cm
                const auto& pj = ptcl[j+index_offset_group_];
                ASSERT(j+index_offset_group_<pred_.getSize());
                Float r2 = manager->interaction.calcAccJerkPairGroupCMGroupCM(_fi, _groupi, _pi, groupj, pj);
                ASSERT(r2>0.0);
                nbi.checkAndAddNeighborGroup(r2, groupj, j+index_offset_group_);
            }
            ASSERT(nbi.n_neighbor_group + nbi.n_neighbor_single == nbi.neighbor_address.getSize());
        }        

        //! calculate one resolved group interaction from all singles and groups
        /*! Neighbor information is also updated
          @param[out] _fi: acc and jerk of particle i
          @param[out] _groupi: group i 
          @param[in] _pi: predicted group i cm 
        */
        template <class Tpi, class Tgroupi>
        inline void calcOneGroupMemberAccJerkNB(ForceH4 &_fi, 
                                                Tgroupi &_groupi, 
                                                const Tpi &_pi) {
            auto& nbi = _groupi.perturber;
            nbi.resetNeighbor();
            
            // only get neighbors
            // single list
            const int* single_list = index_dt_sorted_single_.getDataAddress();
            const int n_single = index_dt_sorted_single_.getSize();
            auto* ptcl = pred_.getDataAddress();
            for (int i=0; i<n_single; i++) {
                const int j = single_list[i];
                ASSERT(j<pred_.getSize());
                const auto& pj = ptcl[j];
                if (pj.mass==0.0) continue;
                ASSERT(_pi.id!=pj.id);
                Float r2 = manager->interaction.calcR2Pair(_pi, pj);
                ASSERT(r2>0.0);
                nbi.checkAndAddNeighborSingle(r2, particles[j], neighbors[j], j);
            }

            auto* group_ptr = groups.getDataAddress();

            // resolved group list
            const int n_group_resolve = index_group_resolve_.getSize();
            for (int i=0; i<n_group_resolve; i++) {
                const int j =index_group_resolve_[i];
                // used predicted particle instead of original cm
                const auto& pj = ptcl[j+index_offset_group_];
                auto& groupj = group_ptr[j];
                if(_pi.id==pj.id) continue;
                if (pj.mass==0.0) continue;
                Float r2 = manager->interaction.calcR2Pair(_pi, pj);
                ASSERT(r2>0.0);
                nbi.checkAndAddNeighborGroup(r2, groupj, j+index_offset_group_);
            }

            // cm group list
            const int n_group_cm = index_group_cm_.getSize();
            for (int i=0; i<n_group_cm; i++) {
                const int j = index_group_cm_[i];
                auto& groupj = group_ptr[j];
                // used predicted particle instead of original cm
                const auto& pj = ptcl[j+index_offset_group_];
                if (pj.mass==0.0) continue;
                ASSERT(_pi.id!=pj.id);
                ASSERT(j+index_offset_group_<pred_.getSize());
                Float r2 = manager->interaction.calcR2Pair(_pi, pj);
                ASSERT(r2>0.0);
                nbi.checkAndAddNeighborGroup(r2, groupj, j+index_offset_group_);
            }
            ASSERT(nbi.n_neighbor_group + nbi.n_neighbor_single == nbi.neighbor_address.getSize());
            
            // clear force
            _fi.clear();
            auto* force_ptr = force_.getDataAddress();
            auto* neighbor_ptr = neighbors.getDataAddress();
            auto* member_adr = _groupi.particles.getOriginAddressArray();
            const int n_member = _groupi.particles.getSize();
#ifdef HERMITE_DEBUG
            Float mcm = 0.0;
#endif
            for (int j=0; j<n_member; j++) {
                // get particle index from memory address offset
                //const int kj = int((Tparticle*)member_adr[j]-ptcl);
                const int kj = _groupi.info.particle_index[j];
                ASSERT(kj<particles.getSize());
                ASSERT(kj>=0);
                auto& pkj = *(Tparticle*)member_adr[j];
                ASSERT((void*)member_adr[j]==(void*)&particles[kj]);
                auto& fkj = force_ptr[kj];
                auto& nbkj = neighbor_ptr[kj];
                calcOneSingleAccJerkNB(fkj, nbkj, pkj, _pi.id);
                // replace the cm. force with the summation of members
                _fi.acc0[0] += pkj.mass * fkj.acc0[0]; 
                _fi.acc0[1] += pkj.mass * fkj.acc0[1]; 
                _fi.acc0[2] += pkj.mass * fkj.acc0[2]; 

                _fi.acc1[0] += pkj.mass * fkj.acc1[0]; 
                _fi.acc1[1] += pkj.mass * fkj.acc1[1]; 
                _fi.acc1[2] += pkj.mass * fkj.acc1[2]; 

                _fi.pot += pkj.mass * fkj.pot;
#ifdef HERMITE_DEBUG
                mcm += pkj.mass;
#endif
            }
#ifdef HERMITE_DEBUG
            ASSERT(abs(mcm-_pi.mass)<1e-10);
#endif
            _fi.acc0[0] /= _pi.mass;
            _fi.acc0[1] /= _pi.mass;
            _fi.acc0[2] /= _pi.mass;

            _fi.acc1[0] /= _pi.mass;
            _fi.acc1[1] /= _pi.mass;
            _fi.acc1[2] /= _pi.mass;

            _fi.pot /= _pi.mass;
            
        }
        
        //! Calculate acc and jerk for lists of particles from all particles and update neighbor lists
        /*! 
          @param[in] _index_single: active particle index for singles
          @param[in] _n_single: number of active singles
          @param[in] _index_group: active particle index for groups
          @param[in] _n_group: number of active groups
        */
        inline void calcAccJerkNBList(const int* _index_single,
                                      const int  _n_single,
                                      const int* _index_group,
                                      const int  _n_group) {

            // predictor 
            //auto* ptcl = particles.getDataAddress();
            auto* pred_ptr = pred_.getDataAddress();
            auto* force_ptr = force_.getDataAddress();
            auto* neighbor_ptr = neighbors.getDataAddress();
            for (int k=0; k<_n_single; k++) {
                const int i = _index_single[k];
                auto& pi = pred_ptr[i];
                auto& fi = force_ptr[i];
                auto& nbi = neighbor_ptr[i];
                calcOneSingleAccJerkNB(fi, nbi, pi, pi.id);
            }

            // for group active particles
            auto* group_ptr = groups.getDataAddress();
            
            for (int k=0; k<_n_group; k++) {
                const int i = _index_group[k];
                auto& groupi = group_ptr[i];
                // use predictor of cm
                auto& pi = pred_ptr[i+index_offset_group_];
                auto& fi = force_ptr[i+index_offset_group_];
                if (groupi.particles.cm.mass>0) {
                    if (groupi.perturber.need_resolve_flag) calcOneGroupMemberAccJerkNB(fi, groupi, pi);
                    else calcOneGroupCMAccJerkNB(fi, groupi, pi);
                }
                else calcOneSingleAccJerkNB(fi, groupi.perturber, pi, pi.id);
            }
        }

        //! reduce one group cm step by half and sort index dt group 
        void reduceGroupCMStepByHalfAndSortDtIndex(const int _index) {
            int i = _index;
            int k = index_dt_sorted_group_[i];
            
            auto& pcm = groups[k].particles.cm;
            // reduce step size to get one more step
            pcm.dt *= 0.5;
            time_next_[k+index_offset_group_] = pcm.time + pcm.dt;
            if (_index>0) {
                int kp = index_dt_sorted_group_[i-1];
                while (time_next_[k+index_offset_group_]<time_next_[kp+index_offset_group_]) {
                    // swap index
                    int tmp=index_dt_sorted_group_[i-1];
                    index_dt_sorted_group_[i-1] = index_dt_sorted_group_[i];
                    index_dt_sorted_group_[i] = tmp;

                    // update i to the new position of k
                    i--;
                    if (i==0) break; // if moved to the begining, break
                    else kp = index_dt_sorted_group_[i-1];
                }
            }
            time_next_min_ = std::min(time_next_min_, time_next_[k+index_offset_group_]);

            // in case the particle move to the boundary of act group, increase n_act_group by one
            if (_index>=n_act_group_ && i<=n_act_group_ && n_act_group_<index_dt_sorted_group_.getSize()) n_act_group_++;
        }
        

        //! update time next 
        /*! update time_next array for a list of particles
          @param[in] _index_single: active particle index for singles
          @param[in] _n_single: number of active singles
          @param[in] _index_group: active particle index for groups
          @param[in] _n_group: number of active groups
        */
        void updateTimeNextList(const int* _index_single,
                                const int _n_single,
                                const int* _index_group,
                                const int _n_group) {
            // for single
            for (int i=0; i<_n_single; i++){
                const int k = _index_single[i];
                ASSERT(k<time_next_.getSize());
                time_next_[k] = particles[k].time + particles[k].dt;
            }

            // for group
            for (int i=0; i<_n_group; i++) {
                const int k = _index_group[i];
                const int kf = k + index_offset_group_;
                ASSERT(kf<time_next_.getSize());
                auto& pcm = groups[k].particles.cm;
                time_next_[kf] = pcm.time + pcm.dt;
            }
        }
            
        //! calculate dr dv of a pair
        /*!
          @param[in] _p1: particle 1
          @param[in] _p2: particle 2
         */
        Float calcDrDv(const H4Ptcl& _p1, const H4Ptcl& _p2) {
            Float dx[3],dv[3];
            dx[0] = _p1.pos[0] - _p2.pos[0];
            dx[1] = _p1.pos[1] - _p2.pos[1];
            dx[2] = _p1.pos[2] - _p2.pos[2];

            dv[0] = _p1.vel[0] - _p2.vel[0];
            dv[1] = _p1.vel[1] - _p2.vel[1];
            dv[2] = _p1.vel[2] - _p2.vel[2];
        
            Float drdv= dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];

            return drdv;
        }

        // for get sorted index of single
        class SortIndexDtSingle{
        private:
            Float * time_;
        public:
            SortIndexDtSingle(Float* _time): time_(_time) {}
            bool operator() (const int & left, const int & right) const {
                return time_[left] < time_[right];
            }
        };

        // for get sorted index of group
        class SortIndexDtGroup{
        private:
            Float * time_;
            int offset_;
        public:
            SortIndexDtGroup(Float * _time, int _offset): time_(_time), offset_(_offset) {}
            bool operator() (const int & left, const int & right) const {
                return time_[left+offset_] < time_[right+offset_];
            }
        };
        
        // ! remove index from bool table
        void removeDtIndexFromTable(COMM::List<int>& _index, int& _n_act, int& _n_init, COMM::List<bool>& _table) {
            int i = 0;
            int imove = 0;
            int n_index = _index.getSize();
            while (i<n_index) {
                if (_table[_index[i]]) {
                    if (i<_n_act) _n_act--;
                    if (i<_n_init) _n_init--;
                    imove++;
                    n_index--;
                }
                else i++;
                if (i>=n_index) break;
                _index[i] = _index[i+imove];
            }
            ASSERT(_n_act>=0);
            ASSERT(_n_act<=n_index);
            ASSERT(_n_init>=0);
            ASSERT(_n_init<=n_index);
            _index.resizeNoInitialize(n_index);
        }

        // ! remove group particle address of a list from table_group_mask_
        /*!
          \return remove number 
         */
        int removeNBGroupAddressFromTable(COMM::List<NBAdr<Tparticle>> & _adr) {
            int k = 0;
            int k_last = _adr.getSize()-1;
            int n_rm = 0;
            while (k<=k_last) {
                if (_adr[k].index>=index_offset_group_) {
                    const int kg = _adr[k].index-index_offset_group_;
                    if (table_group_mask_[kg]) {
                        _adr[k] = _adr[k_last];
                        k_last--;
                        n_rm++;
                    }
                    else k++;
                }
                else k++;
            }
            _adr.resizeNoInitialize(k_last+1);
            return n_rm;
        }

        // ! remove group particle address of a list from table_group_mask_
        /*!
          \return remove number 
         */
        int removeNBSingleAddressFromTable(COMM::List<NBAdr<Tparticle>> & _adr) {
            int k = 0;
            int k_last = _adr.getSize()-1;
            int n_rm = 0;
            while (k<=k_last) {
                if (_adr[k].index<index_offset_group_) {
                    if (table_single_mask_[_adr[k].index]) {
                        _adr[k] = _adr[k_last];
                        k_last--;
                        n_rm++;
                    }
                    else k++;
                }
                else k++;
            }
            _adr.resizeNoInitialize(k_last+1);
            return n_rm;
        }
        

    public:

        //! Initial system array
        /*!@param[in] _time_sys: current set time
         */
        void initialSystemSingle(const Float _time_sys) {
            ASSERT(checkParams());
            particles.setModifiedFalse();

            initial_system_flag_ = true;

            // start case, set n_init_single_ and n_init_group_ consistent with index_dt_sorted array
            // set particle numbers
            const int n_particle = particles.getSize();
            ASSERT(n_particle>1);
            
            // use increase since addgroups may already increase the sizes
            pred_.increaseSizeNoInitialize(n_particle);
            force_.increaseSizeNoInitialize(n_particle);
            time_next_.increaseSizeNoInitialize(n_particle);

            // set number of lists
            neighbors.resizeNoInitialize(n_particle);
            index_dt_sorted_single_.resizeNoInitialize(n_particle);
            table_single_mask_.resizeNoInitialize(n_particle);

            for (int i=0; i<n_particle; i++) {
                // initial index_dt_sorted_single_ first to allow initialIntegration work
                index_dt_sorted_single_[i] = i;
                table_single_mask_[i] = false;
            }
            // set initial time
            time_ = _time_sys;

            n_init_single_ = index_dt_sorted_single_.getSize();
            // if initial step, n_act are not given initially
            n_act_single_  = n_init_single_;
        }

        //! add groups from lists of particle index 
        /*! 
          @param[in] _particle_index: particle index
          @param[in] _n_group_offset: offset to separate groups in _particle_index
          @param[in] _n_group: number of new groups
         */
        /* Add particles from single (hermite.particles) to groups.
           Add new groups into index_dt_sorted_group_ and increase n_init_group_.
           set single mask to true for added particles
         */
        void addGroups(const int* _particle_index, const int* _n_group_offset, const int _n_group) {
            ASSERT(checkParams());
            ASSERT(!particles.isModified());
            ASSERT(initial_system_flag_);

            if (_n_group==0) return;

            modify_system_flag_ = true;

            // gether the new index in groups first
            int group_index[_n_group];
            for (int i=0; i<_n_group; i++) {
                // check suppressed group
                int igroup;
                if(index_group_mask_.getSize()>0) {
                    igroup = index_group_mask_.getLastMember();
                    table_group_mask_[igroup] = false;
                    index_group_mask_.decreaseSizeNoInitialize(1);
                }
                else {
                    // add new
                    igroup = groups.getSize();
                    table_group_mask_.addMember(false);
                    groups.increaseSizeNoInitialize(1);
                    // also increase the cm predictor force and time_next
                    pred_.increaseSizeNoInitialize(1);
                    force_.increaseSizeNoInitialize(1);
                    time_next_.increaseSizeNoInitialize(1);
                }
                group_index[i] = igroup;
            }

            // add new group
            for (int i=0; i<_n_group; i++) {
                auto& group_new = groups[group_index[i]];

                // set manager
                ASSERT(ar_manager!=NULL);
                group_new.manager = ar_manager;

                // allocate memory
                const int n_particle = _n_group_offset[i+1] - _n_group_offset[i];
                ASSERT(n_particle>0);
                ASSERT(n_particle<=particles.getSize());
                group_new.particles.setMode(COMM::ListMode::copy);
                group_new.particles.reserveMem(n_particle);
                group_new.reserveIntegratorMem();
                const int nmax_tot = particles.getSizeMax() + groups.getSizeMax();
                group_new.perturber.neighbor_address.setMode(COMM::ListMode::local);
                group_new.perturber.neighbor_address.reserveMem(nmax_tot);
                group_new.info.reserveMem(n_particle);
            
                // Add members to AR 
                for(int j=_n_group_offset[i]; j<_n_group_offset[i+1]; j++) {
                    const int p_index = _particle_index[j];
                    ASSERT(p_index<particles.getSize());
                    group_new.particles.addMemberAndAddress(particles[p_index]);
                    group_new.info.particle_index.addMember(p_index);
                    group_new.info.r_break_crit = std::max(group_new.info.r_break_crit, particles[p_index].getRGroup());
                    Float r_neighbor_crit = particles[p_index].getRNeighbor();
                    group_new.perturber.r_neighbor_crit_sq = std::max(group_new.perturber.r_neighbor_crit_sq, r_neighbor_crit*r_neighbor_crit);
                    // update single mask table 
                    ASSERT(table_single_mask_[p_index]==false);
                    table_single_mask_[p_index] = true;
                }
                
                // calculate the c.m.
                group_new.particles.calcCenterOfMass();
                // set id to group_index + 1
                group_new.particles.cm.id = - (group_index[i]+1);
                // shift to c.m. frame
                group_new.particles.shiftToCenterOfMassFrame();

                // get binarytree
                group_new.info.generateBinaryTree(group_new.particles,ar_manager->interaction.gravitational_constant);

#ifdef ADJUST_GROUP_DEBUG
                std::cerr<<"Add new group, index: "<<group_index[i]<<" Member_index: ";
                for (int k=0; k<group_new.particles.getSize(); k++) 
                    std::cerr<<group_new.info.particle_index[k]<<" ";
                std::cerr<<"r_break_crit: "<<group_new.info.r_break_crit;
                std::cerr<<std::endl;
                COMM::Binary& bin = group_new.info.getBinaryTreeRoot();
                bin.printColumnTitle(std::cerr);
                std::cerr<<std::endl;
                bin.printColumn(std::cerr);
                std::cerr<<std::endl;
#endif

                // initial perturber
                group_new.perturber.need_resolve_flag = true;
            }
                
            // modify the sorted_dt array 
            const int n_group_old = index_dt_sorted_group_.getSize();
            index_dt_sorted_group_.increaseSizeNoInitialize(_n_group);
            for (int i=n_group_old-1; i>=0; i--) {
                index_dt_sorted_group_[i+_n_group] = index_dt_sorted_group_[i];
            }
            for (int i=0; i<_n_group; i++) {
                index_dt_sorted_group_[i] = group_index[i];
            }
            n_init_group_ += _n_group;
            n_act_group_  += _n_group;

            // remove single from index_dt_sorted
            removeDtIndexFromTable(index_dt_sorted_single_, n_act_single_, n_init_single_, table_single_mask_);
            // clear neighbor lists and add group index
            int n_group = index_dt_sorted_group_.getSize();
            for (int i=0; i<n_group; i++) {
                const int k = index_dt_sorted_group_[i];
                auto& nbk = groups[k].perturber;
                // remove single index 
                int n_rm_single= removeNBSingleAddressFromTable(nbk.neighbor_address);
                nbk.n_neighbor_single -= n_rm_single;
                ASSERT(nbk.n_neighbor_single>=0);
                // add new group index
                for (int j=0; j<_n_group; j++) {
                    const int jk = group_index[j];
                    if (k==jk) continue;
                    nbk.neighbor_address.addMember(NBAdr<Tparticle>(&groups[jk].particles, jk+index_offset_group_));
                    nbk.n_neighbor_group++;
                }
            }
        }

        //! Add group based on a configure file
        /*! @param[in] _fin: std::istream IO for read
          File format: N_group, group_offset_index_lst[N_group number], group_member_particle_index[total group member number]
         */
        void readGroupConfigureAscii(std::istream& _fin) {
            int n_group;
            _fin>>n_group;
            ASSERT(!_fin.eof());
            int n_group_offset[n_group+1];

            if (n_group>0) {
                for (int i=0; i<=n_group; i++) {
                    _fin>>n_group_offset[i];
                    ASSERT(!_fin.eof());
                }
                int n_members = n_group_offset[n_group];
                ASSERT(n_members>1);
                int particle_index[n_members];
                for (int i=0; i<n_members; i++) {
                    _fin>>particle_index[i];
                    ASSERT(!_fin.eof());
                }

                addGroups(particle_index, n_group_offset, n_group);
            }
        }

        //! break groups
        /* Split one group to two sub-groups based on the binarytree upper-most pair
           If sub-groups have only one particles, return it to the single list
           Notice if groups are removed, the group index are not changed, only the removed groups are masked and cleared for reusing
           @param[out] _new_group_particle_index_origin: particle index (for hermite particle data) of new group members
           @param[out] _new_n_group_offset:  group member boundary, first value is defined already 
           @param[out] _new_n_group:    number of new group
           @param[in] _break_group_index_with_offset:   break group index in groups
           @param[in] _n_break_no_add: number of break groups without adding new particles
           @param[in] _n_break:   number of break groups
         */
        void breakGroups(int* _new_group_particle_index_origin, 
                         int* _new_n_group_offset, 
                         int& _new_n_group,
                         const int* _break_group_index_with_offset, 
                         const int _n_break_no_add,
                         const int _n_break) {
            ASSERT(checkParams());
            ASSERT(!particles.isModified());
            ASSERT(initial_system_flag_);
            ASSERT(_n_break<=index_dt_sorted_group_.getSize());

            if (_n_break==0 && _n_break_no_add==0) return;

            modify_system_flag_=true;
            int new_index_single[particles.getSize()];
            int n_single_new=0;

            for (int k=0; k<_n_break; k++) {
                const int i = _break_group_index_with_offset[k] - index_offset_group_;
                ASSERT(i>=0&&i<groups.getSize());

                // before break group, first accummulating energy change 
                accumDESlowDownChangeBreakGroup(i);

                auto& groupi = groups[i];
                const int n_member =groupi.particles.getSize();
                int particle_index_origin[n_member];
                int ibreak = groupi.info.getTwoBranchParticleIndexOriginFromBinaryTree(particle_index_origin, groupi.particles.getDataAddress());

#ifdef ADJUST_GROUP_DEBUG
                auto& sd_root = groupi.info.getBinaryTreeRoot().slowdown;
                std::cerr<<"Break Group:  "
                         <<" Time: "<<time_
                         <<" k: "<<std::setw(2)<<i
                         <<" N_member: "<<std::setw(4)<<groupi.particles.getSize()
                         <<" step: "<<std::setw(12)<<groupi.profile.step_count_sum
                         <<" step(tsyn): "<<std::setw(10)<<groupi.profile.step_count_tsyn_sum
//                         <<" step(sum): "<<std::setw(12)<<profile.ar_step_count
//                         <<" step_tsyn(sum): "<<std::setw(12)<<profile.ar_step_count_tsyn
                         <<" Pert_In: "<<std::setw(20)<<sd_root.getPertIn()
                         <<" Pert_Out: "<<std::setw(20)<<sd_root.getPertOut()
                         <<" SD: "<<std::setw(20)<<sd_root.getSlowDownFactor()
                         <<" SD(org): "<<std::setw(20)<<sd_root.getSlowDownFactorOrigin();
                auto& bin = groupi.info.getBinaryTreeRoot();
                std::cerr<<" semi: "<<std::setw(20)<<bin.semi
                         <<" ecc: "<<std::setw(20)<<bin.ecc
                         <<" NB: "<<std::setw(4)<<groupi.perturber.neighbor_address.getSize()
                         <<std::endl;
#endif

#ifdef ADJUST_GROUP_PRINT
                if (manager->adjust_group_write_flag) {
                    groupi.printGroupInfo(1, manager->fgroup, WRITE_WIDTH, &(particles.cm));
                }
#endif

                // clear group
                groupi.particles.shiftToOriginFrame();
                groupi.particles.template writeBackMemberAll<Tparticle>();
                groupi.clear();
                table_group_mask_[i] = true;
                index_group_mask_.addMember(i);
                
                // check single case
                // left side
                if (ibreak==1) {
                    ASSERT(particle_index_origin[0]<particles.getSize());
                    ASSERT(table_single_mask_[particle_index_origin[0]]==true);
                    table_single_mask_[particle_index_origin[0]] = false;
                    if (particles[particle_index_origin[0]].mass > 0.0) // no zero mass particle
                        new_index_single[n_single_new++] = particle_index_origin[0];
                }
                else {
                    // check first how many particles have mass>0
                    int count_mass = 0;
                    int index_has_mass_last = -1;
                    for (int j=0; j<ibreak; j++) {
                        if (particles[particle_index_origin[j]].mass > 0.0) {
                            count_mass++;
                            index_has_mass_last = particle_index_origin[j];
                        }
                    }
                    if (count_mass>=2) { // only form new group when at least two members have mass >0
                        int offset = _new_n_group_offset[_new_n_group];
                        for (int j=0; j<ibreak; j++) {
                            ASSERT(table_single_mask_[particle_index_origin[j]]==true);
                            table_single_mask_[particle_index_origin[j]] = false;
                            _new_group_particle_index_origin[offset++] = particle_index_origin[j];
                        }
                        _new_n_group_offset[++_new_n_group] = offset;
                    }
                    else if (count_mass==1) { // add single if only one particle has mass >0
                        ASSERT(index_has_mass_last<particles.getSize());
                        ASSERT(index_has_mass_last>=0);
                        ASSERT(table_single_mask_[index_has_mass_last]==true);
                        table_single_mask_[index_has_mass_last] = false;
                        new_index_single[n_single_new++] = index_has_mass_last;
                    }
                }
                // right side
                if (n_member-ibreak==1) {
                    ASSERT(particle_index_origin[ibreak]<particles.getSize());
                    ASSERT(table_single_mask_[particle_index_origin[ibreak]]==true);
                    table_single_mask_[particle_index_origin[ibreak]] = false;
                    if (particles[particle_index_origin[ibreak]].mass >0.0) 
                        new_index_single[n_single_new++] = particle_index_origin[ibreak];
                }
                else {
                    // check first how many particles have mass>0
                    int count_mass = 0;
                    int index_has_mass_last = -1;
                    for (int j=ibreak; j<n_member; j++) {
                        if (particles[particle_index_origin[j]].mass > 0.0) {
                            count_mass++;
                            index_has_mass_last = particle_index_origin[j];
                        }
                    }
                    if (count_mass>=2) { // only form new group when at least two members have mass >0
                        int offset = _new_n_group_offset[_new_n_group];
                        for (int j=ibreak; j<n_member; j++) {
                            ASSERT(table_single_mask_[particle_index_origin[j]]==true);
                            table_single_mask_[particle_index_origin[j]] = false;
                            _new_group_particle_index_origin[offset++] = particle_index_origin[j];
                        }
                        _new_n_group_offset[++_new_n_group] = offset;
                    }
                    else if (count_mass==1) { // add single if only one particle has mass >0
                        ASSERT(index_has_mass_last<particles.getSize());
                        ASSERT(index_has_mass_last>=0);
                        ASSERT(table_single_mask_[index_has_mass_last]==true);
                        table_single_mask_[index_has_mass_last] = false;
                        new_index_single[n_single_new++] = index_has_mass_last;
                    }
                }
            }

            // modify the sorted_dt array of single
            const int n_single_old = index_dt_sorted_single_.getSize();
            index_dt_sorted_single_.increaseSizeNoInitialize(n_single_new);
            for (int i=n_single_old-1; i>=0; i--) {
                index_dt_sorted_single_[i+n_single_new] = index_dt_sorted_single_[i];
            }
            for (int i=0; i<n_single_new; i++) {
                index_dt_sorted_single_[i] = new_index_single[i];
            }
            n_init_single_ += n_single_new;
            n_act_single_  += n_single_new;

            // only clear case
            for (int k=_n_break; k<_n_break+_n_break_no_add; k++) {
                const int i = _break_group_index_with_offset[k] - index_offset_group_;

                // before break group, first accummulating energy change
                accumDESlowDownChangeBreakGroup(i);

                auto& groupi = groups[i];
                groupi.particles.shiftToOriginFrame();
                groupi.particles.template writeBackMemberAll<Tparticle>();
                const int n_member = groupi.particles.getSize();
                ASSERT(n_member==groupi.info.particle_index.getSize());
                for (int j=0; j<n_member; j++) {
                    const int kj = groupi.info.particle_index[j];
                    ASSERT(kj<particles.getSize());
                    table_single_mask_[kj] = false;
                }
                groupi.clear();
                table_group_mask_[i] = true;
                index_group_mask_.addMember(i);
            }
            
            // clear index_dt_sorted_
            removeDtIndexFromTable(index_dt_sorted_group_, n_act_group_, n_init_group_, table_group_mask_);
            // clear neighbor lists and add single index
            int n_group = index_dt_sorted_group_.getSize();
            for (int i=0; i<n_group; i++) {
                const int k = index_dt_sorted_group_[i];
                auto& nbk = groups[k].perturber;
                // remove group index
                int n_rm_group = removeNBGroupAddressFromTable(nbk.neighbor_address);
                nbk.n_neighbor_group -= n_rm_group;
                ASSERT(nbk.n_neighbor_group>=0);
                // add single index
                for (int j=0; j<n_single_new; j++) {
                    const int jk = new_index_single[j];
                    nbk.neighbor_address.addMember(NBAdr<Tparticle>(&particles[jk],jk));
                    nbk.n_neighbor_single++;
                }
            }
        }

        //! Check break condition
        /*! Check whether it is necessary to break the chain\n
          1. Inner distance criterion, for outmost pair: \n
             If r>r_break_crit, and ecca>0.0, break. Also predict next step r to ensure break early enough\n
          2. Perturbation criterion for closed orbit:\n
             For two body case, if slowdown factor >1.0 break. 
             For few-body, if no inner AR slowdown, break when inner kappa is large. \n
          @param[out] _break_group_index_with_offset: group index list to break (with index_offset_group_ added)
          @parma[out] _n_break: number of groups need to break
          @param[in] _start_flag: indicate this is the first adjust of the groups in the integration
        */
        void checkBreak(int* _break_group_index_with_offset, int& _n_break, const bool _start_flag) {
            const int n_group_tot = index_dt_sorted_group_.getSize();
            if (n_group_tot==0) return;

            bool merge_mask[n_group_tot];
            for (int k=0; k<n_group_tot; k++) merge_mask[k] = false;

            // check merger case
            const int n_merger = index_group_merger_.getSize();
            for (int i=0; i<n_merger; i++) {
                const int k = index_group_merger_[i];
                ASSERT(k<groups.getSize());
                ASSERT(table_group_mask_[k]==false);
                auto& groupk = groups[k];
#ifdef ADJUST_GROUP_DEBUG
                auto& bin_root = groupk.info.getBinaryTreeRoot();
                std::cerr<<"Break merger group: i_group: "<<k
                         <<" N_member: "<<groupk.particles.getSize()
                         <<" ecca: "<<bin_root.ecca
                         <<" separation : "<<bin_root.r
                         <<" apo: "<<bin_root.semi*(1.0+bin_root.ecc)
                         <<" dt: "<<groupk.particles.cm.dt
                         <<" r_crit: "<<groupk.info.r_break_crit
                         <<std::endl;
#endif
                // generate binary tree in order to move zero mass to the outer most
                groupk.info.generateBinaryTree(groupk.particles,ar_manager->interaction.gravitational_constant);
                _break_group_index_with_offset[_n_break++] = k + index_offset_group_;
                merge_mask[k] = true;
                
            }
            index_group_merger_.resizeNoInitialize(0);

            // kappa_org criterion for break group kappa_org>kappa_org_crit
            const Float kappa_org_crit = 1e-2;
            for (int i=0; i<n_group_tot; i++) {
                const int k = index_dt_sorted_group_[i];
                ASSERT(table_group_mask_[k]==false);
                if (merge_mask[k]) continue;
                auto& groupk = groups[k];

                const int n_member = groupk.particles.getSize();

                // generate binary tree
                groupk.info.generateBinaryTree(groupk.particles,ar_manager->interaction.gravitational_constant);
                
                auto& bin_root = groupk.info.getBinaryTreeRoot();
                bool outgoing_flag = false; // Indicate whether it is a outgoing case or income case

                // check binary case 
                // ecc anomaly indicates outgoing (ecca>0) or income (ecca<0)
                if (bin_root.semi>0.0 && bin_root.ecca>0.0) {
                    outgoing_flag = true;
                    // check whether separation is larger than distance criterion. 
                    if (bin_root.r > groupk.info.r_break_crit) {
#ifdef ADJUST_GROUP_DEBUG
                        std::cerr<<"Break group: binary escape, time: "<<time_<<" i_group: "<<k<<" N_member: "<<n_member<<" ecca: "<<bin_root.ecca<<" separation : "<<bin_root.r<<" r_crit: "<<groupk.info.r_break_crit<<std::endl;
#endif
                        _break_group_index_with_offset[_n_break++] = k + index_offset_group_;
                        continue;
                    }

                    // in case apo is larger than distance criterion
                    Float apo = bin_root.semi * (1.0 + bin_root.ecc);
                    if (apo>groupk.info.r_break_crit) {
                        Float dr2, drdv;

                        groupk.info.getDrDv(dr2, drdv, *bin_root.getLeftMember(), *bin_root.getRightMember());
                        ASSERT(drdv>=0.0);

                        // check whether next step the separation is larger than distance criterion
                        // Not sure whether it can work correctly or not:
                        // rp = v_r * dt + r
                        Float dr = drdv/bin_root.r*groups[k].particles.cm.dt;
                        Float rp =  dr + bin_root.r;
                        if (rp >groupk.info.r_break_crit) {
                            // in case r is too small, avoid too early quit of group
                            Float rph = 0.5*dr + bin_root.r;
                            if ( rph < groupk.info.r_break_crit && bin_root.r <0.2*groupk.info.r_break_crit) {
                                if (getNextTime()>time_+groups[k].particles.cm.dt) {
#ifdef ADJUST_GROUP_DEBUG
                                    std::cerr<<"Binary will escape but dr is too small, reduce cm step by half, time: "<<time_<<" i_group: "<<k<<" N_member: "<<n_member<<" ecca: "<<bin_root.ecca<<" separation : "<<bin_root.r<<" apo: "<<apo<<" r_pred: "<<rp<<" drdv: "<<drdv<<" dt: "<<groups[k].particles.cm.dt<<" r_crit: "<<groupk.info.r_break_crit<<std::endl;
#endif
                                    reduceGroupCMStepByHalfAndSortDtIndex(i);
                                }
                            }
                            else {
#ifdef ADJUST_GROUP_DEBUG
                                std::cerr<<"Break group: binary will escape, time: "<<time_<<" i_group: "<<k<<" N_member: "<<n_member<<" ecca: "<<bin_root.ecca<<" separation : "<<bin_root.r<<" apo: "<<apo<<" r_pred: "<<rp<<" drdv: "<<drdv<<" dt: "<<groups[k].particles.cm.dt<<" r_crit: "<<groupk.info.r_break_crit<<std::endl;
#endif
                                _break_group_index_with_offset[_n_break++] = k + index_offset_group_;
                            }
                            continue;
                        }
                    }

                }

                // check hyperbolic case
                if (bin_root.semi<0.0) {
                    // hyperbolic case, ecca is not correctly calculated
                    Float dr2, drdv;
                    groupk.info.getDrDv(dr2, drdv, *bin_root.getLeftMember(), *bin_root.getRightMember());
                    if (drdv>0.0) {
                        outgoing_flag = true;
                        // check distance criterion
                        if (bin_root.r > groupk.info.r_break_crit) {
#ifdef ADJUST_GROUP_DEBUG
                            std::cerr<<"Break group: hyperbolic escape, time: "<<time_<<" i_group: "<<k<<" N_member: "<<n_member<<" drdv: "<<drdv<<" separation : "<<bin_root.r<<" r_crit: "<<groupk.info.r_break_crit<<std::endl;
#endif
                            _break_group_index_with_offset[_n_break++] = k + index_offset_group_;
                            continue;
                        }
                        // check for next step
                        Float dr = drdv/bin_root.r*groups[k].particles.cm.dt;
                        Float rp = dr  + bin_root.r;
                        if (rp > groupk.info.r_break_crit) {
                            // in case r is too small, avoid too early quit of group
                            Float rph = 0.5*dr + bin_root.r;
                            if ( rph < groupk.info.r_break_crit && bin_root.r <0.2*groupk.info.r_break_crit) {
                                if (getNextTime()>time_+groups[k].particles.cm.dt) {
#ifdef ADJUST_GROUP_DEBUG
                                    std::cerr<<"Hyperbolic will escape but dr is too small, reduce cm step by half first, time: "<<time_<<" i_group: "<<k<<" N_member: "<<n_member<<" drdv: "<<drdv<<" separation : "<<bin_root.r<<" r_pred: "<<rp<<" drdv: "<<drdv<<" dt: "<<groups[k].particles.cm.dt<<" r_crit: "<<groupk.info.r_break_crit<<std::endl;
#endif
                                    reduceGroupCMStepByHalfAndSortDtIndex(i);
                                }
                            }
                            else {
#ifdef ADJUST_GROUP_DEBUG
                                std::cerr<<"Break group: hyperbolic will escape, time: "<<time_<<" i_group: "<<k<<" N_member: "<<n_member<<" drdv: "<<drdv<<" separation : "<<bin_root.r<<" r_pred: "<<rp<<" drdv: "<<drdv<<" dt: "<<groups[k].particles.cm.dt<<" r_crit: "<<groupk.info.r_break_crit<<std::endl;
#endif
                                _break_group_index_with_offset[_n_break++] = k + index_offset_group_;
                            }
                            continue;
                        }
                    }

                }

                // check perturbation
                // only check further if it is outgoing case
                if (outgoing_flag) {

                    AR::SlowDown sd;
                    auto& sd_group = groupk.info.getBinaryTreeRoot().slowdown;
                    sd.initialSlowDownReference(sd_group.getSlowDownFactorReference(),sd_group.getSlowDownFactorMax());
                    sd.timescale = sd_group.timescale;
                    sd.period = sd_group.period;

                    if (n_member==2) {
                        // check strong perturbed binary case 
                        // calculate slowdown in a consistent way like in checknewgroup to avoid switching
                        // fcm may not properly represent the perturbation force (perturber mass is unknown)
                        sd.pert_in = ar_manager->interaction.calcPertFromBinary(bin_root);
                        Float* acc_cm = groupk.particles.cm.acc0;
                        Float fcm[3] = {acc_cm[0]*bin_root.mass, acc_cm[1]*bin_root.mass, acc_cm[2]*bin_root.mass };
                        sd.pert_out= ar_manager->interaction.calcPertFromForce(fcm, bin_root.mass, bin_root.mass);
                        sd.calcSlowDownFactor();
                        Float kappa_org = sd.getSlowDownFactorOrigin();

                        if (kappa_org<kappa_org_crit && !_start_flag) {
                            // in binary case, only break when apo is larger than distance criterion
                            Float apo = bin_root.semi * (1.0 + bin_root.ecc);
                            if (apo>groupk.info.r_break_crit||bin_root.semi<0.0) {
#ifdef ADJUST_GROUP_DEBUG
                                std::cerr<<"Break group: strong perturbed, time: "<<time_<<" i_group: "<<k<<" N_member: "<<n_member;
                                std::cerr<<" index: ";
                                for (int i=0; i<n_member; i++) 
                                    std::cerr<<groupk.info.particle_index[i]<<" ";
                                auto& sd_root = groupk.info.getBinaryTreeRoot().slowdown;
                                std::cerr<<" pert_in: "<<sd_root.pert_in
                                         <<" pert_out: "<<sd_root.pert_out
                                         <<" kappa_org: "<<kappa_org
                                         <<" dr: "<<bin_root.r
                                         <<" semi: "<<bin_root.semi
                                         <<" ecc: "<<bin_root.ecc
                                         <<" r_break: "<<groupk.info.r_break_crit
                                         <<std::endl;
#endif
                                _break_group_index_with_offset[_n_break++] = k + index_offset_group_;
                                continue;
                            }
                        }
                    }
#if (!defined AR_SLOWDOWN_ARRAY) && (!defined AR_SLOWDOWN_TREE)
                    // check few-body inner perturbation (suppress when use slowdown inner AR)
                    else {
                        for (int j=0; j<2; j++) {
                            if (bin_root.isMemberTree(j)) {
                                auto* bin_sub = bin_root.getMemberAsTree(j);
//                            Float semi_db = 2.0*bin_sub->semi;
//                            // inner hyperbolic case
//                            if(semi_db<0.0 && abs(groupk.getEnergyError()/groupk.getEtot())<100.0*groupk.manager->energy_error_relative_max && bin_root->ecca>0.0) {
//#ifdef ADJUST_GROUP_DEBUG
//                                std::cerr<<"Break group: inner member hyperbolic, time: "<<time_<<" i_group: "<<k<<" i_member: "<<j<<" semi: "<<semi_db<<" ecca: "<<bin_sub->ecca<<std::endl;
//#endif
//                                _break_group_index_with_offset[_n_break++] = k + index_offset_group_;
//                                break;
//                            }

                                // check inner binary slowdown factor
                                if (bin_sub->semi>0.0) {

                                    Float apo_in = bin_sub->semi*(1+bin_sub->ecc);
                                    sd.pert_in = ar_manager->interaction.calcPertFromMR(apo_in, bin_sub->m1, bin_sub->m2);

                                    // present slowdown 
                                    sd.pert_out = ar_manager->interaction.calcPertFromMR(bin_root.r, bin_root.m1, bin_root.m2);
                                    sd.calcSlowDownFactor();
                                    Float kappa_in = sd.getSlowDownFactorOrigin();

                                    // in case slowdown >1
                                    if (kappa_in>1.0) {
                                        Float kappa_in_max = NUMERIC_FLOAT_MAX;
                                        // if outer is binary, estimate slowdown max (apo_out)
                                        if (bin_root.semi>0.0) {
                                            Float apo_out = bin_root.semi*(1+bin_root.ecc);
                                            sd.pert_out = ar_manager->interaction.calcPertFromMR(apo_out, bin_root.m1, bin_root.m2);
                                            sd.calcSlowDownFactor();

                                            kappa_in_max = sd.getSlowDownFactorOrigin();
                                        }
                                    
                                        // if slowdown factor is large, break the group
                                        if (kappa_in_max>5.0) {
                                            // avoid quit at high energy error phase
                                            if (abs(groupk.getEnergyError()/groupk.getEtotRef())<100.0/(1-std::min(bin_sub->ecc,bin_root.ecc))*groupk.manager->energy_error_relative_max) {
#ifdef ADJUST_GROUP_DEBUG
                                                std::cerr<<"Break group: inner kappa large, time: "<<time_<<" i_group: "<<k<<" i_member: "<<j<<" kappa_in:"<<kappa_in<<" kappa_in(max):"<<kappa_in_max
                                                         <<" Energy error:"<<groupk.getEnergyError()
                                                         <<" Etot ref:"<<groupk.getEtotRef()
                                                         <<" ecc(in):"<<bin_sub->ecc
                                                         <<" ecc(out):"<<bin_root.ecc
                                                         <<std::endl;
#endif
                                                _break_group_index_with_offset[_n_break++] = k + index_offset_group_;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
#endif
                }
            }
        }

        //! check the pair with distance below r_crit for ptcl in adr_dt_sorted_
        /*! First check nearest neighbor distance r_min
          If r_min<r_crit, check the direction, if income, accept as group
          @param[out] _new_group_particle_index_origin: particle index (for hermite particle data) of new group members
          @param[out] _new_n_group_offset:  group member boundary, first value is defined already 
          @param[out] _new_n_group:    number of new group
          @param[in,out] _break_group_index_with_offset:   break group index in groups
          @param[out] _n_break_no_add: number of break groups without adding new particles
          @param[in] _n_break:   number of break groups
          @param[in] _start_flag: indicate this is the first adjust of the groups in the integration
        */
        void checkNewGroup(int* _new_group_particle_index_origin,
                           int* _new_n_group_offset, 
                           int& _new_n_group, 
                           int* _break_group_index_with_offset,
                           int& _n_break_no_add,
                           const int _n_break, 
                           const bool _start_flag) {
            // kappa_org criterion for new group kappa_org>kappa_org_crit
            const Float kappa_org_crit = 1e-2;

            const int n_particle = particles.getSize();
            const int n_group = groups.getSize();
            ASSERT(index_offset_group_==n_particle);

            int new_n_particle = 0;
            // used_mask store the current group index of particle i if it is already in new group
            int used_mask[n_particle+n_group];
            // -1 means not yet used 
            for (int k=0; k<n_particle; k++) used_mask[k] = table_single_mask_[k]?-2:-1;
            for (int k=0; k<n_group; k++) used_mask[k+index_offset_group_] = table_group_mask_[k]?-2:-1;
            // -2 means break groups/masked groups
            for (int k=0; k<_n_break; k++) used_mask[_break_group_index_with_offset[k]] = -2;
            
            //const Float r_crit_sq = manager->r_break_crit*manager->r_break_crit;

            // gether list together
            //const int n_check = n_act_single_+n_act_group_;
            //int check_index[n_check];
            //for (int k=0; k<n_act_single_; k++) check_index[k] = index_dt_sorted_single_[k];
            //for (int k=0; k<n_act_group_; k++) check_index[k+n_act_single_] = index_dt_sorted_group_[k] + index_offset_group_;

            // check only active particles 
            // single case
            for (int k=0; k<n_act_single_; k++) {
                const int i = index_dt_sorted_single_[k];
                const int j = neighbors[i].r_min_index;
                ASSERT(j<n_particle+n_group);
                if(j<0) continue; 

                const Float dr2 = neighbors[i].r_min_sq;
                ASSERT(dr2>0.0);

                auto& pi = particles[i];

                // avoid zero mass particle
                if (pi.mass==0.0) continue;

                // distance criterion
                Float r_crit = pi.getRGroup();
                if (j<index_offset_group_) 
                    r_crit = std::max(particles[j].getRGroup(), r_crit);
                else
                    r_crit = std::max(groups[j-index_offset_group_].info.r_break_crit, r_crit);
                Float r_crit_sq = r_crit*r_crit;
                
                if (dr2 < r_crit_sq) {

                    // avoid break/masked member
                    if(used_mask[j]==-2) continue;
                    
                    // avoid double count
                    if(used_mask[i]>=0 && used_mask[j]>=0) continue;

                    H4Ptcl* pj;
                    // neighbor is single 
                    if (j<index_offset_group_) {
                        pj = &particles[j];
                    }
                    else {
                        const int jg = j-index_offset_group_;
//#ifndef AR_SLOWDOWN_ARRAY
//                        // in case without AR slowdown inner, avoid form AR when inner kappa is >1.0
//                        Float kappa_org_j = groups[jg].info.getBinaryTreeRoot().slowdown.getSlowDownFactorOrigin();
//                        if (kappa_org_j>1.0) continue;
//#endif
                        pj = &groups[jg].particles.cm;
                    }

                    
                    // avoid zero mass particle
                    if(pj->mass==0.0) continue;

                    // this increase AR total step too much
                    //bool add_flag=false;
                    //Float semi, ecc, dr, drdv;
                    //AR::BinaryTree<Tparticle>::particleToSemiEcc(semi,ecc,dr,drdv, pi, *pj, ar_manager->interaction.gravitational_constant); 
                    //else {
                    //    Float rp = drdv/dr*dt_limit_ + dr;
                    //    if (rp < r_crit) add_flag = true;
                    //}

                    Float drdv = calcDrDv(pi, *pj);
                    // only inwards or first step case
                    if(drdv<0.0||_start_flag) {
                        //Float mcm = pi.mass + pj->mass;
                        Float fcm[3] = {pi.mass*pi.acc0[0] + pj->mass*pj->acc0[0], 
                                        pi.mass*pi.acc0[1] + pj->mass*pj->acc0[1], 
                                        pi.mass*pi.acc0[2] + pj->mass*pj->acc0[2]};

                        AR::SlowDown sd;
                        Float mcm = pi.mass + pj->mass;
#ifdef AR_SLOWDOWN_MASSRATIO
                        const Float mass_ratio = ar_manager->slowdown_mass_ref/mcm;
                        sd.initialSlowDownReference(mass_ratio*ar_manager->slowdown_pert_ratio_ref, ar_manager->slowdown_timescale_max);
#else
                        sd.initialSlowDownReference(ar_manager->slowdown_pert_ratio_ref, ar_manager->slowdown_timescale_max);
#endif
                        sd.pert_in = ar_manager->interaction.calcPertFromMR(sqrt(dr2), pi.mass, pj->mass);
                        sd.pert_out = ar_manager->interaction.calcPertFromForce(fcm, mcm, mcm);

                        sd.calcSlowDownFactor();
                        Float kappa_org = sd.getSlowDownFactorOrigin();

                        // avoid strong perturbed case, estimate perturbation
                        // if kappa_org < criterion, avoid to form new group, should be consistent as checkbreak
                        if(kappa_org<kappa_org_crit) continue;

#ifdef ADJUST_GROUP_DEBUG
                        if (j<index_offset_group_) {
                            std::cerr<<"Find new group: time: "<<time_
                                     <<" index: "<<i<<" "<<j
                                     <<" dr: "<<sqrt(dr2)
                                     <<" kappa_org: "<<kappa_org<<"\n";
                        }
                        else {
                            auto& bin_root = groups[j-index_offset_group_].info.getBinaryTreeRoot();
                            std::cerr<<"Find new group: time: "<<time_
                                     <<" dr: "<<sqrt(dr2)
                                     <<" kappa_org: "<<kappa_org<<"\n"
                                     <<"       index         slowdown          apo \n"
                                     <<"i1 "
                                     <<std::setw(8)<<i
                                     <<std::setw(16)<<0
                                     <<std::setw(16)<<0;
                            std::cerr<<"\ni2 "
                                     <<std::setw(8)<<j
                                     <<std::setw(16)<<bin_root.slowdown.getSlowDownFactorOrigin()
                                     <<std::setw(16)<<bin_root.semi*(1.0+bin_root.ecc);
                            std::cerr<<std::endl;
                        }
#endif
                        insertParticleIndexToGroup(i, j, used_mask, _new_group_particle_index_origin, _new_n_group_offset, new_n_particle, _new_n_group, _break_group_index_with_offset, _n_break_no_add,  _n_break);
                    }
                }
            }

            // group case
            for (int k=0; k<n_act_group_; k++) {
                const int ig = index_dt_sorted_group_[k];
                const int i = ig + index_offset_group_;
                auto& groupi = groups[ig];
                if (used_mask[i]==-2) continue;

                const int j = groupi.perturber.r_min_index;
                ASSERT(j<n_particle+n_group);
                if(j<0) continue; 
                
                const Float dr2 = groupi.perturber.r_min_sq;
                ASSERT(dr2>0.0);

                Float r_crit = groupi.info.r_break_crit;
                if (j<index_offset_group_) 
                    r_crit = std::max(particles[j].getRGroup(), r_crit);
                else
                    r_crit = std::max(groups[j-index_offset_group_].info.r_break_crit, r_crit);
                Float r_crit_sq = r_crit*r_crit;

                // distance criterion
                if (dr2 < r_crit_sq) {

                    // avoid break member
                    if(used_mask[j]==-2) continue;
                    
                    // avoid double count
                    if(used_mask[i]>=0 && used_mask[j]>=0) continue;

                    auto& pi = groupi.particles.cm;

//#ifndef AR_SLOWDOWN_ARRAY
//                    // avoid kappa>1.0
//                    Float kappa_org_i = groupi.info.getBinaryTreeRoot().slowdown.getSlowDownFactorOrigin();
//                    if (kappa_org_i>1.0) continue;
//#endif

                    H4Ptcl* pj;
                    // neighbor is single 
                    if (j<index_offset_group_) {
                        //if (kappa_org_i>1.0) continue;
                        pj = &particles[j];
                    }
                    else {
                        const int jg = j-index_offset_group_;
//#ifndef AR_SLOWDOWN_ARRAY
//                        Float kappa_org_j = groups[jg].getBinaryTreeRoot().slowdown.getSlowDownFactorOrigin();
//                        if (kappa_org_j>1.0) continue;
//#endif

                        //if (kappa_org_i>1.0&&kappa_org_j>1.0) continue;
                        pj = &groups[jg].particles.cm;

                        // unknown, for test
                        //if (kappa_org_i*kappa_org_j>1.0) continue;
                    }
                    // only inwards or first step case
                    Float drdv = calcDrDv(pi, *pj);
                    if(drdv<0.0||_start_flag) {

                        //Float mcm = pi.mass + pj->mass;
                        Float fcm[3] = {pi.mass*pi.acc0[0] + pj->mass*pj->acc0[0], 
                                        pi.mass*pi.acc0[1] + pj->mass*pj->acc0[1], 
                                        pi.mass*pi.acc0[2] + pj->mass*pj->acc0[2]};

                        AR::SlowDown sd;
                        Float mcm = pi.mass + pj->mass;
#ifdef AR_SLOWDOWN_MASSRATIO
                        const Float mass_ratio = ar_manager->slowdown_mass_ref/mcm;
                        sd.initialSlowDownReference(mass_ratio*ar_manager->slowdown_pert_ratio_ref, ar_manager->slowdown_timescale_max);
#else
                        sd.initialSlowDownReference(ar_manager->slowdown_pert_ratio_ref, ar_manager->slowdown_timescale_max);
#endif
                        sd.pert_in = ar_manager->interaction.calcPertFromMR(sqrt(dr2), pi.mass, pj->mass);
                        // fcm may not properly represent the perturbation force (perturber mass is unknown)
                        sd.pert_out = ar_manager->interaction.calcPertFromForce(fcm, mcm, mcm);

                        sd.calcSlowDownFactor();
                        Float kappa_org = sd.getSlowDownFactorOrigin();

                        // avoid strong (outside) perturbed case, estimate perturbation
                        // if fratiosq >1.5, avoid to form new group, should be consistent as checkbreak
                        if(kappa_org<kappa_org_crit) continue;

#ifdef ADJUST_GROUP_DEBUG
                        auto& bini = groupi.info.getBinaryTreeRoot();
                        std::cerr<<"Find new group: time: "<<time_
                                 <<" dr: "<<sqrt(dr2)
                                 <<" kappa_org: "<<kappa_org
                                 <<"\n       index        slowdown         apo  \n"
                                 <<"i1 "
                                 <<std::setw(8)<<i
                                 <<std::setw(16)<<bini.slowdown.getSlowDownFactorOrigin()
                                 <<std::setw(16)<<bini.semi*(1.0+bini.ecc);
                        if(j<index_offset_group_) {
                            std::cerr<<"\ni2 "
                                     <<std::setw(8)<<j
                                     <<std::setw(16)<<0
                                     <<std::setw(16)<<0;
                        }
                        else {
                            auto& binj = groups[j-index_offset_group_].info.getBinaryTreeRoot();
                            Float kappaj = binj.slowdown.getSlowDownFactorOrigin();
                            std::cerr<<"\ni2 "
                                     <<std::setw(8)<<j
                                     <<std::setw(16)<<kappaj
                                     <<std::setw(16)<<binj.semi*(1.0+binj.ecc);
                        }
                        std::cerr<<std::endl;
#endif
                        insertParticleIndexToGroup(i, j, used_mask, _new_group_particle_index_origin, _new_n_group_offset, new_n_particle, _new_n_group, _break_group_index_with_offset, _n_break_no_add,  _n_break);
                    }
                }
            }
            // for total number of members
            _new_n_group_offset[_new_n_group] = new_n_particle;
            ASSERT(new_n_particle<=particles.getSize());
        }
        

        //! adjust groups
        /* check break and new groups and modify the groups
           @param[in] _start_flag: indicate this is the first adjust of the groups in the integration
         */
        void adjustGroups(const bool _start_flag) {
            ASSERT(checkParams());
            ASSERT(!particles.isModified());
            ASSERT(initial_system_flag_);
            ASSERT(ar_manager->interrupt_detection_option!=2||(ar_manager->interrupt_detection_option==2&&interrupt_binary_.status==AR::InterruptStatus::none));
            modify_system_flag_=true;

            // check
            // break index (with index_offset_group)
            int break_group_index_with_offset[groups.getSize()+1];
            int n_break = 0;
            int n_break_no_add = 0;
            // new index
            int new_group_particle_index[particles.getSize()];
            int new_n_group_offset[groups.getSizeMax()+1];
            int new_n_group = 0;

            checkBreak(break_group_index_with_offset, n_break, _start_flag);
            checkNewGroup(new_group_particle_index, new_n_group_offset, new_n_group, break_group_index_with_offset, n_break_no_add, n_break, _start_flag);
            ASSERT(n_break<=groups.getSize());
            ASSERT(new_n_group<=groups.getSizeMax());
            profile.break_group_count += n_break;
            profile.new_group_count += new_n_group;

            // integrate modified single/groups to current time
            integrateToTimeList(time_, new_group_particle_index, new_n_group_offset[new_n_group]);
            integrateToTimeList(time_, break_group_index_with_offset, n_break);

            // break groups
            breakGroups(new_group_particle_index, new_n_group_offset, new_n_group, break_group_index_with_offset, n_break_no_add, n_break);
            addGroups(new_group_particle_index, new_n_group_offset, new_n_group);

            // initial integration (cannot do it here, in the case AR perturber need initialization first)
            // initialIntegration();
        }

        //! Initial Hermite integrator
        /*! Initial f, fdot and step for new add/remove particles. 
          //If start_flag, the whole system are initialized; else the first n_init_single_ and n_init_group_ of sorted index list will be initialized. 
          The active particle number will be set to the total number of particles
          //@param[in] _start_flag: true: the starting step of integration.
        */
        void initialIntegration() {

            ASSERT(!particles.isModified());
            ASSERT(initial_system_flag_);

            modify_system_flag_=false;


            // check table size and index_dt_sorted size
#ifdef HERMITE_DEBUG
            int particle_index_count[particles.getSize()];
            int group_index_count[groups.getSize()];
            for (int i=0; i<particles.getSize(); i++) particle_index_count[i]=0;
            for (int i=0; i<groups.getSize(); i++) group_index_count[i]=0;
            // single list
            for (int i=0; i<index_dt_sorted_single_.getSize(); i++) {
                ASSERT(index_dt_sorted_single_[i]<particles.getSize());
                particle_index_count[index_dt_sorted_single_[i]]++;
            }
            for (int i=0; i<index_dt_sorted_group_.getSize(); i++) {
                ASSERT(index_dt_sorted_group_[i]<groups.getSize());
                group_index_count[index_dt_sorted_group_[i]]++;
            }
            ASSERT(table_single_mask_.getSize()==particles.getSize());
            for (int i=0; i<table_single_mask_.getSize(); i++) {
                if (table_single_mask_[i]) particle_index_count[i]++;
            }
            ASSERT(table_group_mask_.getSize()==groups.getSize());
            for (int i=0; i<table_group_mask_.getSize(); i++) {
                if (table_group_mask_[i]) group_index_count[i]++;
            }
            for (int i=0; i<particles.getSize(); i++) 
                ASSERT(particle_index_count[i]==1||(particle_index_count[i]==0&&particles[i].mass==0.0));
            for (int i=0; i<groups.getSize(); i++) 
                ASSERT(group_index_count[i]==1);

            ASSERT(n_init_group_<=index_dt_sorted_group_.getSize());
            ASSERT(n_act_group_<=index_dt_sorted_group_.getSize());
            ASSERT(n_init_single_<=index_dt_sorted_single_.getSize());
            ASSERT(n_act_single_<=index_dt_sorted_single_.getSize());
#endif            

            if (n_init_single_==0&&n_init_group_==0) return;

            // single
            int* index_single = index_dt_sorted_single_.getDataAddress();

            auto* ptcl = particles.getDataAddress();
            for(int i=0; i<n_init_single_; i++){
                int k = index_single[i];
                ptcl[k].acc0[0] = ptcl[k].acc0[1] = ptcl[k].acc0[2] = 0.0;
                ptcl[k].acc1[0] = ptcl[k].acc1[1] = ptcl[k].acc1[2] = 0.0;
                ptcl[k].time = time_;
                ptcl[k].dt   = 0.0;
                pred_[k] = ptcl[k];
                // set neighbor radius
                neighbors[k].clearNoFreeMem();
                Float r_neighbor_crit = ptcl[k].getRNeighbor();
                neighbors[k].r_neighbor_crit_sq = r_neighbor_crit*r_neighbor_crit;

                // check parameters
                ASSERT(neighbors[k].checkParams());
            }

            //group
            int* index_group = index_dt_sorted_group_.getDataAddress();

            auto* group_ptr = groups.getDataAddress();
            for(int i=0; i<n_init_group_; i++){
                int k = index_group[i];
                ASSERT(table_group_mask_[k]==false);
                // initial cm
                auto& pcm = group_ptr[k].particles.cm;
                pcm.acc0[0] = pcm.acc0[1] = pcm.acc0[2] = 0.0;
                pcm.acc1[0] = pcm.acc1[1] = pcm.acc1[2] = 0.0;
                pcm.time = time_;
                pcm.dt   = 0.0;
                ASSERT(k+index_offset_group_<pred_.getSize());
                pred_[k+index_offset_group_] = pcm;
            }

            dt_limit_ = step.calcNextDtLimit(time_);

            // slowdown not yet initialized, cannot check
            // checkGroupResolve(n_init_group_);
            writeBackResolvedGroupAndCreateJParticleList(false);

            calcAccJerkNBList(index_single, n_init_single_, index_group, n_init_group_);

            // store predicted force
            for(int i=0; i<n_init_single_; i++){
                const int k = index_single[i];
                ptcl[k].acc0[0] = force_[k].acc0[0];
                ptcl[k].acc0[1] = force_[k].acc0[1];
                ptcl[k].acc0[2] = force_[k].acc0[2];
                ptcl[k].acc1[0] = force_[k].acc1[0];
                ptcl[k].acc1[1] = force_[k].acc1[1];
                ptcl[k].acc1[2] = force_[k].acc1[2];
            }

            for(int i=0; i<n_init_group_; i++){
                const int k = index_group[i];
                auto& pcm = group_ptr[k].particles.cm;
                auto& fcm = force_[k+index_offset_group_];
                
                pcm.acc0[0] = fcm.acc0[0];
                pcm.acc0[1] = fcm.acc0[1];
                pcm.acc0[2] = fcm.acc0[2];
                pcm.acc1[0] = fcm.acc1[0];
                pcm.acc1[1] = fcm.acc1[1];
                pcm.acc1[2] = fcm.acc1[2];

                // initial group integration
                group_ptr[k].initialIntegration(time_);

                // if >2 particles, initial slowdown perturbation and period
                //if (group_ptr[k].particles.getSize()>2) {
                //    auto& bin_root = group_ptr[k].info.binarytree.getLastMember();
                //    group_ptr[k].slowdown.calcPertInBinary(bin_root.semi, bin_root.m1, bin_root.m2);
                //    group_ptr[k].slowdown.calcSlowDownFactor();
                //}

                // set time offset
                group_ptr[k].info.time_offset = time_offset_;
                // set hermite time limit
                group_ptr[k].info.dt_limit = step.getDtMax();
                // in the case of wide binary, make sure the next time step not exceed r_in
                auto& bin_root = group_ptr[k].info.getBinaryTreeRoot();
                if (bin_root.semi*(1.0+bin_root.ecc)>group_ptr[k].info.r_break_crit||bin_root.semi<0) {
                    // ecca <0 indicate binary go to peri-center, ecca=0.0 indicate peri-center, 0-t_peri indicate the time to peri-center
                    // In the initial step, ecca can >0.0 
                    // get boundary position ecca
                    Float ecc_anomaly = bin_root.calcEccAnomaly(group_ptr[k].info.r_break_crit);
                    Float mean_anomaly = bin_root.calcMeanAnomaly(ecc_anomaly, bin_root.ecc);
                    Float kappa = bin_root.slowdown.getSlowDownFactor();
                    // get cross boundary position timescale (half the time to peri-center from boundary), notice slowdown factor should be included
                    // the integrated orbital phase is slow-down by kappa
                    Float t_peri = 0.5*abs(mean_anomaly/12.5663706144*bin_root.period)*kappa;
                    Float dt_limit = group_ptr[k].info.dt_limit;
                    while(dt_limit > t_peri) dt_limit *= 0.5;
                    group_ptr[k].info.dt_limit = dt_limit;
                }

                // check parameters
                ASSERT(group_ptr[k].info.checkParams());
                ASSERT(group_ptr[k].perturber.checkParams());

#ifdef ADJUST_GROUP_PRINT
                if (manager->adjust_group_write_flag) {
                    group_ptr[k].printGroupInfo(0, manager->fgroup, WRITE_WIDTH, &(particles.cm));
                }
#endif

            }

            calcDt2ndList(index_single, n_init_single_, index_group, n_init_group_, dt_limit_);
            
            updateTimeNextList(index_single, n_init_single_, index_group, n_init_group_);

            // reset n_init
            n_init_single_ = n_init_group_ = 0;
        }

        //! Integrate groups
        /*! Integrate all groups to time
          \return interrupted binarytree if exist
         */
        AR::InterruptBinary<Tparticle>& integrateGroupsOneStep() {
            ASSERT(checkParams());
            ASSERT(!particles.isModified());
            ASSERT(initial_system_flag_);
            ASSERT(!modify_system_flag_);
            ASSERT(n_init_group_==0&&n_init_single_==0);

            // get next time
            Float time_next = getNextTime();

            if (ar_manager->interrupt_detection_option==1) {
                interrupt_group_dt_sorted_group_index_ = -1;
                interrupt_binary_.clear();
            }

#ifdef HERMITE_DEBUG            
            if (interrupt_binary_.status!=AR::InterruptStatus::none) ASSERT(interrupt_group_dt_sorted_group_index_>=0);
#endif

            // integrate groups loop 
            const int n_group_tot = index_dt_sorted_group_.getSize();
            const int i_start = interrupt_group_dt_sorted_group_index_>=0 ? interrupt_group_dt_sorted_group_index_ : 0;
            int interrupt_index_dt_group_list[n_group_tot];
            int n_interrupt_change_dt=0;

            for (int i=i_start; i<n_group_tot; i++) {
                const int k = index_dt_sorted_group_[i];

#ifdef HERMITE_DEBUG            
                ASSERT(table_group_mask_[k]==false);
                if (i!=interrupt_group_dt_sorted_group_index_) 
                    ASSERT(abs(groups[k].getTime()-time_)<=ar_manager->time_error_max);
#endif
                // get ds estimation
                groups[k].info.calcDsAndStepOption(ar_manager->step.getOrder(), ar_manager->interaction.gravitational_constant, ar_manager->ds_scale);

                // group integration 
                interrupt_binary_ = groups[k].integrateToTime(time_next);

                // profile
                profile.ar_step_count += groups[k].profile.step_count;
                profile.ar_step_count_tsyn += groups[k].profile.step_count_tsyn;

                // record interrupt group and quit
                if (interrupt_binary_.status!=AR::InterruptStatus::none) {
                    // particle cm is the old cm in original frame
                    auto& pcm = groups[k].particles.cm;

                    if (interrupt_binary_.status==AR::InterruptStatus::merge||interrupt_binary_.status==AR::InterruptStatus::destroy) {
                        index_group_merger_.addMember(k);
                    }
                    else {
                        // set initial step flag 
                        groups[k].perturber.initial_step_flag=true;

                        if (pcm.time + pcm.dt >time_next) {
                            ASSERT(i>=n_act_group_);
                            // set cm step to reach time_next
                            pcm.dt = time_next - pcm.time;
                            time_next_[k+index_offset_group_] = time_next;
                            // recored index in index_dt_sort of the interrupt case for moving later
                            interrupt_index_dt_group_list[n_interrupt_change_dt++] = i;
                        }
                    }
                    
                    // correct cm potential energy 
                    //Float dm = correctMassChangePotEnergyBinaryIter(*interrupt_binary_.adr);
                    // notice the bin_root represent new c.m. in rest frame
                    auto& bink = groups[k].info.getBinaryTreeRoot();
                    Float dm = bink.mass - pcm.mass;
                    Float de_pot = force_[k+index_offset_group_].pot*dm;
                    energy_.de_cum += de_pot;
                    energy_.de_binary_interrupt += de_pot;
                    energy_sd_.de_cum += de_pot;
                    energy_sd_.de_binary_interrupt += de_pot;

                    // correct kinetic energy of cm
                    ASSERT(!groups[k].particles.isOriginFrame());
                   // first remove mass loss kinetic energy assuming original cm velocity
                    auto& vcm = pcm.vel;
                    auto& vbin = bink.vel;
                    auto& vbin_bk = groups[k].info.vcm_record;
                    //Float vbcm[3] = {vcm[0] + vbin_bk[0], vcm[1] + vbin_bk[1], vcm[2] + vbin_bk[2]};
                    Float de_kin = 0.5*dm*(vcm[0]*vcm[0]+vcm[1]*vcm[1]+vcm[2]*vcm[2]);
                    // then correct c.m. motion if the c.m. velocity is shifted, 
                    // this can be by adding m_cm(new) * v_cm(new;rest) \dot v_cm(old;origin)
                    Float dvbin[3] = {vbin[0] - vbin_bk[0], vbin[1] - vbin_bk[1], vbin[2] - vbin_bk[2]};
                    de_kin += bink.mass*(dvbin[0]*vcm[0]+dvbin[1]*vcm[1]+dvbin[2]*vcm[2]);
                    energy_.de_cum += de_kin;
                    energy_.de_binary_interrupt += de_kin;
                    energy_sd_.de_cum += de_kin;
                    energy_sd_.de_binary_interrupt += de_kin;

                    // back up vcm for later on perturbation kinetic energy correction
                    vbin_bk[0] = vbin[0];
                    vbin_bk[1] = vbin[1];
                    vbin_bk[2] = vbin[2];


                    // update particle dm, velocity should not change to be consistent with frame
                    pcm.mass += dm;
                    //particles.cm.mass += dm;

                    if (ar_manager->interrupt_detection_option==2) { // quit integration
                        interrupt_group_dt_sorted_group_index_ = i;
                        ASSERT(time_next - interrupt_binary_.time_now + ar_manager->time_error_max >= 0.0);
                        return interrupt_binary_;
                    }
                }

                ASSERT(abs(groups[k].getTime()-time_next)<=ar_manager->time_error_max);
            }

            // update index_dt_sorted_group_ due to the change of dt
            if (n_interrupt_change_dt>0) {
                // shift position starts from the last interrupt index in sorted list (right to left)
                int ishift_start = interrupt_index_dt_group_list[n_interrupt_change_dt-1]; 
                // back up group index
                int interrupt_index_group_list[n_interrupt_change_dt];
                interrupt_index_group_list[0] = index_dt_sorted_group_[ishift_start];
                // this is to record the offset need to shift index in sorted list
                int shift_offset=1;
                // scan from last interrupt index to first 
                for (int i=n_interrupt_change_dt-2; i>=0; i--) {
                    // shift left edge
                    int k = interrupt_index_dt_group_list[i];
                    // back up next index
                    interrupt_index_group_list[shift_offset] = index_dt_sorted_group_[k];
                    
                    for (int j=ishift_start; j>k+shift_offset; j--) {
                        index_dt_sorted_group_[j] = index_dt_sorted_group_[j-shift_offset];
                    }
                    // set new starting position
                    ishift_start = k+shift_offset;
                    // increase offset
                    shift_offset++;
                }
                ASSERT(shift_offset==n_interrupt_change_dt);
                // shift the range from the first k to beginning of sorted list
                for (int j=ishift_start; j>=shift_offset; j--) {
                    index_dt_sorted_group_[j] = index_dt_sorted_group_[j-shift_offset];
                }

                // add all interrupted group indice in front of index_dt_sorted_group_
                for (int i=0; i<n_interrupt_change_dt; i++) {
                    index_dt_sorted_group_[i] = interrupt_index_group_list[i];
                }
                n_act_group_  += n_interrupt_change_dt;
                n_act_group_ = std::min(n_act_group_, n_group_tot);

#ifdef HERMITE_DEBUG
                // check whether the new list is consistent with table_group_mask_ 
                int table_check_mask[n_group_tot];
                for (int i=0; i<n_group_tot; i++) table_check_mask[i]=0;
                for (int i=0; i<index_dt_sorted_group_.getSize(); i++) table_check_mask[index_dt_sorted_group_[i]]++;
                for (int i=0; i<n_group_tot; i++) {
                    ASSERT((table_check_mask[i]==1&&!table_group_mask_[i])||(table_check_mask[i]==0&&table_group_mask_[i]));
                }
                for (int i=0; i<n_act_group_; i++) {
                    ASSERT(time_next_[index_dt_sorted_group_[i]+index_offset_group_] == time_next);
                }
                if (n_act_group_<n_group_tot) {
                    ASSERT(time_next_[index_dt_sorted_group_[n_act_group_]+index_offset_group_] > time_next);
                }
#endif
            }

            // when no interrupt, clear interrupt records
            //interrupt_group_dt_sorted_group_index_ = -1;
            //interrupt_binary_.clear();

            return interrupt_binary_;
        }

        //! modify single particles due to external functions, update energy
        void modifySingleParticles() {
            int mod_index[n_act_single_];
            int n_mod=0;
            for (int i=0; i<n_act_single_; i++) {
                int k = index_dt_sorted_single_[i];
                // use pred particle to back up velocity, it will be update at predictionAll, thus safe to change
                auto& pk = particles[k];
                auto& rbk = pred_[k].pos;
                auto& vbk = pred_[k].vel;
                auto& mbk = pred_[k].mass;
                vbk[0] = pk.vel[0];
                vbk[1] = pk.vel[1];
                vbk[2] = pk.vel[2];
                mbk = pk.mass;

                int modified_flag=ar_manager->interaction.modifyOneParticle(pk, time_+time_offset_, getNextTime()+time_offset_);
                if (modified_flag) {
                    auto& v = pk.vel;
#ifdef HERMITE_DEBUG_PRINT
                    std::cerr<<"Modify one particle: modify_flag: "<<modified_flag
                             <<" time: "<<pk.time
                             <<" dt: "<<pk.dt
                             <<" mass_bk: "<<mbk
                             <<" mass_new: "<<pk.mass
                             <<" vel_bk: "<<vbk[0]<<" "<<vbk[1]<<" "<<vbk[2]
                             <<" vel_new: "<<v[0]<<" "<<v[1]<<" "<<v[2]
                             <<std::endl;
#endif
                    // correct potential energy 
                    Float de_pot = force_[k].pot*(pk.mass-mbk);
                    energy_.de_cum += de_pot;
                    energy_.de_modify_single += de_pot;
                    energy_sd_.de_cum += de_pot;
                    energy_sd_.de_modify_single += de_pot;
                    
                    // correct kinetic energy
                    // first calc original kinetic energy;
                    Float de_kin = -0.5*mbk*(vbk[0]*vbk[0]+vbk[1]*vbk[1]+vbk[2]*vbk[2]);
                    // then new energy
                    de_kin += 0.5*pk.mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
                    energy_.de_cum += de_kin;
                    energy_.de_modify_single += de_kin;
                    energy_sd_.de_cum += de_kin;
                    energy_sd_.de_modify_single += de_kin;

                    neighbors[k].initial_step_flag = true;

                    // use predictor as template particle with mass of dm
                    mbk = pk.mass-mbk;
                    rbk[0] = pk.pos[0];
                    rbk[1] = pk.pos[1];
                    rbk[2] = pk.pos[2];
                    mod_index[n_mod++] = k;

                    //update time next
                    time_next_[k] = pk.time + pk.dt;
                }
            }

            // fix over-corrected potential
            if (n_mod>1) {
                Float de_pot = 0.0;
                for (int i=0; i<n_mod; i++) {
                    for (int j=0; j<i; j++) {
                        de_pot += manager->interaction.calcEnergyPotSingleSingle(pred_[mod_index[i]],pred_[mod_index[j]]);
                    }
                }
                energy_.de_cum += de_pot;
                energy_.de_modify_single += de_pot;
                energy_sd_.de_cum += de_pot;
                energy_sd_.de_modify_single += de_pot;
            }
        }
        

        //! Integration single active particles and update steps 
        /*! Integrated to next time given by minimum step particle
        */
        void integrateSingleOneStepAct() {
            ASSERT(checkParams());
            ASSERT(!particles.isModified());
            ASSERT(initial_system_flag_);
            ASSERT(!modify_system_flag_);
            ASSERT(n_init_group_==0&&n_init_single_==0);
            ASSERT(ar_manager->interrupt_detection_option!=2||(ar_manager->interrupt_detection_option==2&&interrupt_binary_.status==AR::InterruptStatus::none));
            
            // get next time
            Float time_next = getNextTime();

            dt_limit_ = step.calcNextDtLimit(time_next);

            // prediction positions
            predictAll(time_next);

            // check resolve status
            checkGroupResolve(n_act_group_);
            writeBackResolvedGroupAndCreateJParticleList(true);

            int* index_single = index_dt_sorted_single_.getDataAddress();
            int* index_group = index_dt_sorted_group_.getDataAddress();

            calcAccJerkNBList(index_single, n_act_single_, index_group, n_act_group_);

            correctAndCalcDt4thList(index_single, n_act_single_, index_group, n_act_group_, dt_limit_);

            updateTimeNextList(index_single, n_act_single_, index_group, n_act_group_);

            // update time
            time_ = time_next;

            // profile
            profile.hermite_single_step_count += n_act_single_;
            profile.hermite_group_step_count += n_act_group_;
        }

        //! Integration a list of particle to current time (ingore dt)
        /*! Integrate a list of particle to current time.
          @param[in] _time_next: time to integrate (without offset)
          @param[in] _particle_index: particle index to integrate
          @param[in] _n_particle: number of particles
        */
        void integrateToTimeList(const Float _time_next,
                                 const int* _particle_index,
                                 const int _n_particle) {
            ASSERT(checkParams());
            ASSERT(!particles.isModified());
            ASSERT(initial_system_flag_);

            // adjust dt
            int index_single_select[_n_particle];
            int n_single_select =0;
            int index_group_select[groups.getSize()];
            int n_group_select =0;

            auto* ptcl = particles.getDataAddress();
            auto* group_ptr = groups.getDataAddress();
            for (int i=0; i<_n_particle; i++){
                const int k = _particle_index[i];
                if (k<index_offset_group_) {
                    if (table_single_mask_[k]) continue;
                    if (ptcl[k].time>=_time_next) continue;
                    ptcl[k].dt = _time_next - ptcl[k].time;
                    index_single_select[n_single_select++] = k;
                }
                else {
                    const int kg = k - index_offset_group_;
                    if (table_group_mask_[kg]) continue;
                    auto& pcm = group_ptr[kg].particles.cm;
                    if (pcm.time>=_time_next) continue;
                    pcm.dt = _time_next - pcm.time;
                    index_group_select[n_group_select++] = kg;
                }
            }


            if (n_single_select>0 || n_group_select>0) {

                calcAccJerkNBList(index_single_select, n_single_select, index_group_select, n_group_select);

                correctAndCalcDt4thList(index_single_select, n_single_select, index_group_select, n_group_select, dt_limit_);

                updateTimeNextList(index_single_select, n_single_select, index_group_select, n_group_select);
            }
        }

        //! sort time step array and select active particles
        /*! Make sure time_next_ is updated already
         */
        void sortDtAndSelectActParticle() {
            // sort single
            std::sort(index_dt_sorted_single_.getDataAddress(), index_dt_sorted_single_.getDataAddress()+n_act_single_, SortIndexDtSingle(time_next_.getDataAddress()));
            // sort group
            std::sort(index_dt_sorted_group_.getDataAddress(), index_dt_sorted_group_.getDataAddress()+n_act_group_, SortIndexDtGroup(time_next_.getDataAddress(), index_offset_group_));

            // get minimum next time from single and group
            time_next_min_ = NUMERIC_FLOAT_MAX;
            if (index_dt_sorted_single_.getSize()>0) time_next_min_ = time_next_[index_dt_sorted_single_[0]];
            if (index_dt_sorted_group_.getSize()>0) time_next_min_ = std::min(time_next_min_, time_next_[index_dt_sorted_group_[0]+index_offset_group_]);
            ASSERT(time_next_min_>0.0);

            // find active singles
            const int n_singles = index_dt_sorted_single_.getSize();
            for(n_act_single_=0; n_act_single_<n_singles; n_act_single_++){
                if (time_next_min_ < time_next_[index_dt_sorted_single_[n_act_single_]]) break;
            }

#ifdef HERMITE_DEBUG
            for (int i=0; i<n_singles-1; i++) {
                const int k= index_dt_sorted_single_[i];
                const int k1= index_dt_sorted_single_[i+1];
                ASSERT(particles[k].dt<=particles[k1].dt);
                ASSERT(time_next_[k]<=time_next_[k1]);
            }
#endif

            // find active groups
            const int n_groups = index_dt_sorted_group_.getSize();
            for(n_act_group_=0; n_act_group_<n_groups; n_act_group_++){
                if (time_next_min_ < time_next_[index_dt_sorted_group_[n_act_group_] + index_offset_group_]) break;
            }

#ifdef HERMITE_DEBUG
            for (int i=0; i<n_groups-1; i++) {
                const int k= index_dt_sorted_group_[i];
                const int k1= index_dt_sorted_group_[i+1];
                ASSERT(groups[k].particles.cm.dt<=groups[k1].particles.cm.dt);
                ASSERT(k<groups.getSize());
                ASSERT(k1<groups.getSize());
                ASSERT(time_next_[k+index_offset_group_]<=time_next_[k1+index_offset_group_]);
            }
#endif

            ASSERT(!(n_act_single_==0&&n_act_group_==0));

        }

        //! write back group members to particles
        void writeBackGroupMembers() {
            const int n_group = index_dt_sorted_group_.getSize();
            for (int i=0; i<n_group; i++) {
                const int k = index_dt_sorted_group_[i];
                ASSERT(table_group_mask_[k]==false);
                groups[k].template writeBackParticlesOriginFrame<Tparticle>(); 
            }
        }

        //! correct Etot slowdown reference due to the groups change
        /*! @param[in] _igroup: group index to accumulative de_sd_change_cum
         */
        void accumDESlowDownChangeBreakGroup(const int _igroup) {
            auto& groupi = groups[_igroup];
            Float de_binary_interrupt  = groupi.getDEChangeBinaryInterrupt();
            energy_.de_binary_interrupt += de_binary_interrupt;
            energy_.de_cum += de_binary_interrupt;
            energy_sd_.de_cum += groupi.getDESlowDownChangeCum();
            energy_sd_.de_binary_interrupt += groupi.getDESlowDownChangeBinaryInterrupt();
            // add the change due to the shutdown of slowdown 
            Float etot = groupi.getEkin() + groupi.getEpot();
            Float etot_sd = groupi.getEkinSlowDown() + groupi.getEpotSlowDown();
            energy_sd_.de_cum += etot - etot_sd;
            // add kinetic correction due to vcm change
            auto& vcm = groupi.particles.cm.vel;
            auto& bink = groupi.info.getBinaryTreeRoot();
            auto& vbin = bink.vel;
            auto& vbin_bk = groupi.info.vcm_record;
            Float dvbin[3] = {vbin[0] - vbin_bk[0], vbin[1] - vbin_bk[1], vbin[2] - vbin_bk[2]};
            Float de_kin = bink.mass*(dvbin[0]*vcm[0]+dvbin[1]*vcm[1]+dvbin[2]*vcm[2]);
            Float epert = manager->interaction.calcEnergyPertOneGroup(groupi, perturber);
            energy_.de_cum -= epert - de_kin;
            energy_sd_.de_cum -= epert - de_kin ;
        }

        //! correct Etot slowdown reference due to the groups change
        void accumDESlowDownChange() {
            const int n_group = index_dt_sorted_group_.getSize();
            ASSERT(n_group <= groups.getSize());
            for (int k=0; k<n_group; k++) {
                const int i = index_dt_sorted_group_[k];
                ASSERT(i<groups.getSize());
                ASSERT(i>=0);
                auto& gi = groups[i];

                // interrupt energy
                Float de_binary_interrupt  = gi.getDEChangeBinaryInterrupt();
                energy_.de_binary_interrupt += de_binary_interrupt;
                energy_.de_cum += de_binary_interrupt;
                gi.resetDEChangeBinaryInterrupt();

                // slowdown energy
                energy_sd_.de_cum += gi.getDESlowDownChangeCum();
                gi.resetDESlowDownChangeCum();
                // slowdown interrupt energy
                energy_sd_.de_binary_interrupt += gi.getDESlowDownChangeBinaryInterrupt();
                gi.resetDESlowDownChangeBinaryInterrupt();
            } 
        }

        // correct energy due to mass change
        /*
          @param[in] _bin: interrupted binary to correct
          \return total mass change
        Float correctMassChangePotEnergyBinaryIter(AR::BinaryTree<Tparticle>& _bin) {
            Float dm_tot = 0.0;
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) dm_tot += correctMassChangePotEnergyBinaryIter(*_bin.getMemberAsTree(k));
                else {
                    //int i = _bin.getMemberIndex(k);
                    //Float& pot = force_[i].pot;
                    auto* pi = _bin.getMember(k);
                    dm_tot += pi->dm;
                    //Float de_pot = pot*pi->dm;
                    de_sd_change_cum_ += de_pot;
                    de_sd_change_interrupt_ += de_pot;
                    pi->dm = 0.0;
                }
            }
            return dm_tot;
        }
         */

        //! calculate slowdown energy
        /*! @param[in] _initial_flag: if true, set energy reference
         */
        void calcEnergySlowDown(const bool _initial_flag = false) {
            writeBackGroupMembers();
            // use a tempare array instead of modify pred_. As a good design, the function of calc energy should not influence integration.
            Tparticle ptmp[particles.getSize()];
            const auto* ptcl = particles.getDataAddress();
            const auto* pred = pred_.getDataAddress();

            // in single case, if particle time is not the current time, used predicted particle
            const int n_single = index_dt_sorted_single_.getSize();
            for (int k=0; k<n_single; k++){
                const int i = index_dt_sorted_single_[k];
                if (ptcl[i].time == time_) ptmp[i] = ptcl[i];
                else ptmp[i] = pred[i];
            }

            // in binary case, if c.m. is not up to date, correct member pos and vel
            const int n_group = index_dt_sorted_group_.getSize();
            ASSERT(n_group <= groups.getSize());
            auto* group_ptr = groups.getDataAddress();
            for (int k=0; k<n_group; k++) {
                const int i = index_dt_sorted_group_[k];
                const auto& groupi = group_ptr[i];
                const auto& pcm = groupi.particles.cm;
                auto& predcm = pred[i+index_offset_group_];

                const int n_member = groupi.particles.getSize();
                for (int j=0; j<n_member; j++) {
                    const int kj = groupi.info.particle_index[j];
                    ASSERT(kj<particles.getSize());
                    ptmp[kj] = ptcl[kj];
                    // if not up to date, add difference of predicted and current cm pos and vel to predicted members 
                    if (pcm.time < time_ ) {
                        ptmp[kj].pos[0] += predcm.pos[0] - pcm.pos[0];
                        ptmp[kj].pos[1] += predcm.pos[1] - pcm.pos[1];
                        ptmp[kj].pos[2] += predcm.pos[2] - pcm.pos[2];
                        ptmp[kj].vel[0] += predcm.vel[0] - pcm.vel[0];
                        ptmp[kj].vel[1] += predcm.vel[1] - pcm.vel[1];
                        ptmp[kj].vel[2] += predcm.vel[2] - pcm.vel[2];
                    }
                }
            }            

            // since a part of particles are not up to date, used predict particles to calculate energy, 
            // notice pred_.size is larger than particles.size because group c.m. is in pred_.
            manager->interaction.calcEnergy(energy_, ptmp, particles.getSize(), groups.getDataAddress(), index_dt_sorted_group_.getDataAddress(), index_dt_sorted_group_.getSize(), perturber);
            accumDESlowDownChange();
            // slowdown energy
            energy_sd_.ekin = energy_.ekin;
            energy_sd_.epot = energy_.epot;
            energy_sd_.epert= energy_.epert;
            for (int i=0; i<groups.getSize(); i++) {
                auto& gi = groups[i];
                energy_sd_.ekin -= gi.getEkin();
                energy_sd_.epot -= gi.getEpot(); 

                energy_sd_.ekin += gi.getEkinSlowDown();
                energy_sd_.epot += gi.getEpotSlowDown();
            }
            if (_initial_flag) {
                energy_init_ref_ = energy_.getEtot();
            }
            energy_.etot_ref    = energy_init_ref_ + energy_.de_cum;
            energy_sd_.etot_ref = energy_init_ref_ + energy_sd_.de_cum;
        }

        //! get slowdown energy error
        Float getEnergyErrorSlowDown() const {
            return energy_sd_.getEnergyError();
        }

        //! get energy error 
        /*! \return energy error
         */
        Float getEnergyError() const {
            return energy_.getEnergyError();
        }

        //! Get current total energy from ekin and epot
        /*! \return total integrated energy 
         */
        Float getEtot() const {
            return energy_.getEtot();
        }

        //! Get current total integrated energy 
        /*! \return total integrated energy 
         */
        Float getEtotRef() const {
            return energy_.etot_ref;
        }

        //! Get current total energy with slowdown from ekin_sd and epot_sd
        /*! \return total energy with slowdown
         */
        Float getEtotSlowDown() const {
            return energy_sd_.getEtot();
        }

        //! Get current total integrated energy with inner slowdown
        /*! \return total integrated energy with inner slowdown
         */
        Float getEtotSlowDownRef() const {
            return energy_sd_.etot_ref;
        }

        //! Get current kinetic energy
        /*! \return current kinetic energy
         */
        Float getEkin() const {
            return energy_.ekin;
        }

        //! Get current kinetic energy with slowdown
        /*! \return current slowdown kinetic energy
         */
        Float getEkinSlowDown() const {
            return energy_sd_.ekin;
        }

        //! Get current potential energy
        /*! \return current potetnial energy (negative value for bounded systems)
         */
        Float getEpot() const {
            return energy_.epot;
        }

        //! Get current potential energy with slowdown
        /*! \return current potetnial energy with slowdown (negative value for bounded systems)
         */
        Float getEpotSlowDown() const {
            return energy_sd_.epot;
        }

        //! get cumulative slowdown energy change due to slowdown change
        Float getDESlowDownChangeCum() const {
            return energy_sd_.de_cum;
        }

        //! reset cumulative slowdown energy change due to slowdown change
        void resetDESlowDownChangeCum()  {
            energy_sd_.de_cum = 0.0;
        }

        //! get cumulative slowdown energy change due to interruption
        Float getDESlowDownChangeBinaryInterrupt() const {
            return energy_sd_.de_binary_interrupt;
        }

        //! reset cumulative slowdown energy change due to interruption
        void resetDESlowDownChangeBinaryInterrupt()  {
            energy_sd_.de_binary_interrupt = 0.0;
        }

        //! get cumulative slowdown energy change due to modification of one particle
        Float getDESlowDownChangeModifySingle() const {
            return energy_sd_.de_modify_single;
        }

        //! reset cumulative slowdown energy change due to modification of one particle
        void resetDESlowDownChangeModifySingle()  {
            energy_sd_.de_modify_single = 0.0;
        }

        //! get cumulative energy change 
        Float getDEChangeCum() const {
            return energy_.de_cum;
        }

        //! reset cumulative energy change
        void resetDEChangeCum()  {
            energy_.de_cum = 0.0;
        }

        //! get cumulative energy change due to interruption
        Float getDEChangeBinaryInterrupt() const {
            return energy_.de_binary_interrupt;
        }

        //! reset cumulative energy change due to interruption
        void resetDEChangeBinaryInterrupt()  {
            energy_.de_binary_interrupt = 0.0;
        }

        //! get cumulative energy change due to modification of one particle
        Float getDEChangeModifySingle() const {
            return energy_.de_modify_single;
        }

        //! reset cumulative energy change due to modification of one particle
        void resetDEChangeModifySingle()  {
            energy_.de_modify_single = 0.0;
        }

        //! get next time for intergration
        Float getNextTime() const {
            return time_next_min_;
        }

        //! get current integration time 
        Float getTimeInt() const {
            return time_;
        }

        //! get current real time
        Float getTime() const {
            return time_+time_offset_;
        }

        //! set time offset
        void setTimeOffset(const Float& _time_offset) {
            time_offset_ = _time_offset;
        }

        //! get interrupt group index
        /*! if not exist, return -1
         */
        int getInterruptGroupIndex() const {
            if (interrupt_group_dt_sorted_group_index_>=0) 
                return index_dt_sorted_group_[interrupt_group_dt_sorted_group_index_];
            else return -1;
        }

        //! get interrupt binary (tree) address
        /*! if not exist, return NULL
         */
        AR::InterruptBinary<Tparticle>& getInterruptBinary() {
            return interrupt_binary_;
        }

        //! get active number of particles of singles
        int getNActSingle() const{
            return n_act_single_;
        }

        //! get active number of particles of groups
        int getNActGroup() const{
            return n_act_group_;
        }

        //! get initial number of particles of groups
        int getNInitGroup() const{
            return n_init_group_;
        }

        //! get initial number of particles of singles
        int getNInitSingle() const{
            return n_init_single_;
        }

        //! get N groups
        int getNGroup() const {
            return index_dt_sorted_group_.getSize();
        }

        //! get N single
        int getNSingle() const {
            return index_dt_sorted_single_.getSize();
        }

        //! get group index offset (group index + N_single_max)
        int getIndexOffsetGroup() const {
            return index_offset_group_;
        }

        //! get sorted dt index of singles
        int* getSortDtIndexSingle() {
            return index_dt_sorted_single_.getDataAddress();
        }

        //! get sorted dt index of groups
        int* getSortDtIndexGroup() {
            return index_dt_sorted_group_.getDataAddress();
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width
          @param[in] _n_sd_list: AR inner slowdown numbers per group (list)
          @param[in] _n_group: total AR group number
          @param[in] _n_sd_tot: total slowdown numbers
        */
        void printColumnTitle(std::ostream & _fout, const int _width, const int _n_sd_list[], const int _n_group, const int _n_sd_tot) {
            _fout<<std::setw(_width)<<"Time"
                 <<std::setw(_width)<<"dE"
                 <<std::setw(_width)<<"Etot_ref"
                 <<std::setw(_width)<<"Ekin"
                 <<std::setw(_width)<<"Epot"
                 <<std::setw(_width)<<"Epert"
                 <<std::setw(_width)<<"dE_cum"
                 <<std::setw(_width)<<"dE_intr"
                 <<std::setw(_width)<<"dE_mod"
                 <<std::setw(_width)<<"dE_SD"
                 <<std::setw(_width)<<"Etot_SD_ref"
                 <<std::setw(_width)<<"Ekin_SD"
                 <<std::setw(_width)<<"Epot_SD"
                 <<std::setw(_width)<<"Epert_SD"
                 <<std::setw(_width)<<"dE_SD_cum"
                 <<std::setw(_width)<<"dE_SD_intr"
                 <<std::setw(_width)<<"dE_SD_mod"
                 <<std::setw(_width)<<"N_SD";
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            AR::SlowDown sd_empty;
            int n_sd_count = 0;
            for (int i=0; i<_n_group; i++) {
                n_sd_count += _n_sd_list[i];
                for (int j=0; j<_n_sd_list[i]; j++) 
                    sd_empty.printColumnTitle(_fout, _width);
            }
            ASSERT(_n_sd_tot == n_sd_count);
#endif
            perturber.printColumnTitle(_fout, _width);
            info.printColumnTitle(_fout, _width);
            profile.printColumnTitle(_fout, _width);
            particles.printColumnTitle(_fout, _width);
        }    

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width
          @param[in] _n_sd_list: AR inner slowdown numbers per group (list)
          @param[in] _n_group: total AR group number
          @param[in] _n_sd_tot: total slowdown numbers
        */
        void printColumn(std::ostream & _fout, const int _width, const int _n_sd_list[], const int _n_group, const int _n_sd_tot){
            _fout<<std::setw(_width)<<time_;
            energy_.printColumn(_fout, _width);
            energy_sd_.printColumn(_fout, _width);
            _fout<<std::setw(_width)<<_n_sd_tot;
            AR::SlowDown sd_empty;
            int n_group_now = groups.getSize();
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            int n_sd_count = 0;
            for (int i=0; i<_n_group; i++) {
                n_sd_count += _n_sd_list[i];
                if (i<n_group_now) {
                    auto & gi = groups[i];
#ifdef AR_SLOWDOWN_ARRAY
                    int n_sd_in = gi.binary_slowdown.getSize();
                    for (int j=0; j<_n_sd_list[i]; j++) {
                        if (j<n_sd_in) gi.binary_slowdown[j]->slowdown.printColumn(_fout, _width);
                        else sd_empty.printColumn(_fout, _width);
                    }
#else
                    int n_sd_in = gi.info.binarytree.getSize();
                    for (int j=0; j<_n_sd_list[i]; j++) {
                        if (j<n_sd_in) gi.info.binarytree[j].slowdown.printColumn(_fout, _width);
                        else sd_empty.printColumn(_fout, _width);
                    }
#endif
                }
                else {
                    for (int j=0; j<_n_sd_list[i]; j++) sd_empty.printColumn(_fout, _width);
                    sd_empty.printColumn(_fout, _width);
                }
            }
            ASSERT(_n_sd_tot == n_sd_count);
#endif
            perturber.printColumn(_fout, _width);
            info.printColumn(_fout, _width);
            profile.printColumn(_fout, _width);
            particles.printColumn(_fout, _width);
        }        

        //! print step histogram
        void printStepHist(){
            std::map<Float, int> stephist;
            for(int i=0; i<index_dt_sorted_single_.getSize(); i++) {
                int k = index_dt_sorted_single_[i];
                Float dt = particles[k].dt;
                std::map<Float, int>::iterator p = stephist.find(dt);
                if (p==stephist.end()) stephist[dt] = 1;
                else stephist[dt]++;
            }
            for(int i=0; i<index_dt_sorted_group_.getSize(); i++) {
                int k = index_dt_sorted_group_[i];
                Float dt=groups[k].particles.cm.dt;
                std::map<Float, int>::iterator p = stephist.find(dt);
                if (p==stephist.end()) stephist[dt] = 1;
                else stephist[dt]++;
            }
            std::cerr<<"Step hist: time = "<<time_<<"\n";
            for(auto i=stephist.begin(); i!=stephist.end(); i++) {
                std::cerr<<std::setw(24)<<i->first;
            }
            std::cerr<<std::endl;
            for(auto i=stephist.begin(); i!=stephist.end(); i++) {
                std::cerr<<std::setw(24)<<i->second;
            }
            std::cerr<<std::endl;
        }

    };
}
