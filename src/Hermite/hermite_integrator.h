#pragma once

#include "Common/Float.h"
#include "Common/list.h"
#include "AR/symplectic_integrator.h"
#include "Hermite/ar_information.h"
#include "Hermite/hermite_particle.h"
#include "Hermite/block_time_step.h"
#include "Hermite/neighbor.h"
#include <map>

namespace H4{
    //! Hermite manager class
    template <class Tmethod>
    class HermiteManager{
    public:
        Float r_break_crit; ///> the distance criterion to break the groups
        Float r_neighbor_crit; ///> the distance for neighbor search
        Tmethod interaction; ///> class contain interaction function
        BlockTimeStep4th step; ///> time step calculator

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
            interaction.writeBinary(_fp);
            step.writeBinary(_fp);
        }

        //! read class data to file with binary format
        /*! @param[in] _fp: FILE type file for reading
         */
        void readBinary(FILE *_fin) {
            interaction.readBinary(_fin);
            step.readBinary(_fin);
        }

    };

    //!Hermite integrator class
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tarmethod, class Tinfo>
    class HermiteIntegrator{
    private:
        typedef AR::SymplecticIntegrator<ParticleAR<Tparticle>, ParticleH4<Tparticle>, Neighbor<Tparticle>, Tarmethod, ARInformation<Tparticle>> ARSym;

        // time 
        Float time_;   ///< integrated time (not real physical time if slowdown is on)
        Float time_next_min_; ///< the next minimum time 
        Float dt_limit_;  // maximum step size allown for next integration step calculation

        // active particle number
        int n_act_single_;    /// active particle number of singles
        int n_act_group_;    /// active particle number of groups

        // initial particle number
        int n_init_single_;  /// number of singles need to initial
        int n_init_group_;   /// number of groups need to initial

        // group offset
        int index_offset_group_; /// offset of group index in pred_, force_ and time_next_

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
        HermiteManager<Tmethod>* manager; ///< integration manager
        AR::SymplecticManager<Tarmethod>* ar_manager; ///< integration manager
        COMM::ParticleGroup<ParticleH4<Tparticle>, Tpcm> particles; // particles
        COMM::List<ARSym> groups; // integrator for sub-groups
        COMM::List<Neighbor<Tparticle>> neighbors; // neighbor information of particles
        Tpert perturber; // external perturberx
        Tinfo info; ///< information of the system

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
                const Float dt =manager->step.calcBlockDt2nd(acc0, acc1, _dt_limit);
                particles[k].dt = dt;
            }
            // for group
            for (int i=0; i<_n_group; i++) {
                const int k = _index_group[i];
                const int kf = k + index_offset_group_;
                const Float* acc0 = force_[kf].acc0;
                const Float* acc1 = force_[kf].acc1;
                const Float dt =manager->step.calcBlockDt2nd(acc0, acc1, _dt_limit);
                groups[k].particles.cm.dt = dt;
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
        void correctAndCalcDt4thOne(ParticleH4<Tparticle>& _pi, ForceH4& _fi, const Float _dt_limit, const bool _init_step_flag) {
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
            const Float dt_old = _pi.dt;
            if(_init_step_flag) _pi.dt = manager->step.calcBlockDt2nd(_pi.acc0, _pi.acc1, _dt_limit);
            else _pi.dt = manager->step.calcBlockDt4th(_pi.acc0, _pi.acc1, acc2, acc3, _dt_limit);

            ASSERT((dt_old > 0.0 && _pi.dt >0.0));
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
                correctAndCalcDt4thOne(groupi.particles.cm, force[k+index_offset_group_], _dt_limit, groupi.perturber.initial_step_flag);
                groupi.perturber.initial_step_flag = false;
                //groupi.correctCenterOfMassDrift();
            }
        }

        //! check whether group need to resolve
        /*! @param[in] _n_group: number of groups in index_dt_sorted_group_ to check
         */
        void checkGroupResolve(const int _n_group) {
            ASSERT(_n_group<=index_dt_sorted_group_.getSize());
            for (int i=0; i<_n_group; i++) {
                const int k = index_dt_sorted_group_[i];
                auto& groupk = groups[k];
                const Float kappa = groupk.slowdown.getSlowDownFactor();
                groupk.perturber.checkGroupResolve(kappa);
            }
        }

        //! Generate j particle list from groups for force calculation
        /*! check group status, if the components are need to resolve, writeback the data with slowdown velocity, save the index to _index_group_resolve list, otherwisde save to _index_group_cm list
         */
        void writeBackResolvedGroupAndCreateJParticleList() {
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
                    group_k.writeBackSlowDownParticles(ptcl[k+index_offset_group_]);
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
          @param[in,out] _break_group_index:   break group index in groups
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
                                        int* _break_group_index,
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
                    _break_group_index[_n_break+_n_break_no_add] = insert_index - index_offset_group_;
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
                    _break_group_index[_n_break+_n_break_no_add] = _i - index_offset_group_;
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
                    _break_group_index[_n_break+_n_break_no_add] = _j - index_offset_group_;
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
          @param[in] _pid: particle id (used to avoid self-interaction)
        */
        template <class Tpi>
        inline void calcOneSingleAccJerkNB(ForceH4 &_fi, 
                                           Neighbor<Tparticle> &_nbi,
                                           const Tpi &_pi,
                                           const int _pid) {
            // clear force
            _fi.clear();
            _nbi.neighbor_address.resizeNoInitialize(0);

            // single list
            const int* single_list = index_dt_sorted_single_.getDataAddress();
            const int n_single = index_dt_sorted_single_.getSize();
            auto* ptcl = pred_.getDataAddress();
            for (int i=0; i<n_single; i++) {
                const int j = single_list[i];
                ASSERT(j<pred_.getSize());
                const auto& pj = ptcl[j];
                if (_pid==pj.id) continue;
                Float r2 = manager->interaction.calcAccJerkPair(_fi, _pi, pj);
                ASSERT(r2>0.0);
                if (r2<_nbi.r_crit_sq) _nbi.neighbor_address.addMember(NBAdr<Tparticle>(&particles[j],j));
                // mass weighted nearest neigbor
                if (r2*_nbi.r_min_mass < _nbi.r_min_sq*pj.mass) {
                    _nbi.r_min_sq = r2;
                    _nbi.r_min_index = j;
                    _nbi.r_min_mass = pj.mass;
                }
                // minimum mass
                if(_nbi.mass_min>pj.mass) {
                    _nbi.mass_min = pj.mass;
                    _nbi.mass_min_index = j;
                }
                // if neighbor step need update, also set update flag for i particle
                if(neighbors[j].initial_step_flag) _nbi.initial_step_flag = true;
            }

            auto* group_ptr = groups.getDataAddress();

            // resolved group list
            const int n_group_resolve = index_group_resolve_.getSize();
            for (int i=0; i<n_group_resolve; i++) {
                const int j =index_group_resolve_[i];
                auto& groupj = group_ptr[j];
                if (_pid==groupj.particles.cm.id) continue;
                const int n_member = groupj.particles.getSize();
                auto* member_adr = groupj.particles.getOriginAddressArray();
                bool nb_flag = false;
                Float r2_min = NUMERIC_FLOAT_MAX;
                for (int k=0; k<n_member; k++) {
                    ASSERT(_pi.id!=member_adr[k]->id);
                    const auto& pj = *member_adr[k];
                    Float r2 = manager->interaction.calcAccJerkPair(_fi, _pi, pj);
                    ASSERT(r2>0.0);
                    if (r2<_nbi.r_crit_sq) nb_flag = true;
                    r2_min = std::min(r2_min, r2);
                }
                if (nb_flag) _nbi.neighbor_address.addMember(NBAdr<Tparticle>(&groupj.particles.cm, j+index_offset_group_));
                // mass weighted nearest neigbor
                const Float cm_mass = groupj.particles.cm.mass;
                if (r2_min*_nbi.r_min_mass < _nbi.r_min_sq*cm_mass) {
                    _nbi.r_min_sq = r2_min;
                    _nbi.r_min_index = j+index_offset_group_;
                    _nbi.r_min_mass = cm_mass;
                }
                // minimum mass
                if(_nbi.mass_min > cm_mass) {
                    _nbi.mass_min = cm_mass;
                    _nbi.mass_min_index = j+index_offset_group_;
                }
                // if neighbor step need update, also set update flag for i particle
                if(groupj.perturber.initial_step_flag) _nbi.initial_step_flag = true;
            }

            // cm group list
            const int n_group_cm = index_group_cm_.getSize();
            for (int i=0; i<n_group_cm; i++) {
                const int j = index_group_cm_[i];
                auto& groupj = group_ptr[j];
                if (_pid==groupj.particles.cm.id) continue;
                // used predicted particle instead of original cm
                const auto& pj = ptcl[j+index_offset_group_];
                ASSERT(j+index_offset_group_<pred_.getSize());
                ASSERT(_pi.id!=pj.id);
                Float r2 = manager->interaction.calcAccJerkPair(_fi, _pi, pj);
                ASSERT(r2>0.0);
                if (r2<_nbi.r_crit_sq) _nbi.neighbor_address.addMember(NBAdr<Tparticle>(&groupj.particles.cm,j+index_offset_group_));
                // mass weighted nearest neigbor
                if (r2*_nbi.r_min_mass < _nbi.r_min_sq*pj.mass) {
                    _nbi.r_min_sq = r2;
                    _nbi.r_min_index = j+index_offset_group_;
                    _nbi.r_min_mass = pj.mass;
                }
                // minimum mass
                if(_nbi.mass_min>pj.mass) {
                    _nbi.mass_min = pj.mass;
                    _nbi.mass_min_index = j+index_offset_group_;
                }
                // if neighbor step need update, also set update flag for i particle
                if(groupj.perturber.initial_step_flag) _nbi.initial_step_flag = true;
            }
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
                auto& nbi = groupi.perturber;
                calcOneSingleAccJerkNB(fi, nbi, pi, pi.id);

                // for resolved case
                if(groupi.perturber.need_resolve_flag) {
                    auto* member_adr = groupi.particles.getOriginAddressArray();
                    const int n_member = groupi.particles.getSize();
                    fi.clear();
#ifdef HERMITE_DEBUG
                    Float mcm = 0.0;
#endif
                    for (int j=0; j<n_member; j++) {
                        // get particle index from memory address offset
                        //const int kj = int((Tparticle*)member_adr[j]-ptcl);
                        const int kj = groupi.info.particle_index[j];
                        ASSERT(kj<particles.getSize());
                        ASSERT(kj>=0);
                        auto& pkj = *(Tparticle*)member_adr[j];
                        ASSERT((void*)member_adr[j]==(void*)&particles[kj]);
                        auto& fkj = force_ptr[kj];
                        auto& nbkj = neighbor_ptr[kj];
                        calcOneSingleAccJerkNB(fkj, nbkj, pkj, pi.id);
                        // replace the cm. force with the summation of members
                        fi.acc0[0] += pkj.mass * fkj.acc0[0]; 
                        fi.acc0[1] += pkj.mass * fkj.acc0[1]; 
                        fi.acc0[2] += pkj.mass * fkj.acc0[2]; 

                        fi.acc1[0] += pkj.mass * fkj.acc1[0]; 
                        fi.acc1[1] += pkj.mass * fkj.acc1[1]; 
                        fi.acc1[2] += pkj.mass * fkj.acc1[2]; 
#ifdef HERMITE_DEBUG
                        mcm += pkj.mass;
#endif
                    }
#ifdef HERMITE_DEBUG
                    ASSERT(abs(mcm-pi.mass)<1e-10);
#endif
                    fi.acc0[0] /= pi.mass;
                    fi.acc0[1] /= pi.mass;
                    fi.acc0[2] /= pi.mass;

                    fi.acc1[0] /= pi.mass;
                    fi.acc1[1] /= pi.mass;
                    fi.acc1[2] /= pi.mass;
                }
            }
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
        Float calcDrDv(const ParticleH4<Tparticle>& _p1, const ParticleH4<Tparticle>& _p2) {
            Float dx[3],dv[3];
            const Float* pos1 = _p1.pos;
            const Float* pos2 = _p2.pos;
            const Float* vel1 = _p1.vel;
            const Float* vel2 = _p2.vel;
            dx[0] = pos1[0] - pos2[0];
            dx[1] = pos1[1] - pos2[1];
            dx[2] = pos1[2] - pos2[2];

            dv[0] = vel1[0] - vel2[0];
            dv[1] = vel1[1] - vel2[1];
            dv[2] = vel1[2] - vel2[2];
        
            Float drdv= dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];

            return drdv;
        }

        //! estimate perturbation / inner ratio square for a pair
        /*!
          @param[in] _r2: separation square
          @param[in] _p1: particle 1
          @param[in] _p2: particle 2
         */
        Float calcPertInnerRatioSq(const Float _dr2, const ParticleH4<Tparticle>& _p1, const ParticleH4<Tparticle>& _p2) {
            Float fcm[3] = {_p1.mass*_p1.acc0[0] + _p2.mass*_p2.acc0[0], 
                            _p1.mass*_p1.acc0[1] + _p2.mass*_p2.acc0[1], 
                            _p1.mass*_p1.acc0[2] + _p2.mass*_p2.acc0[2]};
            Float fcm2 = fcm[0]*fcm[0] + fcm[1]*fcm[1] + fcm[2]*fcm[2];
            Float fin2 = _p1.mass*_p2.mass/_dr2;
            Float fratiosq = fcm2/fin2;
            return fratiosq;
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
        void removeNBGroupAddressFromTable(COMM::List<NBAdr<Tparticle>> & _adr) {
            int k = 0;
            int k_last = _adr.getSize()-1;
            while (k<=k_last) {
                if (_adr[k].index>=index_offset_group_) {
                    const int kg = _adr[k].index-index_offset_group_;
                    if (table_group_mask_[kg]) {
                        _adr[k] = _adr[k_last];
                        k_last--;
                    }
                }
                else k++;
            }
            _adr.resizeNoInitialize(k_last+1);
        }

        // ! remove group particle address of a list from table_group_mask_
        void removeNBSingleAddressFromTable(COMM::List<NBAdr<Tparticle>> & _adr) {
            int k = 0;
            int k_last = _adr.getSize()-1;
            while (k<=k_last) {
                if (_adr[k].index<index_offset_group_) {
                    if (table_single_mask_[k]) {
                        _adr[k] = _adr[k_last];
                        k_last--;
                    }
                }
                else k++;
            }
            _adr.resizeNoInitialize(k_last+1);
        }
        

    public:
        //! constructor
        HermiteIntegrator(): time_(Float(0.0)), time_next_min_(Float(0.0)), 
                             n_act_single_(0), n_act_group_(0), 
                             n_init_single_(0), n_init_group_(0), 
                             index_offset_group_(0), 
                             initial_system_flag_(false), modify_system_flag_(false),
                             index_dt_sorted_single_(), index_dt_sorted_group_(), 
                             index_group_resolve_(), index_group_cm_(), 
                             pred_(), force_(), time_next_(), 
                             index_group_mask_(), table_group_mask_(), table_single_mask_(), 
                             manager(NULL), ar_manager(NULL), particles(), groups(), neighbors(), perturber(), info() {}

        //! clear function
        void clear() {
            time_ = 0.0;
            time_next_min_ = 0.0;
            n_act_single_ = n_act_group_ = n_init_single_ = n_init_group_ = 0;
            index_offset_group_ = 0;
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
            index_group_resolve_.setMode(COMM::ListMode::local);
            index_group_cm_.setMode(COMM::ListMode::local);
            index_group_mask_.setMode(COMM::ListMode::local);
            table_group_mask_.setMode(COMM::ListMode::local);
            table_single_mask_.setMode(COMM::ListMode::local);
            pred_.setMode(COMM::ListMode::local);
            force_.setMode(COMM::ListMode::local);
            time_next_.setMode(COMM::ListMode::local);
            neighbors.setMode(COMM::ListMode::local);
            
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

        //! Initial system array
        /*!@param[in] _time_sys: current set time
         */
        void initialSystemSingle(const Float _time_sys) {
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
                // set neighbor radius
                neighbors[i].r_crit_sq = manager->r_neighbor_crit*manager->r_neighbor_crit;
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
                group_new.perturber.r_crit_sq = manager->r_neighbor_crit*manager->r_neighbor_crit;
                group_new.info.reserveMem(n_particle);
            
                // Add members to AR 
                for(int j=_n_group_offset[i]; j<_n_group_offset[i+1]; j++) {
                    const int p_index = _particle_index[j];
                    ASSERT(p_index<particles.getSize());
                    group_new.particles.addMemberAndAddress(particles[p_index]);
                    group_new.info.particle_index.addMember(p_index);
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
                group_new.info.generateBinaryTree(group_new.particles);

                // initial perturber
                group_new.perturber.need_resolve_flag = false;
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
            // clear neighbor lists
            int n_group = index_dt_sorted_group_.getSize();
            for (int i=0; i<n_group; i++) {
                const int k = index_dt_sorted_group_[i];
                removeNBSingleAddressFromTable(groups[k].perturber.neighbor_address);
            }
        }

        //! Add group based on a configure file
        /*! @param[in] _fin: std::istream IO for read
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
           @param[in] _break_group_index:   break group index in groups
           @param[in] _n_break_no_add: number of break groups without adding new particles
           @param[in] _n_break:   number of break groups
         */
        void breakGroups(int* _new_group_particle_index_origin, 
                         int* _new_n_group_offset, 
                         int& _new_n_group,
                         const int* _break_group_index, 
                         const int _n_break_no_add,
                         const int _n_break) {
            ASSERT(!particles.isModified());
            ASSERT(initial_system_flag_);
            ASSERT(_n_break<=index_dt_sorted_group_.getSize());

            if (_n_break==0 && _n_break_no_add==0) return;

            modify_system_flag_=true;
            int new_index_single[particles.getSize()];
            int n_single_new=0;

            for (int k=0; k<_n_break; k++) {
                const int i = _break_group_index[k];
                ASSERT(i>=0&&i<groups.getSize());
                auto& groupi = groups[i];
                const int n_member =groupi.particles.getSize();
                int particle_index_origin[n_member];
                int ibreak = groupi.info.getTwoBranchParticleIndexOriginFromBinaryTree(particle_index_origin, groupi.particles.getDataAddress());
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
                    new_index_single[n_single_new++] = particle_index_origin[0];
                }
                else {
                    int offset = _new_n_group_offset[_new_n_group];
                    for (int j=0; j<ibreak; j++) {
                        ASSERT(table_single_mask_[particle_index_origin[j]]==true);
                        table_single_mask_[particle_index_origin[j]] = false;
                        _new_group_particle_index_origin[offset++] = particle_index_origin[j];
                    }
                    _new_n_group_offset[++_new_n_group] = offset;
                }
                // right side
                if (n_member-ibreak==1) {
                    ASSERT(particle_index_origin[ibreak]<particles.getSize());
                    ASSERT(table_single_mask_[particle_index_origin[ibreak]]==true);
                    table_single_mask_[particle_index_origin[ibreak]] = false;
                    new_index_single[n_single_new++] = particle_index_origin[ibreak];
                }
                else {
                    int offset = _new_n_group_offset[_new_n_group];
                    for (int j=ibreak; j<n_member; j++) {
                        ASSERT(table_single_mask_[particle_index_origin[j]]==true);
                        table_single_mask_[particle_index_origin[j]] = false;
                        _new_group_particle_index_origin[offset++] = particle_index_origin[j];
                    }
                    _new_n_group_offset[++_new_n_group] = offset;
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
            for (int i=_n_break; i<_n_break+_n_break_no_add; i++) {
                auto& groupi = groups[i];
                groupi.particles.shiftToOriginFrame();
                groupi.particles.template writeBackMemberAll<Tparticle>();
                const int n_member = groupi.particles.getSize();
                ASSERT(n_member==groupi.info.particle_index.getSize());
                for (int j=0; j<n_member; j++) {
                    const int k = groupi.info.particle_index[j];
                    ASSERT(k<particles.getSize());
                    table_single_mask_[k] = false;
                }
                groupi.clear();
                table_group_mask_[i] = true;
            }
            
            // clear index_dt_sorted_
            removeDtIndexFromTable(index_dt_sorted_group_, n_act_group_, n_init_group_, table_group_mask_);
            // clear neighbor lists
            int n_group = index_dt_sorted_group_.getSize();
            for (int i=0; i<n_group; i++) {
                const int k = index_dt_sorted_group_[i];
                removeNBGroupAddressFromTable(groups[k].perturber.neighbor_address);
            }
        }

        //! Check break condition
        /*! Check whether it is necessary to break the chain\n
          1. Inner distance criterion:
          If r>r_break_crit, and ecca>0.0, break.\n
          2. Perturbation criterion for closed orbit:
          If inner slowdown factor >1.0 break.\n
          @param[out] _break_group_list: group index list to break
          @parma[out] _n_break: number of groups need to break
        */
        void checkBreak(int* _break_group_index, int& _n_break) {
            const int n_group_tot = index_dt_sorted_group_.getSize();
            for (int i=0; i<n_group_tot; i++) {
                const int k = index_dt_sorted_group_[i];
                ASSERT(table_group_mask_[k]==false);
                auto& groupk = groups[k];

                const int n_member = groupk.particles.getSize();
                // generate binary tree
                groupk.info.generateBinaryTree(groupk.particles);
                
                auto& bin_root = groupk.info.binarytree.getLastMember();

                // check distance criterion and outcome (ecca>0) or income (ecca<0)
                if (bin_root.r > manager->r_break_crit && bin_root.ecca>0.0) {
#ifdef ADJUST_GROUP_DEBUG
                    std::cerr<<"Break group: escape, i_group: "<<k<<" N_member: "<<n_member<<" ecca: "<<bin_root.ecca<<" separation : "<<bin_root.r<<" r_crit: "<<manager->r_break_crit<<std::endl;
#endif
                    _break_group_index[_n_break++] = k;
                    continue;
                }
                    
                // check strong perturbed binary case
                Float kappa_org = groupk.slowdown.getSlowDownFactorOrigin();
                if (kappa_org<0.01 && bin_root.semi>0) {
#ifdef ADJUST_GROUP_DEBUG
                    std::cerr<<"Break group: strong perturbed, i_group: "<<k<<" N_member: "<<n_member<<" kappa_org: "<<kappa_org<<" semi: "<<bin_root.semi<<std::endl;
#endif
                    _break_group_index[_n_break++] = k;
                    continue;
                }

                // check few-body inner perturbation
                if (n_member>2 && bin_root.ecca>0.0) {
                    AR::SlowDown sd;
                    sd.initialSlowDownReference(groupk.slowdown.getSlowDownFactorReference(),groupk.slowdown.getSlowDownFactorMax());
                    for (int j=0; j<2; j++) {
                        if (bin_root.getMember(j)->id<0) {
                            auto* bin_sub = (COMM::BinaryTree<Tparticle>*) bin_root.getMember(j);
                            Float semi_db = 2.0*bin_sub->semi;
                            // inner hyperbolic case
                            if(semi_db<0.0) {
#ifdef ADJUST_GROUP_DEBUG
                                std::cerr<<"Break group: inner member hyperbolic, i_group: "<<k<<" i_member: "<<j<<" semi: "<<semi_db<<std::endl;
#endif
                                _break_group_index[_n_break++] = k;
                                break;
                            }
                            // if slowdown factor is large, break the group
                            Float f_in = bin_sub->m1*bin_sub->m2/(semi_db*semi_db*semi_db);
                            Float f_out = bin_sub->mass*bin_root.getMember(1-j)->mass/(bin_root.r*bin_root.r*bin_root.r);
                            Float kappa_in = sd.calcSlowDownFactor(f_in, f_out);
                            if (kappa_in>1.0) {
#ifdef ADJUST_GROUP_DEBUG
                                std::cerr<<"Break group: inner kappa large, i_group: "<<k<<" i_member: "<<j<<" kappa_in:"<<kappa_in<<std::endl;
#endif
                                _break_group_index[_n_break++] = k;
                                break;
                            }
                        }
                    }
                }
            }
        }

        //! check the pair with distance below r_crit for ptcl in adr_dt_sorted_
        /*! First check nearest neighbor distance r_min
          If r_min<r_crit, check the direction, if income, accept as group
          @param[out] _new_group_particle_index_origin: particle index (for hermite particle data) of new group members
          @param[out] _new_n_group_offset:  group member boundary, first value is defined already 
          @param[out] _new_n_group:    number of new group
          @param[in,out] _break_group_index:   break group index in groups
          @param[out] _n_break_no_add: number of break groups without adding new particles
          @param[in] _n_break:   number of break groups
          @param[in] _start_flag: indicate this is the first adjust of the groups in the integration
        */
        void checkNewGroup(int* _new_group_particle_index_origin,
                           int* _new_n_group_offset, 
                           int& _new_n_group, 
                           int* _break_group_index,
                           int& _n_break_no_add,
                           const int _n_break, 
                           const bool _start_flag) {
            const int n_particle = particles.getSize();
            const int n_group = groups.getSize();
            ASSERT(index_offset_group_==n_particle);

            int new_n_particle = 0;
            // used_mask store the current group index of particle i if it is already in new group
            int used_mask[n_particle+n_group];
            // -1 means not yet used 
            for (int k=0; k<n_particle; k++) used_mask[k] = -1;
            for (int k=index_offset_group_; k<n_group+index_offset_group_; k++) used_mask[k] = -1;
            // -2 means break groups
            for (int k=0; k<_n_break; k++) used_mask[_break_group_index[k]+index_offset_group_] = -2;
            
            const Float r_crit_sq = manager->r_break_crit*manager->r_break_crit;

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

                // distance criterion
                if (dr2 < r_crit_sq) {

                    // avoid break member
                    if(used_mask[j]==-2) continue;
                    
                    // avoid double count
                    if(used_mask[i]>=0 && used_mask[j]>=0) continue;

                    auto& pi = particles[i];
                    ParticleH4<Tparticle>* pj;
                    // neighbor is single 
                    if (j<index_offset_group_) {
                        pj = &particles[j];
                    }
                    else {
                        const int jg = j-index_offset_group_;
                        Float kappa_org_j = groups[jg].slowdown.getSlowDownFactorOrigin();
                        if (kappa_org_j>1.0) continue;
                        pj = &groups[jg].particles.cm;
                    }
                    // only inwards or first step case
                    Float drdv = calcDrDv(pi, *pj);
                    if(drdv<0.0||_start_flag) {
                            
                        Float fratiosq = calcPertInnerRatioSq(dr2, pi, *pj);

                        // avoid strong perturbed case, estimate perturbation
                        // if fratiosq >1.5, avoid to form new group, should be consistent as checkbreak
                        if(fratiosq>1.5) continue;

#ifdef ADJUST_GROUP_DEBUG
                        if (j<index_offset_group_) {
                            std::cerr<<"Find new group   index: "<<i<<" "<<j<<"  ftid_sq: "<<fratiosq<<"\n";
                        }
                        else {
                            auto& bin_root = groups[j-index_offset_group_].info.binarytree.getLastMember();
                            std::cerr<<"Find new group      index      slowdown      apo      ftid_sq \n"
                                     <<"i1              "
                                     <<std::setw(8)<<i
                                     <<std::setw(16)<<0
                                     <<std::setw(16)<<0
                                     <<std::setw(16)<<fratiosq;
                            std::cerr<<"\ni2              "
                                     <<std::setw(8)<<j
                                     <<std::setw(16)<<groups[j-index_offset_group_].slowdown.getSlowDownFactorOrigin()
                                     <<std::setw(16)<<bin_root.semi*(1.0+bin_root.ecc)
                                     <<std::setw(16)<<fratiosq;
                            std::cerr<<std::endl;
                        }
#endif
                        insertParticleIndexToGroup(i, j, used_mask, _new_group_particle_index_origin, _new_n_group_offset, new_n_particle, _new_n_group, _break_group_index, _n_break_no_add,  _n_break);
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

                // distance criterion
                if (dr2 < r_crit_sq) {

                    // avoid break member
                    if(used_mask[j]==-2) continue;
                    
                    // avoid double count
                    if(used_mask[i]>=0 && used_mask[j]>=0) continue;

                    auto& pi = groupi.particles.cm;

                    // avoid kappa>1.0
                    Float kappa_org_i = groupi.slowdown.getSlowDownFactorOrigin();
                    if (kappa_org_i>1.0) continue;

                    ParticleH4<Tparticle>* pj;
                    // neighbor is single 
                    if (j<index_offset_group_) {
                        //if (kappa_org_i>1.0) continue;
                        pj = &particles[j];
                    }
                    else {
                        const int jg = j-index_offset_group_;
                        Float kappa_org_j = groups[jg].slowdown.getSlowDownFactorOrigin();
                        if (kappa_org_j>1.0) continue;

                        //if (kappa_org_i>1.0&&kappa_org_j>1.0) continue;
                        pj = &groups[jg].particles.cm;

                        // unknown, for test
                        //if (kappa_org_i*kappa_org_j>1.0) continue;
                    }
                    // only inwards or first step case
                    Float drdv = calcDrDv(pi, *pj);
                    if(drdv<0.0||_start_flag) {

                        Float fratiosq = calcPertInnerRatioSq(dr2, pi, *pj);

                        // avoid strong (outside) perturbed case, estimate perturbation
                        // if fratiosq >1.5, avoid to form new group, should be consistent as checkbreak
                        if(fratiosq>1.5) continue;

#ifdef ADJUST_GROUP_DEBUG
                        auto& bini = groupi.info.binarytree.getLastMember();
                        std::cerr<<"Find new group      index      slowdown       apo      ftid_sq \n"
                                 <<"i1              "
                                 <<std::setw(8)<<i
                                 <<std::setw(16)<<kappa_org_i
                                 <<std::setw(16)<<bini.semi*(1.0+bini.ecc)
                                 <<std::setw(16)<<fratiosq;
                        if(j<index_offset_group_) {
                            std::cerr<<"\ni2              "
                                     <<std::setw(8)<<j
                                     <<std::setw(16)<<0
                                     <<std::setw(16)<<0
                                     <<std::setw(16)<<fratiosq;
                        }
                        else {
                            auto& binj = groups[j-index_offset_group_].info.binarytree.getLastMember();
                            Float kappaj = groups[j-index_offset_group_].slowdown.getSlowDownFactorOrigin();
                            std::cerr<<"\ni2              "
                                     <<std::setw(8)<<j
                                     <<std::setw(16)<<kappaj
                                     <<std::setw(16)<<binj.semi*(1.0+binj.ecc)
                                     <<std::setw(16)<<fratiosq;
                        }
                        std::cerr<<std::endl;
#endif
                        insertParticleIndexToGroup(i, j, used_mask, _new_group_particle_index_origin, _new_n_group_offset, new_n_particle, _new_n_group, _break_group_index, _n_break_no_add,  _n_break);
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
            ASSERT(!particles.isModified());
            ASSERT(initial_system_flag_);
            modify_system_flag_=true;

            // check
            // break index
            int break_group_index[groups.getSize()+1];
            int n_break = 0;
            int n_break_no_add = 0;
            // new index
            int new_group_particle_index[particles.getSize()];
            int new_n_group_offset[groups.getSize()+1];
            int new_n_group = 0;

            checkBreak(break_group_index, n_break);
            checkNewGroup(new_group_particle_index, new_n_group_offset, new_n_group, break_group_index, n_break_no_add, n_break, _start_flag);

            // integrate modified single/groups to current time
            integrateToTimeList(time_, new_group_particle_index, new_n_group_offset[new_n_group]);
            integrateToTimeList(time_, break_group_index, n_break);

            // break groups
            breakGroups(new_group_particle_index, new_n_group_offset, new_n_group, break_group_index, n_break_no_add, n_break);
            addGroups(new_group_particle_index, new_n_group_offset, new_n_group);

            // initial integration
            initialIntegration();
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
                ASSERT(particle_index_count[i]==1);
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

            dt_limit_ = manager->step.calcNextDtLimit(time_);

            // check the resolved cases
            checkGroupResolve(n_init_group_);
            writeBackResolvedGroupAndCreateJParticleList();

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

                // get ds estimation
                group_ptr[k].info.calcDsAndStepOption(group_ptr[k].slowdown.getSlowDownFactorOrigin(), ar_manager->step.getOrder());
            }

            calcDt2ndList(index_single, n_init_single_, index_group, n_init_group_, dt_limit_);
            
            updateTimeNextList(index_single, n_init_single_, index_group, n_init_group_);

            // reset n_init
            n_init_single_ = n_init_group_ = 0;
        }


        //! Integration active particles
        /*! Integrated to time_sys
        */
        void integrateOneStepAct() {
            ASSERT(!particles.isModified());
            ASSERT(initial_system_flag_);
            ASSERT(!modify_system_flag_);
            ASSERT(n_init_group_==0&&n_init_single_==0);
            
            // get next time
            Float time_next = getNextTime();

            dt_limit_ = manager->step.calcNextDtLimit(time_next);

            // integrate groups first
            const int n_group_tot = index_dt_sorted_group_.getSize();
            for (int i=0; i<n_group_tot; i++) {
                const int k = index_dt_sorted_group_[i];
                ASSERT(table_group_mask_[k]==false);
                const Float ds = groups[k].info.ds;
                groups[k].integrateToTime(ds, time_next, groups[k].info.fix_step_option);
            }

            // prediction positions
            predictAll(time_next);

            // check resolve status
            checkGroupResolve(n_group_tot);
            writeBackResolvedGroupAndCreateJParticleList();

            int* index_single = index_dt_sorted_single_.getDataAddress();
            int* index_group = index_dt_sorted_group_.getDataAddress();

            calcAccJerkNBList(index_single, n_act_single_, index_group, n_act_group_);

            correctAndCalcDt4thList(index_single, n_act_single_, index_group, n_act_group_, dt_limit_);

            updateTimeNextList(index_single, n_act_single_, index_group, n_act_group_);

            // update time
            time_ = time_next;
        }

        //! Integration a list of particle to current time (ingore dt)
        /*! Integrate a list of particle to current time.
          @param[in] _time_next: time to integrate
          @param[in] _particle_index: particle index to integrate
          @param[in] _n_particle: number of particles
        */
        void integrateToTimeList(const Float _time_next,
                                 const int* _particle_index,
                                 const int _n_particle) {
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
                ASSERT(k+index_offset_group_<n_singles + n_groups);
                ASSERT(k1+index_offset_group_<n_singles + n_groups);
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

        //! get next time for intergration
        Float getNextTime() const {
            return time_next_min_;
        }

        //! get current time
        Float getTime() const {
            return time_;
        }

        //! get active number of particles of singles
        int getNActSingle() const{
            return n_act_single_;
        }

        //! get active number of particles of groups
        int getNActGroup() const{
            return n_act_group_;
        }

        //! get sorted dt index of singles
        int* getSortDtIndexSingle() {
            return index_dt_sorted_single_.getDataAddress();
        }

        //! get sorted dt index of groups
        int* getSortDtIndexGroup() {
            return index_dt_sorted_group_.getDataAddress();
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
            std::cerr<<"Step hist:\n";
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
