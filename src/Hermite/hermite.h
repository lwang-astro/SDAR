#pragma once

#include "Common/Float.h"
#include "Hermite/particle.h"
#include "Hermite/block_time_step.h"
#include "Hermite/neighbor.h"

namespace H4{
    //! Hermite manager class
    template <class Tmethod>
    class HermiteManager{
    public:
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
            interaciton.writeBinary(_fp);
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
    template <class Tparticle, class Tpcm, class Tpert, class Tmethod, class Tgint, class Tinfo>
    class HermiteIntegrator{
    private:
        typedef ParticleH4<Tparticle> PtclH4;
        // for get sorted index of single
        class SortIndexDtSingle{
        private:
            Float * time_;
        public:
            SortIndexDt(Float * _time): time_(_time) {}
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
            SortIndexDt(Float * _time, int _offset): time_(_time), offset_(_offset) {}
            bool operator() (const int & left, const int & right) const {
                return time_[left+_offset] < time_[right+_offset];
            }
        };

        Float time_;   ///< integrated time (not real physical time if slowdown is on)
        Float time_next_min_; ///< the next minimum time 
        int n_act_single_;    /// active particle number of singles
        int n_act_group_;    /// active particle number of groups
        int n_init_single_;  /// number of singles need to initial
        int n_init_group_;   /// number of groups need to initial
        int index_offset_group_; /// offset of group index in pred_, force_ and time_next_
        bool initial_system_flag_; /// flag to indicate whether the system is initialized

        COMM::List<int> index_dt_sorted_single_; // index of particles with time next sorted order (small to large)
        COMM::List<int> index_dt_sorted_group_; // index list of active groups
        COMM::List<int> index_group_resolve_;   // index of resolved group for integration
        COMM::List<int> index_group_cm_;        // index of group with cm for integration
        COMM::List<Tparticle> pred_; // predictor
        COMM::List<ForceH4> force_;  // force
        COMM::List<Float> time_next_; // next integrated time of particles
        COMM::List<bool> index_single_mask_; // bool mask to indicate whether the particle of index (single) is masked (true) or used (false)
        COMM::List<bool> index_group_mask_; // bool mask to indicate whether the particle of index (group) is masked (true) or used (false)
    public:
        HermiteManager<Tmethod>* manager; ///< integration manager
        COMM::ParticleGroup<ParticleH4<Tparticle>, Tpcm> particles; // particles
        COMM::List<Tgint> groups; // integrator for sub-groups
        COMM::List<Neighbor> neighbors; // neighbor information of particles
        Tinfo info; ///< information of the system

    private:
        //! Calculate 2nd order time step for lists particles 
        /*! Calculate 2nd order time step 
          @param[in] _index_single: active particle index for singles
          @param[in] _n_single: number of active singles
          @param[in] _index_group: active particle index for groups
          @param[in] _n_group: number of active groups
          \return fail_flag
        */
        template <class Tptcl>
        bool calcDt2ndList(const int* _index_single,
                           const int _n_single,
                           const int* _index_group,
                           const int _n_group) {
            assert(manager!=NULL);
            // for single
            for (int i=0; i<_n_single; i++){
                const int k = _index_single(i);
                const Float* acc0 = force_[k].acc0;
                const Float* acc1 = force_[k].acc1;
                const dt =manager->step.calcDt2nd(acc0, acc1);
                particles[k].dt = dt;
                if(dt<0.0) {
                    std::cerr<<"Particle: index="<<k<<" id="<<_ptcl[k].id<<std::endl;
                    return true;
                }
            }
            // for group
            for (int i=0; i<_n_group; i++) {
                const int k = _index_group[i];
                const int kf = k + index_offset_group_;
                const Float* acc0 = force_[kf].acc0;
                const Float* acc1 = force_[kf].acc1;
                const dt =manager->step.calcDt2nd(acc0, acc1);
                groups[k].particles.cm.dt = dt;
                if(dt<0.0) {
                    std::cerr<<"Group: index="<<k<<" id="<<groups[k].particles.cm.id<<std::endl;
                    return true;
                }
            }
        }



        //! predict particles to the time
        /*! @param[in] _time_pred: time for prediction
         */
        void predictAll(const Float _time_pred) {
            static thread_local const Float inv3 = 1.0 / 3.0;
            // single
            const int n_single = index_dt_sorted_single_.getSize();
            assert(n_single <= particles.getSize());
            const auto* ptcl = particles.getDataAddress();
            auto* pred = pred_.getDataAddress();
            for (int k=0; k<n_single; k++){
                const int i = index_dt_sorted_single_[k];
                const Float dt = _time_pred - ptcl[i].time;
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
            assert(n_group <= groups.getSize());
            auto* group_ptr = groups.getDataAddress();
            for (int k=0; k<n_group; k++) {
                const int i = index_dt_sorted_group_[k];
                const auto& pcm = group_ptr[i].particles.cm;
                const Float dt = _time_pred - cm.dt;
                // group predictor is after single predictor with offset n_single
                auto& predcm = pred[i+n_single];
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
         */
        bool correctAndCalcDt4thOne(ParticleH4<Tparticle>& pi, ForceH4& fi) {
            const Float dt = pi.dt;
            const Float h = 0.5 * dt;
            const Float hinv = 2.0 / dt;
            const Float A0p[3] = {fi.acc0[0] + pi.acc0[0], fi.acc0[1] + pi.acc0[1], fi.acc0[2] + pi.acc0[2]};
            const Float A0m[3] = {fi.acc0[0] - pi.acc0[0], fi.acc0[1] - pi.acc0[1], fi.acc0[2] - pi.acc0[2]};
            const Float A1p[3] = {(fi.acc1[0] + pi.acc1[0])*h, (fi.acc1[1] + pi.acc1[1])*h, (fi.acc1[2] + pi.acc1[2])*h};
            const Float A1m[3] = {(fi.acc1[0] - pi.acc1[0])*h, (fi.acc1[1] - pi.acc1[1])*h, (fi.acc1[2] - pi.acc1[2])*h};

            const Float vel_new[3] = {pi.vel + h*( A0p[0] - inv3*A1m[0] ), 
                                      pi.vel + h*( A0p[1] - inv3*A1m[1] ), 
                                      pi.vel + h*( A0p[2] - inv3*A1m[2] )};
            pi.pos[0] += h*( (pi.vel[0] + vel_new[0]) + h*(-inv3*A0m[0]));
            pi.pos[1] += h*( (pi.vel[1] + vel_new[1]) + h*(-inv3*A0m[1]));
            pi.pos[2] += h*( (pi.vel[2] + vel_new[2]) + h*(-inv3*A0m[2]));

            pi.vel[0] = vel_new[0];
            pi.vel[1] = vel_new[1];
            pi.vel[2] = vel_new[2];

            pi.acc0[0] = fi.acc0[0];
            pi.acc0[1] = fi.acc0[1];
            pi.acc0[2] = fi.acc0[2];

            pi.acc1[0] = fi.acc1[0];
            pi.acc1[1] = fi.acc1[1];
            pi.acc1[2] = fi.acc1[2];

            pi.time += dt;

            const Float acc3[3] = {(1.5*hinv*hinv*hinv) * (A1p[0] - A0m[0]),
                                   (1.5*hinv*hinv*hinv) * (A1p[1] - A0m[1]),
                                   (1.5*hinv*hinv*hinv) * (A1p[2] - A0m[2])};
            const Float acc2[3] = {(0.5*hinv*hinv) * A1m[0] + h*acc3[0], 
                                   (0.5*hinv*hinv) * A1m[1] + h*acc3[1], 
                                   (0.5*hinv*hinv) * A1m[2] + h*acc3[2]};

#ifdef HERMITE_DEBUG
            // for debug
            assert(!std::isnan(pi.pos[0]));
            assert(!std::isnan(pi.pos[1]));
            assert(!std::isnan(pi.pos[2]));
            assert(!std::isnan(pi.vel[0]));
            assert(!std::isnan(pi.vel[1]));
            assert(!std::isnan(pi.vel[2]));
#ifdef HERMITE_DEBUG_ACC
            pi.acc2[0] = acc2[0];
            pi.acc2[1] = acc2[1];
            pi.acc2[2] = acc2[2];

            pi.acc3[0] = acc3[0];
            pi.acc3[1] = acc3[1];
            pi.acc3[2] = acc3[2];
#endif
#endif
            const Float dt_old = pi.dt;
            pi.dt = manager->step.calcDt4th(pi.acc0, pi.acc1, acc2, acc3);

#ifdef HERMITE_DEBUG
#ifdef HERMITE_DEBUG_DUMP
            if (dt_old! <= 0.0|| pi.dt <= 0.0) {
                pi.print(std::cerr);
                return true;
            }
#else
            assert(dt_old != 0.0);
            assert(pi.dt != 0.0);
#endif
#endif
        }

        //! correct particle and calculate step 
        /*! Correct particle and calculate next time step
          @param[in] _index_single: active particle index for singles
          @param[in] _n_single: number of active singles
          @param[in] _index_group: active particle index for groups
          @param[in] _n_group: number of active groups
          \return fail_flag: if the time step < dt_min, return true (failure)
        */
        bool correctAndCalcDt4thList(const int* _index_single,
                                     const int _n_single,
                                     const int* _index_group,
                                     const int _n_group) {
            static thread_local const Float inv3 = 1.0 / 3.0;

            assert(_n_single<=index_dt_sorted_single_.getSize());
            assert(_n_group<=index_dt_sorted_group_.getSize());

            // single
            auto* ptcl = particles.getDataAddress();
            ForceH4* force = force_.getDataAddress();
            for(int i=0; i<_n_single; i++){
                const int k = _index_single[i];
                auto&    pi = ptcl[k];
                ForceH4& fi = force[k];
                bool fail_flag = correctAndCalcDt4thOne(pi, fi);
                if(fail_flag) return true;
            }
            // group
            auto* group_ptr = groups.getDataAddress();
            for(int i=0; i<_n_group; i++) {
                const int k = _index_group[i];
                auto&  pi = group_ptr[k].particles.cm;
                ForceH4& fi = force[k+index_offset_group_];
                bool fail_flag = correctAndCalcDt4thOne(pi, fi);
                if(fail_flag) return true;
            }

            return false;
        }

        //! check whether group need to resolve
        /*! @param[in] _n_group: number of groups in index_dt_sorted_group_ to check
         */
        void checkGroupResolve(const int _n_group) {
            assert(_n_group<=index_dt_sorted_group_.getSize());
            for (int i=0; i<_n_group; i++) {
                const int k = index_dt_sorted_group_[i];
                auto& groupk = groups[k];
                const Float kappa = groupk.slowdown.getSlowDownFactor();
                auto& pert = groupk.perturber;
                if (kappa==1.0 && !pert.need_resolve_flag) {
                    pert.need_resolve_flag = true;
                    pert.neighbors.initial_step_flag = true;
                }
                if (kappa>3.0 && pert.need_resolve_flag) {
                    pert.need_resolve_flag = false;
                    pert.neighbors.initial_step_flag = true;
                }
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
                const int k = index_dt_sorted_group_[i];
                auto& group_k = group_ptr[k];
                if (group_k.perturber.need_resolve_flag) {
                    group_k.writeBackSlowDownParticles(ptcl[k+index_offset_group_]);
                    index_group_resolve_.addMember(k);
                }
                else index_group_cm_.addMember(k);
            }
        }

        //! calculate one particle interaction from all singles and groups
        /*! Neighbor information is also updated
          @param[out] _fi: acc and jerk of particle i
          @param[out] _nbi: neighbor information of i
          @param[in] _pi: i particle
        */
        template <class Tpi>
        inline void calcOneSingleAccJerkNB(ForceH4 &_fi, 
                                           Neighbor &_nbi,
                                           const Tpi &_pi) {
            // clear force
            _fi.clear();

            // single list
            _nbi.single_list.resizeNoInitialize(0);
            const int* single_list = index_dt_sorted_single_.getDataAddress();
            const int n_single = index_dt_sorted_single_.getSize();
            auto* ptcl = pred_.getDataAddress();
            for (int i=0; i<n_single; i++) {
                const int j = single_list[i];
                assert(j<pred_.getSize());
                const auto& pj = ptcl[j];
                if (_pi.id==pj.id) continue;
                Float r2 = manager->interaction.calcAccJerkPair(_fi, _pi, pj);
                assert(r2>0.0);
                if (r2<_nbi.r_crit_sq) _nbi.single_list.addMember(j);
                // mass weighted nearest neigbor
                if (r2*_nbi.r_min_mass < _nbi.r_min_sq*pj.mass) {
                    _nbi.r_min_sq = r2;
                    _nbi.r_min_index = j;
                    _nbi.r_min_mass = _pj.mass;
                }
                // minimum mass
                if(_nbi.min_mass<_pj.mass) {
                    _nbi.min_mass = _pj.mass;
                    _nbi.min_mass_index = j;
                }
                // if neighbor step need update, also set update flag for i particle
                if(neighbors[j].initial_step_flag) _nbi.initial_step_flag = true;
            }

            _nbi.group_list.resizeNoInitialize(0);
            auto* group_ptr = groups.getDataAddress();

            // resolved group list
            cont int n_group_resolve = index_group_resolve_.getSize();
            for (int i=0; i<n_group_resolve; i++) {
                const int j =index_group_resolve_[i];
                auto& groupj = group_ptr[j];
                assert(_pi.id!=groupj.particles.cm.id);
                const int n_member = groupj.particles.getSize();
                auto* member_adr = groupj.particles.getMemberOriginAddress();
                bool nb_flag = false;
                Float r2_min = NUMERIC_FLOAT_MAX;
                for (int k=0; k<n_member; k++) {
                    assert(_pi.id!=member_adr[k]->id);
                    const auto& pj = *member_adr[k];
                    Float r2 = manager->interaction.calcAccJerkPair(_fi, _pi, pj);
                    assert(r2>0.0);
                    if (r2<_nbi.r_crit_sq) nb_flag = true;
                    r2_min = std::min(r2_min, r2);
                }
                if (nb_flag) _nbi.group_list.addMember(j);
                // mass weighted nearest neigbor
                const Float cm_mass = groupj.particles.cm.mass;
                if (r2_min*_nbi.r_min_mass < _nbi.r_min_sq*cm_mass) {
                    _nbi.r_min_sq = r2_min;
                    _nbi.r_min_index = j;
                    _nbi.r_min_mass = cm_mass;
                }
                // minimum mass
                if(_nbi.min_mass < cm_mass) {
                    _nbi.min_mass = cm_mass;
                    _nbi.min_mass_index = j;
                }
                // if neighbor step need update, also set update flag for i particle
                if(groupj.perturber.neighbors.initial_step_flag) _nbi.initial_step_flag = true;
            }

            // cm group list
            cont int n_group_cm = index_group_cm.getSize();
            for (int i=0; i<n_group_cm; i++) {
                const int j = index_group_cm_[i];
                auto& groupj = group_ptr[j];
                // used predicted particle instead of original cm
                const auto& pj = ptcl[j+index_offset_group_];
                assert(j+index_offset_group_<pred_.getSize());
                assert(_pi.id!=pj.id);
                Float r2 = manager->interaction.calcAccJerkPair(_fi, _pi, pj);
                assert(r2>0.0);
                if (r2<_nbi.r_crit_sq) _nbi.group_list.addMember(j);
                // mass weighted nearest neigbor
                if (r2*_nbi.r_min_mass < _nbi.r_min_sq*pj.mass) {
                    _nbi.r_min_sq = r2;
                    _nbi.r_min_index = j;
                    _nbi.r_min_mass = _pj.mass;
                }
                // minimum mass
                if(_nbi.min_mass<_pj.mass) {
                    _nbi.min_mass = _pj.mass;
                    _nbi.min_mass_index = j;
                }
                // if neighbor step need update, also set update flag for i particle
                if(groupj.perturber.neighbors.initial_step_flag) _nbi.initial_step_flag = true;
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
            auto* ptcl = pred_.getDataAddress();
            auto* force_ptr = force_.getDataAddress();
            auto* neighbor_ptr = neighbors.getDataAddress();
            for (int k=0; k<_n_single; k++) {
                const int i = _index_single[k];
                auto& pi = ptcl[i];
                auto& fi = force_ptr[i];
                auto& nbi = neigbor_ptr[i];
                calcOneSingleAccJerkNB(fi, nbi, pi);
            }

            // for group active particles
            auto* group_ptr = groups.getDataAddress();
            
            for (int k=0; k<_n_group; k++) {
                const int i = _index_group[k];
                auto& groupi = group_ptr[i];
                // use predictor of cm
                auto& pi = ptcl[i+index_offset_group_];
                auto& nbi = groupi.perturber.neighbors;
                calcOneSingleAccJerkNB(fi, nbi, pi);

                // for resolved case
                if(groupi.pertuber.need_resolve_flag) {
                    auot* member_adr = groupi.particles.getMemberOriginAddress();
                    const int n_member = groupi.particles.getSize();
                    for (int j=0; j<n_member; j++) {
                        // get particle index from memory address offset
                        const int kj = int(member_adr[j]-ptcl);
                        auto& pkj = *member_adr[j];
                        auto& fkj = force_ptr[kj];
                        auto& nbkj = neighbor_ptr[kj];
                        calcOneSingleAccJerkNB(fkj, nbkj, pkj);
                    }
                    // replace the cm. force with the summation of members
                    pi.calcAverageForce(member_adr, n_member);
                }
            }
        }
        
    public:
        //! constructor
        HermiteIntegrator(): time_(Float(0.0)), time_next_min_(Float(0.0)), 
                             n_act_single_(0), n_act_group_(0), 
                             n_init_single_(0), n_init_group_(0), 
                             index_offset_group_(0), 
                             initial_system_flag_(false), 
                             index_dt_sorted_single_(), index_dt_sorted_group_(), 
                             index_group_resolve_(), index_group_cm_(), 
                             pred_(), force_(), time_next_(), 
                             index_single_mask_(), index_group_mask_(), 
                             manager(NULL), particles(), groups(), neighbors(), info() {}
        //! reserve memory 
        /*! The size of memory depends on the particle data size. Thus particles should be added or reserve first before call this function
          @param[in] _nmax_group: maximum number of groups
        */
        void reserveMem(const int _nmax_group) {
            assert(_nmax_group>0);
            const int nmax = particles.getSizeMax();
            assert(nmax>0);

            index_dt_sorted_single_.setMode(ListMode::local);
            index_dt_sorted_group_.setMode(ListMode::local);
            index_group_resolve_.setMode(ListMode::local);
            index_group_cm_.setMode(ListMode::local);
            index_single_mask_.setMode(ListMode::local);
            index_group_mask_.setMode(ListMode::local);
            pred_.setMode(ListMode::local);
            force_.setMode(ListMode::local);
            time_next_.setMode(ListMode::local);
            neighbors.setMode(ListMode::local);
            groups.setMode(ListMode::local);
            
            index_dt_sorted_single_.reserveMem(nmax);
            index_dt_sorted_group_.reserveMem(_nmax_group);

            index_group_resolve_.reserveMem(_nmax_group);
            index_group_cm_.reserveMem(_nmax_group);

            index_single_mask_.reserveMem(nmax);
            index_group_mask_.reserveMem(_nmax_group);

            // total allocation size includes both single and group c.m. size
            index_offset_group_ = nmax;
            const int nmax_tot = nmax + _nmax_group;
            pred_.reserveMem(nmax_tot);
            force_.reserveMem(nmax_tot);
            time_next_.reserveMem(nmax_tot);

            groups.reserveMem(_nmax_group);

            // reserve for neighbor list
            neighbors.reserveMem(nmax);
            for (int i=0; i<nmax_tot; i++) neighbors.reserveMem(nmax, _nmax_group);
            for (int i=0; i<_nmax_group; i++) groups[i].perturber.neigbors.reserveMem(nmax, _nmax_group);
        }

        //! clear function
        void clear() {
            particles.clear();
            index_dt_sorted_single_.clear();
            index_dt_sorted_group_.clear();
            index_group_resolve_.clear();
            index_group_cm_.clear();
            index_single_mask_.clear();
            index_group_mask_.clear();
            pred_.clear();
            force_.clear();
            time_next_.clear();
            neighbors.clear();
            groups.clear();
            initial_system_flag = false;
        }

        //! Initial system after adding particles
        /*! set number of members for different arrayies. Add group should be done after this funciton
         */
        void initialSystem() {
            particles.setModifiedFalse();
            initial_system_flag = true;

            const int n_particle = particles.getSize();
            assert(n_particle>1);
            
            pred_.resizeNoInitialize(n_particle);
            force_.resizeNoInitialize(n_particle);
            time_next_.resizeNoInitialize(n_particle);
            neighbors.resizeNoInitialize(n_particle);
            index_dt_sorted_single_.resizeNoInitialize(n_particle);
            index_single_mask_.resizeNoInitialize(n_particle);
            // initial index_dt_sorted_single_ first to allow initialIntegration work
            for (int i=0; i<n_particle; i++) {
                index_dt_sorted_single_[i] = i;
                index_single_mask_[i] = false;
            }
        }

        //! Initial Hermite integrator
        /*! Initial f, fdot and step for new add/remove particles, if start_flag, the whole system are initialized
          @param[in] _time_sys: current set time
          @param[in] _start_flag: indicate whether it is the first step 
        */
        void initialIntegration(const Float _time_sys,
                                const bool _start_flag) {

            assert(!particles.isModified());
            assert(initial_system_flag);

            // start case, set n_init_single_ and n_init_group_ consistent with index_dt_sorted array
            if(_start_flag) {
                n_init_single_ = index_dt_sorted_single_.getSize();
                n_init_group_  = index_dt_sorted_group_.getSize();
                // if initial step, n_act are not given initially
                n_act_single_  = n_init_single_;
                n_act_group_   = n_init_group_;
            }

            // single
            int* index_single = index_dt_sorted_single_.getDataAddress();

            auto* ptcl = particles.getDataAddress();
            for(int i=0; i<n_init_single_; i++){
                int k = index_single[i];
                ptcl[k].acc0[0] = ptcl[k].acc0[1] = ptcl[k].acc0[2] = 0.0;
                ptcl[k].acc1[0] = ptcl[k].acc1[1] = ptcl[k].acc1[2] = 0.0;
                ptcl[k].time = _time_sys;
                ptcl[k].dt   = 0.0;
                pred_[k] = ptcl[k];
            }

            //group
            int* index_group = index_dt_sorted_group_.getDataAddress();

            auto* group_ptr = groups.getDataAddress();
            for(int i=0; i<n_init_group_; i++){
                int k = index_group[i];
                // initial group integrator
                group_ptr[k].initial(_time_sys);
                // initial cm
                auto& pcm = group_ptr[k].particles.cm;
                pcm.acc0[0] = pcm.acc0[1] = pcm.acc0[2] = 0.0;
                pcm.acc1[0] = pcm.acc1[1] = pcm.acc1[2] = 0.0;
                pcm.time = _time_sys;
                pcm.dt   = 0.0;
                assert(k+index_offset_group_<pred_.getSize());
                pred_[k+index_offset_group_] = pcm[k];
            }

            Float dt_limit = mamager->step.calcNextDtLimit(_time_sys);

            if(!_start_flag) predictAll(_time_sys);

            // check the resolved cases
            const int n_group_tot = index_dt_sorted_group_.getSize();
            int index_group_resolve[n_group_tot], index_group_cm[n_group_tot];
            int n_group_resolve, n_group_cm;

            checkGroupResolve(_n_group);
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
                aut& fcm = force_[k+index_offset_group_];
                
                pcm.acc0[0] = fcm.acc0[0];
                pcm.acc0[1] = fcm.acc0[1];
                pcm.acc0[2] = fcm.acc0[2];
                pcm.acc1[0] = fcm.acc1[0];
                pcm.acc1[1] = fcm.acc1[1];
                pcm.acc1[2] = fcm.acc1[2];
            }

            bool fail_flag = calcDt2ndAct(index_single, n_init_single_, index_group, n_init_group_);
            if (fail_flag) return true;

            updateTimeNextList(index_single, n_init_single_, index_group, n_init_group_);
            
            return false;
        }


        //! Integration active particles
        /*! Integrated to time_sys
        */
        void integrateOneStepAct() {

            // get next time
            Float time_sys = getNextTime();

            Float dt_limit = mamager->step.calcNextDtLimit(_time_sys);

            // integrate groups first
            const int n_group_tot = index_dt_sorted_group_.getSize();
            for (int i=0; i<n_group_tot; i++) {
                const int k = index_dt_sorted_group_[i];
                const Float ds = groups[k].info.getStepSize();
                groups[k].integrateToTime(ds, _time_sys, groups[k].info.getStepOption());
            }

            // prediction positions
            predictAll(time_sys);

            // check resolve status
            checkGroupResolve(_n_group_tot);

            int index_group_resolve[n_group_tot], index_group_cm[n_group_tot];
            int n_group_resolve, n_group_cm;
            writeBackResolvedGroupAndCreateJParticleList(index_group_resolve, n_group_resolve, index_group_cm, n_group_cm);

            int* index_single = index_dt_sorted_single_.getDataAddress();
            int* index_group = index_dt_sorted_group_.getDataAddress();

            calcAccJerkNBList(index_single, n_act_single_, index_group, n_act_group_);

            const int n_single_tot = index_dt_sorted_single_.getSize();
            correctAndCalcDt4thList(index_single, n_act_single_, index_group, n_act_group_);

            updatedTimeNextList(index_single, n_act_single_, index_group, n_act_group_);
        }

        //! Integration a list of particle to current time (ingore dt)
        /*! Integrate a list of particle to current time.
          @param[in] _time_sys: time to integrate
          @param[in] _index_single: active particle index for singles
          @param[in] _n_single: number of active singles
          @param[in] _index_group: active particle index for groups
          @param[in] _n_group: number of active groups
          \return fail_flag: If step size < dt_min, return true
        */
        template <class ARCint, class Tpars>
        bool integrateToTimeList(const Float _time_sys,
                                 const int* _index_single,
                                 const int _n_single,
                                 const int* _index_group,
                                 const int _n_group) {

            // adjust dt
            // single
            int index_single_select[_n_single];
            int n_single_select =0;

            auto* ptcl = particles.getDataAddress();
            for (int i=0; i<_n_single; i++){
                const int k = _index_single(i);
                if (ptcl[k].time>=time_sys) continue;
                ptcl[k].dt = time_sys - ptcl[k].time;
                index_single_select[n_single_select++] = k;
            }

            // group
            int index_group_select[_n_group];
            int n_group_select =0;
            auto* group_ptr = groups.getDataAddress();
            for (int i=0; i<_n_group; i++) {
                const int k = _index_group[i];
                auto& pcm = group_ptr[k].particles.cm;
                if (pcm[k].time>=time_sys) continue;
                pcm[k].dt = time_sys - pcm[k].time;
                index_group_select[n_group_select++] = k;
            }            

            if (n_single_select==0 && n_group_select==0) return false;

            calcAccJerkNBList(index_single_select, n_single_select, index_group_select, n_group_select);

            bool fail_flag=correctAndCalcDt4thAct(index_single_select, n_single_select, index_group_select, n_group_select);


#ifdef HERMITE_DEBUG_DUMP
            if (fail_flag) return true;
#endif

            updateTimeNextList(index_single_select, n_single_select, index_group_select, n_group_select);

            return false;
        }

        //! sort time step array and select active particles
        /*! Make sure time_next_ is updated already
         */
        void sortDtAndSelectActParticle() {
            // sort single
            std::sort(index_dt_sorted_single_, index_dt_sorted_single_+n_act_single_, SortIndexDtSingle(time_next_));
            // sort group
            std::sort(index_dt_sorted_group_, index_dt_sorted_group_+n_act_group_, SortIndexDtGroup(time_next_, index_offset_group_));

            // get minimum next time from single and group
            time_next_min_ = std::min(time_next_[index_dt_sorted_single_[0]], time_next_[index_dt_sorted_group_[0]));
            assert(time_ref>0.0);

            // find active singles
            for(n_act_single_=0; n_act_single_<index_offset_group_; n_act_single_++){
                if (time_next_min_ < time_next_[index_dt_sorted_single_[0]]) break;
            }
            // find active groups
            const int n_groups = index_dt_sorted_group_.getSize();
            for(n_act_group_=0; n_act_group_<n_groups; n_act_group_++){
                if (time_next_min_ < time_next_[index_dt_sorted_group_[0] + index_offset_group]) break;
            }
            assert(n_act_single_==0&&n_act_group_==0);
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
        const int* getSortDtIndexSingle() const{
            return index_dt_sorted_single_.getDataAddress();
        }

        //! get sorted dt index of groups
        const int* getSortDtIndexGroup() const{
            return index_dt_sorted_group_.getDataAddress();
        }

        //! print step histogram
        void printStepHist(){
            std::map<Float, int> stephist;
            for(int i=0; i<index_dt_sorted_.size(); i++) {
                int k = index_dt_sorted_[i];
                std::map<Float, int>::iterator p = stephist.find(ptcl_[k].dt);
                if (p==stephist.end()) stephist[ptcl_[k].dt]=1;
                else stephist[ptcl_[k].dt]++;
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
