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
        int n_act_single_;    /// active particle number
        int n_act_group_;    /// active particle number
        COMM::List<int> index_dt_sorted_single_; // index of particles with time next sorted order (small to large)
        COMM::List<int> index_dt_sorted_group_; // index list of active groups
        COMM::List<Tparticle> pred_; // predictor
        COMM::List<ForceH4> force_;  // force
        COMM::List<Float> time_next_; // next integrated time of particles
    public:
        HermiteManager<Tmethod>* manager; ///< integration manager
        COMM::ParticleGroup<ParticleH4<Tparticle>, Tpcm> particles; // particles
        COMM::List<Neighbor> neighbors; // neighbor information of particles
        COMM::List<Tgint> groups; // integrator for sub-groups
        Tinfo info; ///< information of the system

    private:
        //! Calculate 2nd order time step for active particles 
        /*! Calculate 2nd order time step 
          \return fail_flag: if the new dt < dt_min, return true (failure)
        */
        template <class Tptcl>
        bool calcDt2ndAct() {
            assert(manager!=NULL);
            // for single
            for (int i=0; i<n_act_single_; i++){
                const int k = index_dt_sorted_single_[i];
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
            const int index_offset_group = index_dt_sorted_single_.getSize();
            for (int i=0; i<n_act_group_; i++) {
                const int k = index_dt_sorted_group_[i];
                const int kf = k + index_offset_group;
                const Float* acc0 = force_[kf].acc0;
                const Float* acc1 = force_[kf].acc1;
                const dt =manager->step.calcDt2nd(acc0, acc1);
                groups[k].particles.cm.dt = dt;
                if(dt<0.0) {
                    std::cerr<<"Group: index="<<k<<" id="<<groups[k].particles.cm.id<<std::endl;
                    return true;
                }
            }
            return false;
        }


        //! sort time step array and select active particles
        /*! Make sure time_next_ is updated already
         */
        void sortDtAndSelectActParticle() {
            // sort single
            std::sort(index_dt_sorted_single_, index_dt_sorted_single_+n_act_single_, SortIndexDtSingle(time_next_));
            // sort group
            const int index_offset_group = index_dt_sorted_single_.getSize();
            std::sort(index_dt_sorted_group_, index_dt_sorted_group_+n_act_group_, SortIndexDtGroup(time_next_, index_offset_group));

            // get minimum next time from single and group
            const Float time_ref = std::min(time_next_[index_dt_sorted_single_[0]], time_next_[index_dt_sorted_group_[0]));
            assert(time_ref>0.0);

            // find active singles
            for(n_act_single_=0; n_act_single_<index_offset_group; n_act_single_++){
                if (time_ref < time_next_[index_dt_sorted_single_[0]]) break;
            }
            // find active groups
            const int n_groups = index_dt_sorted_group_.getSize();
            for(n_act_group_=0; n_act_group_<n_groups; n_act_group_++){
                if (time_ref < time_next_[index_dt_sorted_group_[0] + index_offset_group]) break;
            }
            assert(n_act_single_==0&&n_act_group_==0);
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
            for (int i=0; i<n_single; i++){
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
            const int n_group = groups.getSize();
            assert(n_group == index_dt_sorted_group_.getSize());
            auto* group_ptr = groups.getDataAddress();
            for (int i=0; i<n_group; i++) {
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
            const int index_offset_group = index_dt_sorted_single_.getSize();
            for(int i=0; i<_n_group; i++) {
                const int k = _index_group[i];
                auto&  pi = group_ptr[k].particles.cm;
                ForceH4& fi = force[k+index_offset_group];
                bool fail_flag = correctAndCalcDt4thOne(pi, fi);
                if(fail_flag) return true;
            }

            return false;
        }

        //! check the pair with distance below r_crit for ptcl in adr_dt_sorted_
        /*! First check nearest neighbor distance r_min
          If r_min<r_crit, check the direction, if income, accept as group
          @param[out] _new_group_member_index: new group member index in ptcl_
          @param[out] _new_group_member_adr: new group member address, echo two neighbor is one group
          @param[out] _new_group_offset: new group offset, last one show total number of members.
          @param[in] _r_crit2: group distance criterion square
          @param[in] _mask_list: a list contain the index that would not be in the new group
          @param[in] _n_mask: number of masked indices
          @param[in] _Aint: ARC integrator class
          @param[in] _first_step_flag: if it is first step, even the out going case will be included for the chaotic situation
          \return number of new groups
        */
        template <class ARCint>
        PS::S32 checkNewGroup(PS::S32 _new_group_member_index[], Tphard* _new_group_member_adr[], PS::S32 _new_group_offset[], const PS::F64 _r_crit2, const PS::S32* _mask_list, const PS::S32 _n_mask, const ARCint* _Aint, const bool _first_step_flag) {
#ifdef HARD_DEBUG
            assert(nb_info_.size()==ptcl_.size());
#endif
            PS::S32 n_group = 0;
            if (_Aint!=NULL) n_group = _Aint->getNGroups();
            PS::S32 n_new_group=0, offset=0;

            // here ptcl_.size() is used since any ptcl index can be checked!
            PS::S32 used_mask[ptcl_.size()];
            for (PS::S32 k=0; k<ptcl_.size(); k++) used_mask[k] = -1;
            // for the index not allowed to be grouped
            for (PS::S32 k=0; k<_n_mask; k++) used_mask[_mask_list[k]] = -2; 
        

            // check adr_dt_sorted_ list (avoid suppressed ptcl)
            PS::S32 n_check = adr_dt_sorted_.size();
            for (PS::S32 k=0; k<n_check; k++) {
                const PS::S32 i = adr_dt_sorted_[k];
#ifdef ADJUST_GROUP_DEBUG
                if(i<n_group) {
                    PS::F64 mass_ratio = nb_info_[i].r_min_mass>ptcl_[i].mass ? (nb_info_[i].r_min_mass/ptcl_[i].mass) : (ptcl_[i].mass/nb_info_[i].r_min_mass);
                    std::cerr<<"SD= "<<_Aint->getSlowDown(i)<<", i="<<i<<" nearest neighbor: j="<<nb_info_[i].r_min_index<<" r_min2="<<nb_info_[i].r_min2<<" mass_ratio="<<mass_ratio<<" resolve="<<nb_info_[i].resolve_flag<<" init"<<nb_info_[i].init_flag<<" r_crit2="<<_r_crit2<<std::endl;
                }
#endif

                // get original slowdown factor
                PS::F64 sdi=0.0, frsi=0.0;
                if(i<n_group) {
                    sdi = _Aint->getSlowDown(i);
                    frsi = _Aint->getFratioSq(i);
                }

                // if radius criterion satisify or slowdown factor <1.0
                //if(nb_info_[i].r_min2 < _r_crit2 || (sdi>0.0&&sdi<1.0)) {
                if(nb_info_[i].r_min2 < _r_crit2) {
                    const PS::S32 j = nb_info_[i].r_min_index;
#ifdef HARD_DEBUG
                    assert(j<ptcl_.size());
                    assert(ptcl_[i].mass>0.0);
#endif                
                    if(j<0) continue; 

                    // avoid masked member
                    if(used_mask[i]==-2||used_mask[j]==-2) continue;

                    if(!(used_mask[i]>=0 && used_mask[j]>=0)) { // avoid double count
                        bool out_flag=getDirection(i, j);
                        if(!out_flag||_first_step_flag) {
                            PS::F64 sdj=0.0, frsj=0.0;
                            if(j<n_group) {
                                sdj = _Aint->getSlowDown(j);
                                frsj = _Aint->getFratioSq(j);
                            }
                            if(sdi>1.0&&sdj>1.0) continue;
                            if(sdi>1.0&&sdj==0.0) continue;
                            if(sdi==0.0&&sdj>1.0) continue;
                            // to avoid extremely long time integration, for very large slowdown factor, no merge group
                            // current slowdown factor 100.0 and fratioSq 1e-8 is experimental value
                            //if ((sdi>100.0&&fpj<1e-2)||(sdj>100.0&&fpi<1e-2)) continue;
                            //if ((sdi>100.0&&sdi*sdi*fpj<1.0)||(sdj>100.0&&sdj*sdj*fpi<1.0)) continue;

                            // check the tidal effect of strong bound binary
                            PS::S32 i_bin_strong = -1;  // strong bound binary index (-1 means no slowdown exist)
                            PS::F64 fratio_weak_sq;     // weak binary force ratio square
                            if (sdi>1.0) {
                                i_bin_strong = i;
                                fratio_weak_sq = frsj;
                            }
                            if (sdj>1.0) {
                                i_bin_strong = j;
                                fratio_weak_sq = frsi;
                            }
                            if(i_bin_strong>=0) {
                                PS::F64 apo = _Aint->bininfo[i_bin_strong].semi*(1.0 + _Aint->bininfo[i_bin_strong].ecc);
                                // tidal effect estimated by apo (strong) / rij
                                PS::F64 ftid_strong_sq = apo*apo/nb_info_[i].r_min2;
                                if (ftid_strong_sq*fratio_weak_sq<1e-6) continue;
                            }

                            // avoid strong perturbed case, estimate perturbation
                            PS::F64vec dr = ptcl_[j].pos - ptcl_[i].pos;
                            PS::F64 dr2 = dr*dr;
                            PS::F64 invr = 1/std::sqrt(dr2);
                            PS::F64 invr3 = invr*invr*invr;
                            PS::F64vec daccin = (ptcl_[j].mass + ptcl_[i].mass)*invr3*dr;
                            //PS::F64vec fpi = ptcl_[i].acc0 - ptcl_[j].mass*invr3*dr;
                            //PS::F64vec fpj = ptcl_[j].acc0 + ptcl_[i].mass*invr3*dr;
                            PS::F64vec daccp = ptcl_[i].acc0 - ptcl_[j].acc0 - daccin;
                            PS::F64 daccin2 = daccin*daccin;
                            PS::F64 daccp2 = daccp*daccp;
                            PS::F64 fratiosq = daccp2/daccin2;
                            // if mass ratio >1.5, avoid to form new group, should be consistent as checkbreak
                            if(fratiosq>1.5) continue;

#ifdef ADJUST_GROUP_DEBUG
                            std::cout<<"Find new group      index      slowdown      fratio_sq       apo      ftid_sq \n"
                                     <<"i1              "
                                     <<std::setw(8)<<i
                                     <<std::setw(16)<<sdi
                                     <<std::setw(16)<<frsi;
                            if(i<n_group)  {
                                PS::F64 apo= _Aint->bininfo[i].semi*(1.0+_Aint->bininfo[i].ecc);
                                std::cout<<std::setw(16)<<apo
                                         <<std::setw(16)<<apo*apo/nb_info_[i].r_min2;
                            }
                            std::cout<<"\ni2              "
                                     <<std::setw(8)<<j
                                     <<std::setw(16)<<sdj
                                     <<std::setw(16)<<frsj;
                            if(j<n_group)  {
                                PS::F64 apo= _Aint->bininfo[j].semi*(1.0+_Aint->bininfo[j].ecc);
                                std::cout<<std::setw(16)<<apo
                                         <<std::setw(16)<<apo*apo/nb_info_[i].r_min2;
                            }
                            std::cout<<std::endl;
                            std::cout<<"I_binary_strong: "<<i_bin_strong<<std::endl;
#endif
                            PS::S32 insert_group=-1, insert_index=-1;
                            if(used_mask[i]>=0) {
                                insert_group = used_mask[i];
                                insert_index = j;
                            }
                            else if(used_mask[j]>=0) {
                                insert_group = used_mask[j];
                                insert_index = i;
                            }
                            // the case of merging group
                            if(insert_group>=0) {
                                // shift the first index in last group to last 
                                PS::S32 last_group_offset = _new_group_offset[n_new_group-1];
                                _new_group_member_index[offset] = _new_group_member_index[last_group_offset];
                                _new_group_member_adr  [offset] = _new_group_member_adr  [last_group_offset];
                                offset++;

                                // in the case insert_group is not the last group
                                if(insert_group<n_new_group-1) {
                                    // shift first index in each group to the end to allow to insert j in i_group
                                    for (PS::S32 k=n_new_group-1; k>insert_group; k--) {
                                        PS::S32 k_group_offset = _new_group_offset[k];
                                        PS::S32 k0_group_offset = _new_group_offset[k-1];
                                        _new_group_member_index[k_group_offset] = _new_group_member_index[k0_group_offset];
                                        _new_group_member_adr  [k_group_offset] = _new_group_member_adr  [k0_group_offset];
                                        _new_group_offset[k]++;
                                    }
                                }
                                // replace the first position of insert_group with insert ptcl
                                PS::S32 insert_group_offset=_new_group_offset[insert_group];
                                _new_group_member_index[insert_group_offset] = insert_index;
                                _new_group_member_adr[insert_group_offset] = ptcl_ptr_[insert_index];
                                used_mask[insert_index] = insert_group;
                            }
                            else {   // new group case
                                _new_group_offset[n_new_group] = offset;
                                _new_group_member_index[offset] = i;
                                _new_group_member_adr[offset] = ptcl_ptr_[i];
                                used_mask[i] = n_new_group;
                                offset++;

                                _new_group_member_index[offset] = j;
                                _new_group_member_adr[offset] = ptcl_ptr_[j];
                                used_mask[j] = n_new_group;
                                offset++;

                                n_new_group++;
                            }
                        }
                    }
                }
            }
            // for total number of members
            _new_group_offset[n_new_group] = offset;
#ifdef HARD_DEBUG
            assert(offset<=ptcl_.size());
#endif
            return n_new_group;
        }

        //! Add a list of single particles at the end of particles
        /*! Add a list of particles at the end of Hint.ptcl_, notice all new particles have index in front of index_dt_sorted_.
          @param[in] _ptcl_origin: particle array to be added
          @param[in] _list_origin: adding particle index in _ptcl_origin
          @param[in] _n_list: number of adding particles
          @param[in] _n_nb_correct: number of ptcl (count from first): need to correct the neighbor list
          @param[in] _time_sys: new ptcl time
          @param[in] _adr_dt_front_flag: add new ptcl index in front of index_dt_sorted_ (true) or end (false)
        */
        template <class Tptcl>
        void addSingleIndexList(Tptcl* _ptcl_origin,
                                const int* _list_origin,
                                const int _n_list,
                                const int _n_nb_correct,
                                const Float _time_sys,
                                const bool _adr_dt_front_flag) {

            // quite if no new
            if(_n_list==0) return;
#ifdef HERMITE_DEBUG
            assert(_n_list>0);
#endif

            // original size
            const int n_org = ptcl_.size();
            const int n_adr_dt_org = index_dt_sorted_.size();

#ifdef HARD_DEBUG
            assert(ptcl_.capacity()>=n_org+_n_list);
#endif        
            // increase all data size by _n_list
            ptcl_.increaseSize(_n_list);
            ptcl_ptr_.increaseSize(_n_list);
            pred_.increaseSize(_n_list);
            force_.increaseSize(_n_list);
            index_dt_sorted_.increaseSize(_n_list);
            time_next_.increaseSize(_n_list);
            nb_info_.increaseSize(_n_list);
            for (int i=0; i<_n_list; i++) {
                nb_info_[i+n_org].reserveNBFlag(ptcl_.capacity());
                nb_info_[i+n_org].reserveNBList(ptcl_.capacity());
            }

#ifdef HERMITE_DEBUG
            assert(ptcl_.size()<ARRAY_ALLOW_LIMIT);
            assert(ptcl_.size()==ptcl_ptr_.size());
            assert(ptcl_.size()==pred_.size());
            assert(ptcl_.size()==force_.size());
            assert(ptcl_.size()==time_next_.size());
            assert(ptcl_.size()==nb_info_.size());
#endif

            // add ptcl in order
            if (_list_origin==NULL) {
                for (int i=0; i<_n_list; i++) {
                    const int inew= n_org+i;
                    ptcl_[inew].DataCopy(_ptcl_origin[i]);
                    ptcl_[inew].time = _time_sys;
                    ptcl_[inew].acc0 = ptcl_[inew].acc1 = force_[inew].acc0 = force_[inew].acc1 = 0.0;
                    ptcl_[inew].dt = 0.0;
                    ptcl_ptr_[inew] = (Tptcl*)&_ptcl_origin[i];
                    //_single_index_origin[_n_single+i] = i;
                }
            }
            // add ptcl from list
            else {
                for (int i=0; i<_n_list; i++) {
                    const int iadr= _list_origin[i];
                    const int inew= n_org+i;
                    ptcl_[inew].DataCopy(_ptcl_origin[iadr]);
                    ptcl_[inew].time = _time_sys;
                    ptcl_[inew].acc0 = ptcl_[inew].acc1 = force_[inew].acc0 = force_[inew].acc1 = 0.0;
                    ptcl_[inew].dt = 0.0;
                    ptcl_ptr_[inew] = (Tptcl*)&_ptcl_origin[iadr];
                    //_single_index_origin[_n_single+i] = iadr;
                }
            }
            // nb_list disp

            //_n_single += _n_list;
        
            // add new ptcl in front of adr_dt_sorted
            if(_adr_dt_front_flag) {
                // shift adr_dt_sorted
                for (int i=n_adr_dt_org-1; i>=0; i--) 
                    index_dt_sorted_[i+_n_list] = index_dt_sorted_[i];
                // add new ptcl to adr_dt_sorted
                for (int i=0; i<_n_list; i++) 
                    index_dt_sorted_[i] = n_org+i;
            }
            // add at then end of adr_dt_sorted
            else {
                for (int i=n_adr_dt_org; i<n_adr_dt_org+_n_list; i++) 
                    index_dt_sorted_[i] = i;
            }
            n_act_ += _n_list;

            // add new ptcl to neighbor list of c.m.
            for (int i=0; i<_n_nb_correct; i++) { 
                // escape the suppressed case
                if(ptcl_[i].mass==0&&nb_info_[i].getNBN()==0) continue;
                for (int j=n_org; j<n_org+_n_list; j++) 
                    nb_info_[i].addNB(j);
            }
        }
                    
        //! Remove a list of particle
        /*! Remove particles from _list, keep first _n_correct in ptcl_ undeleted (but mass to zero) and update their neighbor lists
          @param[in] _list: removing particle index for ptcl_
          @param[in] _n_list: number of removing particles
          @param[in] _n_group: number of c.m. ptcl (count from first): ptcl will not be delected but only remove from index_dt_sorted_, also the offset of the index between Hint.ptcl_ _single_index_origin
          @param[in] _n_nb_correct: number of ptcl (count from first) to correct neighbor list
        */
        void removeSingleIndexList(const int* _list,
                                   const int _n_list,
                                   const int _n_group,
                                   const int _n_nb_correct) {
            // quit if no deleted ones
            if(_n_list==0) return;
#ifdef HERMITE_DEBUG
            assert(_n_list>0);
#endif

            const int n_org=ptcl_.size();
            const int n_new = n_org - _n_list; // new ptcl number
            // moving trace (initial -1)
            int trace[n_org];
            for (int i=0; i<n_org; i++) trace[i] = -1;
            int ilast = n_org-1;

            // create moving table (trace)
            // record the new position if the ptcl is moved, if deleted, set to ptcl_.size()
            int n_decrease=0;
            for (int i=0; i<_n_list; i++) {
                int k=_list[i];
                // check no delete case
                if(k<_n_group) {
                    trace[k]=n_org;
                    // set suppressed group mass to zero
                    ptcl_[k].mass = 0.0;
                    // set status to -20 to identify the suppressed ptcl
                    ptcl_[k].status = -20;
                    // set neighbor to zero to avoid issue
                    nb_list_n_[k] = 0;
                    // set r_min_index to -1
                    nb_info_[k].r_min_index = -1;
                    continue;
                }

                int idel = k;
                // set ilast >= idel for safety.
                ilast = std::max(ilast, idel);
            
                // check whether position is already moved, if so, check the moved new position and set current to delete
                if(trace[k]>=0) {
                    idel=trace[k];
                    trace[k]=n_org;
                }
#ifdef HERMITE_DEBUG
                assert(k<n_org);
                assert(idel<n_org);
#endif

                // check the last avaiable particle that can be moved to idel
                // If the ilast is already moved, check whether the new moved position trace[ilast] is before the current idel.
                // If trace[ilast] is after the current idel, update the trace[ilast] to idel and set trace[ilast] as new idel to check, until trace[ilast]=-1 or idel >= ilast
                while (idel>=0) {
                    int itrlast = -1;
                    //cond: last is moved && ( moved pos before new n_ptcl || last is del ) &&  ilast > idel 
                    while(trace[ilast]>=0 && (trace[ilast]<n_new || trace[ilast]==n_org) && ilast>idel) ilast--;

                    // if ilast is not yet moved or (new pos of ilast > idel and ilast is not del, move ilast to idel)
                    if(trace[ilast]<0||(idel<trace[ilast]&&trace[ilast]<n_org)) {
                        itrlast=trace[ilast];
                        trace[ilast]=idel;
                    }
                    // idel is already at last, remove it
                    trace[idel]=n_org; 
                    idel = itrlast;
                }

                n_decrease++;
            }

            // move data
            for (int i=0; i<n_org; i++) {
                // no change case
                if(trace[i]<0) continue;
                // remove case
                else if(trace[i]<n_org) {
                    const int inew=trace[i];
                    ptcl_[inew] = ptcl_[i];
                    ptcl_ptr_[inew] = ptcl_ptr_[i];
                    pred_[inew] = pred_[i];
                    force_[inew] = force_[i];
                    time_next_[inew] = time_next_[i];
                    nb_info_[inew] = nb_info_[i];
                    //_single_index_origin[inew-_n_group] = _single_index_origin[i-_n_group];

                    // shift nb_list
                    nb_list_n_[inew] = nb_list_n_[i];
                    nb_list_n_[i] = 0;
                    int disp_i = nb_list_disp_[i];
                    int disp_inew = nb_list_disp_[inew];
                    for (int i=0; i<nb_list_n_[inew]; i++) 
                        nb_list_[disp_inew+i] = nb_list_[disp_i+i];
                }
            }
            ptcl_.decreaseSize(n_decrease);
            ptcl_ptr_.decreaseSize(n_decrease);
            pred_.decreaseSize(n_decrease);
            force_.decreaseSize(n_decrease);
            time_next_.decreaseSize(n_decrease);
            nb_info_.decreaseSize(n_decrease);
            nb_list_n_.decreaseSize(n_decrease);
            nb_list_disp_.decreaseSize(n_decrease);
            nb_list_.decreaseSize(n_decrease*n_nb_off_);
            //_n_single -= n_decrease;

            // check nb_info
            for (int i=0; i<ptcl_.size(); i++) {
                const int i_min=nb_info_[i].r_min_index;
                if(i_min>=0) {
                    const int jtr=trace[i_min];
                    if(jtr>=0) {
                        if(jtr<n_org) nb_info_[i].r_min_index=jtr;
                        else nb_info_[i].r_min_index=-1;
                    }
                }
            }

            // update adr_dt_sorted
            int index_dt_sorted_new_size=modifyList(index_dt_sorted_.getPointer(), index_dt_sorted_.size(), trace, n_org);
            index_dt_sorted_.decreaseSize(_n_list);

            assert(index_dt_sorted_new_size==index_dt_sorted_.size());

        
            // correct neighbor list
            for (int i=0; i<_n_nb_correct; i++) {
                // escape suppressed one
                if (trace[i]==n_org) continue;
            
                int nb_list_i_new_size=modifyList(nb_list_.getPointer(nb_list_disp_[i]), nb_list_n_[i], trace, n_org);
                nb_list_n_[i] = nb_list_i_new_size;
#ifdef HERMITE_DEBUG
                assert(nb_list_n_[i]>=0);
#endif
            }
            
        }                       

    public:

        //! reserve memory 
        /*! The size of memory depends on the particle data size. Thus particles should be added or reserve first before call this function
          @param[in] _nmax_group: maximum number of groups
        */
        void reserveMem(const int _nmax_group) {
            int nmax = particles.getSizeMax();
            assert(nmax>0);

            index_dt_sorted_single_.setMode(ListMode::local);
            index_dt_sorted_group_.setMode(ListMode::local);
            pred_.setMode(ListMode::local);
            force_.setMode(ListMode::local);
            time_next_.setMode(ListMode::local);
            neighbors.setMode(ListMode::local);
            

            index_dt_sorted_single_.reserveMem(nmax);
            index_dt_sorted_group_.reserveMem(nmax);
            pred_.reserveMem(nmax);
            force_.reserveMem(nmax);
            time_next_.reserveMem(nmax);
            neighbors.reserveMem(nmax);

            groups.setMode(ListMode::local);
            groups.reserveMem(_nmax_group);
        }

        //! clear function
        void clear() {
            particles.clear();
            index_dt_sorted_single_.clear();
            index_dt_sorted_group_.clear();
            pred_.clear();
            force_.clear();
            time_next_.clear();
            neighbors.clear();
            groups.clear();
        }


        //! Write back particle data to ptcl
        /* @param[out] _ptcl: ptcl array to write
           @param[in] _ptcl_list: particle address in _ptcl
           @param[in] _n_ptcl: number of particles need for copy
           @param[in] _i_start: start index in local particle array for copy
        */
        template <class Tptcl>
        void writeBackPtcl(Tptcl * _ptcl, 
                           const int* _ptcl_list,
                           const int _n_ptcl, 
                           const int _i_start) {
#ifdef HERMITE_DEBUG
            assert(_i_start+_n_ptcl<=ptcl_.size());
#endif
            for (int i=0; i<_n_ptcl; i++) {
#ifdef HERMITE_DEBUG
                assert(_ptcl[_ptcl_list[i]].id==ptcl_[i+_i_start].id);
                assert(!std::isnan(ptcl_[i+_i_start].pos[0]));
                assert(!std::isnan(ptcl_[i+_i_start].pos[1]));
                assert(!std::isnan(ptcl_[i+_i_start].pos[2]));
                assert(!std::isnan(ptcl_[i+_i_start].vel[0]));
                assert(!std::isnan(ptcl_[i+_i_start].vel[1]));
                assert(!std::isnan(ptcl_[i+_i_start].vel[2]));
#endif
                _ptcl[_ptcl_list[i]].DataCopy(ptcl_[i+_i_start]);
            }
        }

        //! Write back particle data to original array
        /*!
          @param[in] _i_start: starting index in ptcl_ to copy
        */
        void writeBackPtcl(const int _i_start) {
            for (int i=_i_start; i<ptcl_.size(); i++) {
#ifdef HERMITE_DEBUG
                assert(!std::isnan(ptcl_[i].pos[0]));
                assert(!std::isnan(ptcl_[i].pos[1]));
                assert(!std::isnan(ptcl_[i].pos[2]));
                assert(!std::isnan(ptcl_[i].vel[0]));
                assert(!std::isnan(ptcl_[i].vel[1]));
                assert(!std::isnan(ptcl_[i].vel[2]));
//            assert(ptcl_[i].pos[0]==ptcl_[i].pos[0]);
//            assert(ptcl_[i].pos[1]==ptcl_[i].pos[1]);
//            assert(ptcl_[i].pos[2]==ptcl_[i].pos[2]);
//            assert(ptcl_[i].vel[0]==ptcl_[i].vel[0]);
//            assert(ptcl_[i].vel[1]==ptcl_[i].vel[1]);
//            assert(ptcl_[i].vel[2]==ptcl_[i].vel[2]);
                assert(ptcl_ptr_[i]->id==ptcl_[i].id);
#endif
                ptcl_ptr_[i]->DataCopy(ptcl_[i]);
            }
        }

        //! Write back a list of particle data to original array 
        /*!
          @param[in] _ptcl_list: particle index in ptcl_
          @param[in] _n_ptcl: number of particles need for copy
          @param[in] _n_avoid: if particle index is below _n_avoid, no write (to avoid write group c.m.)
        */
        void writeBackPtcl(const int* _ptcl_list,
                           const int _n_ptcl,
                           const int _n_avoid) {
            for (int i=0; i<_n_ptcl; i++) {
                const int k = _ptcl_list[i];
                if(k<_n_avoid) continue;
#ifdef HERMITE_DEBUG
                assert(k<ptcl_.size()&&k>=0);
                assert(!std::isnan(ptcl_[k].pos[0]));
                assert(!std::isnan(ptcl_[k].pos[1]));
                assert(!std::isnan(ptcl_[k].pos[2]));
                assert(!std::isnan(ptcl_[k].vel[0]));
                assert(!std::isnan(ptcl_[k].vel[1]));
                assert(!std::isnan(ptcl_[k].vel[2]));
//            assert(ptcl_[k].pos[0]==ptcl_[k].pos[0]);
//            assert(ptcl_[k].pos[1]==ptcl_[k].pos[1]);
//            assert(ptcl_[k].pos[2]==ptcl_[k].pos[2]);
//            assert(ptcl_[k].vel[0]==ptcl_[k].vel[0]);
//            assert(ptcl_[k].vel[1]==ptcl_[k].vel[1]);
//            assert(ptcl_[k].vel[2]==ptcl_[k].vel[2]);
                assert(ptcl_ptr_[k]->id==ptcl_[k].id);
#endif            
                ptcl_ptr_[k]->DataCopy(ptcl_[k]);
            }
        }

        //! Write back one particle data to original array 
        /*!
          @param[in] _i: particle index in ptcl_
        */
        void writeBackOnePtcl(const int _i) {
#ifdef HERMITE_DEBUG
            assert(_i<ptcl_.size()&&_i>=0);
            assert(!std::isnan(ptcl_[_i].pos[0]));
            assert(!std::isnan(ptcl_[_i].pos[1]));
            assert(!std::isnan(ptcl_[_i].pos[2]));
            assert(!std::isnan(ptcl_[_i].vel[0]));
            assert(!std::isnan(ptcl_[_i].vel[1]));
            assert(!std::isnan(ptcl_[_i].vel[2]));
//        assert(ptcl_[_i].pos[0]==ptcl_[_i].pos[0]);
//        assert(ptcl_[_i].pos[1]==ptcl_[_i].pos[1]);
//        assert(ptcl_[_i].pos[2]==ptcl_[_i].pos[2]);
//        assert(ptcl_[_i].vel[0]==ptcl_[_i].vel[0]);
//        assert(ptcl_[_i].vel[1]==ptcl_[_i].vel[1]);
//        assert(ptcl_[_i].vel[2]==ptcl_[_i].vel[2]);
            assert(ptcl_ptr_[_i]->id==ptcl_[_i].id);
#endif            
            ptcl_ptr_[_i]->DataCopy(ptcl_[_i]);
        }

        //! Search perturber and neighbor
        /*! Assume the suppressed ptcl has zero mass
          @param[in] _apo_bin: apo-center distance of binary 
          @param[in] _n_bin: number of binaries (locate at begining of ptcl_)
          @param[in] _n_search: number of particles to search perturber
        */
        void searchPerturberBin(Float _apo_bin[], const int _n_bin, const int _n_search) {
            int n = ptcl_.size();
        
            // find perturber
            for(int i=0; i<_n_search; i++) {
                int disp=nb_list_disp_[i];
                int n_pert=0;
                nb_info_[i].r_min2 = PS::LARGE_FLOAT;
                nb_info_[i].min_mass = PS::LARGE_FLOAT;
                for(int j=0; j<_n_bin; j++) {
                    Floatvec dr = ptcl_[i].pos-ptcl_[j].pos;
                    Float r2 = dr*dr;
                    Float r_search = std::max(ptcl_[i].r_search,ptcl_[j].r_search) + _apo_bin[j];
                    Float r_search2 = r_search*r_search;
                    if (r2<r_search2&&i!=j&&ptcl_[j].mass>0.0) {
                        nb_list_[disp+n_pert]=j;
                        n_pert++;
                        if(r2<nb_info_[i].r_min2) {
                            nb_info_[i].r_min2 = r2;
                            nb_info_[i].r_min_index = j;
                        }
                        nb_info_[i].min_mass = std::min(nb_info_[i].min_mass, ptcl_[j].mass);
                    }
                }
                for(int j=_n_bin; j<n; j++) {
                    Floatvec dr = ptcl_[i].pos-ptcl_[j].pos;
                    Float r2 = dr*dr;
                    Float r_search = std::max(ptcl_[i].r_search,ptcl_[j].r_search);
                    Float r_search2 = r_search*r_search;
                    if (r2<r_search2&&i!=j&&ptcl_[j].mass>0.0) {
                        nb_list_[disp+n_pert]=j;
                        n_pert++;
                        if(r2<nb_info_[i].r_min2) {
                            nb_info_[i].r_min2 = r2;
                            nb_info_[i].r_min_index = j;
                        }
                        nb_info_[i].min_mass = std::min(nb_info_[i].min_mass, ptcl_[j].mass);
                    }
                }
//            nb_list_disp_[i] = disp;
                nb_list_n_[i] = n_pert;
#ifdef HERMITE_DEBUG
                assert(n_pert<=n_nb_off_);
#endif
            }
        }

        //! Search perturber and find nearest particles for the first N particles in time step sorted list
        /*!
          @param[in] _n_search: number of particles to search perturber
        */
        void searchNeighborWithInfoN(const int _n_search) {
            for(int i=0; i<_n_search; i++) searchNeighborWithInfoOne(i);
        }

        //! Search neighbor and find nearest particles for one particle
        /*! Search perturber and neighbor for one particle from adr_dt_sorted list.
          In this case, the suppressed particle can be avoid
          @param[in] _i: index of particle index in adr_dt_sorted to search perturber
        */
        inline void searchNeighborWithInfoOne(const int _i) {
            int n = index_dt_sorted_.size();
            // find perturber
            int disp=nb_list_disp_[_i];
            int n_pert=0;
            nb_info_[_i].r_min2 = PS::LARGE_FLOAT;
            nb_info_[_i].min_mass = PS::LARGE_FLOAT;
            for(int k=0; k<n; k++) {
                int j = index_dt_sorted_[k];
                Floatvec dr = ptcl_[_i].pos-ptcl_[j].pos;
                Float r2 = dr*dr;
                Float r_search = std::max(ptcl_[_i].r_search,ptcl_[j].r_search);
                Float r_search2 = r_search*r_search;
                if (r2<r_search2&&_i!=j) {
#ifdef HERMITE_DEBUG
                    assert(ptcl_[j].mass>0);
#endif
                    nb_list_[disp+n_pert]=j;
                    n_pert++;
                    if(r2<nb_info_[_i].r_min2) {
                        nb_info_[_i].r_min2 = r2;
                        nb_info_[_i].r_min_index = j;
                    }
                    nb_info_[_i].min_mass = std::min(nb_info_[_i].min_mass, ptcl_[j].mass);
                }
            }
//            nb_list_disp_[i] = disp;
            nb_list_n_[_i] = n_pert;
#ifdef HERMITE_DEBUG
            assert(n_pert<=n_nb_off_);
#endif
        }

        //! update rsearch of ptcl
        /*! 
          @param[in] _i_start: start index to calculate
          @param[in] _dt_tree: tree time step
          @param[in] _v_max: maximum velocity used to calcualte r_search
          \return the maximum research
        */
        Float updateRSearch(const int _i_start, const Float _dt_tree, const Float _v_max) {
            Float dt_reduce_factor=1.0;
            for(int i=_i_start; i<ptcl_.size(); i++) {
                Float dt_reduce_fi = ptcl_[i].calcRSearch(_dt_tree, _v_max);
                dt_reduce_factor = std::max(dt_reduce_fi, dt_reduce_factor);
            }
            return dt_reduce_factor;
        }

        const int* getPertList(const int i) {
            return &nb_list_[nb_list_disp_[i]];
        }

        int getPertN(const int i) const {
            return nb_list_n_[i];
        }

        int getPertListSize() const {
            return nb_list_.size();
        }

        const PtclH4* getPtcl() const {
            return ptcl_.getPointer();
        }
    
        const PtclForce* getForce() const {
            return force_.getPointer();
        }

        int getPtclN() const {
            return ptcl_.size();
        }

        const Tgptcl** getPtclAdr() const {
            return ptcl_ptr_.getPointer();
        }

        const Tneighbor& getNbInfo(const int i) const {
            return nb_info_[i];
        }

#ifdef HERMITE_DEBUG_PRINT
        void writePtcl(FILE* _fout, const int _i_start) const{
            for (int i=_i_start; i<ptcl_.size(); i++) {
                ptcl_[i].ParticleBase::writeAscii(_fout);
            }
        }
#endif

        //! Set neighbor list offset (maximum number)
        /*!
          @param[in] _n_nb_off: neighbor list offset (maximum number)
        */
        void setNBOff(const int _n_nb_off) {
            n_nb_off_ = _n_nb_off;
        }

        //! Get Minimum mass of ptcl
        /*!
          \return minimum mass
        */
        const Float getMassMin(){
            Float mass_min = PS::LARGE_FLOAT;
            for(int i=0; i<ptcl_.size(); i++){
                if(mass_min > ptcl_[i].mass)  mass_min = ptcl_[i].mass;
            }
            return mass_min;
        }

        //! Initial Hermite ptcl force and step 
        /*! Initial f, fdot and step
          @param[in] _list: ptcl index to initialize
          @param[in] _n_list: number of ptcl to initialize
          @param[in] _time_sys: current set time
          @param[in] _int_pars: integration parameters (include changeover functions)
          @param[in,out] _Aint: ARC integrator class (resolve and shift is used)
          @param[in] _start_flag: indicate whether it is the initial step; true: set Nact to zero; false: predict particles for force calculation
        */
        template <class ARCint, class Tpars>
        bool initial(const int* _list,
                     const int _n_list,
                     const Float _time_sys,
                     const Tpars * _int_pars,
                     ARCint* _Aint,
                     const bool _start_flag = true) {

            // if no particle need initial, quit
            if(_n_list==0) return false;
#ifdef HERMITE_DEBUG
            assert(_n_list>0);
#endif

            const int* ptcl_list=_list;
            if(ptcl_list==NULL) ptcl_list = index_dt_sorted_.getPointer();

            Float mass_min = PS::LARGE_FLOAT;
            for(int i=0; i<_n_list; i++){
                int iadr = ptcl_list[i];
                pred_[iadr].r_search = ptcl_[iadr].r_search;
                ptcl_[iadr].time = _time_sys;
                ptcl_[iadr].dt   = _time_sys;
                ptcl_[iadr].acc0 = ptcl_[iadr].acc1 = 0.0;
                mass_min = std::min(mass_min, ptcl_[iadr].mass);
            }

            time_step_.setTimeAndCalcNextDtMax(_time_sys);

            if(_start_flag) {
                n_act_ = 0;
                time_step_.calcAcc0OffsetSq(mass_min, _int_pars->changeover.getRout());
            }
            else PredictAll(pred_.getPointer(), ptcl_.getPointer(), ptcl_.size(), _time_sys);

            if(_Aint!=NULL) {
                _Aint->updateCM(ptcl_.getPointer());
                _Aint->resolve();
            }

            // force::acc0,acc1, neighbor list updated
//        if(_calc_full_flag) 
            bool fail_flag=CalcActAcc0Acc1(force_.getPointer(), 
                                           nb_info_.getPointer(),
                                           nb_list_.getPointer(), nb_list_n_.getPointer(), 
                                           nb_list_disp_.getPointer(), 
                                           ptcl_.getPointer(), ptcl_.size(), 
                                           ptcl_list, _n_list, _int_pars,
                                           _Aint);
#ifdef HERMITE_DEBUG_DUMP
            if(fail_flag) return true;
#endif
//        // only neighbor force calculation
//        else CalcAcc0Acc1ActNb(force_.getPointer(), 
//                               ptcl_.getPointer(), 
//                               ptcl_list, _n_list,
//                               nb_list_.getPointer(), nb_list_disp_.getPointer(), nb_list_n_.getPointer(), 
//                               r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, _Aint, _n_group);
    
            // store predicted force
            for(int i=0; i<_n_list; i++){
                int iadr = ptcl_list[i];
                ptcl_[iadr].acc0 = force_[iadr].acc0;
                ptcl_[iadr].acc1 = force_[iadr].acc1;
            }

            if(_Aint!=NULL) _Aint->shift2CM();

            fail_flag=calcDt2ndAct(ptcl_.getPointer(), force_.getPointer(), index_dt_sorted_.getPointer(), _n_list, time_step_);

#ifdef HERMITE_DEBUG_DUMP
            if(fail_flag) return true;
#endif
            for(int i=0; i<_n_list; i++){
                int iadr = ptcl_list[i];
                time_next_[iadr] = ptcl_[iadr].time + ptcl_[iadr].dt;
            }

            return fail_flag;
        }
    
//    template<class Energy>
//    void CalcEnergy(Energy & eng) {
//        CalcEnergy(ptcl_.getPointer(), ptcl_.size(), eng, r_in_, r_out_, eps_sq_);
//    }

        Float getNextTime() const {
            return time_next_[index_dt_sorted_[0]];
        }

        Float getTime() const {
            return time_step_.getTime();
        }

        //Float* getNextTimeList() const {
        //    return time_next_.getPointer();
        //}
    
        Float getOneTime(const std::size_t i) const {
            return ptcl_[i].time;
        }

        Float getOneDt(const std::size_t i) const {
            return ptcl_[i].dt;
        }

        //! Integration active particles
        /*! Integrated to time_sys
          @param[in] _int_pars: integration parameters (include changeover functions)
          @param[in,out] _Aint: ARC integrator class (resolve and shift is used)
          \return fail_flag: If step size < dt_min, return true
        */
        template <class ARCint, class Tpars>
        bool integrateOneStepAct(const Tpars* _int_pars,
                                 ARCint* _Aint) {
            int n_group = 0;

            // get next time
            Float time_sys = getNextTime();
            time_step_.setTimeAndCalcNextDtMax(time_sys);

            // pred::mass,pos,vel updated
            PredictAll(pred_.getPointer(), ptcl_.getPointer(), ptcl_.size(), time_sys);
            if(_Aint!=NULL) {
                n_group = _Aint->getNGroups();
                _Aint->updateCM(pred_.getPointer());
                _Aint->resolve();
            }

            // force::acc0,acc1, neighbor list updated
//        if(_calc_full_flag) 
            bool fail_flag=CalcAcc0Acc1Act(force_.getPointer(), 
                                           nb_info_.getPointer(),
                                           nb_list_.getPointer(), nb_list_n_.getPointer(), 
                                           nb_list_disp_.getPointer(), 
                                           pred_.getPointer(), ptcl_.size(), 
                                           index_dt_sorted_.getPointer(), n_act_, 
                                           _int_pars,
                                           _Aint);
#ifdef HERMITE_DEBUG_DUMP
            if (fail_flag) return true;
#endif
//        // only neighbor force calculation
//        else CalcAcc0Acc1ActNb(force_.getPointer(), 
//                               pred_.getPointer(), 
//                               index_dt_sorted_.getPointer(), n_act_, 
//                               nb_list_.getPointer(), nb_list_disp_.getPointer(), nb_list_n_.getPointer(), 
//                               r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, _Aint, _n_group);

            // ptcl_org::pos,vel; pred::time,dt,acc0,acc1,acc2,acc3 updated
            fail_flag=CorrectAndCalcDt4thAct(ptcl_.getPointer(), force_.getPointer(), index_dt_sorted_.getPointer(), n_act_, time_step_);
#ifdef HERMITE_DEBUG_DUMP
            if (fail_flag) return true;
#endif

            for(int i=0; i<n_act_; i++){
                int adr = index_dt_sorted_[i];
                time_next_[adr] = ptcl_[adr].time + ptcl_[adr].dt;
                // update new perturber list
                if(adr<n_group) _Aint->updatePertOneGroup(adr, ptcl_.getPointer(), force_.getPointer(), getPertList(adr), getPertN(adr));
            }

            if(_Aint!=NULL) {
                // shift member to c.m. frame
                _Aint->shift2CM();
                // go back to original time
                _Aint->updateCM(ptcl_.getPointer());
                //Aint.updateCM(Hint.getPtcl(), group_act_list.getPointer(), group_act_n);
            }
        
            return fail_flag;
        }

        //! Integration a list of particle to current time (ingore dt)
        /*! Integrate a list of particle to current time.
          @param[in] _ptcl_list: particle index list to integrate
          @param[in] _n_list: number of particles
          @param[in] _int_pars: integrator paramters
          @param[in,out] _Aint: ARC integrator class (resolve and shift is used)
          \return fail_flag: If step size < dt_min, return true
        */
        template <class ARCint, class Tpars>
        bool integrateOneStepListNoPred(const int* _ptcl_list,
                                        const int _n_list,
                                        const Tpars* _int_pars,
                                        ARCint* _Aint) {

            if(_n_list==0) return false;
#ifdef HERMITE_DEBUG
#ifdef HERMITE_DEBUG_DUMP
            if(_n_list<0) {
                std::cerr<<"Error: number of particle in list is negative!\n";
                return true;
            }
#else
            assert(_n_list>0);
#endif
#endif
            int n_group = 0;

            // get next time
            Float time_sys = time_step_.getTime();

            if(_Aint!=NULL) n_group = _Aint->getNGroups();

            // adjust dt
            int int_list[_n_list];
            int group_int_list[n_group+1];
            int n_group_int=0;
            int n_int=0;
            for (int i=0; i<_n_list; i++) {
                int iadr = _ptcl_list[i];
                if(ptcl_[iadr].time>=time_sys) continue;
                ptcl_[iadr].dt = time_sys - ptcl_[iadr].time;
                int_list[n_int++] = iadr;
                if(iadr<n_group) group_int_list[n_group_int++] = iadr;
            }

            if (n_int==0) return false;

            // update next time step limit

            if(_Aint!=NULL) _Aint->resolve();

            // force::acc0,acc1, neighbor list updated
//        if(_calc_full_flag) 
            bool fail_flag=CalcActAcc0Acc1(force_.getPointer(), 
                                           nb_info_.getPointer(),
                                           nb_list_.getPointer(), nb_list_n_.getPointer(), 
                                           nb_list_disp_.getPointer(), 
                                           pred_.getPointer(), ptcl_.size(), 
                                           int_list, n_int,
                                           _int_pars,
                                           _Aint);
//        // only neighbor force calculation
//        else CalcAcc0Acc1ActNb(force_.getPointer(), 
//                               pred_.getPointer(), 
//                               int_list, n_int,
//                               nb_list_.getPointer(), nb_list_disp_.getPointer(), nb_list_n_.getPointer(), 
//                               r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, _Aint, _n_group);

#ifdef HERMITE_DEBUG_DUMP
            if (fail_flag) return true;
#endif

            fail_flag=CorrectAndCalcDt4thAct(ptcl_.getPointer(), 
                                             force_.getPointer(), 
                                             int_list, n_int, time_step_);


#ifdef HERMITE_DEBUG_DUMP
            if (fail_flag) return true;
#endif

            // shift member to c.m. frame
            if(_Aint!=NULL) {
                _Aint->shift2CM();
                // update c.m. of _i
                _Aint->updateCM(ptcl_.getPointer(), group_int_list, n_group_int);
            }

            for(int i=0; i<n_int; i++){
                int iadr=int_list[i];
                // update time table
                time_next_[iadr] = ptcl_[iadr].time + ptcl_[iadr].dt;
            }

            // update new perturber list
            for(int i=0; i<n_group_int; i++) {
                int iadr=group_int_list[i];
                _Aint->updatePertOneGroup(iadr, ptcl_.getPointer(), force_.getPointer(), getPertList(iadr), getPertN(iadr));
            }

            return false;
        }

        void SortAndSelectIp(int group_act_list[],
                             int &group_act_n,
                             const int n_groups) {
            SortAndSelectIp(index_dt_sorted_.getPointer(), time_next_.getPointer(), n_act_, time_next_.size(), group_act_list, group_act_n, n_groups);
        
        }

        void SortAndSelectIp() {
            SortAndSelectIp(index_dt_sorted_.getPointer(), time_next_.getPointer(), n_act_, index_dt_sorted_.size());
        }
    
        int getNact() const{
            return n_act_;
        }

        const int* getActList() const{
            return index_dt_sorted_.getPointer();
        }

#ifdef HERMITE_DEBUG
        template <class ARCint>
        bool checkAdrList(const ARCint& _Aint) {
            bool fail_flag = false;
            const int n_ptcl = ptcl_.size();
            const int n_adr = index_dt_sorted_.size();
            int ncheck[n_ptcl];
            for (int i=0; i<n_ptcl; i++) ncheck[i]=0;
            for (int i=0; i<n_adr; i++) {
                int k =index_dt_sorted_[i];
#ifdef HERMITE_DEBUG_DUMP
                if (time_next_[k] != ptcl_[k].time + ptcl_[k].dt) fail_flag = true;
#else
                assert(time_next_[k]==ptcl_[k].time + ptcl_[k].dt);
#endif
                Float tr=int(ptcl_[k].time/ptcl_[k].dt);
#ifdef HERMITE_DEBUG_DUMP
                if (tr*ptcl_[k].dt != ptcl_[k].time) fail_flag = true;
#else
                assert(tr*ptcl_[k].dt==ptcl_[k].time);
#endif
                ncheck[index_dt_sorted_[i]]++;
            }
            const int n_group = _Aint.getNGroups();
            for (int i=0; i<n_group; i++) {
#ifdef HERMITE_DEBUG_DUMP
                if(_Aint.getMask(i)) { 
                    if (ncheck[i]!=0) fail_flag = true;
                }
                else if (ncheck[i]!=1) 
                    fail_flag = true;
#else
                if(_Aint.getMask(i)) assert(ncheck[i]==0);
                else assert(ncheck[i]==1);
#endif
            }
            for (int i=n_group; i<n_ptcl; i++) {
#ifdef HERMITE_DEBUG_DUMP
                if (ncheck[i]!=1) fail_flag = true;
#else
                assert(ncheck[i]==1);
#endif
            }
            for (int i=0; i<n_adr-1; i++) {
#ifdef HERMITE_DEBUG_DUMP
                if (ptcl_[index_dt_sorted_[i]].dt>ptcl_[index_dt_sorted_[i+1]].dt) fail_flag = true;
                if (time_next_[index_dt_sorted_[i]]>time_next_[index_dt_sorted_[i+1]]) fail_flag = true;
#else
                assert(ptcl_[index_dt_sorted_[i]].dt<=ptcl_[index_dt_sorted_[i+1]].dt);
                assert(time_next_[index_dt_sorted_[i]]<=time_next_[index_dt_sorted_[i+1]]);
#endif

            }
            return fail_flag;
        }

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
#endif

    };

}
