#pragma once

#include "Common/Float.h"
#include "Common/list.h"
#include "Common/binary_tree.h"
#include "Hermite/hermite_particle.h"

namespace H4{
    //! contain group information
    template <class Tparticle>
    class ARInformation{
    private:
        typedef ParticleAR<Tparticle> ARPtcl;
        //! calculate minimum kepler ds for ARC
        static Float calcDsKeplerIter(const Float& _ds, COMM::BinaryTree<ARPtcl>& _bin) {
            Float ds;
            //kepler orbit, step ds=dt*m1*m2/r estimation (1/64 orbit): pi/4*sqrt(semi/(m1+m2))*m1*m2 
            if (_bin.semi>0) ds = 0.09817477042*std::sqrt(_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
            //hyperbolic orbit, step ds=dt*m1*m2/r estimation (1/256 orbit): pi/128*sqrt(semi/(m1+m2))*m1*m2         
            else ds = 0.0245436926*std::sqrt(-_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);   
            return std::min(ds, _ds);
        }

        //! add one particle and relink the particle pointer to the local one in _particles
        static void addOneParticleAndReLinkPointer(COMM::ParticleGroup<ARPtcl, ParticleH4<Tparticle>>& _particles, ARPtcl*& _ptcl) {
            _particles.addMember(*_ptcl);
            _ptcl = &_particles.getLastMember();
        }

        //! struture to obtain the original particle index from binary_tree leaf iterator
        struct ParticleIndexAdr{
            int n_particle;
            ARPtcl* first_particle_index_local_address; // the local particle data first address to calculate local index
            int* particle_index_origin_output; // to store the particle index 
            const int* particle_index_origin_all;    // to obtain the original index from local index

            ParticleIndexAdr(ARPtcl* _first_adr, int* _index_output, const int* _index_org): n_particle(0), first_particle_index_local_address(_first_adr), particle_index_origin_output(_index_output), particle_index_origin_all(_index_org) {}
        };

        //! calculate particle index
        static void calcParticleIndex(ParticleIndexAdr& _index, ARPtcl*& _ptcl) {
            _index.particle_index_origin_output[_index.n_particle++] = _index.particle_index_origin_all[int(_ptcl - _index.first_particle_index_local_address)];
        }

    public:
        Float ds;  ///> estimated step size for AR integration
        Float dt_limit; ///> hermite time step limit for this group
        Float r_break_crit;    // group break radius criterion
        AR::FixStepOption fix_step_option; ///> fixt step option for integration
        bool need_resolve_flag; // indicate whether the members need to be resolved for outside
        COMM::List<int> particle_index; // particle index in original array (Hermite particles)
        COMM::List<COMM::BinaryTree<ARPtcl>> binarytree;

        ARInformation(): ds(Float(0.0)), dt_limit(NUMERIC_FLOAT_MAX), r_break_crit(-1.0), fix_step_option(AR::FixStepOption::none), need_resolve_flag(false), particle_index(), binarytree() {}

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            ASSERT(r_break_crit>=0.0);
            return true;
        }        
    
        //! reserve memory
        void reserveMem(const int _nmax) {
            particle_index.setMode(COMM::ListMode::local);
            binarytree.setMode(COMM::ListMode::local);
            binarytree.reserveMem(_nmax);
            particle_index.reserveMem(_nmax);
        }
    
        void clear() {
            ds=0.0;
            dt_limit = NUMERIC_FLOAT_MAX;
            r_break_crit=-1.0;
            fix_step_option = AR::FixStepOption::none;
            need_resolve_flag=false;
            particle_index.clear();
            binarytree.clear();
        }

        COMM::BinaryTree<ARPtcl>& getBinaryTreeRoot() const {
            int n = binarytree.getSize();
            ASSERT(n>0);
            return binarytree[n-1];
        }

        //! calculate ds from binary tree information and adjust by slowdown, determine the fix step option
        /*! Estimate ds first from binary orbit (eccentricity anomaly), then adjust based on original slowdown factor and order of integrator
          @param[in] _sd_org: original slowdown factor
          @param[in] _int_order: integrator order 
         */
        void calcDsAndStepOption(const Float _sd_org, const int _int_order) {
            ds = NUMERIC_FLOAT_MAX;
            ds = getBinaryTreeRoot().processRootIter(ds, calcDsKeplerIter);
            if (_sd_org<1.0) ds *= 1.0/8.0*pow(_sd_org, 1.0/Float(_int_order));
            auto& bin_root = getBinaryTreeRoot();
            const int n_particle = bin_root.getMemberN();

            // determine the fix step option
            fix_step_option = AR::FixStepOption::none;
            // for two-body case, determine the step at begining then fix
            if (n_particle==2) fix_step_option = AR::FixStepOption::later;
            // for multiple case, check whether outer peri-center is close to inner apo-center, if not, use fix step
            if (n_particle>2) {
                Float apo_in_max = 0;
                for (int j=0; j<2; j++) {
                    if (bin_root.getMember(j)->id<0) {
                        auto* bin_sub = (COMM::BinaryTree<ARPtcl>*) bin_root.getMember(j);
                        apo_in_max = std::max(apo_in_max,bin_sub->semi*(1+bin_sub->ecc));
                    }
                }
                Float peri_out = bin_root.semi*(1-bin_root.ecc);
                if (peri_out>3*apo_in_max) fix_step_option = AR::FixStepOption::later;
            }
        }

        //! generate binary tree for a particle group
        /*! @param[in] _particles: particle group 
         */
        void generateBinaryTree(COMM::ParticleGroup<ARPtcl, ParticleH4<Tparticle>>& _particles) {
            const int n_particle = _particles.getSize();
            ASSERT(n_particle>1);
            binarytree.resizeNoInitialize(n_particle-1);
            int particle_index_local[n_particle];
            for (int i=0; i<n_particle; i++) particle_index_local[i] = i;
            COMM::BinaryTree<ARPtcl>::generateBinaryTree(binarytree.getDataAddress(), particle_index_local, n_particle, _particles.getDataAddress());
        }

        //! Initialize group of particles from a binarytree
        /*! Add particles to local, copy Keplertree to local and relink the leaf particle address to local _particles data
          Make sure the original frame is used for particles linked in _bin
          @param[in] _particles: particle group to add particles
          @param[in] _bin: Binary tree root contain member particles
        */
        void addParticlesAndCopyBinaryTree(COMM::ParticleGroup<ARPtcl, ParticleH4<Tparticle>>& _particles, COMM::BinaryTree<ARPtcl>& _bin) {
            int n_members = _bin.getMemberN();
            // copy KeplerTree first
            binarytree.resizeNoInitialize(n_members-1);
            _bin.getherBinaryTreeIter(binarytree.getDataAddress());
            // add particle and relink the address in bin_ leaf
            ASSERT(binarytree.getLastMember().getMemberN() == n_members);
            binarytree.getLastMember().processLeafIter(_particles, addOneParticleAndReLinkPointer);
            // relink the original address based on the Particle local copy (adr_org) in ParticleAR type
            auto* padr_org = _particles.getMemberOriginAddress();
            for (int i=0; i<n_members; i++) {
                padr_org[i] = _particles[i].ARPtcl::adr;
            }
        }

        //! get two branch particle index
        /*!
          @param[out] _particle_index_origin_output: origin index of particles for output
          @param[in] _origin_particle_address: particle begining address to calculate index
          \return the split index to separate two branches, index is the right boundary 
         */
        int getTwoBranchParticleIndexOriginFromBinaryTree(int* _particle_index_origin_output, ARPtcl* _first_particle_address) {
            ParticleIndexAdr p_index(_first_particle_address, _particle_index_origin_output, particle_index.getDataAddress());
            auto& bin_root = binarytree.getLastMember();
            bin_root.processLeafIter(p_index, calcParticleIndex);
            if (bin_root.getMemberN()==2) return 1;
            else {
                auto* member1 = bin_root.getMember(0);
                if (member1->id<0) return ((COMM::BinaryTree<ARPtcl>*)member1)->getMemberN();
                else return 1;
            }
        }

        //! get dr * dv for two particles
        /*!
          @param[out] _dr2: dr*dr
          @param[out] _drdr: dr*dv
          @param[in] _p1: particle 1
          @param[in] _p2: particle 2
        */
        void getDrDv(Float& _dr2, Float& _drdv, ARPtcl& _p1, ARPtcl& _p2) {
            Float dx[3],dv[3];
            const Float* pos1 = _p1.getPos();
            const Float* pos2 = _p2.getPos();
            const Float* vel1 = _p1.getVel();
            const Float* vel2 = _p2.getVel();
            dx[0] = pos1[0] - pos2[0];
            dx[1] = pos1[1] - pos2[1];
            dx[2] = pos1[2] - pos2[2];

            dv[0] = vel1[0] - vel2[0];
            dv[1] = vel1[1] - vel2[1];
            dv[2] = vel1[2] - vel2[2];
        
            _dr2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
            _drdv= dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
        }
    };
}
