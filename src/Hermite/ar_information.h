#pragma once

#include "Common/Float.h"
#include "Common/list.h"
#include "Common/binary_tree.h"
#include "Hermite/hermite_particle.h"
#include "AR/information.h"

namespace H4{
    //! contain group information

    template <class Tparticle>
    class ARInformation: public AR::Information<ParticleAR<Tparticle>, ParticleH4<Tparticle>> {
    private:
        typedef ParticleAR<Tparticle> ARPtcl;
        typedef AR::Information<ParticleAR<Tparticle>, ParticleH4<Tparticle>> ARInfoBase;

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
        Float dt_limit;       ///> hermite time step limit for this group
        Float r_break_crit;    // group break radius criterion
        COMM::List<int> particle_index; // particle index in original array (Hermite particles)

        ARInformation(): ARInfoBase(), dt_limit(NUMERIC_FLOAT_MAX), r_break_crit(-1.0), particle_index() {}

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            ASSERT(r_break_crit>=0.0);
            return true;
        }        

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            ARInfoBase::printColumnTitle(_fout, _width);
            _fout<<std::setw(_width)<<"dt_limit"
                 <<std::setw(_width)<<"r_break_crit";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            ARInfoBase::printColumn(_fout, _width);
            _fout<<std::setw(_width)<<dt_limit
                 <<std::setw(_width)<<r_break_crit;
        }

    
        //! reserve memory
        void reserveMem(const int _nmax) {
            ARInfoBase::reserveMem(_nmax);
            particle_index.setMode(COMM::ListMode::local);
            particle_index.reserveMem(_nmax);
        }
    
        void clear() {
            ARInfoBase::clear();
            dt_limit = NUMERIC_FLOAT_MAX;
            r_break_crit=-1.0;
            particle_index.clear();
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
            ARInfoBase::binarytree.resizeNoInitialize(n_members-1);
            _bin.getherBinaryTreeIter(ARInfoBase::binarytree.getDataAddress());
            // add particle and relink the address in bin_ leaf
            ASSERT(ARInfoBase::binarytree.getLastMember().getMemberN() == n_members);
            ARInfoBase::binarytree.getLastMember().processLeafIter(_particles, addOneParticleAndReLinkPointer);
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
            auto& bin_root = ARInfoBase::binarytree.getLastMember();
            bin_root.processLeafIter(p_index, calcParticleIndex);
            if (bin_root.getMemberN()==2) return 1;
            else {
                if (bin_root.isMemberTree(0)) return bin_root.getLeftMemberAsTree()->getMemberN();
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
