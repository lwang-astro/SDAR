#pragma once

#include <cassert>
#include "Float.h"
#include "list.h"

namespace AR {

//! Particle group class to store and manage a group of particle
/*! A list that storing particle memory addresses and their copy (based on template class particle_)
 */
    template <class Tparticle>
    class ParticleGroup: public List<Tparticle>{
    private:
        typedef List<Tparticle> TList;
        bool origin_frame_flag; //!< true: particles are in original frame; false: particles are in c.m. frame

    public:
        Tparticle cm;  //! center of mass particle for the group

        //! Constructor 
        /*! Set particle number to zero, clear pointers
         */
        ParticleGroup(): TList(), origin_frame_flag(true) {};
  
        //! Clear function
        /*! Free dynamical memory space allocated
         */
        void clear() {
            TList::clear();
            origin_frame_flag = true;
        }

        //! operator = is copy
        /* Copy function will remove the local data and also copy the particle data or the link
         */
        ParticleGroup& operator = (const ParticleGroup& _particle_group) {
            TList::clear();
            *(List<Tparticle>*)this = *(List<Tparticle>*)&_particle_group;
            origin_frame_flag = _particle_group.origin_frame_flag;
            cm = _particle_group.cm;

            return *this;
        }


        //! destructor
        ~ParticleGroup() {
            clear();
        }

        //! Backup member particle position and velocity
        /*!
          \return backup array size
        */
        int backupParticlePosVel(Float* _bk) {
            for (int i=0; i<TList::num_; i++) {
                const int k=6*i;
                Float* pos = TList::data_[i].pos;
                Float* vel = TList::data_[i].vel;
                _bk[k  ] = pos[0];
                _bk[k+1] = pos[1];
                _bk[k+2] = pos[2];
                _bk[k+3] = vel[0];
                _bk[k+4] = vel[1];
                _bk[k+5] = vel[2];
            }
            return 6*TList::num_;
        }

        //! restore member particle position and velocity
        /*!
          \return backup array size
        */
        int restoreParticlePosVel(Float* _bk) {
            for (int i=0; i<TList::num_; i++) {
                const int k=6*i;
                Float* pos = TList::data_[i].pos;
                Float* vel = TList::data_[i].vel;
                pos[0] = _bk[k  ];
                pos[1] = _bk[k+1];
                pos[2] = _bk[k+2];
                vel[0] = _bk[k+3];
                vel[1] = _bk[k+4];
                vel[2] = _bk[k+5];
            }
            return 6*TList::num_;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"N";
            for (int i=0; i<TList::num_; i++) TList::data_[i].printColumnTitle(_fout, _width);
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<TList::num_;
            for (int i=0; i<TList::num_; i++) TList::data_[i].printColumn(_fout, _width);
        }

        //! write particle data to files (notice original address is lost)
        /*! write particle data into file with BINARY format. Number of particles is written first, then the data of particles 
          @param [in] _fout: FILE IO for writing
        */
        void writeBinary(FILE* _fout) {
            fwrite(&TList::num_,sizeof(int), 1, _fout);
            for (int i=0; i<TList::num_; i++) TList::data_[i].writeBinary(_fout);
        }

        ////! write particle data to files (notice original address is lost)
        ///*! write particle data into file with ASCII format. Number of particles is written first, then the data of particles 
        //  @param [in] _fout: FILE IO for writing
        //*/
        //void writeAscii(FILE* _fout) const{
        //    fprintf(_fout, "%d ", TList::num_);
        //    for (int i=0; i<TList::num_; i++) TList::data_[i].writeAscii(_fout);
        //}

        //! write particle data to files (notice original address is lost)
        /*! write particle data into file with ASCII format. Number of particles is written first, then the data of particles 
          @param [in] _fout: std::ostream IO for writing
        */
        void writeAscii(std::ostream& _fout) const{
            _fout<<TList::num_<<" ";
            for (int i=0; i<TList::num_; i++) TList::data_[i].writeAscii(_fout);
        }

        //! Read particle data from file 
        /*! Read particle data from file with BINARY format. Number of particles should be first variable, then the data of particles.
          Only work for #ListMode::local case
          @param [in] _fin: FILE IO for reading.
        */
        void readBinary(FILE* _fin) {
            assert(TList::mode_==ListMode::local);
            assert(TList::num_==0);
            assert(TList::nmax_==0);
            int n_new;
            int rn = fread(&n_new, sizeof(int),1, _fin);
            if(rn<1) {
                std::cerr<<"Error: cannot read particle number!\n";
                abort();
            }
            if(n_new<=0) {
                std::cerr<<"Error: reading particle number "<<n_new<<"<=0!\n";
                abort();
            }
            TList::reserveMem(n_new);
            for (int i=0; i<n_new; i++) TList::data_[i].readBinary(_fin);
            TList::num_ = n_new;
        }

        //! Read particle data from file 
        /*! Read particle data from file with BINARY format. Number of particles should be first variable, then the data of particles.
          Notice the memory should be allocated first, and the free space is enough to save the reading particles
          @param [in] _fin: std::istream IO for reading.
        */
        void readAscii(std::istream& _fin) {
            assert(TList::mode_==ListMode::local);
            assert(TList::num_==0);
            assert(TList::nmax_==0);
            int n_new;
            _fin>>n_new;
            if(n_new<=0) {
                std::cerr<<"Error: reading particle number "<<n_new<<"<=0!\n";
                abort();
            }
            TList::reserveMem(n_new);
            for (int i=0; i<n_new; i++) TList::data_[i].readAscii(_fin);
            TList::num_ = n_new;
        }

        //! shift particle to their c.m. frame
        /*! Shift positions and velocities of particles from original frame to their center-of-mass frame\n
          Notice the center-of-mass position and velocity use values from #cm
        */
        void shiftToCM() {
            if (origin_frame_flag) {
                const Float *rc = cm.pos;
                const Float *vc = cm.vel;
                for (int i=0;i<TList::num_;i++) {
                    Float *ri = TList::data_[i].pos;
                    Float *vi = TList::data_[i].vel;
                    ri[0] -= rc[0];
                    ri[1] -= rc[1];
                    ri[2] -= rc[2];
                    vi[0] -= vc[0];
                    vi[1] -= vc[1];
                    vi[2] -= vc[2];
                }
                origin_frame_flag = false;
            }
            else {
                std::cerr<<"Warning: particles are already in the center-of-mass frame!\n";
            }      
        }
        
        //! shift particle to their original frame
        /*! Shift positions and velocities of particles from center-of-mass frame to original frame\n
          Notice the center-of-mass position and velocity use values from #cm
        */
        void shiftToOrigin() {
            if (origin_frame_flag) {
                std::cerr<<"Warning: particles are already in original frame!\n";
            }
            else {
                const Float *rc = cm.pos;
                const Float *vc = cm.vel;
                for (int i=0;i<TList::num_;i++) {
                    Float *ri = TList::data_[i].pos;
                    Float *vi = TList::data_[i].vel;
                    ri[0] += rc[0];
                    ri[1] += rc[1];
                    ri[2] += rc[2];
                    vi[0] += vc[0];
                    vi[1] += vc[1];
                    vi[2] += vc[2];
                }
                origin_frame_flag = true;
            }
        }

        //! calculate center-of-mass
        void calcCenterOfMass() {
            cm.pos[0] = cm.pos[1] = cm.pos[2] = 0.0;
            cm.vel[0] = cm.vel[1] = cm.vel[2] = 0.0;
            cm.mass = 0.0;
            for (int i=0;i<TList::num_;i++) {
                const Float *ri = TList::data_[i].pos;
                const Float *vi = TList::data_[i].vel;
                const Float mi = TList::data_[i].mass;
                cm.pos[0] += ri[0] * mi;
                cm.pos[1] += ri[1] * mi;
                cm.pos[2] += ri[2] * mi;

                cm.vel[0] += vi[0] * mi;
                cm.vel[1] += vi[1] * mi;
                cm.vel[2] += vi[2] * mi;

                cm.mass += mi;
            }
            cm.pos[0] /= cm.mass; 
            cm.pos[1] /= cm.mass; 
            cm.pos[2] /= cm.mass; 
            cm.vel[0] /= cm.mass; 
            cm.vel[1] /= cm.mass; 
            cm.vel[2] /= cm.mass;
        }

    };
    
}
