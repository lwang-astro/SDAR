#pragma once

#include <iomanip>
#include "Float.h"
#include "list.h"

namespace COMM {

//! Particle group class to store and manage a group of particle
/*! A list that storing particle memory addresses and their copy (based on template class particle_)
 */
    template <class Tparticle, class Tpcm>
    class ParticleGroup: public List<Tparticle>{
    private:
        typedef List<Tparticle> TList;
        bool origin_frame_flag; //!< true: particles are in original frame; false: particles are in c.m. frame

    public:
        Tpcm cm;  //! center of mass particle for the group

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

        //! get backup data size 
        /*! \return the data array size for backupParticlePosVel()
         */
        int getBackupDataSize() const {
            return 6*TList::num_;
        }

        //! Backup member particle position and velocity
        /*!
          \return backup array size
        */
        int backupParticlePosVel(Float* _bk) {
            for (int i=0; i<TList::num_; i++) {
                const int k=6*i;
                Float* pos = TList::data_[i].getPos();
                Float* vel = TList::data_[i].getVel();
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
                Float* pos = TList::data_[i].getPos();
                Float* vel = TList::data_[i].getVel();
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
            int num = TList::num_;
            fwrite(&num,sizeof(int), 1, _fout);
            for (int i=0; i<TList::num_; i++) TList::data_[i].writeBinary(_fout);
            fwrite(&origin_frame_flag, sizeof(bool), 1, _fout);
            cm.writeBinary(_fout);
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
        void writeMemberAscii(std::ostream& _fout) const{
            _fout<<TList::num_<<" ";
            for (int i=0; i<TList::num_; i++) TList::data_[i].writeAscii(_fout);
        }

        //! Read particle data from file 
        /*! Read particle data from file with BINARY format. Number of particles should be first variable, then the data of particles.
          Only work for #ListMode::local case
          @param [in] _fin: FILE IO for reading.
        */
        void readBinary(FILE* _fin) {
            ASSERT(TList::mode_==ListMode::local);
            ASSERT(TList::num_==0);
            ASSERT(TList::nmax_==0);
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
            rn = fread(&origin_frame_flag, sizeof(bool), 1, _fin);
            if(rn<1) {
                std::cerr<<"Error: cannot read origin_frame_flag!\n";
                abort();
            }
            cm.readBinary(_fin);
        }

        //! Read particle data from file 
        /*! Read particle data from file with BINARY format. Number of particles should be first variable, then the data of particles.
          Notice the memory should be allocated first, and the free space is enough to save the reading particles
          @param [in] _fin: std::istream IO for reading.
        */
        void readMemberAscii(std::istream& _fin) {
            ASSERT(TList::mode_==ListMode::local);
            ASSERT(TList::num_==0);
            ASSERT(TList::nmax_==0);
            int n_new;
            _fin>>n_new;
            ASSERT(!_fin.eof());
            if(n_new<=0) {
                std::cerr<<"Error: reading particle number "<<n_new<<"<=0!\n";
                abort();
            }
            TList::reserveMem(n_new);
            for (int i=0; i<n_new; i++) {
                TList::data_[i].readAscii(_fin);
                ASSERT(!_fin.eof());
            }
            TList::num_ = n_new;
        }

        //! shift particle to their c.m. frame
        /*! Shift positions and velocities of particles from original frame to their center-of-mass frame\n
          Notice the center-of-mass position and velocity use values from #cm
        */
        void shiftToCenterOfMassFrame() {
            if (origin_frame_flag) {
                const Float *rc = cm.getPos();
                const Float *vc = cm.getVel();
                for (int i=0;i<TList::num_;i++) {
                    Float *ri = TList::data_[i].getPos();
                    Float *vi = TList::data_[i].getVel();
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
        void shiftToOriginFrame() {
            if (origin_frame_flag) {
                std::cerr<<"Warning: particles are already in original frame!\n";
            }
            else {
                const Float *rc = cm.getPos();
                const Float *vc = cm.getVel();
                for (int i=0;i<TList::num_;i++) {
                    Float *ri = TList::data_[i].getPos();
                    Float *vi = TList::data_[i].getVel();
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
            Float *rc = cm.getPos();
            Float *vc = cm.getVel();
            rc[0] = rc[1] = rc[2] = 0.0;
            vc[0] = vc[1] = vc[2] = 0.0;
            cm.mass = 0.0;
            for (int i=0;i<TList::num_;i++) {
                const Float *ri = TList::data_[i].getPos();
                const Float *vi = TList::data_[i].getVel();
                const Float mi = TList::data_[i].mass;
                rc[0] += ri[0] * mi;
                rc[1] += ri[1] * mi;
                rc[2] += ri[2] * mi;

                vc[0] += vi[0] * mi;
                vc[1] += vi[1] * mi;
                vc[2] += vi[2] * mi;

                cm.mass += mi;
            }
            rc[0] /= cm.mass; 
            rc[1] /= cm.mass; 
            rc[2] /= cm.mass; 
            vc[0] /= cm.mass; 
            vc[1] /= cm.mass; 
            vc[2] /= cm.mass;
        }
        
        //! return true if the system is the in their origin frame
        bool isOriginFrame() const {
            return origin_frame_flag;
        }
    };

}
