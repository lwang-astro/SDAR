#pragma once

#include <cassert>
#include "Float.h"

namespace AR {

//! Particle group class to store and manage a group of particle
/*! A list that storing particle memory addresses and their copy (based on template class particle_)
 */
    template <class Tparticle>
    class ParticleGroup{
    private:
        int num_;          //!< number of current particles in the list 
        int nmax_;         //!< maximum number of particles allocated in memory
        Tparticle* data_;      //!< particle array to store the data, or a link to existed particle array (not allocated)
        Tparticle** adr_; //!< original particle address of particle members
        
        //flag
        bool origin_frame_flag; //!< true: particles are in original frame; false: particles are in c.m. frame
        bool modified_flag;     //!< true: particle list is modified; used for safety checking

    public:
        Tparticle cm;  //! center of mass particle for the group

        //! Constructor 
        /*! Set particle number to zero, clear pointers
         */
        ParticleGroup(): num_(0), nmax_(0), data_(NULL), adr_(NULL), origin_frame_flag(true), modified_flag(false) {};
  
        //! Constructor with memory allocation for particle data
        /*! @param[in] _nmax: maximum number of particles for memory allocation
         */
        ParticleGroup(const int _nmax) {
            reserveMem(_nmax); 
        }

        //! Memory allocation for storing particles
        /*! allocate memory for storing particles and their original addresses.
          @param[in] _nmax: maximum number of particles for memory allocation
         */
        void reserveMem(const int _nmax) {
            assert(nmax_==0);
            assert(_nmax>0);
            assert(adr_==NULL);
            assert(data_==NULL);

            adr_=new Tparticle*[_nmax];
            data_=new Tparticle[_nmax];
            nmax_ = _nmax;
        }

        //! Clear function
        /*! Free dynamical memory space allocated
         */
        void clear() {
            if (nmax_>0) {
                assert(data_!=NULL);
                assert(adr_!=NULL);
                delete[] data_;
                delete[] adr_;
                adr_ = NULL;
            }
            data_ = NULL; // in link case also set to NULL
            num_ = 0;
            nmax_ = 0;
            origin_frame_flag = true;
            modified_flag = false;
        }

        //! operator = is copy
        /* Copy function will remove the local data and also copy the particle data or the link
         */
        ParticleGroup& operator = (const ParticleGroup& _particle_group) {
            clear();
            // in allocated case
            if (_particle_group.nmax_>0) {
                reserveMem(_particle_group.nmax_); 
                num_ = _particle_group.num_; 
                for(int i=0; i<num_; i++)  adr_[i] = _particle_group.adr_[i];
                for(int i=0; i<num_; i++) data_[i] = _particle_group.data_[i];
            }
            else { // in link case
                num_ = _particle_group.num_; 
                data_ = _particle_group.data_;
            }
            origin_frame_flag = _particle_group.origin_frame_flag;
            modified_flag     = _particle_group.modified_flag;
            cm = _particle_group.cm;

            return *this;
        }


        //! destructor
        ~ParticleGroup() {
            clear();
        }

        //! get current particle number
        /*! \return Current particle number in particle address list #p
         */
        int getParticleNumber() const {
            return num_;
        }

        //! return one particle data reference
        /*!
          @param [in] _index: index of particle in group
         */
        Tparticle& getParticle(const int _index) {
            assert(_index<num_);
            assert(_index>=0);
            return data_[_index];
        }

        //! return particle data array address
        Tparticle* getParticleDataAddress() {
            return data_;
        }

        //! return one particle data reference operator
        /*! overload []
          @param [in] _index: index of particle in group
        */
        Tparticle& operator [](const int _index) const {
            assert(_index<num_);
            assert(_index>=0);
            return data_[_index];
        }

        //! return one particle original address 
        /*!
          \return: one particle original address array 
        */
        Tparticle* getParticleOriginAddress(const int _index) const{
            assert(_index<num_);
            assert(_index<nmax_);
            assert(_index>=0);
            return adr_[_index];
        }

        //! Get maximum particle number allow to store 
        /*! \return allocated number of particles in memory
         */
        int getParticleMemNumberMax() const {
            return nmax_;
        }

        //! store one particle
        /*! add a particle by copying particle into local allocated memory and storing the original memory address of it
          @param [in] a: a particle to store (push back at the end of local array)
        */
        void addParticle(Tparticle &_particle) {
            assert(num_<nmax_);
            adr_[num_]  = &_particle;
            data_[num_] = _particle;
            num_++;
            modified_flag = true;
        }

        //! link a particle array load a list of particles
        /*! point the local particle data to the begining of the input \a _particle, no data copy, any modification later will be directly on original particle array
          @param [in] _particle: array of particles
          @param [in] _n_particle: number of particles
        */
        void linkParticleList(Tparticle _particle[], const int _n_particle) {
            assert(nmax_==0);
            assert(data_==NULL);
            data_= _particle;
            num_ = _n_particle;
            modified_flag = true;
        }

        //! create remove table for removing a list of particles
        /*!
          @param[in] _index: particle index list for remove
          @param[in] _n_list: number of removing particles
          @param[in] _n_particle: number of total particles in the table
          @param[out] _remove_table: remove_table generated, the size should be at least number of particles to avoid overflow. for index i in table, value <0: no remove; <\a _number_particle: new position; = \a _number_particle: delete
         */
        static void createRemoveTable(const int* _index, const int _n_index, const int _n_particle, int* _remove_table) {
            assert(_n_index>0);
            assert(_n_particle>=_n_index);
            
            // last particle index 
            int ilast = _n_particle-1;
            const int n_new = _n_particle - _n_index; // new ptcl number
            // initial table
            for (int i=0; i<_n_particle; i++) _remove_table[i] = -1;

            for (int i=0; i<_n_index; i++) {
                int k=_index[i];
                
                int idel = k;
                
                // set ilast >= idel for safety.
                ilast = std::max(ilast, idel);
            
                // check whether position is already moved, if so, check the moved new position and set current to delete
                if(_remove_table[k]>=0) {
                    idel=_remove_table[k];
                    _remove_table[k]=_n_particle;
                }

                assert(k<_n_particle);
                assert(idel<_n_particle);
                
                // check the last avaiable particle that can be moved to idel
                // If the ilast is already moved, check whether the new moved position _remove_table[ilast] is before the current idel.
                // If _remove_table[ilast] is after the current idel, update the _remove_table[ilast] to idel and set _remove_table[ilast] as new idel to check, until _remove_table[ilast]=-1 or idel >= ilast
                while (idel>=0) {
                    int itrlast = -1;
                    //cond: last is moved && ( moved pos before new n_ptcl || last is del ) &&  ilast > idel 
                    while(_remove_table[ilast]>=0 && (_remove_table[ilast]<n_new || _remove_table[ilast]==_n_particle) && ilast>idel) ilast--;

                    // if ilast is not yet moved or (new pos of ilast > idel and ilast is not del, move ilast to idel)
                    if(_remove_table[ilast]<0||(idel<_remove_table[ilast]&&_remove_table[ilast]<_n_particle)) {
                        itrlast=_remove_table[ilast];
                        _remove_table[ilast]=idel;
                    }
                    // idel is already at last, remove it
                    _remove_table[idel]=_n_particle; 
                    idel = itrlast;
                }
            }
        }
    
        //! remove one particle
        /*! remove one particle from local particle array
          @param [in] _index: particle index to be removed
          @param [in] _shift_last_only_flag: true: shift last particle to current position (defaulted); false: shift all right particle to left by one
        */
        void removeParticle(const int _index, const bool _shift_last_only_flag) {
            assert(_index>=0);
            assert(_index<num_);

            int ilast = num_-1;
            if (_index<ilast) {
                if (_shift_last_only_flag) {
                    data_[_index] = data_[ilast];
                    adr_[_index] = adr_[ilast];
                }
                else {
                    for (int j=_index; j<ilast; j++) {
                        data_[j] = data_[j+1];
                        adr_[j]  = adr_[j+1];
                    }
                }
            }
            num_--;
            modified_flag = true;
        }


        //! remove a list of particle based on remove table
        /*!
          @param[in] _remove_table: table map the particle new position and deleted ones, generated from #createRemoveTable
         */
        void removeParticleTable(const int* remove_table) {
            const int num_org = num_;
            for (int i=0; i<num_org; i++) {
                assert(remove_table[i]<=num_org);
                // no change case
                if(remove_table[i]<0) continue;
                // move case
                else if(remove_table[i]<num_org) {
                    const int inew = remove_table[i];
                    data_[inew] = data_[i];
                    if(nmax_>0) adr_[inew] = adr_[i];
                }
                // remove count
                else num_--;
            }
            modified_flag = true;
        }

        //! remove a list of particle based on an index list
        /*!
          @param[in] _index: particle index list for remove
          @param[in] _n_list: number of removing particles
         */
        void removeParticleList(const int* _index, const int _n_index) {
            assert(_n_index<=num_);
            const int num_org = num_;
            int remove_table[num_];
            createRemoveTable(_index, _n_index, num_, remove_table);
            removeParticleTable(remove_table);
            assert(num_+_n_index==num_org);
        }

        //! Backup member particle position and velocity
        /*!
          \return backup array size
        */
        int backupParticlePosVel(Float* _bk) {
            for (int i=0; i<num_; i++) {
                const int k=6*i;
                Float* pos = data_[i].pos;
                Float* vel = data_[i].vel;
                _bk[k  ] = pos[0];
                _bk[k+1] = pos[1];
                _bk[k+2] = pos[2];
                _bk[k+3] = vel[0];
                _bk[k+4] = vel[1];
                _bk[k+5] = vel[2];
            }
            return 6*num_;
        }

        //! restore member particle position and velocity
        /*!
          \return backup array size
        */
        int restoreParticlePosVel(Float* _bk) {
            for (int i=0; i<num_; i++) {
                const int k=6*i;
                Float* pos = data_[i].pos;
                Float* vel = data_[i].vel;
                pos[0] = _bk[k  ];
                pos[1] = _bk[k+1];
                pos[2] = _bk[k+2];
                vel[0] = _bk[k+3];
                vel[1] = _bk[k+4];
                vel[2] = _bk[k+5];
            }
            return 6*num_;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"N";
            for (int i=0; i<num_; i++) data_[i].printColumnTitle(_fout, _width);
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<num_;
            for (int i=0; i<num_; i++) data_[i].printColumn(_fout, _width);
        }

        //! write particle data to files (notice original address is lost)
        /*! write particle data into file with BINARY format. Number of particles is written first, then the data of particles 
          @param [in] _fout: FILE IO for writing
        */
        void writeBinary(FILE* _fout) {
            fwrite(&num_,sizeof(int), 1, _fout);
            for (int i=0; i<num_; i++) data_[i].writeBinary(_fout);
        }

        ////! write particle data to files (notice original address is lost)
        ///*! write particle data into file with ASCII format. Number of particles is written first, then the data of particles 
        //  @param [in] _fout: FILE IO for writing
        //*/
        //void writeAscii(FILE* _fout) const{
        //    fprintf(_fout, "%d ", num_);
        //    for (int i=0; i<num_; i++) data_[i].writeAscii(_fout);
        //}

        //! write particle data to files (notice original address is lost)
        /*! write particle data into file with ASCII format. Number of particles is written first, then the data of particles 
          @param [in] _fout: std::ostream IO for writing
        */
        void writeAscii(std::ostream& _fout) const{
            _fout<<num_<<" ";
            for (int i=0; i<num_; i++) data_[i].writeAscii(_fout);
        }

        //! Read particle data from file 
        /*! Read particle data from file with BINARY format. Number of particles should be first variable, then the data of particles.
          Notice the memory should be allocated first, and the free space is enough to save the reading particles
          @param [in] _fin: FILE IO for reading.
        */
        void readBinary(FILE* _fin) {
            int n_new;
            int rn = fread(n_new,sizeof(int),1, _fin);
            if(rn<1) {
                std::cerr<<"Error: cannot read particle number!\n";
                abort();
            }
            if(n_new<=0) {
                std::cerr<<"Error: reading particle number "<<n_new<<"<=0!\n";
                abort();
            }
            if(nmax_>0) {
                if (n_new+num_>nmax_) {
                    std::cerr<<"Error: reading overflow, local particle number is "<<num_<<" maximum space is "<<nmax_<<" reading number is "<<n_new<<std::endl;
                    abort();
                }
            }
            else {
                if(num_>0) {
                    std::cerr<<"Error: link particle data exist, cannot read more particles\n";
                    abort();
                }
                reserveMem(n_new);
            }
            for (int i=num_; i<num_+n_new; i++) data_[i].readBinary(_fin);
            for (int i=num_; i<num_+n_new; i++) adr_[i]=NULL;
            num_ += n_new;
        }

        ////! Read particle data from file 
        ///*! Read particle data from file with BINARY format. Number of particles should be first variable, then the data of particles.
        //  Notice the memory should be allocated first, and the free space is enough to save the reading particles
        //  @param [in] _fin: FILE IO for reading.
        //*/
        //void readAscii(FILE* _fin) {
        //    int n_new;
        //    int rn = fscanf(_fin, "%d ", &n_new);
        //    if(rn<1) {
        //        std::cerr<<"Error: cannot read particle number!\n";
        //        abort();
        //    }
        //    if(n_new<=0) {
        //        std::cerr<<"Error: reading particle number <0!\n";
        //        abort();
        //    }
        //    if(n_new+num_>nmax_) {
        //        std::cerr<<"Error: reading overflow, local particle number is "<<num_<<" maximum space is "<<nmax_<<" reading number is "<<n_new<<std::endl;
        //        abort();
        //    }
        //    for (int i=num_; i<num_+n_new; i++) data_[i].readAscii(_fin);
        //    for (int i=num_; i<num_+n_new; i++) adr_[i]=NULL;
        //    num_ += n_new;
        //}

        //! Read particle data from file 
        /*! Read particle data from file with BINARY format. Number of particles should be first variable, then the data of particles.
          Notice the memory should be allocated first, and the free space is enough to save the reading particles
          @param [in] _fin: std::istream IO for reading.
        */
        void readAscii(std::istream& _fin) {
            int n_new;
            _fin>>n_new;
            if(n_new<=0) {
                std::cerr<<"Error: reading particle number "<<n_new<<"<=0!\n";
                abort();
            }
            if(nmax_>0) {
                if (n_new+num_>nmax_) {
                    std::cerr<<"Error: reading overflow, local particle number is "<<num_<<" maximum space is "<<nmax_<<" reading number is "<<n_new<<std::endl;
                    abort();
                }
            }
            else {
                if(num_>0) {
                    std::cerr<<"Error: link particle data exist, cannot read more particles\n";
                    abort();
                }
                reserveMem(n_new);
            }
            for (int i=num_; i<num_+n_new; i++) data_[i].readAscii(_fin);
            for (int i=num_; i<num_+n_new; i++) adr_[i]=NULL;
            num_ += n_new;
        }

        //! copy all particle data back to original address
        void writeBackParticleAll() {
            assert(nmax_>0);
            for (int i=0; i<num_; i++) {
                if(adr_[i]!=NULL) *adr_[i] = data_[i];
            }
        }

        //! copy a list of particle data back to original address
        /*!
          @param[in] _index: list of particle index
          @param[in] _n:  number of particles
         */
        void writeBackParticleList(const int* _index, const int _n) {
            assert(nmax_>0);
            for (int i=0; i<_n; i++) {
                int k=_index[i];
                if(adr_[k]!=NULL) *adr_[k] = data_[k];
            }
        }

        //! shift particle to their c.m. frame
        /*! Shift positions and velocities of particles from original frame to their center-of-mass frame\n
          Notice the center-of-mass position and velocity use values from #cm
        */
        void shiftToCM() {
            if (origin_frame_flag) {
                const Float *rc = cm.pos;
                const Float *vc = cm.vel;
                for (int i=0;i<num_;i++) {
                    Float *ri = data_[i].pos;
                    Float *vi = data_[i].vel;
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
                for (int i=0;i<num_;i++) {
                    Float *ri = data_[i].pos;
                    Float *vi = data_[i].vel;
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
            for (int i=0;i<num_;i++) {
                const Float *ri = data_[i].pos;
                const Float *vi = data_[i].vel;
                const Float mi = data_[i].mass;
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

        //! Get modified status
        bool isModified() const {
            return modified_flag;
        }

        //! Reset modified status to false
        void setModifiedFalse() {
            modified_flag = false;
        }
    };
    
}
