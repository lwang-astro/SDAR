#pragma once

//! Algorithmic regularization chain (ARC) namespace
/*!
  All major ARC classes and related acceleration functions (typedef) are defined
*/
namespace AR {

//! Particle group class to store and manage a group of particle
/*! A list that storing particle memory addresses and their copy (based on template class particle_)
 */
    template <class Tparticle>
    class ParticleGroup{
    private:
        int num_;          //!< number of current particles in the list #p
        int nmax_;         //!< maximum number of particles allocated in memory
        Tparticle* data_;      //!< particle array to store the data, or a link to existed particle array (not allocated)
        Tparticle** adr_; //!< original particle address of particle members
        
        //flag
        bool origin_frame_flag; //!< true: particles are in original frame; false: particles are in c.m. frame

    public:
        Tparticle cm;  //! center of mass particle for the group

        //! Constructor 
        /*! Set particle number to zero, clear pointers
         */
        chainlist(): num_(0), nmax_(0), data_(NULL), adr_(NULL), origin_frame_flag(true) {};
  
        //! Constructor with memory allocation for particle data
        /*! @param[in] _nmax: maximum number of particles for memory allocation
         */
        chainlist(const int _nmax) {
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
                data_ = NULL;
                adr_ = NULL;
            }
            nmax_ = 0;
        }

        // copy operator
        ParticleGroup& operator = (const ParticleGroup& _particle_group) {
            clear();
            // in allocated case
            if (_particle_group.nmax_>0) {
                reserveMem(_particle_group.nmax_); 
                num_ = _particle_group.num_; 
                for(int i=0; i<num; i++)  adr_[i] = _particle_group.adr_[i];
                for(int i=0; i<num; i++) data_[i] = _particle_group.data_[i];
            }
            else { // in link case
                num_ = _particle_group.num_; 
                data_ = _particle_group.data_;
            }
            return *this;
        }


        //! destructor
        ~chainlist() {
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
        Tparticle& getParticle(const int _index) const {
            assert(_index<num_);
            assert(_index>=0);
            return data_[_index];
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
        particle** getParticleOriginAddress(const int _index) const{
            assert(_index<num_);
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
            num++;
        }

        //! link a particle array load a list of particles
        /*! point the local particle data to the begining of the input \a _particle, no data copy, any modification later will be directly on original particle array
          @param [in] _particle: array of particles
          @param [in] _n: number of particles
        */
        void linkParticleList(Tparticle _particle[], const int _n) {
            assert(nmax_==0);
            assert(data_==NULL);
            data_= _particle;
            num_ = _number;
        }

        //! create remove table for removing a list of particles
        /*!
          @param[in] _index: particle index list for remove
          @param[in] _n_list: number of removing particles
          @param[in] _n_particle: number of total particles in the table
          @param[out] _remove_table: remove_table generated, the size should be at least number of particles to avoid overflow. for index i in table, value <0: no remove; <\a _number_particle: new position; = \a _number_particle: delete
         */
        static void createRemoveTable(const int* _index, const int _n_index; const int _n_particle, int* _remove_table) {
            assert(_n_index>0);
            assert(_n_particle>=_n_index);
            
            // last particle index 
            int ilast = _n_particle-1;
            const int n_new = _n_particle - _n_index; // new ptcl number

            for (int i=0; i<_number_index; i++) {
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
                    PS::S32 itrlast = -1;
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
                    for (int j=i; j<ilast; j++) {
                        data_[j] = data_[j+1];
                        adr_[j]  = adr_[j+1];
                    }
                }
            }
            num--;
        }


        //! remove a list of particle based on remove table
        /*!
          @param[in] _remove_table: table map the particle new position and deleted ones, generated from #createRemoveTable
         */
        void removeParticleList(const int* remove_table) {
            for (int i=0; i<num_; i++) {
                assert(remove_table[i]<=num_);
                // no change case
                if(remove_table[i]<0) continue;
                // remove case
                else if(remove_table[i]<num_) {
                    const int inew = remove_table[i];
                    data_[inew] = data_[i];
                    adr_[inew] = adr_[i];
                    num_--;
                }
            }
        }

        //! write particle data to files (notice original address is lost)
        /*! write particle data into file with BINARY format. Number of particles is written first, then the data of particles 
          @param [in] _fout: FILE IO for writing
        */
        void writeBinary(FILE* _fout) {
            fwrite(&num_,sizeof(int), 1, _fout);
            fwrite(data_, sizeof(Tparticle), num_, _fout);
        }

        //! Read particle data from file 
        /*! Read particle data from file with BINARY format. Number of particles should be first variable, then the data of particles.
          Notice the memory should be allocated first, and the free space is enough to save the reading particles
          @param [in] _fin: FILE IO for reading.
        */
        void readBINARY(FILE* _fin) {
            int n_new;
            int rn = fread(n_new,sizeof(int),1, _fin);
            if(rn<1) {
                std::cerr<<"Error: cannot read particle number!\n";
                abort();
            }
            if(n_new<=0) {
                std::cerr<<"Error: reading particle number <0!\n";
                abort();
            }
            if(n_new+num_>nmax_) {
                std::cerr<<"Error: reading overflow, local particle number is "<<num_<<" maximum space is "<<nmax_<<" reading number is "<<n_new<<std::endl;
                abort();
            }
            rn = fread(&data[i+num_], sizeof(Tparticle), n_new , _fin);
            if (rn<n_new) {
                std::cerr<<"Error: cannot read particle "<<rn<<", total particle number is "<<n_new<<"!\n";
                abort();
            }
            for (int i=0; i<n_new; i++) adr_[i+num_]=NULL;
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

        
        
    };
    
}
