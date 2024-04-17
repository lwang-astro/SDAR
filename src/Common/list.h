#pragma once

#include "Common/Float.h"

namespace COMM {

    //! list mode identification
    /*! copy: copy data and link original address (data_, adr_ allocated); \n
       local: allocate data_ memory locally, no allocation of adr_ (no original address); \n
       link: no allocation, link data from outside; \n
       none: initial case.
     */
    enum class ListMode {copy, local, link ,none};

    //! list class to store and manage a group of member
    /*! A list that storing member memory addresses and their copy (based on template class Ttype)
     */
    template <class Ttype>
    class List{
    protected:
        int num_;          //!< number of current members in the list 
        int nmax_;         //!< maximum number of members allocated in memory
        Ttype* data_;      //!< member array to store the data, or a link to existed member array (not allocated)
        Ttype** adr_;      //!< original member address of member members
        ListMode mode_;    //!< mode indicator

        //flag
        bool modified_flag_;     //!< true: member list is modified; used for safety checking
    public:

        //! Constructor 
        /*! Set member number to zero, clear pointers
         */
        List(): num_(0), nmax_(0), data_(NULL), adr_(NULL), mode_(ListMode::none), modified_flag_(false) {};

        //! set mode
        /*! Only set once, to set new mode, clear must be done first
         */
        void setMode(const ListMode _mode) {
            ASSERT(mode_==ListMode::none);
            mode_ = _mode;
        }
        
        //! Memory allocation for storing members
        /*! allocate memory for storing members and their original addresses.
          @param[in] _nmax: maximum number of members for memory allocation
         */
        void reserveMem(const int _nmax) {
            ASSERT(nmax_==0);
            ASSERT(_nmax>0);
            ASSERT(adr_==NULL);
            ASSERT(data_==NULL);
            ASSERT(mode_!=ListMode::none);

            switch(mode_) {
            case ListMode::copy: 
                data_=new Ttype[_nmax];
                adr_=new Ttype*[_nmax];
                nmax_ = _nmax;
                break;
            case ListMode::local:
                data_=new Ttype[_nmax];
                nmax_ = _nmax;
                break;
            case ListMode::link:
                break;
            case ListMode::none:
                break;
            }
        }


        //! Clear function
        /*! Free dynamical memory space allocated
         */
        void clear() {
            switch(mode_) {
            case ListMode::copy:
                ASSERT(data_!=NULL);
                ASSERT(adr_!=NULL);
                delete[] data_;
                delete[] adr_;
                num_ = 0;
                nmax_ =0;
                data_= NULL;
                adr_ = NULL;
                break;
            case ListMode::local:
                ASSERT(data_!=NULL);
                ASSERT(adr_==NULL);
                delete[] data_;
                num_ = 0;
                nmax_ =0;
                data_=NULL;
                break;
            case ListMode::link:
                ASSERT(adr_==NULL);
                num_ = 0;
                nmax_ = 0;
                data_=NULL;
                break;
            case ListMode::none:
                ASSERT(num_==0);
                ASSERT(nmax_==0);
                ASSERT(data_==NULL);
                ASSERT(adr_==NULL);
                break;
            }
            mode_ = ListMode::none;
            modified_flag_ = false;
        }
    
        //! operator = is copy
        /* Copy function will remove the local data and also copy the member data or the link
         */
        List& operator = (const List& _list) {
            clear();
            mode_ = _list.mode_;
            switch(mode_) {
            case ListMode::copy:
                reserveMem(_list.nmax_, mode_); 
                num_ = _list.num_; 
                for(int i=0; i<num_; i++) data_[i] = _list.data_[i];
                for(int i=0; i<num_; i++)  adr_[i] = _list.adr_[i];
                break;
            case ListMode::local:
                reserveMem(_list.nmax_, mode_);
                num_ = _list.num_; 
                for(int i=0; i<num_; i++) data_[i] = _list.data_[i];
                break;
            case ListMode::link:
                num_ = _list.num_; 
                nmax_ = _list.nmax_;
                data_ = _list.data_;
                break;
            case ListMode::none:
                break;
            }
            modified_flag_     = _list.modified_flag_;
            return *this;
        }

        //! destructor
        ~List() {
            clear();
        }

        //! get current member number
        /*! \return Current member number in member address list #data_
         */
        int getSize() const {
            return num_;
        }

        ListMode getMode() const {
            return mode_;
        }
    
        //! return one member data reference
        /*!
          @param [in] _index: index of member in group
         */
        Ttype& getMember(const int _index) {
            ASSERT(_index<num_);
            ASSERT(_index>=0);
            return data_[_index];
        }

        //! return member data array address
        Ttype* getDataAddress() const {
            return data_;
        }

        //! return member original address array
        Ttype** getOriginAddressArray() const {
            return adr_;
        }

        //! return last member
        Ttype& getLastMember() {
            return data_[num_-1];
        }

        //! return one member data reference operator
        /*! overload []
          @param [in] _index: index of member in group
        */
        Ttype& operator [](const int _index) const {
            ASSERT(_index<num_);
            ASSERT(_index>=0);
            return data_[_index];
        }

        //! return one member original address 
        /*! Only work when mode is #ListMode::copy
          \return: one member original address array 
        */
        Ttype* getMemberOriginAddress(const int _index) const{
            ASSERT(mode_==ListMode::copy);
            ASSERT(_index<num_);
            ASSERT(_index<nmax_);
            ASSERT(_index>=0);
            return adr_[_index];
        }

        //! Get maximum member number allow to store 
        /*! \return allocated number of members in memory
         */
        int getSizeMax() const {
            return nmax_;
        }
    
        //! copy one member and its address
        /*! add a member by copying member into local allocated memory and storing the original memory address of it.
          Only work for #ListMode::copy
          @param [in] _member: a member to store (push back at the end of local array)
        */
        template <class T>
        void addMemberAndAddress(T &_member) {
            ASSERT(mode_==ListMode::copy);
            ASSERT(num_<nmax_);
            data_[num_] = _member;
            adr_[num_]  = (Ttype*)&_member;
            num_++;
            modified_flag_ = true;
        }

        //! copy one member and its address
        /*! add a member by copying it into local allocated memory.
          Only work for #ListMode::local
          @param [in] _member: a member to store (push back at the end of local array)
        */
        template <class T>
        void addMember(const T &_member) {
            ASSERT(mode_==ListMode::local);
            ASSERT(num_<nmax_);
            data_[num_] = _member;
            num_++;
            modified_flag_ = true;
        }

        //! increase size without initialization 
        /*! Work for #ListMode::local case. 
          Require the memory size is big enough to contain the new number of members.
         */
        void increaseSizeNoInitialize(const int _n) {
            ASSERT(mode_==ListMode::local||mode_==ListMode::link);
            ASSERT(num_+_n<=nmax_);
            num_ += _n;
            modified_flag_ = true;
        }

        //! increase size without initialization 
        /*! Work for #ListMode::local and link case. 
          Require the memory size is big enough to contain the new number of members.
         */
        void decreaseSizeNoInitialize(const int _n) {
            ASSERT(mode_==ListMode::local||mode_==ListMode::link);
            ASSERT(num_-_n>=0);
            num_ -= _n;
            modified_flag_ = true;
        }

        //! increase size without initialization 
        /*! Work for #ListMode::local case. 
          Require the memory size is big enough to contain the new number of members.
         */
        void resizeNoInitialize(const int _n) {
            ASSERT(mode_==ListMode::local||mode_==ListMode::link);
            ASSERT(_n<=nmax_);
            num_ = _n;
            modified_flag_ = true;
        }
    
        //! link a member array 
        /*! point the local member data to the begining of the input \a _member, no data copy, any modification later will be directly on original member array. Work for #ListMode::link
          @param [in] _member: array of members
          @param [in] _n_member: number of members
        */
        void linkMemberArray(Ttype _member[], const int _n_member) {
            ASSERT(mode_==ListMode::link);
            ASSERT(nmax_==0);
            ASSERT(num_==0);
            ASSERT(data_==NULL);
            data_= _member;
            num_ = _n_member;
            nmax_= _n_member;
            modified_flag_ = true;
        }

        //! create remove table for removing a list of members
        /*!
          @param[in] _index: member index list for remove
          @param[in] _n_index: number of removing member indice
          @param[in] _n_member: number of total members in the table
          @param[out] _remove_table: remove_table generated, the size should be at least number of members to avoid overflow. for index i in table, value <0: no remove; <\a _number_member: new position; = \a _number_member: delete
         */
        static void createRemoveTable(const int* _index, const int _n_index, const int _n_member, int* _remove_table) {
            ASSERT(_n_index>0);
            ASSERT(_n_member>=_n_index);
            
            // last member index 
            int ilast = _n_member-1;
            const int n_new = _n_member - _n_index; // new ptcl number
            // initial table
            for (int i=0; i<_n_member; i++) _remove_table[i] = -1;

            for (int i=0; i<_n_index; i++) {
                int k=_index[i];
                
                int idel = k;
                
                // set ilast >= idel for safety.
                ilast = std::max(ilast, idel);
            
                // check whether position is already moved, if so, check the moved new position and set current to delete
                if(_remove_table[k]>=0) {
                    idel=_remove_table[k];
                    _remove_table[k]=_n_member;
                }

                ASSERT(k<_n_member);
                ASSERT(idel<_n_member);
                
                // check the last avaiable member that can be moved to idel
                // If the ilast is already moved, check whether the new moved position _remove_table[ilast] is before the current idel.
                // If _remove_table[ilast] is after the current idel, update the _remove_table[ilast] to idel and set _remove_table[ilast] as new idel to check, until _remove_table[ilast]=-1 or idel >= ilast
                while (idel>=0) {
                    int itrlast = -1;
                    //cond: last is moved && ( moved pos before new n_ptcl || last is del ) &&  ilast > idel 
                    while(_remove_table[ilast]>=0 && (_remove_table[ilast]<n_new || _remove_table[ilast]==_n_member) && ilast>idel) ilast--;

                    // if ilast is not yet moved or (new pos of ilast > idel and ilast is not del, move ilast to idel)
                    if(_remove_table[ilast]<0||(idel<_remove_table[ilast]&&_remove_table[ilast]<_n_member)) {
                        itrlast=_remove_table[ilast];
                        _remove_table[ilast]=idel;
                    }
                    // idel is already at last, remove it
                    _remove_table[idel]=_n_member; 
                    idel = itrlast;
                }
            }
        }
    
        //! remove one member
        /*! remove one member from local member array
          @param [in] _index: member index to be removed
          @param [in] _shift_last_only_flag: true: shift last member to current position (defaulted); false: shift all right member to left by one
        */
        void removeMember(const int _index, const bool _shift_last_only_flag) {
            ASSERT(_index>=0);
            ASSERT(_index<num_);

            int ilast = num_-1;
            if (_index<ilast) {
                if (_shift_last_only_flag) {
                    data_[_index] = data_[ilast];
                    if(mode_==ListMode::copy) adr_[_index] = adr_[ilast];
                }
                else {
                    if(mode_==ListMode::copy) {
                        for (int j=_index; j<ilast; j++) {
                            data_[j] = data_[j+1];
                            adr_[j]  = adr_[j+1];
                        }
                    }
                    else {
                        for (int j=_index; j<ilast; j++) data_[j] = data_[j+1];
                    }
                }
            }
            num_--;
            modified_flag_ = true;
        }


        //! remove a list of member based on remove table
        /*!
          @param[in] remove_table: table map the member new position and deleted ones, generated from #createRemoveTable
         */
        void removeMemberTable(const int* remove_table) {
            const int num_org = num_;
            if (mode_==ListMode::copy) {
                for (int i=0; i<num_org; i++) {
                    ASSERT(remove_table[i]<=num_org);
                    // no change case
                    if(remove_table[i]<0) continue;
                    // move case
                    else if(remove_table[i]<num_org) {
                        const int inew = remove_table[i];
                        data_[inew] = data_[i];
                        adr_[inew] = adr_[i];
                    }
                    // remove count
                    else num_--;
                }
            }
            else {
                for (int i=0; i<num_org; i++) {
                    ASSERT(remove_table[i]<=num_org);
                    // no change case
                    if(remove_table[i]<0) continue;
                    // move case
                    else if(remove_table[i]<num_org) {
                        const int inew = remove_table[i];
                        data_[inew] = data_[i];
                    }
                    // remove count
                    else num_--;
                }
            }
            modified_flag_ = true;
        }

        //! remove a list of member based on an index list
        /*!
          @param[in] _index: member index list for remove
          @param[in] _n_index: number of removing member indices
         */
        void removeMemberList(const int* _index, const int _n_index) {
            ASSERT(_n_index<=num_);
            const int num_org = num_;
            int remove_table[num_];
            createRemoveTable(_index, _n_index, num_, remove_table);
            removeMemberTable(remove_table);
            ASSERT(num_+_n_index==num_org);
        }

        //! copy all member data back to original address
        /*! work for #ListMode::copy
         */
        template <class T>
        void writeBackMemberAll() {
            ASSERT(mode_==ListMode::copy);
            ASSERT(nmax_>0);
            for (int i=0; i<num_; i++) {
                if(adr_[i]!=NULL) *(T*)adr_[i] = data_[i];
            }
        }

        //! copy a list of member data back to original address
        /*! work for #ListMode::copy
          @param[in] _index: list of member index
          @param[in] _n:  number of members
         */
        template <class T>
        void writeBackMemberList(const int* _index, const int _n) {
            ASSERT(mode_==ListMode::copy);
            ASSERT(nmax_>0);
            for (int i=0; i<_n; i++) {
                int k=_index[i];
                if(adr_[k]!=NULL) *(T*)adr_[k] = data_[k];
            }
        }

        //! Get modified status
        bool isModified() const {
            return modified_flag_;
        }

        //! Reset modified status to false
        void setModifiedFalse() {
            modified_flag_ = false;
        }

    };
}
