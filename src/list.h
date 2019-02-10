#pragma once

#include <cassert>
#include "Float.h"

namespace AR {

    //! list mode indicator
    /* copy: copy data and link original address (data_, adr_ allocated); \n
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
            assert(mode_==ListMode::none);
            mode_ = _mode;
        }
        
        //! Memory allocation for storing members
        /*! allocate memory for storing members and their original addresses.
          @param[in] _nmax: maximum number of members for memory allocation
         */
        void reserveMem(const int _nmax) {
            assert(nmax_==0);
            assert(_nmax>0);
            assert(adr_==NULL);
            assert(data_==NULL);

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
                assert(data_!=NULL);
                assert(adr_!=NULL);
                delete[] data_;
                delete[] adr_;
                data_= NULL;
                adr_ = NULL;
                nmax_ =0;
                num_ = 0;
                break;
            case ListMode::local:
                assert(data_!=NULL);
                assert(adr_==NULL);
                delete[] data_;
                data_=NULL;
                nmax_ =0;
                num_ = 0;
                break;
            case ListMode::link:
                assert(nmax_==0);
                assert(adr_==NULL);
                data_=NULL;
                num_ = 0;
                break;
            case ListMode::none:
                assert(data_==NULL);
                assert(adr_==NULL);
                assert(nmax_==0);
                assert(num_==0);
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
        /*! \return Current member number in member address list #p
         */
        int getSize() const {
            return num_;
        }
    
        //! return one member data reference
        /*!
          @param [in] _index: index of member in group
         */
        Ttype& getMember(const int _index) {
            assert(_index<num_);
            assert(_index>=0);
            return data_[_index];
        }

        //! return member data array address
        Ttype* getDataAddress() {
            return data_;
        }

        //! return one member data reference operator
        /*! overload []
          @param [in] _index: index of member in group
        */
        Ttype& operator [](const int _index) const {
            assert(_index<num_);
            assert(_index>=0);
            return data_[_index];
        }

        //! return one member original address 
        /*! Only work when mode is #ListMode::copy
          \return: one member original address array 
        */
        Ttype* getMemberOriginAddress(const int _index) const{
            assert(mode_==ListMode::copy);
            assert(_index<num_);
            assert(_index<nmax_);
            assert(_index>=0);
            return adr_[_index];
        }

        //! Get maximum member number allow to store 
        /*! \return allocated number of members in memory
         */
        int getSizeMax() const {
            return nmax_;
        }
    
        //! copy one member 
        /*! add a member by copying member into local allocated memory and storing the original memory address of it.
          Only work for #ListMode::copy
          @param [in] a: a member to store (push back at the end of local array)
        */
        void addMember(Ttype &_member) {
            assert(mode_==ListMode::copy);
            assert(num_<nmax_);
            data_[num_] = _member;
            adr_[num_]  = &_member;
            num_++;
            modified_flag_ = true;
        }

        //! increase size without initialization 
        /*! Work for #ListMode::local case. 
          Require the memory size is big enough to contain the new number of members.
         */
        void increaseSizeNoInitial(const int _n) {
            assert(mode_==ListMode::local);
            assert(num_+_n<=nmax_);
            num_ += _n;
            modified_flag_ = true;
        }
    
        //! link a member array 
        /*! point the local member data to the begining of the input \a _member, no data copy, any modification later will be directly on original member array. Work for #ListMode::link
          @param [in] _member: array of members
          @param [in] _n_member: number of members
        */
        void linkMemberList(Ttype _member[], const int _n_member) {
            assert(mode_==ListMode::link);
            assert(nmax_==0);
            assert(data_==NULL);
            data_= _member;
            num_ = _n_member;
            modified_flag_ = true;
        }

        //! create remove table for removing a list of members
        /*!
          @param[in] _index: member index list for remove
          @param[in] _n_list: number of removing members
          @param[in] _n_member: number of total members in the table
          @param[out] _remove_table: remove_table generated, the size should be at least number of members to avoid overflow. for index i in table, value <0: no remove; <\a _number_member: new position; = \a _number_member: delete
         */
        static void createRemoveTable(const int* _index, const int _n_index, const int _n_member, int* _remove_table) {
            assert(_n_index>0);
            assert(_n_member>=_n_index);
            
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

                assert(k<_n_member);
                assert(idel<_n_member);
                
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
            assert(_index>=0);
            assert(_index<num_);

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
          @param[in] _remove_table: table map the member new position and deleted ones, generated from #createRemoveTable
         */
        void removeMemberTable(const int* remove_table) {
            const int num_org = num_;
            if (mode_==ListMode::copy) {
                for (int i=0; i<num_org; i++) {
                    assert(remove_table[i]<=num_org);
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
                    assert(remove_table[i]<=num_org);
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
          @param[in] _n_list: number of removing members
         */
        void removeMemberList(const int* _index, const int _n_index) {
            assert(_n_index<=num_);
            const int num_org = num_;
            int remove_table[num_];
            createRemoveTable(_index, _n_index, num_, remove_table);
            removeMemberTable(remove_table);
            assert(num_+_n_index==num_org);
        }

        //! copy all member data back to original address
        /*! work for #ListMode::copy
         */
        void writeBackMemberAll() {
            assert(mode_==ListMode::copy);
            assert(nmax_>0);
            for (int i=0; i<num_; i++) {
                if(adr_[i]!=NULL) *adr_[i] = data_[i];
            }
        }

        //! copy a list of member data back to original address
        /*! work for #ListMode::copy
          @param[in] _index: list of member index
          @param[in] _n:  number of members
         */
        void writeBackMemberList(const int* _index, const int _n) {
            assert(mode_==ListMode::copy);
            assert(nmax_>0);
            for (int i=0; i<_n; i++) {
                int k=_index[i];
                if(adr_[k]!=NULL) *adr_[k] = data_[k];
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
