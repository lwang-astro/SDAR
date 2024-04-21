#pragma once

#include "AR/information.h"

namespace AR {

    //! binary interrupt status
    enum class InterruptStatus{change=1, merge=2, destroy=3, none=0};

    //! Binary interrupt recoreder
    template <class Tparticle>
    struct InterruptBinary {
    private:
        bool is_local_backup; // indicator whether adr has local memory allocation
        AR::BinaryTree<Tparticle>* adr; // binary tree address
    public:
        Float time_now;  // current time
        Float time_end;  // finishing time 
        InterruptStatus status; // binary status

        InterruptBinary(): is_local_backup(false), adr(NULL), time_now(0.0), time_end(0.0), status(InterruptStatus::none) {}

        //InterruptBinary(AR::BinaryTree<Tparticle>* _adr, Tparticle* _p, const Float& _time_now, const Float& _time_end, const InterruptStatus& _status): adr(_adr), time_now(_time_now), time_end(_time_end), status(_status) {}

        //InterruptBinary(InterruptBinary<Tparticle>& _bin) {
        //    adr = _bin.adr;
        //    time_now= _bin.time_now;
        //    time_end= _bin.time_end;
        //    status = _bin.status;
        //}

        //! set binary tree address
        /*! 
          @param[in] clear_backup: if true and there is local backup of binarytree, clear up the backup first
         */
        void setBinaryTreeAddress(AR::BinaryTree<Tparticle>* _adr) {
            if (is_local_backup) {
                delete adr;
                is_local_backup = false;
            }
            ASSERT(!is_local_backup);
            adr = _adr;
        }

        //! get binary tree address
        AR::BinaryTree<Tparticle>* getBinaryTreeAddress() const {
            return adr;
        }

        //! allocate binary tree memory locally and backup the binarytree information
        /*! Only support a binary, no high-order multiple system
         */
        void backupBinaryTreeLocal() {
            ASSERT(adr->getMemberN()==2);
            ASSERT(!is_local_backup);
            ASSERT(adr!=NULL);
            auto backup_adr = adr;
            adr = new AR::BinaryTree<Tparticle>;
            *adr = *backup_adr;
            if (!backup_adr->isAllocatedMembers()) {
                adr->allocateParticleMember();
                for (int i=0; i<2; i++) {
                    auto pi = adr->getMember(i);
                    *pi = *(backup_adr->getMember(i));
                }
            }
            is_local_backup = true;
        }

        //! operator = 
        /*! Copy function will remove the local allocated data and backup copied data.
         */
        InterruptBinary & operator = (const InterruptBinary& _bin) {
            clear();
            adr = _bin.adr;
            time_now = _bin.time_now;
            time_end = _bin.time_end;
            status = _bin.status;
            if (_bin.is_local_backup) backupBinaryTreeLocal();
            return *this;
        }

        void clear() { 
            if (is_local_backup) {
                delete adr;
                is_local_backup = false;
            }
            adr=NULL; 
            time_now=0.0; 
            time_end=0.0; 
            status=InterruptStatus::none;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted: 20)
          @param[in] print_member_flag: if true, print two members (defaulted: false)
        */
        static void printColumnTitle(std::ostream & _fout, const int _width=20, const bool print_member_flag=false) {
            _fout<<std::setw(_width)<<"time_now"
                 <<std::setw(_width)<<"time_end"
                 <<std::setw(_width)<<"status";
            AR::BinaryTree<Tparticle>::printColumnTitle(_fout, _width);
            if (print_member_flag) {
                Tparticle::printColumnTitle(_fout, _width);
                Tparticle::printColumnTitle(_fout, _width);
            }
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
          @param[in] print_member_flag: if true, print two members (defaulted: false)
        */
        void printColumn(std::ostream & _fout, const int _width=20, const bool print_member_flag=false) {
            _fout<<std::setw(_width)<<time_now
                 <<std::setw(_width)<<time_end
                 <<std::setw(_width)<<int(status);
            adr->printColumn(_fout, _width);
            if (print_member_flag) {
                for (int j=0; j<2; j++) 
                    adr->getMember(j)->printColumn(_fout, _width);
            }
        }

        bool checkParams() {
            ASSERT((status!=InterruptStatus::none&&adr!=NULL)||status==InterruptStatus::none);
            return true;
        }
    };
}
