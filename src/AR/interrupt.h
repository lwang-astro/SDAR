#pragma once

#include "AR/information.h"

namespace AR {

    //! binary interrupt status
    enum class InterruptStatus{change=1, merge=2, destroy=3, none=0};

    //! Binary interrupt recoreder
    template <class Tparticle>
    struct InterruptBinary {
        AR::BinaryTree<Tparticle>* adr; // binary tree address
        Float time_now;  // current time
        Float time_end;  // finishing time 
        InterruptStatus status; // binary status

        InterruptBinary(): adr(NULL), time_now(0.0), time_end(0.0), status(InterruptStatus::none) {}

        //InterruptBinary(AR::BinaryTree<Tparticle>* _adr, Tparticle* _p, const Float& _time_now, const Float& _time_end, const InterruptStatus& _status): adr(_adr), time_now(_time_now), time_end(_time_end), status(_status) {}

        //InterruptBinary(InterruptBinary<Tparticle>& _bin) {
        //    adr = _bin.adr;
        //    time_now= _bin.time_now;
        //    time_end= _bin.time_end;
        //    status = _bin.status;
        //}

        void clear() { 
            adr=NULL; 
            time_now=0.0; 
            time_end=0.0; 
            status=InterruptStatus::none;
        }

        bool checkParams() {
            ASSERT((status!=InterruptStatus::none&&adr!=NULL)||status==InterruptStatus::none);
            return true;
        }
    };
}
