#pragma once

#include "AR/information.h"

namespace AR {

    //! binary interrupt status
    enum class InterruptStatus{change=1, merge=2, none=0};

    //! Binary interrupt recoreder
    template <class Tparticle>
    struct InterruptBinary {
        AR::BinaryTree<Tparticle>* adr; // binary tree address
        Float dm; // mass change
        Float time_now;  // current time
        Float time_end;  // finishing time 
        InterruptStatus status; // binary status

        InterruptBinary(): adr(NULL), time_now(0.0), time_end(0.0), status(InterruptStatus::none) {}

        InterruptBinary(AR::BinaryTree<Tparticle>* _adr, Float _time_now, Float _time_end, InterruptStatus _status): adr(_adr), time_now(_time_now), time_end(_time_end), status(_status) {}

        void clear() { 
            adr=NULL; 
            dm =0.0;
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
