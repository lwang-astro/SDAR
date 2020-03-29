#pragma once

namespace AR {

    //! binary interrupt status
    enum class InterruptStatus{change=1, merge=2, none=0};

    //! Binary interrupt recoreder
    template <class Tparticle>
    struct InterruptBinary {
        COMM::BinaryTree<Tparticle>* adr;
        Float time_now;
        Float time_end;
        InterruptStatus status;

        InterruptBinary(): adr(NULL), time_now(0.0), time_end(0.0), status(InterruptStatus::none) {}

        InterruptBinary(COMM::BinaryTree<Tparticle>* _adr, Float _time_now, Float _time_end, InterruptStatus _status): adr(_adr), time_now(_time_now), time_end(_time_end), status(_status) {}

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
