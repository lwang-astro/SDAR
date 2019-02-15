#pragma once

#include "Hermite/hermite_particle.h"

namespace H4{
    //! AR information used for Hermite integrator
    template <class Tparticle>
    class ARPerturber{
    public:
        Neighbor neighbors<Tparticle>;  // neighbor information
        bool need_resolve_flag; // indicate whether the members need to be resolved for outside

        ARPerturber(): neighbors(), need_resolve_flag(false) {}

        //! clear function
        void clear() {
            neigbors.clear();
            need_resolve_flag = false;
        }
    };
};
