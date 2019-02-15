#pragma once

#include "Hermite/hermite_particle.h"

namespace H4{
    //! AR information used for Hermite integrator
    template <class Tparticle, class Tarint>
    class ARPerturber{
    public:
        PtclH4<Tparticle>* single_particles; // single particle data address
        Tarint* group_particles; // group data address
        Neighbor neighbors;  // neighbor information
        bool need_resolve_flag; // indicate whether the members need to be resolved for outside

        ARPerturber(): single_particles(NULL), group_particles(NULL), neighbors(), need_resolve_flag(false) {}

        //! clear function
        void clear() {
            single_particles = NULL;
            group_particles = NULL;
            neigbors.clear();
            need_resolve_flag = false;
        }
    };
};
