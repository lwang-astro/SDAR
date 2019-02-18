#pragma once

#include "Common/Float.h"
#include "Hermite/hermite_particle.h"
#include <limits>

namespace H4 {
    // neighbor address 
    template <class Tparticle>
    struct NBAdr{ 
        ParticleH4<Tparticle>* adr;
        int index;
            
        NBAdr(): adr(NULL), index(-1) {}
        NBAdr(ParticleH4<Tparticle>* _adr, const int _index): adr(_adr), index(_index) {}

        NBAdr& operator = (const NBAdr& _nb) {
            adr = _nb.adr;
            index =_nb.index;
            return *this;
        }

    };

    
//! Neighbor information collector
    template <class Tparticle>
    class Neighbor{
    public:

        int r_min_index;   // nearest neighbor index for each ptcl
        int mass_min_index;    // mimimum mass in neighbors
        Float r_min_sq;    // nearest neighbor distance square
        Float r_min_mass;   // nearest neighbor index for each ptcl
        Float mass_min;    // mimimum mass in neighbors
        Float r_crit_sq;   // neighbor radius criterion
        bool initial_step_flag; // indicate whether the time step need to be initialized due to the change of neighbors
        COMM::List<NBAdr<Tparticle>> neighbor_address; // neighbor perturber address

        //! constructor
        Neighbor(): r_min_index(-1), mass_min_index(-1), r_min_sq(NUMERIC_FLOAT_MAX), mass_min(NUMERIC_FLOAT_MAX), r_crit_sq(0.0), initial_step_flag(false), neighbor_address() {}

        //! reserve memory for neighbor lists
        /*! 
          @param[in] _nmax: maximum number of neighbors
        */
        void reserveMem(const int _nmax) {
            neighbor_address.setMode(COMM::ListMode::local);
            neighbor_address.reserveMem(_nmax);
        }

        //! clear function
        void clear() {
            neighbor_address.clear();
        }
    };
}
