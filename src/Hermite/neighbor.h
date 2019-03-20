#pragma once

#include "Common/Float.h"
#include "Common/particle_group.h"
#include "Hermite/hermite_particle.h"
#include <limits>

namespace H4 {
    //! type of neighbor address
    enum class NBType {single, group, none};

    // neighbor address 
    template <class Tparticle>
    struct NBAdr{ 
        typedef ParticleH4<Tparticle> Single;
        typedef COMM::ParticleGroup<ParticleAR<Tparticle>, ParticleH4<Tparticle>> Group;
        // adr;
        void* adr;
        int index;
        NBType type; // -1: undefine; 0: H4particle; 1: ParticleGroup
            
        NBAdr(): adr(NULL), index(-1), type(NBType::none) {}

        NBAdr(Single* _adr, const int _index): adr((void*)_adr), index(_index), type(NBType::single) {}

        NBAdr(Group* _adr, const int _index): adr((void*)_adr), index(_index), type(NBType::group) {}

        NBAdr& operator = (const NBAdr& _nb) {
            adr = _nb.adr;
            index =_nb.index;
            type = _nb.type;
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
        bool need_resolve_flag; // for group
        bool initial_step_flag; // indicate whether the time step need to be initialized due to the change of neighbors
        int n_neighbor_group; // number of group neighbor
        int n_neighbor_single; // number of single neighbor
        COMM::List<NBAdr<Tparticle>> neighbor_address; // neighbor perturber address

        //! constructor
        Neighbor(): r_min_index(-1), mass_min_index(-1), r_min_sq(NUMERIC_FLOAT_MAX), mass_min(NUMERIC_FLOAT_MAX), r_crit_sq(-1.0), need_resolve_flag(false), initial_step_flag(false), n_neighbor_group(0), n_neighbor_single(0), neighbor_address() {}

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            ASSERT(r_crit_sq>0.0);
            return true;
        }        

        //! reserve memory for neighbor lists
        /*! 
          @param[in] _nmax: maximum number of neighbors
        */
        void reserveMem(const int _nmax) {
            neighbor_address.setMode(COMM::ListMode::local);
            neighbor_address.reserveMem(_nmax);
        }

        //! clear function
        void clearNoFreeMemNoResizeNeighborAdress() {
            r_min_index = -1;
            mass_min_index = -1;
            r_min_sq = NUMERIC_FLOAT_MAX;
            mass_min = NUMERIC_FLOAT_MAX;
            r_crit_sq = -1.0;
            need_resolve_flag = false;
            initial_step_flag = false;
        }

        //! clear function
        void clear() {
            clearNoFreeMemNoResizeNeighborAdress();
            n_neighbor_group = 0;
            n_neighbor_single = 0;
            neighbor_address.clear();
        }

        //! clear no release memory
        void clearNoFreeMem() {
            clearNoFreeMemNoResizeNeighborAdress();
            n_neighbor_group = 0;
            n_neighbor_single = 0;
            neighbor_address.resizeNoInitialize(0);            
        }

        //! check whether members should be resolved for outside
        /*! @param[in] _kappa: slowdown factor
         */
        void checkGroupResolve(const Float _kappa) {
            if (_kappa==1.0 && !need_resolve_flag) {
                need_resolve_flag = true;
                initial_step_flag = true;
            }
            // _kappa>3.0 seems cause oscillation 
            if (_kappa>10.0 && need_resolve_flag) {
                need_resolve_flag = false;
                initial_step_flag = true;
            }
        }

    };
}
