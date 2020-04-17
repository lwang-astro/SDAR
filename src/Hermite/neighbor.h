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
        Float r_neighbor_crit_sq; // neighbor radius criterion
        bool need_resolve_flag; // indicate whether the members need to be resolved for outside
        bool initial_step_flag; // indicate whether the time step need to be initialized due to the change of neighbors
        int n_neighbor_group; // number of group neighbor
        int n_neighbor_single; // number of single neighbor
        COMM::List<NBAdr<Tparticle>> neighbor_address; // neighbor perturber address

        //! constructor
        Neighbor(): r_min_index(-1), mass_min_index(-1), r_min_sq(NUMERIC_FLOAT_MAX), mass_min(NUMERIC_FLOAT_MAX), r_neighbor_crit_sq(-1.0), need_resolve_flag(false), initial_step_flag(false), n_neighbor_group(0), n_neighbor_single(0), neighbor_address() {}

        //! check whether parameters values are correct
        /*! \return true: all correct
         */
        bool checkParams() {
            ASSERT(r_neighbor_crit_sq>0.0);
            return true;
        }        


        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"r_min_index"
                 <<std::setw(_width)<<"mass_min_index"
                 <<std::setw(_width)<<"r_min_sq"
                 <<std::setw(_width)<<"r_min_mass"
                 <<std::setw(_width)<<"mass_min"
                 <<std::setw(_width)<<"r_neighbor_crit_sq"
                 <<std::setw(_width)<<"need_resolve_flag"
                 <<std::setw(_width)<<"initial_step_flag"
                 <<std::setw(_width)<<"n_neighbor_group"
                 <<std::setw(_width)<<"n_neighbor_single";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<r_min_index
                 <<std::setw(_width)<<mass_min_index
                 <<std::setw(_width)<<r_min_sq
                 <<std::setw(_width)<<r_min_mass
                 <<std::setw(_width)<<mass_min
                 <<std::setw(_width)<<r_neighbor_crit_sq
                 <<std::setw(_width)<<need_resolve_flag
                 <<std::setw(_width)<<initial_step_flag
                 <<std::setw(_width)<<n_neighbor_group
                 <<std::setw(_width)<<n_neighbor_single;
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
            need_resolve_flag = false;
            initial_step_flag = false;
        }

        //! reset neighbor information
        void resetNeighbor() {
            r_min_index = -1;
            mass_min_index = -1;
            r_min_sq = NUMERIC_FLOAT_MAX;
            mass_min = NUMERIC_FLOAT_MAX;
            n_neighbor_group = 0;
            n_neighbor_single = 0;
            neighbor_address.resizeNoInitialize(0);            
        }

        //! clear function
        void clear() {
            clearNoFreeMemNoResizeNeighborAdress();
            r_neighbor_crit_sq = -1.0;
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

        //! check and add neighbor of single
        template <class Tp>
        void checkAndAddNeighborSingle(const Float _r2, Tp& _particle, const Neighbor<Tparticle>& _nbp, const int _index) {
            if (_r2<std::max(r_neighbor_crit_sq,_nbp.r_neighbor_crit_sq)) {
                neighbor_address.addMember(NBAdr<Tparticle>(&_particle, _index));
                n_neighbor_single++;
            }
            // mass weighted nearest neigbor
            Float mass = _particle.mass;
            if (_r2*r_min_mass < r_min_sq*mass) {
                r_min_sq = _r2;
                r_min_index = _index;
                r_min_mass = mass;
            }
            // minimum mass
            if(mass_min>mass) {
                mass_min = mass;
                mass_min_index = _index;
            }
            // if neighbor step need update, also set update flag for i particle
            if(_nbp.initial_step_flag) initial_step_flag = true;
        }

        //! check and add neighbor of group
        template <class Tgroup>
        void checkAndAddNeighborGroup(const Float _r2, Tgroup& _group, const int _index) {
            ASSERT(_r2>0.0);
            if (_r2<std::max(r_neighbor_crit_sq, _group.perturber.r_neighbor_crit_sq)) {
                neighbor_address.addMember(NBAdr<Tparticle>(&_group.particles,_index));
                n_neighbor_group++;
            }
            // mass weighted nearest neigbor
            Float mass_cm = _group.particles.cm.mass;
            if (_r2*r_min_mass < r_min_sq*mass_cm) {
                r_min_sq = _r2;
                r_min_index = _index;
                r_min_mass = mass_cm;
            }
            // minimum mass
            if(mass_min>mass_cm) {
                mass_min = mass_cm;
                mass_min_index = _index;
            }
            // if neighbor step need update, also set update flag for i particle
            if(_group.perturber.initial_step_flag) initial_step_flag = true;
        }

        //! check whether members should be resolved for outside
        /*! @param[in] _kappa: slowdown factor
         */
        void checkGroupResolve() {
            // always resolve
            if (!need_resolve_flag) {
                need_resolve_flag = true;
                initial_step_flag = true;
            }
            // _kappa>3.0 seems cause oscillation 
            /* suppress unresolved flag 
            if (_kappa>11.0 && need_resolve_flag) {
                need_resolve_flag = false;
                initial_step_flag = true;
            }
            */
        }

    };
}
