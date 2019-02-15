#pragma once

#include "Common/Float.h"
#include <limits>

namespace H4 {
    
//! Neighbor information collector
    class Neighbor{
    public:
        int r_min_index;   // nearest neighbor index for each ptcl
        int mass_min_index;    // mimimum mass in neighbors
        Float r_min_sq;    // nearest neighbor distance square
        Float r_min_mass;   // nearest neighbor index for each ptcl
        Float mass_min;    // mimimum mass in neighbors
        Float r_crit_sq;   // neighbor radius criterion
        bool initial_step_flag; // indicate whether the time step need to be initialized due to the change of neighbors
        COMM::List<int> single_list; // neighbor list index of singles
        COMM::List<int> group_list;  // neighbor list index of groups

        //! constructor
        Neighbor(): r_min_index(-1), mass_min_index(-1), r_min_sq(NUMERIC_FLOAT_MAX), mass_min(NUMERIC_FLOAT_MAX), r_crit_sq(0.0), initial_step_flag(false), nsingle_list(), group_list() {}

        //! reserve memory for neighbor lists
        /*! 
          @param[in] _nmax_single: maximum number of neighbors of singles
          @param[in] _nmax_group: maximum number of neighbors of groups
        */
        void reserveMem(const int _nmax_single, const int _nmax_group) {
            single_list.setMode(ListMode::local);
            single_list.reserveMem(_nmax_single);
            group_list.setMode(ListMode::local);
            group_list.reserveMem(_nmax_group);
        }

        //! clear function
        void clear() {
            single_list.clear();
            group_list.clear();
        }

        //! replace one index for nearest particle and neighbor lists
        /*!
          @parma[in] _i_org: original index to be replaced
          @param[in] _i_new: new index
          @param[in] _nb_modify_flag: if true, also modify the neighbor list index
        */
        void replaceIndex(const int _i_org, const int _i_new, const bool nb_modify_flag) {
            if (r_min_index_==_i_org) r_min_index_ = _i_new;
            if (nb_modify_flag) {
                const int n = nb_list_.size();
                for (int j=0; j<n; j++) {
                    if(nb_list_[j]==_i_org) nb_list_[j]=_i_new;
                }
            }
        }
    
    };
}
