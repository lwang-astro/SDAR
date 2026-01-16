#pragma once

#include "Common/profile.h"

namespace H4{
    class Profile{
    public:
        typedef long long unsigned int UInt64;
        UInt64 hermite_single_step_count; // number of integration steps of hermite single
        UInt64 hermite_group_step_count; // number of integration steps of hermite groups
        UInt64 hermite_single_interact_count; // number of interactions of hermite single
        UInt64 hermite_group_interact_count; // number of interactions of hermite groups
        UInt64 ar_step_count; // number of integration steps of ar
        UInt64 ar_step_count_tsyn; // number of integration steps of ar
        UInt64 break_group_count; // times of break groups
        UInt64 new_group_count; // times of new groups
        UInt64 merge_group_count; // times of merged groups
#ifdef SDAR_TIME_MEASURE        
        COMM::TimeMeasure prof_tot; // time measure of total time
        COMM::TimeMeasure prof_hermite_single; // time measure of hermite single
        COMM::TimeMeasure prof_hermite_group; // time measure of hermite group
        COMM::TimeMeasure prof_adjust; // time measure of adjust groups
        COMM::TimeMeasure prof_init; // time measure of initialization
        COMM::TimeMeasure prof_modify_single; // time measure of modify single
        COMM::TimeMeasure prof_select_act; // time measure of select active particles
        COMM::TimeMeasure prof_ar; // time measure of ar integration in groups
#endif

        Profile() {clear();} 
    
        void clear() {
            hermite_single_step_count = hermite_group_step_count = 0;
            hermite_single_interact_count = hermite_group_interact_count = 0;
            ar_step_count = ar_step_count_tsyn = 0;
            break_group_count = 0;
            new_group_count = 0;
            merge_group_count = 0;
#ifdef SDAR_TIME_MEASURE
            prof_tot.time = 0.0;
            prof_hermite_single.time = 0.0;
            prof_hermite_group.time = 0.0;
            prof_adjust.time = 0.0;
            prof_init.time = 0.0;
            prof_modify_single.time = 0.0;
            prof_select_act.time = 0.0;
            prof_ar.time = 0.0;
#endif
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"H4_step_single"
                 <<std::setw(_width)<<"H4_step_group"
                 <<std::setw(_width)<<"H4_force_single"
                 <<std::setw(_width)<<"H4_force_group"
                 <<std::setw(_width)<<"AR_step"
                 <<std::setw(_width)<<"AR_step_tsyn"
                 <<std::setw(_width)<<"break_group"
                 <<std::setw(_width)<<"new_group"
                 <<std::setw(_width)<<"merge_group";
#ifdef SDAR_TIME_MEASURE
            _fout<<std::setw(_width)<<"prof_tot[s]"
                 <<std::setw(_width)<<"prof_H4_single[s]"
                 <<std::setw(_width)<<"prof_H4_group[s]"
                 <<std::setw(_width)<<"prof_adjust[s]"
                 <<std::setw(_width)<<"prof_init[s]"
                 <<std::setw(_width)<<"prof_modify[s]"
                 <<std::setw(_width)<<"prof_select[s]"
                 <<std::setw(_width)<<"prof_AR[s]";
#endif
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<hermite_single_step_count
                 <<std::setw(_width)<<hermite_group_step_count
                 <<std::setw(_width)<<hermite_single_interact_count
                 <<std::setw(_width)<<hermite_group_interact_count
                 <<std::setw(_width)<<ar_step_count
                 <<std::setw(_width)<<ar_step_count_tsyn
                 <<std::setw(_width)<<break_group_count
                 <<std::setw(_width)<<new_group_count
                 <<std::setw(_width)<<merge_group_count;
#ifdef SDAR_TIME_MEASURE
            _fout<<std::setw(_width)<<prof_tot.time
                 <<std::setw(_width)<<prof_hermite_single.time
                 <<std::setw(_width)<<prof_hermite_group.time
                 <<std::setw(_width)<<prof_adjust.time
                 <<std::setw(_width)<<prof_init.time
                 <<std::setw(_width)<<prof_modify_single.time
                 <<std::setw(_width)<<prof_select_act.time
                 <<std::setw(_width)<<prof_ar.time;
#endif
        }

    };
}
