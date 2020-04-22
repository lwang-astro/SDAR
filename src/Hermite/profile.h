#pragma once

namespace H4{
    class Profile{
    public:
        typedef long long unsigned int UInt64;
        UInt64 hermite_single_step_count; // number of integration steps of hermite single
        UInt64 hermite_group_step_count; // number of integration steps of hermite groups
        UInt64 ar_step_count; // number of integration steps of ar
        UInt64 ar_step_count_tsyn; // number of integration steps of ar
        UInt64 break_group_count; // times of break groups
        UInt64 new_group_count; // times of new groups

        Profile() {clear();} 
    
        void clear() {
            hermite_single_step_count = hermite_group_step_count = 0;
            ar_step_count = ar_step_count_tsyn = 0;
            break_group_count = 0;
            new_group_count = 0;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"H4_step_single"
                 <<std::setw(_width)<<"H4_step_group"
                 <<std::setw(_width)<<"AR_step"
                 <<std::setw(_width)<<"AR_step_tsyn"
                 <<std::setw(_width)<<"break_group"
                 <<std::setw(_width)<<"new_group";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<hermite_single_step_count
                 <<std::setw(_width)<<hermite_group_step_count
                 <<std::setw(_width)<<ar_step_count
                 <<std::setw(_width)<<ar_step_count_tsyn
                 <<std::setw(_width)<<break_group_count
                 <<std::setw(_width)<<new_group_count;
        }

    };
}
