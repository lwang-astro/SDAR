#pragma once

#include "Common/profile.h"

namespace AR{
    //! profiling class for AR integrator
    class Profile {
    public:
        typedef long long unsigned int UInt64;
        UInt64 step_count_sum; // number of integration steps summation
        UInt64 step_count_tsyn_sum; // number of integration steps during time synchronization summation
        UInt64 step_count; // number of integration steps from last step 
        UInt64 step_count_tsyn; // number of integration steps during time synchronization from last step
#ifdef SDAR_TIME_MEASURE
        COMM::TimeMeasure prof_tot; // total wallclock time
        COMM::TimeMeasure prof_int; // only integration wallclock time
        COMM::TimeMeasure prof_int_tsyn; // only integration wallclock time for time synchronization
#endif

        // constructor
        Profile(): step_count_sum(0), step_count_tsyn_sum(0), step_count(0), step_count_tsyn(0)
#ifdef SDAR_TIME_MEASURE
                  , prof_tot(), prof_int(), prof_int_tsyn() 
#endif
                  {}

        // clear function
        void clear() {
            step_count = step_count_tsyn = 0;
            step_count_sum = step_count_tsyn_sum = 0;
#ifdef SDAR_TIME_MEASURE
            prof_int.time = prof_int_tsyn.time = prof_tot.time = 0.0; 
#endif
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"Nstep(sum)"
                 <<std::setw(_width)<<"Nstep_tsyn(sum)"
                 <<std::setw(_width)<<"Nstep"
                 <<std::setw(_width)<<"Nstep_tsyn";
#ifdef SDAR_TIME_MEASURE
            _fout<<std::setw(_width)<<"Total(s)"
                 <<std::setw(_width)<<"Int(s)"
                 <<std::setw(_width)<<"Int_tsyn(s)";
#endif
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<step_count_sum
                 <<std::setw(_width)<<step_count_tsyn_sum
                 <<std::setw(_width)<<step_count
                 <<std::setw(_width)<<step_count_tsyn;
#ifdef SDAR_TIME_MEASURE
            _fout<<std::setw(_width)<<prof_tot.time
                 <<std::setw(_width)<<prof_int.time
                 <<std::setw(_width)<<prof_int_tsyn.time;
#endif
        }

        //! write class data with BINARY format
        /*! @param[in] _fout: file IO for write
         */
        void writeBinary(FILE *_fout) {
            fwrite(this, sizeof(*this),1,_fout);
        }

        //! read class data with BINARY format and initial the array
        /*! @param[in] _fin: file IO for read
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(this, sizeof(*this), 1, _fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
    };
}
