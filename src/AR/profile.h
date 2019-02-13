#pragma once

#include <sys/time.h>

namespace AR{
    //! Profile class to measure the performance
    struct TimeMeasure{
        // time measure function
        static double get_wtime(){
            struct timeval tv;
            gettimeofday(&tv, NULL);
            return tv.tv_sec + 1.e-6 * tv.tv_usec;
        }
        double time;

        // time measure start
        void start() {
            time -= get_wtime();
        }

        // time measure end
        void end() {
            time += get_wtime();
        }
    };

    //! profiling class for AR integrator
    class Profile {
    public:
        typedef long long unsigned int UInt64;
        UInt64 step_count; // number of integration steps
        UInt64 step_count_tsyn; // number of integration steps during time synchronization

        // constructor
        Profile(): step_count(0), step_count_tsyn(0) {}

        // clear function
        void clear() {
            step_count = step_count_tsyn = 0;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"Nstep"
                 <<std::setw(_width)<<"Nstep_tsyn";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<step_count
                 <<std::setw(_width)<<step_count_tsyn;
        }
    };
}
