#pragma once

#include <sys/time.h>

namespace COMM {
    //! Profile class to measure the performance
    struct TimeMeasure{
        double time;

        TimeMeasure(): time(0.0) {}

        // time measure function
        static double get_wtime() {
            struct timespec ts;
            clock_gettime(CLOCK_MONOTONIC, &ts);
            return ts.tv_sec + 1.e-9 * ts.tv_nsec;
        }

        // time measure start
        void start() {
            time -= get_wtime();
        }

        // time measure end
        void end() {
            time += get_wtime();
        }
    };
}
