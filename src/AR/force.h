#pragma once

#include "Common/Float.h"

namespace AR {
    //! class to store the acceleration, perturbation and time transformation function gradient for AR integration
    struct Force {
    public:
        Float acc_in[3];      // total acceleration for one particle
        Float acc_pert[3]; // perturbation 
#ifdef AR_TTL
        Float gtgrad[3];   // time transformation function gradient

        Force(): acc_in{Float(0.0),Float(0.0),Float(0.0)}, 
                 acc_pert{Float(0.0),Float(0.0),Float(0.0)},
                 gtgrad{Float(0.0),Float(0.0),Float(0.0)} {}
#else
        Force(): acc_in{Float(0.0),Float(0.0),Float(0.0)}, 
                 acc_pert{Float(0.0),Float(0.0),Float(0.0)} {}
#endif

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
