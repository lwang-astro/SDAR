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
    };
}
