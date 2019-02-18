#pragma once

#include <algorithm>
#include "Common/Float.h"

namespace AR {
    //! class to manager symplectic KD steps
    /*! The Yoshida (1990) type symplectic methods
     */
    class SymplecticStep {
    private:
        // cd pair type
        typedef std::array<Float,2> CDPair;
        //! Structure to store cumsum of c_k and the index of step k
        struct CumSumCkIndex{
            Float cck;  ///> cumsum of c_k
            int index;   ///> index of step k
        };

        // members
        int sym_order_;                       // symplectic order 
        int sym_type_;                        // symplectic method (Yoshida) solution 1 or 2
        int cd_pair_array_size_;              // cd_pair array size 
        CDPair* cd_pair_;                     // Ck Dk pair for KDKD... or DKDK...
        CumSumCkIndex* sorted_cumsum_ck_index_;  // increasing order of cumsum ck with index of cd_pair

        //! Recursive Symplectic cofficients generator 
        /*!
          With Halmitionian Splitting \f$ H = T(p) + U(q) \f$
          The corresponding Symplectic mapping with drift (D) for q and Kick (K) for p is
          \f$ \dot z = \prod_{i=0}^k (exp[c_{k} t D] exp[c_{k+1} t K] exp[c_{k+2} t D]) z \f$
          where \f$ z=(p,q), k = 3^(n-1) \f$
      
          @param[in] coff_array: array to store cofficients \f$ c_k \f$ for steps, size of n^3
          @param[in] n: half of the order of integrator (order is 2n)
        */
        void recursiveSymplecticCofficientsSplit(Float* _ck, const int _n) {
            if (_n==1) {
                _ck[0] = 0.5;
                _ck[1] = 1.0;
                _ck[2] = 0.5;
            }
            else if (_n>1) {
                // order 2(n-1) 
                recursiveSymplecticCofficientsSplit(_ck, _n-1);

                // last order array size
                const int n_size=std::pow(3, _n-1);

                // backup last order cofficients
                Float ck_last[n_size];
                for (int i=0; i<n_size; i++) ck_last[i] = _ck[i];

                // 2^(1/(2n-1))
                Float ap = pow(2.0, 1.0/(2*_n-1));

                Float w[3];
                w[0] = 1.0/(2.0-ap);
                w[1] = - ap*w[0];
                w[2] = w[0];

                for (int i=0; i<3; i++) {
                    for (int j=0; j<n_size; j++) {
                        _ck[i*n_size+j] = w[i]*ck_last[j];
                    }
                }
            }
            else {
                std::cerr<<"Error: integrator order should be positive, given value "<<_n<<std::endl;
                abort();
            }
        }

        //! Symplectic cofficients generator
        /*!
          With Halmitionian Splitting \f$ H = T(p) + U(q) \f$
          The corresponding Symplectic mapping with drift (D) for q and Kick (K) for p is
          \f$ \dot z = \prod_{i=0}^k (exp[c_{k} t D] exp[d_{k} t K]) z \f$
          where \f$ z=(p,q), k = 3^(n-1)+1 \f$
      
          call #recursiveSymplecticCofficientsSplit first, then gether the \f$ c_k, d_k \f$ to cd_pair (KDKD... or DKDK...)
          @params[out] _cd_pair: cofficients pair for KDKD... or DKDK... [1:k][c_k,d_k]
          @params[out] _sorted_cumsum_ck_index: two dimensional array [1:k][index,cumsum(c_k)]; index is the increasing order index for cumsum(c_k)
          @params[in] _n: half of the order of integrator, if positive, provide solution one of Yoshida (1990), if negative, provide solution two (only first groups in Table 1/2 of Yoshida 1990), notice only -3 (6th order) and -4 (8th order) works
          @params[in] _cd_pair_size: _cd_pair array size (for safety check), should be >=k (3^(n-1)+1)
        */
        void symplecticCofficientsGenerator(CDPair _cd_pair[], CumSumCkIndex _sorted_cumsum_ck_index[], const int _n, const int _cd_pair_size) {
            if (_n>0) {
                const int n_size = std::pow(3,_n);
                Float coff_array[n_size];

                recursiveSymplecticCofficientsSplit(coff_array, _n);

                // DKD group number
                const int n_group= n_size/3;

                // _cd_pair array size
                const int k = n_group+1;
        
                if(k>_cd_pair_size) {
                    std::cerr<<"Error: _cd_pair array size "<<_cd_pair_size<<" not big enough for order "<<_n<<", required size "<<k<<std::endl;
                    abort();
                }

                _cd_pair[0][0] = coff_array[0];
                _sorted_cumsum_ck_index[0].cck = coff_array[0];
                _sorted_cumsum_ck_index[0].index = 0;
                for (int i=0; i<n_group-1; i++) {
                    _cd_pair[i  ][1] = coff_array[3*i+1];
                    _cd_pair[i+1][0] = coff_array[3*i+2]+coff_array[3*i+3];
                    _sorted_cumsum_ck_index[i+1].cck = _sorted_cumsum_ck_index[i].cck + _cd_pair[i+1][0];
                    _sorted_cumsum_ck_index[i+1].index = i+1;
                }
                _cd_pair[k-2][1] = coff_array[n_size-2];
                _cd_pair[k-1][0] = coff_array[n_size-1];
                _cd_pair[k-1][1] = Float(0.0);
                _sorted_cumsum_ck_index[k-1].cck = _sorted_cumsum_ck_index[k-2].cck + _cd_pair[k-1][0];
                _sorted_cumsum_ck_index[k-1].index = k-1;

                std::sort(_sorted_cumsum_ck_index, _sorted_cumsum_ck_index+k, [](const CumSumCkIndex &a, const CumSumCkIndex &b){ return a.cck<b.cck;} );
            
            }
            else if (_n==-3) {
                if(_cd_pair_size<8) {
                    std::cerr<<"Error: _cd_pair array size "<<_cd_pair_size<<" not big enough for order "<<_n<<", required size "<<8<<std::endl;
                    abort();
                }
                // Solution A
                _cd_pair[0][0] =  0.3922568052387800;
                _cd_pair[0][1] =  0.7845136104775600;
                _cd_pair[1][0] =  0.5100434119184585;
                _cd_pair[1][1] =  0.2355732133593570;
                _cd_pair[2][0] = -0.4710533854097566;
                _cd_pair[2][1] = -1.1776799841788701;
                _cd_pair[3][0] =  0.0687531682525181;
                _cd_pair[3][1] =  1.3151863206839063;
                _cd_pair[4][0] =  0.0687531682525181;
                _cd_pair[4][1] = -1.1776799841788701;
                _cd_pair[5][0] = -0.4710533854097566;
                _cd_pair[5][1] =  0.2355732133593570;
                _cd_pair[6][0] =  0.5100434119184585;
                _cd_pair[6][1] =  0.7845136104775600;
                _cd_pair[7][0] =  0.3922568052387800;
                _cd_pair[7][1] =  0.0000000000000000;
                _sorted_cumsum_ck_index[0].cck =   0.0976997828427615;
                _sorted_cumsum_ck_index[0].index = 5;
                _sorted_cumsum_ck_index[1].cck =   0.3922568052387800;
                _sorted_cumsum_ck_index[1].index = 0;
                _sorted_cumsum_ck_index[2].cck =   0.4312468317474820;
                _sorted_cumsum_ck_index[2].index = 2;
                _sorted_cumsum_ck_index[3].cck =   0.5000000000000000;
                _sorted_cumsum_ck_index[3].index = 3;
                _sorted_cumsum_ck_index[4].cck =   0.5687531682525181;
                _sorted_cumsum_ck_index[4].index = 4;
                _sorted_cumsum_ck_index[5].cck =   0.6077431947612200;
                _sorted_cumsum_ck_index[5].index = 6;
                _sorted_cumsum_ck_index[6].cck =   0.9023002171572385;
                _sorted_cumsum_ck_index[6].index = 1;
                _sorted_cumsum_ck_index[7].cck =   1.0000000000000000;
                _sorted_cumsum_ck_index[7].index = 7;
            }
            else if (_n==-4) {
                if(_cd_pair_size<16) {
                    std::cerr<<"Error: _cd_pair array size "<<_cd_pair_size<<" not big enough for order "<<_n<<", required size "<<16<<std::endl;
                    abort();
                }
                //Solution B 
                _cd_pair[0][0] =  0.4574221231148700;
                _cd_pair[0][1] =  0.9148442462297400;
                _cd_pair[1][0] =  0.5842687913979845;
                _cd_pair[1][1] =  0.2536933365662290;
                _cd_pair[2][0] = -0.5955794501471254;
                _cd_pair[2][1] = -1.4448522368604799;
                _cd_pair[3][0] = -0.8015464361143615;
                _cd_pair[3][1] = -0.1582406353682430;
                _cd_pair[4][0] =  0.8899492511272584;
                _cd_pair[4][1] =  1.9381391376227599;
                _cd_pair[5][0] = -0.0112355476763650;
                _cd_pair[5][1] = -1.9606102329754900;
                _cd_pair[6][0] = -0.9289051917917525;
                _cd_pair[6][1] =  0.1027998493919850;
                _cd_pair[7][0] =  0.9056264600894914;
                _cd_pair[7][1] =  1.7084530707869978;
                _cd_pair[8][0] =  0.9056264600894914;
                _cd_pair[8][1] =  0.1027998493919850;
                _cd_pair[9][0] = -0.9289051917917525;
                _cd_pair[9][1] = -1.9606102329754900;
                _cd_pair[10][0] = -0.0112355476763650;
                _cd_pair[10][1] =  1.9381391376227599;
                _cd_pair[11][0] =  0.8899492511272584;
                _cd_pair[11][1] = -0.1582406353682430;
                _cd_pair[12][0] = -0.8015464361143615;
                _cd_pair[12][1] = -1.4448522368604799;
                _cd_pair[13][0] = -0.5955794501471254;
                _cd_pair[13][1] =  0.2536933365662290;
                _cd_pair[14][0] =  0.5842687913979845;
                _cd_pair[14][1] =  0.9148442462297400;
                _cd_pair[15][0] =  0.4574221231148700;
                _cd_pair[15][1] =  0.0000000000000000;
                _sorted_cumsum_ck_index[0].cck =  -0.4056264600894914;
                _sorted_cumsum_ck_index[0].index = 6;
                _sorted_cumsum_ck_index[1].cck =  -0.3554349717486324;
                _sorted_cumsum_ck_index[1].index = 3;
                _sorted_cumsum_ck_index[2].cck =  -0.0416909145128544;
                _sorted_cumsum_ck_index[2].index = 13;
                _sorted_cumsum_ck_index[3].cck =   0.4461114643657291;
                _sorted_cumsum_ck_index[3].index = 2;
                _sorted_cumsum_ck_index[4].cck =   0.4574221231148700;
                _sorted_cumsum_ck_index[4].index = 0;
                _sorted_cumsum_ck_index[5].cck =   0.4654857206213739;
                _sorted_cumsum_ck_index[5].index = 10;
                _sorted_cumsum_ck_index[6].cck =   0.4767212682977390;
                _sorted_cumsum_ck_index[6].index = 9;
                _sorted_cumsum_ck_index[7].cck =   0.5000000000000000;
                _sorted_cumsum_ck_index[7].index = 7;
                _sorted_cumsum_ck_index[8].cck =   0.5232787317022610;
                _sorted_cumsum_ck_index[8].index = 5;
                _sorted_cumsum_ck_index[9].cck =   0.5345142793786261;
                _sorted_cumsum_ck_index[9].index = 4;
                _sorted_cumsum_ck_index[10].cck =   0.5425778768851300;
                _sorted_cumsum_ck_index[10].index = 14;
                _sorted_cumsum_ck_index[11].cck =   0.5538885356342710;
                _sorted_cumsum_ck_index[11].index = 12;
                _sorted_cumsum_ck_index[12].cck =   1.0000000000000000;
                _sorted_cumsum_ck_index[12].index = 15;
                _sorted_cumsum_ck_index[13].cck =   1.0416909145128546;
                _sorted_cumsum_ck_index[13].index = 1;
                _sorted_cumsum_ck_index[14].cck =   1.3554349717486325;
                _sorted_cumsum_ck_index[14].index = 11;
                _sorted_cumsum_ck_index[15].cck =   1.4056264600894914;
                _sorted_cumsum_ck_index[15].index = 8;
            }
        }
    public:

        //!constructor
        SymplecticStep(): sym_order_(0), sym_type_(0), cd_pair_array_size_(0),  cd_pair_(NULL), sorted_cumsum_ck_index_(NULL) {}

        //! clear function
        void clear() {
            if (cd_pair_array_size_>0) {
                ASSERT(cd_pair_!=NULL);
                delete [] cd_pair_;
                cd_pair_ = NULL;
                ASSERT(sorted_cumsum_ck_index_!=NULL);
                delete [] sorted_cumsum_ck_index_;
                sorted_cumsum_ck_index_ = NULL;
                cd_pair_array_size_ = 0;
            }
            sym_order_ = 0;
            sym_type_ = 0;
        }
        
        //! Symplectic coefficients generation for input order
        /*!
          @param [in] n: symplectic integrator order, should be even, otherwise reduce to closest even number; if 0, not set; if negative, use second solution from Yoshida (1990), but only -6 and -8 works
         */
        void initialSymplecticCofficients(const int _n) {
            // if already initialized, clear
            if(cd_pair_array_size_>0) clear();

            // calculate cd_pair array size
            int k = _n/2;
            if (_n>0) cd_pair_array_size_ = std::pow(3,k-1)+1;
            else if(_n==-6) cd_pair_array_size_ = 8;
            else if(_n==-8) cd_pair_array_size_ = 16;
            else {
                std::cerr<<"Error: input order "<<_n<<" cannot be generated !\n";
                abort();
            }

            // create new array
            cd_pair_ = new CDPair[cd_pair_array_size_];
            sorted_cumsum_ck_index_ = new CumSumCkIndex[cd_pair_array_size_];
            // generate coefficients
            symplecticCofficientsGenerator(cd_pair_, sorted_cumsum_ck_index_, k, cd_pair_array_size_);

            // store the order
            sym_order_ = _n<0?-_n:_n;
            
            // store the symplectic method type
            sym_type_ = _n<0?2:1;
        }

        //! get coefficient c_k
        Float getCK(const int _k) const {
            ASSERT(_k<cd_pair_array_size_);
            return cd_pair_[_k][0];
        }

        //! get coefficient d_k
        Float getDK(const int _k) const {
            ASSERT(_k<cd_pair_array_size_);
            return cd_pair_[_k][1];
        }

        //! get cd_pair array size
        int getCDPairSize() const {
            return cd_pair_array_size_;
        }

        //! get sorted cumsum CK index
        int getSortCumSumCKIndex(const int _k) const {
            return sorted_cumsum_ck_index_[_k].index;
        }

        //! get sorted cumsum CK 
        Float getSortCumSumCK(const int _k) const {
            return sorted_cumsum_ck_index_[_k].cck;
        }

        //! calculate the step modify factor based on the order scaling and input energy error ratio
        /*! 
          @param[in] _error_new_over_old: integration error change ratio expected to reach after modify step size (expected new error / old error)
         */
        Float calcStepModifyFactorFromErrorRatio(const Float _error_new_over_old) {
            return pow(_error_new_over_old, Float(1.0/sym_order_));
        }

        //! calculate the enstep modify factor based on the order scaling and input energy error ratio
        /*! 
          @param[in] _step_new_over_old: step size modify ratio (new step size / old)
         */
        Float calcErrorRatioFromStepModifyFactor(const Float _step_new_over_old) {
            return pow(_step_new_over_old, Float(sym_order_));
        }

        //! get Symplectic order
        int getOrder() const {
            return sym_order_;
        }
    };

}
