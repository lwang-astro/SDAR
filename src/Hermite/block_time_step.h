#pragma once

#include "Common/Float.h"

namespace H4{
    // forth order Time step calculator class
    class TimeStep4th{
    public:
        Float eta_4th;  ///> time step coefficient (outside sqrt) for forth order
        Float eta_2nd;  ///> time step coefficient (outside sqrt) for second order

        TimeStep4th(): eta_4th(0.02), eta_2nd(0.005) {}

        //! calculate 2nd order time step 
        /*! calculate time step based on Acc and its derivatives
          @param[in] acc0: acceleration
          @param[in] acc1: first order derivative of acc0
          \return time step
        */
        Float calcDt2nd(const Float* acc0, 
                        const Float* acc1) const {
            const Float s0 = acc0[0] * acc0[0] + acc0[1] * acc0[1] + acc0[2] * acc0[2];
            const Float s1 = acc1[0] * acc1[0] + acc1[1] * acc1[1] + acc1[2] * acc1[2];
            if(s0 == 0.0 || s1 == 0.0)
                return NUMERIC_FLOAT_MAX;
            else 
                return eta_2nd * sqrt( s0 / s1 );
        }

        //! calculate 4th order time step 
        /*! calculate time step based on Acc and its derivatives
          @param[in] acc0: acceleration
          @param[in] acc1: first order derivative of acc0
          @param[in] acc2: second order derivative of acc0
          @param[in] acc3: thrid  order derivative of acc0
          \return time step
        */
        Float calcDt4th(const Float* acc0,
                        const Float* acc1,
                        const Float* acc2,
                        const Float* acc3) const {
            const Float s0 = acc0[0] * acc0[0] + acc0[1] * acc0[1] + acc0[2] * acc0[2];
            const Float s1 = acc1[0] * acc1[0] + acc1[1] * acc1[1] + acc1[2] * acc1[2];
            const Float s2 = acc2[0] * acc2[0] + acc2[1] * acc2[1] + acc2[2] * acc2[2];
            const Float s3 = acc3[0] * acc3[0] + acc3[1] * acc3[1] + acc3[2] * acc3[2];
            if(!(s1==0.0||s3==0.0)&&s2==0.0)
                return NUMERIC_FLOAT_MAX;
            else 
                return eta_4th * sqrt( (sqrt(s0*s2) + s1) / (sqrt(s1*s3) + s2) );
        }

        void print(std::ostream & _fout) const{
            _fout<<" eta_4th="<<eta_4th
                 <<" eta_2nd="<<eta_2nd;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"Eta(4th)"
                 <<std::setw(_width)<<"Eta(2nd)";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<eta_4th
                 <<std::setw(_width)<<eta_2nd;
        }

        //! write class data to file with binary format
        /*! @param[in] _fp: FILE type file for output
         */
        void writeBinary(FILE *_fp) const {
            fwrite(this, sizeof(*this),1,_fp);
        }

        //! read class data to file with binary format
        /*! @param[in] _fp: FILE type file for reading
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(this, sizeof(*this),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
    };

    // Block time step class
    class BlockTimeStep4th: public TimeStep4th{
    private:
        Float dt_max_; // maximum time step
        Float dt_min_; // minimum time step
    public:
        //! contructor
        BlockTimeStep4th(): TimeStep4th(), dt_max_(-1.0), dt_min_(-1.0) {}

        //! set dt limit (max and min)
        /*! 
          @param[in] _dt_max: maximum step size
          @param[in] _pow_index_min: the power index (k) to set the minimum step size: _dt_max*(0.5^k)
         */
        void setDtRange(const Float _dt_max, const int _pow_index_min) {
            dt_max_ = _dt_max;
            dt_min_ = _dt_max;
            for (int i=0; i<_pow_index_min; i++) {
                dt_min_ *= Float(0.5);
            }
        }

        //! get maximum time step
        /*! \return maximum time step
         */
        Float getDtMax() const {
            return dt_max_;
        }

        //! get minimum time step
        /*! \return miniimum time step
         */
        Float getDtMin() const {
            return dt_min_;
        }

        //! calculate the maximum time step limit for next block step based on the input (current) time
        /*! 
          Basic algorithm: the integer of time/dt_min is the binary tree for block step, counting from the minimum digital, the last zero indicate the maximum block step level allown for next step
          @param[in] _time: current time
        */
        Float calcNextDtLimit(const Float _time) {
            ASSERT(dt_max_>dt_min_);
            ASSERT(dt_min_>0.0);
            // for first step, the maximum time step is OK
            if(_time==0.0) return dt_max_;
            else {
                // the binary tree for current time position in block step 
                unsigned long long int bitmap = to_int(_time/dt_min_);
                //#ifdef __GNUC__ 
                //        PS::S64 dts = __builtin_ctz(bitmap) ;
                //        PS::U64 c = (1<<dts);
                ////        std::cerr<<"time = "<<_time<<"  dt_min = "<<_dt_min<<"  bitmap = "<<bitmap<<"  dts = "<<dts<<std::endl;
                //#else

                // block step multiply factor 
                unsigned long long int c=1;
                // find the last zero in the binary tree to obtain the current block step level
                while((bitmap&1)==0) {
                    bitmap = (bitmap>>1);
                    c = (c<<1);
                }
                //#endif
            
                // return the maximum step allown
                return std::min(c*dt_min_,dt_max_);
            }
        }

        //! Calculate 2nd order block time step
        /*! calculate block time step based on Acc and its derivatives
          @param[in] _acc0: acceleration
          @param[in] _acc1: first order derivative of acc0
          @param[in] _dt_limit: maximum step limit
          \return step size, if dt<dt_min, return -dt;
        */
        Float calcBlockDt2nd(const Float* _acc0, 
                             const Float* _acc1,
                             const Float _dt_limit) const{
            ASSERT(dt_max_>dt_min_);
            ASSERT(_dt_limit<=dt_max_);
            ASSERT(_dt_limit>=dt_min_);

            const Float dt_ref = TimeStep4th::calcDt2nd(_acc0, _acc1);
            Float dt = _dt_limit;
            while(dt > dt_ref) dt *= 0.5;

            if(dt<dt_min_) {
                std::cerr<<"Error: time step size too small: ("<<dt<<") < dt_min ("<<dt_min_<<")!"<<std::endl;
                DATADUMP();
                abort();
            }
            else return dt;
        }

        //! Calculate 4th order block time step
        /*! calculate 4th order block time step based on Acc and its derivatives
          @param[in] acc0: acceleration
          @param[in] acc1: first order derivative of acc0
          @param[in] acc2: second order derivative of acc0
          @param[in] acc3: thrid  order derivative of acc0
          \return step size, if dt<dt_min, return -dt;
        */
        Float calcBlockDt4th(const Float* acc0, 
                             const Float* acc1,
                             const Float* acc2,
                             const Float* acc3,
                             const Float _dt_limit) const {
            ASSERT(dt_max_>dt_min_);
            ASSERT(_dt_limit<=dt_max_);
            ASSERT(_dt_limit>=dt_min_);

            const Float dt_ref = TimeStep4th::calcDt4th(acc0, acc1, acc2, acc3);
            Float dt = _dt_limit;
            while(dt > dt_ref) dt *= 0.5;

            if(dt<dt_min_) {
                std::cerr<<"Error: time step size too small: ("<<dt<<") < dt_min ("<<dt_min_<<")!"<<std::endl;
                DATADUMP();
                abort();
            }
            else return dt;
        }

        void print(std::ostream & _fout) const{
            TimeStep4th::print(_fout);
            _fout<<" dt_max="<<dt_max_
                 <<" dt_min="<<dt_min_;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            TimeStep4th::printColumnTitle(_fout, _width);
            _fout<<std::setw(_width)<<"Dt_max"
                 <<std::setw(_width)<<"Dt_min";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            TimeStep4th::printColumn(_fout, _width);
            _fout<<std::setw(_width)<<dt_max_
                 <<std::setw(_width)<<dt_min_;
        }

        //! write class data to file with binary format
        /*! @param[in] _fp: FILE type file for output
         */
        void writeBinary(FILE *_fp) const {
            fwrite(this, sizeof(*this),1,_fp);
        }

        //! read class data to file with binary format
        /*! @param[in] _fp: FILE type file for reading
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(this, sizeof(*this),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
    };

}
