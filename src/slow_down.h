#pragma once

#include <cassert>
#include "Float.h"

namespace AR{
    
    //! Slow-down parameter control class
    /*! Determine the slow-down factor due to the perturbation and internal force
      \f$ \kappa = k_0 / [F_{pert,max}/F_{inner}] \f$
    */
    class SlowDown{
    private:
        Float time_real_;      // current time
        Float kappa_;          // slow-down factor
        Float kappa_org_;      // original slow-down factor without kappa limit (1.0, kappa_max)
        Float kappa_max_;      // maximum kappa factor
        Float kappa_ref_;      ///< reference kappa factor; slow-down factor kappa = max(1,kref/perturbation_factor)

    public:
        //! defaulted constructor
        SlowDown(): time_real_(Float(0.0)), kappa_(Float(1.0)), kappa_org_(Float(1.0)), kappa_max_(Float(1.0)), kappa_ref_(Float(1.0e-6)) {}
    
        //! clear function
        void clear(){
            time_real_ = Float(0.0);
            kappa_ = kappa_org_ = kappa_max_ = Float(1.0);
            kappa_ref_ = Float(1.0e-6);
        }

        //! initialize slow-down parameters
        /*! Set slow-down parameters, slow-down method will be switched on
          @param[in] _kappa_ref: reference kappa factor; slow-down factor kappa = max(1,kref/perturbation_factor)
          @param[in] _kappa_max: maximum slow-down factor
        */
        void initialSlowDownReference(const Float _kappa_ref, const Float _kappa_max) {
            assert(_kappa_ref>0.0);
            assert(_kappa_max>1.0);
            kappa_ref_ = _kappa_ref;
            kappa_max_ = _kappa_max;
        }

        //! set real time value
        void setRealTime(const Float _time_real) {
            time_real_ = _time_real;
        }

        //! get real time value
        Float getRealTime() const {
            return time_real_;
        }

        //! drift real time
        /*! drift real time with integrated time step x kappa
          @param[in] _dt_int: integrated time step
        */
        void driftRealTime(const Float _dt_int) {
            time_real_ += _dt_int * kappa_;
        }

        //! manually set kappa
        /* auto-adjust in the range (1.0, kappa_max)
         */
        void setSlowDownFactor(const Float _kappa) {
            kappa_ = _kappa;
            kappa_org_ = _kappa;
            kappa_ = std::min(kappa_, kappa_max_);
            kappa_ = std::max(Float(1.0), kappa_);
        }

        //! calculate slowdown factor based on perturbation force and inner binding energy for a binary
        /* if it is a hyperbolic encounter, (ebin_>0), set slowdown factor to 1.0
          @param[in] _ebin: binary binding energy
          @param[in] _fpert: perturbation force (Notice this is scaler force not accelertion), if zero, slowdown factor is 1.0
          \return slowdown factor
         */
        Float calcSlowDownFactorBinary(const Float _ebin, const Float _fpert) {
            // hyberbolic case
            if(_ebin>0.0||_fpert==0.0) {
                kappa_org_ = kappa_ = Float(1.0);
            }
            else { // binary case
                kappa_org_ = - kappa_ref_*_ebin/_fpert;
                kappa_ = std::min(kappa_, kappa_max_);
                kappa_ = std::max(Float(1.0), kappa_);
            }
            return kappa_;
        }

        // Get slow-down factor
        /*!
          \return get adjusted kappa by keeping phase corrected
        */
        Float getSlowDownFactor() const {
            return kappa_;
        }

        // Get original slow-down factor
        /*!
          \return kappa_origin
        */
        Float getSlowDownFactorOrigin() const {
            return kappa_org_;
        }

        // Get sd reference factor
        Float getSlowDownFactorReference() const {
            return kappa_ref_;
        }

        //! write class data with BINARY format
        /*! @param[in] _fout: file IO for write
         */
        void writeBINARY(FILE *_fout) {
            fwrite(this, sizeof(*this),1,_fout);
        }

        //! read class data with BINARY format
        /*! @param[in] _fin: file IO for read
         */
        void read(FILE *_fin) {
            size_t rcount = fread(this, sizeof(*this), 1, _fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
        
        //! get backup data size 
        /*! \return the data array size for backupSlowDownFactorAndTimeReal()
         */
        int getBackupDataSize() const {
            return 2;
        }

        //! backup real time and force ratio
        /*! @param[in] _bk: backup data array, should be size of getBackupDataSize() (2)
          \return backup array size
         */
        int backupSlowDownFactorAndTimeReal(Float* _bk) {
            _bk[0] = kappa_;
            _bk[1] = time_real_;
            return 2;
        }

        //! restore real time and force ratio
        /*! @param[in] _bk: restore data array[3]
          \return backup array size
         */
        int restoreSlowDownFactorAndTimeReal(Float* _bk) {
            kappa_     =   _bk[0];
            time_real_ =   _bk[1];
            return 2;
        }
    
        //! print slowdown data
        /*! Print slowdown data 
          @param[in] fout: ofstream for printing
          @param[in] precision: printed precision for one variable
          @param[in] width: printing width for one variable
        */
        void print(std::ostream & fout, const int precision=15, const int width=23) {
            assert(width>0);
            assert(precision>0);
            fout<<"time_real= "<<std::setw(width)<<time_real_
                <<"kappa= "<<std::setw(width)<<kappa_
                <<"kappa_org= "<<std::setw(width)<<kappa_org_
                <<"kappa_max= "<<std::setw(width)<<kappa_max_
                <<"kappa_ref= "<<std::setw(width)<<kappa_ref_;
        }
    
        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"Time_real"
                 <<std::setw(_width)<<"SD_factor"
                 <<std::setw(_width)<<"SD_factor_org";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<time_real_
                 <<std::setw(_width)<<kappa_
                 <<std::setw(_width)<<kappa_org_;
        }
    };
}
