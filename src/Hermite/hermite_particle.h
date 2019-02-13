#pragma once

#include "Float.h"

namespace H4{

    //! Particle type for hermite integrator
    template <Tparticle>
    class ParticleH4: public Tparticle{
    public:
        Float dt;
        Float time;
        Float acc0[3];
        Float acc1[3];
#ifdef HERMITE_DEBUG_ACC
        Float acc2[3]; // for debug
        Float acc3[3]; // for debug
#endif

        //! print function for one line
        void print(std::ostream & _fout) const{
            Ptcl::print(_fout);
            _fout<<" dt="<<dt
                 <<" time="<<time
                 <<" acc0="<<acc0
                 <<" acc1="<<acc0;
#ifdef HERMITE_DEBUG_ACC
            _fout<<" acc2="<<acc2
                 <<" acc3="<<acc3;
#endif
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            Ptcl::printColumnTitle(_fout, _width);
            _fout<<std::setw(_width)<<"dt"
                 <<std::setw(_width)<<"time"
                 <<std::setw(_width)<<"acc0.x"
                 <<std::setw(_width)<<"acc0.y"
                 <<std::setw(_width)<<"acc0.z"
                 <<std::setw(_width)<<"acc1.x"
                 <<std::setw(_width)<<"acc1.y"
                 <<std::setw(_width)<<"acc1.z";
#ifdef HERMITE_DEBUG_ACC
            _fout<<std::setw(_width)<<"acc2.x"
                 <<std::setw(_width)<<"acc2.y"
                 <<std::setw(_width)<<"acc2.z"
                 <<std::setw(_width)<<"acc3.x"
                 <<std::setw(_width)<<"acc3.y"
                 <<std::setw(_width)<<"acc3.z";
#endif
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<dt
                 <<std::setw(_width)<<time
                 <<std::setw(_width)<<acc0[0]
                 <<std::setw(_width)<<acc0[1]
                 <<std::setw(_width)<<acc0[2]
                 <<std::setw(_width)<<acc1[0]
                 <<std::setw(_width)<<acc1[1]
                 <<std::setw(_width)<<acc1[2];
#ifdef HERMITE_DEBUG_ACC
            _fout<<std::setw(_width)<<acc2[0]
                 <<std::setw(_width)<<acc2[1]
                 <<std::setw(_width)<<acc2[2]
                 <<std::setw(_width)<<acc3[0]
                 <<std::setw(_width)<<acc3[1]
                 <<std::setw(_width)<<acc3[2];
#endif
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

    //! Particle force for hermite integrator
    class ForceH4{
    public:
        Float acc0[3]; //
        Float acc1[3]; //
    };

}
