#pragma once

#include "Common/Float.h"

namespace H4{
    template <class Tparticle> class ParticleH4;

    //! particle type for AR integrator
    template <class Tparticle>
    class ParticleAR: public Tparticle{
    public:
        Tparticle* adr; // original particle address

        ParticleAR& operator = (Tparticle & _p) {
            *(Tparticle*)this = _p;
            adr = &_p;
            return *this;
        }

        ParticleAR& operator = (ParticleH4<Tparticle> & _p) {
            *(Tparticle*)this = *(Tparticle*)&_p;
            adr = &_p;
            return *this;
        }
        
        ParticleAR& operator = (ParticleAR & _p) {
            *(Tparticle*)this = *(Tparticle*)&_p;
            adr = _p.adr;
            return *this;
        }
    };

    //! Particle type for hermite integrator
    template <class Tparticle>
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
        Float pot;

        ParticleH4() {}

        ParticleH4(const Tparticle & _p) {
            *(Tparticle*)this = *(Tparticle*)&_p;
        }

        ParticleH4(const ParticleH4<Tparticle> & _p) {
            *(Tparticle*)this = *(Tparticle*)&_p;
            dt = _p.dt;
            time = _p.time;
            acc0[0] = _p.acc0[0];
            acc0[1] = _p.acc0[1];
            acc0[2] = _p.acc0[2];
            acc1[0] = _p.acc1[0];
            acc1[1] = _p.acc1[1];
            acc1[2] = _p.acc1[2];
#ifdef HERMITE_DEBUG_ACC
            acc2[0] = _p.acc2[0];
            acc2[1] = _p.acc2[1];
            acc2[2] = _p.acc2[2];
            acc3[0] = _p.acc3[0];
            acc3[1] = _p.acc3[1];
            acc3[2] = _p.acc3[2];
#endif
            pot = _p.pot;
        }

        ParticleH4& operator = (const Tparticle & _p) {
            *(Tparticle*)this = *(Tparticle*)&_p;
            return *this;
        }

        ParticleH4& operator = (const ParticleH4<Tparticle> & _p) {
            *(Tparticle*)this = *(Tparticle*)&_p;
            dt = _p.dt;
            time = _p.time;
            acc0[0] = _p.acc0[0];
            acc0[1] = _p.acc0[1];
            acc0[2] = _p.acc0[2];
            acc1[0] = _p.acc1[0];
            acc1[1] = _p.acc1[1];
            acc1[2] = _p.acc1[2];
#ifdef HERMITE_DEBUG_ACC
            acc2[0] = _p.acc2[0];
            acc2[1] = _p.acc2[1];
            acc2[2] = _p.acc2[2];
            acc3[0] = _p.acc3[0];
            acc3[1] = _p.acc3[1];
            acc3[2] = _p.acc3[2];
#endif
            pot = _p.pot;
            return *this;
        }

        //! print function for one line
        void print(std::ostream & _fout) const{
            Tparticle::print(_fout);
            _fout<<" dt="<<dt
                 <<" time="<<time
                 <<" acc0="<<acc0
                 <<" acc1="<<acc0;
#ifdef HERMITE_DEBUG_ACC
            _fout<<" acc2="<<acc2
                 <<" acc3="<<acc3;
#endif
            _fout<<" pot="<<pot;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            Tparticle::printColumnTitle(_fout, _width);
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
            _fout<<std::setw(_width)<<"pot";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            Tparticle::printColumn(_fout, _width);
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
            _fout<<std::setw(_width)<<pot;
        }

        ////! write class data to file with binary format
        ///*! @param[in] _fp: FILE type file for output
        // */
        //void writeBinary(FILE *_fp) const {
        //    fwrite(this, sizeof(*this),1,_fp);
        //}
        // 
        ////! read class data to file with binary format
        ///*! @param[in] _fp: FILE type file for reading
        // */
        //void readBinary(FILE *_fin) {
        //    size_t rcount = fread(this, sizeof(*this),1,_fin);
        //    if (rcount<1) {
        //        std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
        //        abort();
        //    }
        //}

    };

    //! Particle force for hermite integrator
    class ForceH4{
    public:
        Float acc0[3]; //
        Float acc1[3]; //
        Float pot;

        //! clear function
        void clear() {
            acc0[0] = acc0[1] = acc0[2] = Float(0.0);
            acc1[0] = acc1[1] = acc1[2] = Float(0.0);
            pot = 0.0;
        }
    };

}
