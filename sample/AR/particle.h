#pragma once
#include <iomanip>
#include "Common/Float.h"

#ifndef NAN_CHECK
#define NAN_CHECK(val) ASSERT((val) == (val));
#endif

enum class Status{single=1, premerge=2, unused=0};

//! A sample particle class
/*! A particle class should contain public members:
  Float mass, Float pos[3], Float vel[3], 
*/
class Particle{
public:
    long long int id;
    Float mass;
    Float pos[3];
    Float vel[3];
    Float radius;
    Float time_check; // time to check next interrupt
    Status status;

    Particle(): id(-1), mass(0.0), pos{0,0,0}, vel{0,0,0}, radius(0.0), time_check(0.0), status(Status::single) {}

    //! Get position (required)
    /*! \return position vector (Float[3])
     */
    Float* getPos() {
        return pos;
    }

    //! Get velocity (required)
    /*! \return velocity vector (Float[3])
     */
    Float* getVel() {
        return vel;
    }

    //! write class data to file with binary format (required)
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fout) const {
        fwrite(this, sizeof(*this),1,_fout);
    }


    //! read class data to file with binary format (required)
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! write class data to file with ASCII format (required)
    /*! @param[in] _fout: std:osteram file for output
     */
    void writeAscii(std::ostream& _fout) const {
        _fout<<mass<<" "
             <<pos[0]<<" "
             <<pos[1]<<" " 
             <<pos[2]<<" " 
             <<vel[0]<<" " 
             <<vel[1]<<" " 
             <<vel[2]<<" "
             <<radius<<" ";
    }

    //! read class data to file with ASCII format (required)
    /*! @param[in] _fin: std::istream file for input
     */
    void readAscii(std::istream&  _fin) {
        _fin>>mass>>pos[0]>>pos[1]>>pos[2]>>vel[0]>>vel[1]>>vel[2]>>radius;
    }
    
    //! print titles of class members using column style (required)
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"mass"
             <<std::setw(_width)<<"pos.x"
             <<std::setw(_width)<<"pos.y"
             <<std::setw(_width)<<"pos.z"
             <<std::setw(_width)<<"vel.x"
             <<std::setw(_width)<<"vel.y"
             <<std::setw(_width)<<"vel.z"
             <<std::setw(_width)<<"radius"
             <<std::setw(_width)<<"id";
    }

    //! print data of class members using column style (required)
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<mass
             <<std::setw(_width)<<pos[0]
             <<std::setw(_width)<<pos[1]
             <<std::setw(_width)<<pos[2]
             <<std::setw(_width)<<vel[0]
             <<std::setw(_width)<<vel[1]
             <<std::setw(_width)<<vel[2]
             <<std::setw(_width)<<radius
             <<std::setw(_width)<<id;
    }
    
};
