#pragma once
#include <cassert>
#include <iomanip>
#include "Common/Float.h"

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

enum class BinaryInterruptState:int {none = 0, form = 1, exchange = 2, collision = 3};
#define BINARY_STATE_ID_SHIFT 4
#define BINARY_INTERRUPT_STATE_MASKER 0xF

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
    Float dm;
    Float time_check; // time to check next interrupt
    long long int binary_state; // contain two parts, low bits (first BINARY_STATE_ID_SHIFT bits) is binary interrupt state and high bits are pair ID
    static Float r_neighbor_crit;
    static Float r_break_crit;

    Particle(): id(-1), mass(0.0), pos{0,0,0}, vel{0,0,0}, radius(0.0), dm(0.0), time_check(NUMERIC_FLOAT_MAX), binary_state(0) {}

    //! save pair id in binary_state with shift bit size of BINARY_STATE_ID_SHIFT
    void setBinaryPairID(const int _id) {
        binary_state = (binary_state&BINARY_INTERRUPT_STATE_MASKER) | (_id<<BINARY_STATE_ID_SHIFT);
    }

    //! save binary interrupt state in the first  BINARY_STATE_ID_SHIFT bit in binary_state
    void setBinaryInterruptState(const BinaryInterruptState _state) {
        binary_state = ((binary_state>>BINARY_STATE_ID_SHIFT)<<BINARY_STATE_ID_SHIFT) | int(_state);
    }

    //! get binary interrupt state from binary_state
    BinaryInterruptState getBinaryInterruptState() const {
        return static_cast<BinaryInterruptState>(binary_state&BINARY_INTERRUPT_STATE_MASKER);
    }

    //! get pair ID from binary_state 
    int getBinaryPairID() const {
        return (binary_state>>BINARY_STATE_ID_SHIFT);
    }

    //! Get position 
    /*! \return position vector (Float[3])
     */
    Float* getPos() {
        return pos;
    }

    //! Get velocity 
    /*! \return velocity vector (Float[3])
     */
    Float* getVel() {
        return vel;
    }

    //! Get neighbor distance criterion 
    Float getRNeighbor() {
        return r_neighbor_crit;
    }

    //! Get Group distance criterion 
    Float getRGroup() {
        return r_break_crit;
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fout) const {
        fwrite(this, sizeof(*this),1,_fout);
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

    ////! write class data to file with ASCII format
    ///*! @param[in] _fp: FILE type file for output
    // */
    //void writeAscii(FILE *_fout) const {
    //    fprintf(_fout, "%26.17e %26.17e %26.17e %26.17e %26.17e %26.17e %26.17e ",
    //            this->mass, 
    //            this->pos[0], this->pos[1], this->pos[2],  
    //            this->vel[0], this->vel[1], this->vel[2]);
    //}
    // 
    ////! read class data to file with ASCII format
    ///*! @param[in] _fin: FILE type file for input
    // */
    //void readAscii(FILE* _fin) {
    //    int rcount=fscanf(_fin, "%lf %lf %lf %lf %lf %lf %lf ",
    //                      &this->mass, 
    //                      &this->pos[0], &this->pos[1], &this->pos[2],
    //                      &this->vel[0], &this->vel[1], &this->vel[2]);
    //    if(rcount<7) {
    //        std::cerr<<"Error: Data reading fails! requiring data number is 7, only obtain "<<rcount<<".\n";
    //        abort();
    //    }
    //}

    //! write class data to file with ASCII format
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

    //! read class data to file with ASCII format
    /*! @param[in] _fin: std::istream file for input
     */
    void readAscii(std::istream&  _fin) {
        _fin>>mass>>pos[0]>>pos[1]>>pos[2]>>vel[0]>>vel[1]>>vel[2]>>radius;
    }
    
    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"id"
             <<std::setw(_width)<<"mass"
             <<std::setw(_width)<<"pos.x"
             <<std::setw(_width)<<"pos.y"
             <<std::setw(_width)<<"pos.z"
             <<std::setw(_width)<<"vel.x"
             <<std::setw(_width)<<"vel.y"
             <<std::setw(_width)<<"vel.z"
             <<std::setw(_width)<<"radius";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<id
             <<std::setw(_width)<<mass
             <<std::setw(_width)<<pos[0]
             <<std::setw(_width)<<pos[1]
             <<std::setw(_width)<<pos[2]
             <<std::setw(_width)<<vel[0]
             <<std::setw(_width)<<vel[1]
             <<std::setw(_width)<<vel[2]
             <<std::setw(_width)<<radius;
    }
    
};

Float Particle::r_neighbor_crit = -1.0;
Float Particle::r_break_crit = -1.0;


