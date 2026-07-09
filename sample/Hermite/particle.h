#pragma once
#include <cassert>
#include <iomanip>
#include <sstream>
#include "Common/Float.h"

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

enum class BinaryInterruptState:int {none = 0, form = 1, exchange = 2, collision = 3, delaycollision = 4};
#define BINARY_STATE_ID_SHIFT 4
#define BINARY_INTERRUPT_STATE_MASKER 0xF

//! A sample particle class
/*! A particle class should contain public members:
  Float mass, Float pos[3], Float vel[3], 
*/
class Particle{
public:
    Float mass;
    Float pos[3];
    Float vel[3];
    Float radius;
    long long int id;
    Float dm;
    Float time_check; // time to check next interrupt
    long long int binary_state; // contain two parts, low bits (first BINARY_STATE_ID_SHIFT bits) is binary interrupt state and high bits are pair ID
    Float r_group_crit;      // per-particle group radius criterion
    Float r_neighbor_crit;   // per-particle neighbor radius criterion

    Particle(): mass(0.0), pos{0,0,0}, vel{0,0,0}, radius(0.0), id(-1), dm(0.0), time_check(NUMERIC_FLOAT_MAX), binary_state(0), r_group_crit(-1.0), r_neighbor_crit(-1.0) {}

    //! Initialize group and neighbor radii with mass-dependent weighting
    /*! Similar to PeTar's ChangeOver::setR() logic.
      @param[in] _r_break: base group radius (r_break input)
      @param[in] _r_neighbor_over_group: coefficient to compute neighbor radius from group radius
      @param[in] _mass_ref: reference mass (typically average mass)
    */
    void setRGroupAndNeighbor(const Float _r_break,
                              const Float _r_neighbor_over_group,
                              const Float _mass_ref) {
        Float mass_factor = std::max(pow(mass / _mass_ref, Float(1.0/3.0)), Float(1.0));
        r_group_crit    = _r_break * mass_factor;
        r_neighbor_crit = r_group_crit * _r_neighbor_over_group;
    }

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
    Float getRNeighbor() const {
        return r_neighbor_crit;
    }

    //! set neighbor distance criterion
    void setRNeighbor(const Float _r_neighbor) {
        r_neighbor_crit = _r_neighbor;
    }

    //! Get Group distance criterion 
    Float getRGroup() const {
        return r_group_crit;
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
             <<radius<<" "
             <<r_group_crit<<" "
             <<r_neighbor_crit<<" ";
    }

    //! read class data to file with ASCII format
    /*! @param[in] _fin: std::istream file for input
        Supports both new format (with r_group_crit, r_neighbor_crit) and
        old format (without them) via line-based parsing.
     */
    void readAscii(std::istream&  _fin) {
        std::string line;
        std::getline(_fin, line);
        std::istringstream iss(line);
        iss>>mass>>pos[0]>>pos[1]>>pos[2]>>vel[0]>>vel[1]>>vel[2]>>radius;
        if (iss>>r_group_crit) {
            iss>>r_neighbor_crit;
        } else {
            r_group_crit = -1.0;
            r_neighbor_crit = -1.0;
        }
    }
    
    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    static void printColumnTitleAscii(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"mass"
             <<std::setw(_width)<<"pos.x"
             <<std::setw(_width)<<"pos.y"
             <<std::setw(_width)<<"pos.z"
             <<std::setw(_width)<<"vel.x"
             <<std::setw(_width)<<"vel.y"
             <<std::setw(_width)<<"vel.z"
             <<std::setw(_width)<<"radius"
             <<std::setw(_width)<<"id"
             <<std::setw(_width)<<"r_group"
             <<std::setw(_width)<<"r_neighbor";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnAscii(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<mass
             <<std::setw(_width)<<pos[0]
             <<std::setw(_width)<<pos[1]
             <<std::setw(_width)<<pos[2]
             <<std::setw(_width)<<vel[0]
             <<std::setw(_width)<<vel[1]
             <<std::setw(_width)<<vel[2]
             <<std::setw(_width)<<radius
             <<std::setw(_width)<<id
             <<std::setw(_width)<<r_group_crit
             <<std::setw(_width)<<r_neighbor_crit;
    }
    
};


