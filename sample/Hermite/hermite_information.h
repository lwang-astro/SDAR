#pragma once

#include "Common/Float.h"
#include "hermite_interaction.h"

class HermiteInformation{
public:
    Float time; // current time
    Float de;    // energy difference
    Float etot_init; // initial total energy;
    Float etot; // total energy
    Float ekin; // kinetic energy
    Float epot; // potential energy
    Float de_sd; // slowdown energy difference
    Float etot_sd_ref; // reference slowdown energy for calculating de
    Float etot_sd; // slowdown energy record
    Float ekin_sd; // slowdown kinetic energy record
    Float epot_sd; // slowdown potential energy record

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        return true;
    }        

    //! calculate energy of particle group
    template <class TList>
    void calcEnergy(TList& _particles, HermiteInteraction& _interaction, const bool _initial_flag) {
        ekin = epot = etot = 0.0;
        const int n = _particles.getSize();
        for (int i=0; i<n; i++) {
            auto& pi = _particles[i];
            ekin += pi.mass* (pi.vel[0]*pi.vel[0] + pi.vel[1]*pi.vel[1] + pi.vel[2]*pi.vel[2]);
            Float poti = 0.0;
            for (int j=0; j<i; j++) {
                poti +=_interaction.calcPotPair(pi, _particles[j]);
            }
            epot += poti*pi.mass;
        }
        ekin *= 0.5;
        etot = ekin + epot;

        if (_initial_flag) {
            etot_init = etot;
            de = 0.0;
        }
        else de = etot - etot_init;
    }

    //! correct Etot slowdown reference due to the groups change
    template <class TGroupList>
    void correctEtotSlowDownRef(TGroupList& _groups) {
        for (int i=0; i<_groups.getSize(); i++) {
            auto& gi = _groups[i];
            etot_sd_ref += gi.getDESlowDownCum();
            gi.resetDESlowDownCum();
        } 
    }

    //! calculate energy of particle group
    template <class TList, class TGroupList>
    void calcEnergySlowDown(TList& _particles, TGroupList& _groups, HermiteInteraction& _interaction, const bool _initial_flag) {
        ekin = epot = etot = 0.0;
        const int n = _particles.getSize();
        for (int i=0; i<n; i++) {
            auto& pi = _particles[i];
            ekin += pi.mass* (pi.vel[0]*pi.vel[0] + pi.vel[1]*pi.vel[1] + pi.vel[2]*pi.vel[2]);
            Float poti = 0.0;
            for (int j=0; j<i; j++) {
                poti +=_interaction.calcPotPair(pi, _particles[j]);
            }
            epot += poti*pi.mass;
        }
        ekin *= 0.5;
        etot = ekin + epot;

        ekin_sd = ekin;
        epot_sd = epot;
        for (int i=0; i<_groups.getSize(); i++) {
            auto& gi = _groups[i];
            Float kappa = gi.slowdown.getSlowDownFactor();
            Float kappa_inv = 1.0/kappa;
            ekin_sd -= gi.getEkin();
            epot_sd -= gi.getEpot(); 
#ifdef AR_TTL_SLOWDOWN_INNER
            ekin_sd += kappa_inv * gi.getEkinSlowDownInner();
            epot_sd += kappa_inv * gi.getEpotSlowDownInner();
#else
            ekin_sd += kappa_inv * gi.getEkin();
            epot_sd += kappa_inv * gi.getEpot(); 
#endif
        }
        etot_sd = ekin_sd + epot_sd;

        if (_initial_flag) {
            etot_init = etot;
            etot_sd_ref = etot_sd;
            de = 0.0;
            de_sd = 0.0;
        }
        else {
            de = etot - etot_init;
            de_sd = etot_sd - etot_sd_ref;
        }
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"Time_int"
             <<std::setw(_width)<<"dE"
             <<std::setw(_width)<<"Etot"
             <<std::setw(_width)<<"Ekin"
             <<std::setw(_width)<<"Epot"
             <<std::setw(_width)<<"dE_SD"
             <<std::setw(_width)<<"Etot_SD_ref"
             <<std::setw(_width)<<"Etot_SD"
             <<std::setw(_width)<<"Ekin_SD"
             <<std::setw(_width)<<"Epot_SD";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<time
             <<std::setw(_width)<<de
             <<std::setw(_width)<<etot
             <<std::setw(_width)<<ekin
             <<std::setw(_width)<<epot
             <<std::setw(_width)<<de_sd
             <<std::setw(_width)<<etot_sd_ref
             <<std::setw(_width)<<etot_sd
             <<std::setw(_width)<<ekin_sd
             <<std::setw(_width)<<epot_sd;
    }
};
