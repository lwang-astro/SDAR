#pragma once
#include "ptcl.hpp"
#include "neighbor.hpp"
#include "changeover.hpp"
#include "timestep.hpp"

#ifndef ARRAY_ALLOW_LIMIT
#define ARRAY_ALLOW_LIMIT 1000000000
#endif

class PtclH4: public Ptcl{
public:
    PS::F64 dt;
    PS::F64 time;
    PS::F64vec acc0;
    PS::F64vec acc1;
#ifdef HERMITE_DEBUG_ACC
    PS::F64vec acc2; // for debug
    PS::F64vec acc3; // for debug
#endif

    PtclH4(): dt(0.0), time(0.0) {}
    
    template<class Tptcl>
    PtclH4(const Tptcl &_p): Ptcl(_p), dt(0.0), time(0.0) {}

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
             <<std::setw(_width)<<acc0.x
             <<std::setw(_width)<<acc0.y
             <<std::setw(_width)<<acc0.z
             <<std::setw(_width)<<acc1.x
             <<std::setw(_width)<<acc1.y
             <<std::setw(_width)<<acc1.z;
#ifdef HERMITE_DEBUG_ACC
        _fout<<std::setw(_width)<<acc2.x
             <<std::setw(_width)<<acc2.y
             <<std::setw(_width)<<acc2.z
             <<std::setw(_width)<<acc3.x
             <<std::setw(_width)<<acc3.y
             <<std::setw(_width)<<acc3.z;
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

class PtclPred{
public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 mass;
    PS::F64 r_search;

    void print(std::ostream & _fout) const{
        _fout<<" pos="<<pos
             <<" vel="<<vel
             <<" mass="<<mass
             <<" r_search="<<r_search;
    }
};

class PtclForce{
public:
    PS::F64vec acc0; //
    PS::F64vec acc1; //
    void clear(){
        acc0 = acc1 = 0.0;
    }
};

//! Hermite paramter class
class HermiteParams{
public:
    PS::F64 eps_sq;
    ChangeOver changeover;

    void print(std::ostream & _fout) const{
        _fout<<" eps_sq="<<r_in_;
        changeover.print(_fout);
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"Eps^2";
        changeover.printColumnTitle(_fout, _width);

    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<eps_sq_;
        changeover.printColumn(_fout, _width);
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


//!Hermite integrator class
template <class Tgptcl>
class HermiteIntegrator{
private:
    PS::ReallocatableArray<PtclPred> pred_;  // predictor ptcl
    PS::ReallocatableArray<PtclForce> force_;  // force for ptcl
    PS::ReallocatableArray<PS::S32> adr_dt_sorted_;  // sorting list of time_next
    PS::ReallocatableArray<PS::F64> time_next_;      // next integrated time list

    PS::ReallocatableArray<PtclH4> ptcl_;   // first c.m.; second single
    PS::ReallocatableArray<Tgptcl*> ptcl_ptr_; // original ptcl address
    PS::ReallocatableArray<Tneighbor> nb_; // neighbor information 

    PS::S32 n_act_;        /// active particle number
    Tgptcl pcm_; // c.m. of the cluster

public:
    BlockTimeStep time_step;  // time step manager class

private:
    //! Center-of-mass frame shift for particles
    /*! Shift positions and velocities of particles in #p from their original frame to their center-of-mass frame\n
      Notice the center-of-mass position and velocity use values from #cm
    */
    void calcAndShiftToCM(Tgptcl &pcm, PtclH4* ptcl, const PS::S32 n) {
        pcm.mass = 0.0;
        pcm.pos = pcm.vel = 0.0;
        for (PS::S32 i=0; i<n; i++) {
            PS::F64 mi=ptcl[i].mass;
            pcm.mass += mi;
            pcm.pos += mi*ptcl[i].pos;
            pcm.vel += mi*ptcl[i].vel;
        }
        PS::F64 invm = 1.0/pcm.mass;
        pcm.pos *= invm;
        pcm.vel *= invm;
        
        for (PS::S32 i=0;i<n;i++) {
            // std::cerr<<std::setprecision(14)<<"cm sbore i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
            ptcl[i].pos -= pcm.pos;
            ptcl[i].vel -= pcm.vel;
            //#ifdef HERMITE_DEBUG
            // std::cerr<<"cm shift i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
//#endif
        }
    }

    //! shift ptcl from c.m. frame to original frame
    /*!
     */
    void shiftToOrigin(const ParticleBase &pcm, PtclH4* ptcl, const PS::S32 n) {
#ifdef HERMITE_DEBUG
        PS::F64 m=0;
        for (PS::S32 i=0;i<n;i++) m+= ptcl[i].mass;
        assert(abs(m-pcm.mass)<1e-10);
        //PS::F64vec pos=0,vel=0;
        //for (PS::S32 i=0;i<n;i++) {
        //    pos += ptcl[i].mass*ptcl[i].pos;
        //    vel += ptcl[i].mass*ptcl[i].vel;
        //}
        //pos /=m;
        //vel /=m;
        //std::cerr<<std::setprecision(14)<<"cm pos "<<pos<<" vel "<<vel<<std::endl;
#endif
        for (PS::S32 i=0;i<n;i++) {
            //std::cerr<<"cm bbck i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
            ptcl[i].pos += pcm.pos;
            ptcl[i].vel += pcm.vel;
            //std::cerr<<"cm back i "<<i<<" pos "<<ptcl[i].pos<<" vel "<<ptcl[i].vel<<std::endl;
        }
    }

    //! Calculate 2nd order time step for list of particles
    /*! Calculate 2nd order time step for list of particles
      @param[in]: _ptcl: particle array contain dt
      @param[in]: _force: force array contain acc0, acc1
      @parma[in]: _list: active particle index list to calc time step
      @param[in]: _n_act: active particle number
      @param[in]: _step: time step manager class
      \return fail_flag: if the new dt < dt_min, return true (failure)
     */
    template <class Tptcl>
    bool calcDt2ndAct(Tptcl _ptcl[],
                      const PtclForce _force[],
                      const PS::S32 _list[],
                      const PS::S32 _n_act,
                      const BlockTimeStep& _step) {
        for(PS::S32 i=0; i<_n_act; i++){
            const PS::S32 k = _list[i];
            const PS::F64vec acc0 = _force[k].acc0;
            const PS::F64vec acc1 = _force[k].acc1;
            _ptcl[k].dt = _step.calcDt2nd(acc0, acc1);
            if(_ptcl[k].dt<0) {
                std::cerr<<"Ptcl: index="<<k<<" id="<<_ptcl[k].id<<std::endl;
                return true;
            }
        }
        return false;
    }

    template<cclass Tptcl>
    inline bool CalcAcc0Acc1R2Cutoff(PtclForce & _fi,
                                     PS::F64 & _r2,
                                     const Tptcl& _pi,
                                     const Tptcl& _pj,
                                     const HermiteParams& _pars) {
        const PS::F64vec dR = _pi.pos - _pj.pos;
        r2 = dR*dR;
        const PS::F64 dr2_eps = r2 + pars.eps_sq;
        const PS::F64 r_out = _pars.changeover.getRout();
        if(dr2_eps <= r_out*r_out){
            const PS::F64vec dV = _pi.vel - _pj.vel;
            const PS::F64 drdv = dR * dV;
            const PS::F64 dr_eps = sqrt(dr2_eps);
            const PS::F64 dri = 1.0/dr_eps;
#ifdef HARD_DEBUG
#ifdef HARD_DEBUG_DUMP
            if (std::isnan(dri)) {
                std::cerr<<"Error: relative position is zero!\n";
                return true;
            }
#else
            assert(!std::isnan(dri));
#endif        
#endif
            //const PS::F64 R = 1.0 / sqrt(r2_eps);
            const PS::F64 dri2 = dri*dri;
            const PS::F64 dri3 = dri2*dri;
            const PS::F64 k = pars.changeover.calcAcc0W(r_eps);
            const PS::F64 dk = pars.changeover.calcAcc1W(r_eps) * pars.changeover.getNorm();
            const PS::F64 mri3 = _pj.mass*dri3;
            const PS::F64vec F0 = -mri3*k*dR;
            const PS::F64vec F1 = -mri3*k*dV - 3.0*drdv*(dri2*F0 + mri3*dk*dri*dR);
            fi.acc0 += F0;
            fi.acc1 += F1;

        }
        return false;
    }

    //! calculate one particle f and fdot from all j particles 
    /*! Calculate f and fdot from all j particles, and return neighbor map
       @param[out] _fi:  acc and jerk for i particle
       @param[out] _nbi: neighbor information collector for i particle
       @param[in] _pi: i particle
       @param[in] _ptcl: j particle array
       @param[in] _n_tot: total j particle number
       @param[in] _n_group: number of groups
       @param[in] _int_pars: integration parameters (include changeover functions)
       @param[in] _Aint: ARC integrator class
       \return fail flag (true for failure)
     */
    template <class Tpi, class Tptcl, class ARCint>
    inline bool calcAcc0Acc1One(PtclForce& _fi,
                                Tneighbor& _nbi,
                                const Tpi&  _pi,
                                const Tptcl _ptcl[],
                                const PS::S32 _n_tot,
                                const HermiteParams& _int_pars,
                                const ARCint* _Aint = NULL) {
        if(_Aint!=NULL) n_group = _Aint->getNGroups();
        for(PS::S32 j=0; j<n_group; j++) {
            if(pi.id==_ptcl[j].id) continue;
            if(_Aint->getMask(j)) continue;
#ifdef HERMITE_DEBUG
            PS::F64 mcmcheck =0.0;
#endif
            const auto* pj = _Aint->getGroupPtcl(j);
            PS::F64 sdj = 1.0/_Aint->getSlowDown(j);
            PS::F64vec vcmsdj = (1.0-sdj)*_ptcl[j].vel;

            for(PS::S32 k=0; k<_Aint->getGroupN(j); k++) {
                PS::F64 r2 = 0.0;

                // slowdown
                ParticleBase pksd;
                pksd.mass = pj[k].mass;
                pksd.pos = pj[k].pos;
                pksd.vel = pj[k].vel*sdj + vcmsdj;

                bool fail_flag = CalcAcc0Acc1R2Cutoff(_force, r2, _pi, pksd, _int_pars);

#ifdef HERMITE_DEBUG
                if (fail_flag) {
                    std::cerr<<"Particle i:\n";
                    _pi.print(std::cerr);
                    std::cerr<<"Particle j:\n";
                    pj[k].print(std::cerr);
                    return true;
                }
                mcmcheck += pj[k].mass;
#endif

                // neighbor information
                PS::F64 rs=std::max(_pi.r_search, _pj[k].r_search);
                if(r2<=rs*rs) _nbi.setNBFlag(j);
                _nbi.updateRmin(j, r2);
                _nbi.updateMassMin(pksd.mass);
            }

#ifdef HERMITE_DEBUG
#ifdef HERMITE_DEBUG_DUMP
            if(abs(mcmcheck-_ptcl[j].mass)>1e-10) {
                std::cerr<<"Error: c.m. mass ("<<_ptcl[j].mass<<") not match group membe mass sum ("<<mcmcheck<<") diff="<<abs(mcmcheck-_ptcl[j].mass)<<std::endl;
                std::cerr<<"Particle i:\n";
                _pi.print(std::cerr);
                std::cerr<<"Particle j:\n";
                _ptcl[j].print(std::cerr);
                return true;
            }
            if(mcmcheck<=0.0) {
                std::cerr<<"Error: group total mass is zero"<<std::endl;
                std::cerr<<"Particle i:\n";
                _pi.print(std::cerr);
                std::cerr<<"Particle j:\n";
                _ptcl[j].print(std::cerr);
                return true;
            }
#else
            assert(abs(mcmcheck-_ptcl[j].mass)<1e-10);
            assert(mcmcheck>0.0);
#endif
#endif                    
        }
        for(PS::S32 j=n_group; j<_n_tot; j++){
            if(_pi.id==_ptcl[j].id) continue;

            PS::F64 r2 = 0.0;
            bool fail_flag = CalcAcc0Acc1R2Cutoff(_fi, r2, _pi, _ptcl[j]. _int_pars);
#ifdef HERMITE_DEBUG
            if (fail_flag) {
                std::cerr<<"Particle i:\n";
                _pi.print(std::cerr);
                std::cerr<<"Particle j:\n";
                _ptcl[j].print(std::cerr);
                return true;
            }
#endif            
            // if(r2 < ((ptcl[adr].r_merge + ptcl[j].r_merge)*(ptcl[adr].r_merge + ptcl[j].r_merge)) && ptcl[j].mass > 0.0){
            //     merge_pair.push_back( std::make_pair(adr, j) );
            // }

            // determine neighbor
            PS::F64 rs = std::max(_pi.r_search, _ptcl[j].r_search);
            if(r2<=rs*rs) _nbi.setNBFlag(j);
            _nbi.updateRmin(j, r2);
            _nbi.updateMassMin(_ptcl[j].mass);
        }
        return false;
    }

    //! Calculate acc and jerk for active particles from all particles and update neighbor lists
    /*! calculate force and update neighbor lists
       @param[out] _force: acc and jerk array
       @param[out] _nb_info: neighbor information collector
       @param[in] _ptcl: particle list
       @param[in] _n_tot: total number of particles
       @param[in] _list_act: active i particle index in ptcl_
       @param[in] _n_act: active particle number
       @param[in] _int_pars: integration parameters (include changeover functions)
       @param[in] _Aint: ARC integrator class
       \return fail flag (true for failure)
     */
    template <class Tptcl, class ARCint>
    inline bool CalcAcc0Acc1Act(PtclForce _force[],
                                Tneighbor _nb_info[],
                                const Tptcl _ptcl[],
                                const PS::S32 _n_tot,
                                const PS::S32 _list_act[],
                                const PS::S32 _n_act,
                                const HermiteParams& _int_pars,
                                const ARCint* _Aint) {
        bool fail_flag = false;
        PS::S32 n_group = 0;
        if(_Aint!=NULL) n_group = _Aint->getNGroups();

        const PS::F64vec vzero = PS::F64vec(0.0);
        // PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > & merge_pair ){

        // active iparticles loop
        for(PS::S32 i=0; i<_n_act; i++){
            const PS::S32 iadr = _list_act[i];
            _force[iadr].clear();
            
            _nb_info[iadr].resetRmin();
            _nb_info[iadr].resetMassMin();
            _nb_info[iadr].resetNBFlag();
            
            // for group particle
            if (iadr<n_group) {
#ifdef HERMITE_DEBUG
                PS::F64 mcmcheck =0.0;
#endif
                const PS::S32 ni = _Aint->getGroupN(iadr);             // number of members
#ifdef HERMITE_DEBUG
                assert(ni>0);
#endif
                const auto* pi = _Aint->getGroupPtcl(iadr);
                const PS::F64 sdi = 1.0/_Aint->getSlowDown(iadr);      // slowdown factor
                const PS::F64vec vcmsdi = (1.0-sdi)*_ptcl[iadr].vel;   // slowdown velocity
//                PS::F64 r2min[_n_tot]={PS::LARGE_FLOAT};

                for (PS::S32 j=0; j<ni; j++) {
                    PtclForce fpj;
                    fpj.clear();

                    // slowdown particle
                    Ptcl pisd; 
                    pisd.mass = pi[j].mass;
                    pisd.pos = pi[j].pos;
                    pisd.vel = pi[j].vel*sdi + vcmsdi;
                    pisd.id = _ptcl[iadr].id;
                    pisd.r_search = pi[j].r_search;

                    fail_flag = calcAcc0Acc1One(fpj, _nb_info[iadr], pisd, _ptcl, _n_tot, _int_pars, _Aint);
                    // c.m. force
                    _force[iadr].acc0 += pi[j].mass*fp[j].acc0;
                    _force[iadr].acc1 += pi[j].mass*fp[j].acc1;
                    
#ifdef HERMITE_DEBUG
                    if (fail_flag) return true;
                    mcmcheck += pi[j].mass;
#endif
                }

#ifdef HERMITE_DEBUG
#ifdef HERMITE_DEBUG_DUMP
                if(abs(mcmcheck-_ptcl[j].mass)>1e-10) {
                    std::cerr<<"Error: c.m. mass ("<<_ptcl[iadr].mass<<") not match group membe mass sum ("<<mcmcheck<<") diff="<<abs(mcmcheck-_ptcl[iadr].mass)<<std::endl;
                    std::cerr<<"Particle:\n";
                    _ptcl[iadr].print(std::cerr);
                    return true;
                }
                if(mcmcheck<=0.0) {
                    std::cerr<<"Error: group total mass is zero"<<std::endl;
                    std::cerr<<"Particle:\n";
                    _ptcl[iadr].print(std::cerr);
                    return true;
                }
#else
                assert(abs(mcmcheck-_ptcl[iadr].mass)<1e-10);
                assert(mcmcheck>0.0);
                assert(_ptcl[iadr].mass>0);
#endif
#endif                    
                // c.m. force
                _force[iadr].acc0 /= _ptcl[iadr].mass;
                _force[iadr].acc1 /= _ptcl[iadr].mass;

            }
            else {
                fail_flag=calcAcc0Acc1One(_force[iadr], _nb_info[iadr], _ptcl[iadr],  _ptcl, _n_tot, _int_pars, _Aint);

#ifdef HERMITE_DEBUG
                if (fail_flag) return true;
#endif
            }

            // Update neighbors
            _nb_info[iadr].replaceNBListFromFlag();

#ifdef HERMITE_DEBUG
#ifdef HERMITE_DEBUG_DUMP
            if(_nb_info[iadr].getNBN()>_n_tot) {
                std::cerr<<"Error: neighbor list overflow (nb="<<nb<<" >n_tot="<<_n_tot<<"\n";
                std::cerr<<"Particle:\n"
                _ptcl[iadr].print(std::cerr);
                return true;
            }
#else
            assert(_nb_info[iadr].getNBN()<=_n_tot);
#endif
#endif
        }
        return fail_flag;
    }
    
    class SortAdr{
    public:
        PS::F64 * time;
        SortAdr(PS::F64 * _time): time(_time){}
        bool operator() (const PS::S32 & left, const PS::S32 & right) const {
            return time[left] < time[right];
        }
    };

    //! sort time step array and select active particles
    /*! also return the sub list of active group index in adr_dt_sorted
     */
    void SortAndSelectIp(PS::S32 adr_dt_sorted[],
                         PS::F64 time_next[],
                         PS::S32 & n_act,
                         const PS::S32 n_tot,
                         PS::S32 list_group_act[],
                         PS::S32 & n_group_act,
                         const PS::F64 n_group){
        // const PS::S32 n_tot = time_next.size();
        //std::cerr<<"before sort"<<std::endl;
        /*
          for(PS::S32 ip=0; ip<ni_old; ip++){
          const PS::S32 adr = adr_dt_sorted[ip];
          time_next[adr] += ptcl[adr].dt; // n_act only
          }
        */
        std::sort(adr_dt_sorted, adr_dt_sorted+n_act, SortAdr(time_next));

        const PS::F64 time_ref = time_next[adr_dt_sorted[0]];
        n_group_act = 0;
        for(n_act=1; n_act<n_tot; n_act++){
            if(time_ref < time_next[adr_dt_sorted[n_act]]) {
                break;
            }
            if(n_act<n_group) list_group_act[n_group_act++] = n_act;
        }
    }

    void SortAndSelectIp(PS::S32 adr_dt_sorted[],
                         PS::F64 time_next[],
                         PS::S32 & n_act,
                         const PS::S32 n_tot){
        // const PS::S32 n_tot = time_next.size();
        //std::cerr<<"before sort"<<std::endl;
        /*
          for(PS::S32 ip=0; ip<ni_old; ip++){
          const PS::S32 adr = adr_dt_sorted[ip];
          time_next[adr] += ptcl[adr].dt; // n_act only
          }
        */
        std::sort(adr_dt_sorted, adr_dt_sorted+std::min(n_act,n_tot), SortAdr(time_next));

        const PS::F64 time_ref = time_next[adr_dt_sorted[0]];
        for(n_act=1; n_act<n_tot; n_act++){
            if(time_ref < time_next[adr_dt_sorted[n_act]]) {
                break;
            }
        }
    }
    
    void PredictAll(PtclPred pred[],
                    const PtclH4 ptcl[],
                    const PS::S32 n_tot,
                    const PS::F64 time_next){
        static thread_local const PS::F64 inv3 = 1.0 / 3.0;
        for(PS::S32 i=0; i<n_tot; i++){
            const PS::F64 dt = time_next - ptcl[i].time;
            pred[i].pos = ptcl[i].pos + dt*(ptcl[i].vel  + 0.5*dt*(ptcl[i].acc0 + inv3*dt*ptcl[i].acc1));
            pred[i].vel = ptcl[i].vel + dt*(ptcl[i].acc0 + 0.5*dt*ptcl[i].acc1);
            // pred[i].r_out = ptcl[i].r_out;
            pred[i].mass = ptcl[i].mass;
        }
    /*
      if(PS::Comm::getRank() == 0){
      for(PS::S32 i=0; i<n_tot; i++){
      std::cerr<<"pred[i].pos= "<<pred[i].pos
      <<" ptcl [i].pos="<<ptcl [i].pos<<std::endl;
      std::cerr<<"pred[i].vel= "<<pred[i].vel
      <<" ptcl [i].vel="<<ptcl [i].vel<<std::endl;
      }
      }
    */
    }

    //! correct ptcl and calculate step 
    /*! Correct ptcl and calculate next time step
      @param[in,out] _ptcl: particle array store the dt, time, previous acc0, acc1, the acc0, acc1 are updated with new ones
      @param[in] _force: predicted forces
      @param[in] _list: active particle index array
      @param[in] _n_act: number of active particles
      @param[in] _step: time step manager class
      \return fail_flag: if the time step < dt_min, return true (failure)
     */
    bool CorrectAndCalcDt4thAct(PtclH4 _ptcl[],
                                const PtclForce _force[],
                                const PS::S32 _list[], 
                                const PS::S32 _n_act,
                                const BlockTimeStep& _step) {

        static thread_local const PS::F64 inv3 = 1.0 / 3.0;
        for(PS::S32 i=0; i<_n_act; i++){
            const PS::S32 k = _list[i];
            PtclH4*          pti = &_ptcl[k];
            const PtclForce* fpi = &_force[k];

            const PS::F64 dt = pti->dt;
            const PS::F64 h = 0.5 * dt;
            const PS::F64 hinv = 2.0 / dt;
            const PS::F64vec A0p = (fpi->acc0 + pti->acc0);
            const PS::F64vec A0m = (fpi->acc0 - pti->acc0);
            const PS::F64vec A1p = (fpi->acc1 + pti->acc1)*h;
            const PS::F64vec A1m = (fpi->acc1 - pti->acc1)*h;

            const PS::F64vec vel_new = pti->vel + h*( A0p - inv3*A1m );
            pti->pos += h*( (pti->vel + vel_new) + h*(-inv3*A0m));
            pti->vel = vel_new;

            pti->acc0 = fpi->acc0;
            pti->acc1 = fpi->acc1;
            pti->time += dt;

            const PS::F64vec acc3 = (1.5*hinv*hinv*hinv) * (A1p - A0m);
            const PS::F64vec acc2 = (0.5*hinv*hinv) * A1m + h*acc3;
#ifdef HERMITE_DEBUG
            // for debug
            assert(!std::isnan(pti->pos[0]));
            assert(!std::isnan(pti->pos[1]));
            assert(!std::isnan(pti->pos[2]));
            assert(!std::isnan(pti->vel[0]));
            assert(!std::isnan(pti->vel[1]));
            assert(!std::isnan(pti->vel[2]));
            //assert(pti->pos[0]==pti->pos[0]);
            //assert(pti->pos[1]==pti->pos[1]);
            //assert(pti->pos[2]==pti->pos[2]);
            //assert(pti->vel[0]==pti->vel[0]);
            //assert(pti->vel[1]==pti->vel[1]);
            //assert(pti->vel[2]==pti->vel[2]);
#ifdef HERMITE_DEBUG_ACC
            pti->acc2 = acc2;
            pti->acc3 = acc3;
#endif
#endif
            const PS::F64 dt_old = pti->dt;
            pti->dt = _step.calcDt4th(pti->acc0, pti->acc1, acc2, acc3);
//            pti->dt = dt_old*2 < pti->dt ?  dt_old*2 : pti->dt;
//            if(int(pti->time/pti->dt)*pti->dt<pti->time) pti->dt *= 0.5;

#ifdef HERMITE_DEBUG
#ifdef HERMITE_DEBUG_DUMP
            if (dt_old! <= 0.0|| pti->dt <= 0.0) {
                pti->print(std::cerr);
                return true;
            }
#else
            assert(dt_old != 0.0);
            assert(pti->dt != 0.0);
#endif
#endif
        }
        return false;
    }

    //! modify list based on the trace table
    /*! modify list based on the trace table
      @param[out] _list: list to modify
      @param[in] _n_list: list size
      @param[in] _trace: moving table for _list, <0: no change, ==_del_key: remove, others: new position in _list
      @param[in] _del_key: indicator for removing
      \return new list size
     */
    PS::S32 modifyList(PS::S32* _list, 
                       const PS::S32 _n_list,
                       const PS::S32* _trace,
                       const PS::S32 _del_key) {
        if(_n_list==0) return 0;
#ifdef HERMITE_DEBUG
        assert(_n_list>0);
#endif
        PS::S32 list_size=_n_list;
        // moving offset
        PS::S32 imove=0;
        // index 
        PS::S32 i=0;
        while(i<list_size) {
            PS::S32 iadr=_list[i];
            PS::S32 itrace=_trace[iadr];
            // delete case
            if(itrace==_del_key) {
                imove++;
                list_size--;
            }
            else {
                // replace case
                if(itrace>=0) _list[i] = itrace;
                // move to next
                i++;
            }
            if(i>=list_size) break;
#ifdef HERMITE_DEBUG
            assert(i+imove<_n_list);
#endif
            _list[i] = _list[i+imove];
        }

        return list_size;
    }

    //! get whether two pair is leaving each other or coming
    /*!
      @param[in] _i: first index in ptcl_
      @param[in] _j: second index in ptcl_
      \return true: go away, otherwise income
    */
    bool getDirection(const PS::S32 _i, const PS::S32 _j) {
#ifdef HERMITE_DEBUG
        assert(_i>=0&&_i<ptcl_.size());
        assert(_j>=0&&_j<ptcl_.size());
#endif                
        PS::F64 xv2=(ptcl_[_i].pos-ptcl_[_j].pos)*(ptcl_[_i].vel-ptcl_[_j].vel);
        return (xv2>0);
    }


    //! check the pair with distance below r_crit for ptcl in adr_dt_sorted_
    /*! First check nearest neighbor distance r_min
        If r_min<r_crit, check the direction, if income, accept as group
      @param[out] _new_group_member_index: new group member index in ptcl_
      @param[out] _new_group_member_adr: new group member address, echo two neighbor is one group
      @param[out] _new_group_offset: new group offset, last one show total number of members.
      @param[in] _r_crit2: group distance criterion square
      @param[in] _mask_list: a list contain the index that would not be in the new group
      @param[in] _n_mask: number of masked indices
      @param[in] _Aint: ARC integrator class
      \return number of new groups
     */
    template <class ARCint>
    PS::S32 checkNewGroup(PS::S32 _new_group_member_index[], 
                          Tgptcl* _new_group_member_adr[], 
                          PS::S32 _new_group_offset[], 
                          const PS::F64 _r_crit2, 
                          const PS::S32* _mask_list,
                          const PS::S32 _n_mask,
                          const ARCint* _Aint) {
#ifdef HERMITE_DEBUG
        assert(nb_info_.size()==ptcl_.size());
#endif
        PS::S32 n_group = 0;
        if (_Aint!=NULL) n_group = _Aint->getNGroups();
        PS::S32 n_new_group=0, offset=0;

        // here ptcl_.size() is used since any ptcl index can be checked!
        PS::S32 used_mask[ptcl_.size()];
        for (PS::S32 k=0; k<ptcl_.size(); k++) used_mask[k] = -1;
        // for the index not allowed to be grouped
        for (PS::S32 k=0; k<_n_mask; k++) used_mask[_mask_list[k]] = -2; 
        

        // check adr_dt_sorted_ list (avoid suppressed ptcl)
        PS::S32 n_check = adr_dt_sorted_.size();
        for (PS::S32 k=0; k<n_check; k++) {
            const PS::S32 i = adr_dt_sorted_[k];
            if(nb_info_[i].r_min2<_r_crit2) {
                const PS::S32 j = nb_info_[i].r_min_index;
#ifdef HERMITE_DEBUG
                assert(j<ptcl_.size());
                assert(ptcl_[i].mass>0.0);
#endif                
                if(j<0) continue; 

                // avoid masked member
                if(used_mask[i]==-2||used_mask[j]==-2) continue;

                if(!(used_mask[i]>=0 && used_mask[j]>=0)) { // avoid double count
                    bool out_flag=getDirection(i, j);
                    if(!out_flag) {
                        PS::F64 sdi=0.0, sdj=0.0, frsi=0.0, frsj=0.0;
                        if(i<n_group) {
                            sdi = _Aint->getSlowDown(i);
                            frsi = _Aint->getFratioSq(i);
                        }
                        if(j<n_group) {
                            sdj = _Aint->getSlowDown(j);
                            frsj = _Aint->getFratioSq(j);
                        }
                        if(sdi>1.0&&sdj>1.0) continue;
                        if(sdi>1.0&&sdj==0.0) continue;
                        if(sdi==0.0&&sdj>1.0) continue;
                        // to avoid extremely long time integration, for very large slowdown factor, no merge group
                        // current slowdown factor 100.0 and fratioSq 1e-8 is experimental value
                        //if ((sdi>100.0&&fpj<1e-2)||(sdj>100.0&&fpi<1e-2)) continue;
                        //if ((sdi>100.0&&sdi*sdi*fpj<1.0)||(sdj>100.0&&sdj*sdj*fpi<1.0)) continue;

                        // check the tidal effect of strong bound binary
                        PS::S32 i_bin_strong = -1;  // strong bound binary index (-1 means no slowdown exist)
                        PS::F64 fratio_weak_sq;     // weak binary force ratio square
                        if (sdi>1.0) {
                            i_bin_strong = i;
                            fratio_weak_sq = frsj;
                        }
                        if (sdj>1.0) {
                            i_bin_strong = j;
                            fratio_weak_sq = frsi;
                        }
                        if(i_bin_strong>=0) {
                            PS::F64 apo = _Aint->bininfo[i_bin_strong].semi*(1.0 + _Aint->bininfo[i_bin_strong].ecc);
                            // tidal effect estimated by apo (strong) / rij
                            PS::F64 ftid_strong_sq = apo*apo/nb_info_[i].r_min2;
                            if (ftid_strong_sq*fratio_weak_sq<1e-8) continue;
                        }

#ifdef ADJUST_GROUP_DEBUG
                        std::cout<<"Find new group      index      slowdown      fratio_sq       apo      ftid_sq \n"
                                 <<"i1              "
                                 <<std::setw(8)<<i
                                 <<std::setw(16)<<sdi
                                 <<std::setw(16)<<frsi;
                        if(i<n_group)  {
                            PS::F64 apo= _Aint->bininfo[i].semi*(1.0+_Aint->bininfo[i].ecc);
                            std::cout<<std::setw(16)<<apo
                                     <<std::setw(16)<<apo*apo/nb_info_[i].r_min2;
                        }
                        std::cout<<"\ni2              "
                                 <<std::setw(8)<<j
                                 <<std::setw(16)<<sdj
                                 <<std::setw(16)<<frsj;
                        if(j<n_group)  {
                            PS::F64 apo= _Aint->bininfo[j].semi*(1.0+_Aint->bininfo[j].ecc);
                            std::cout<<std::setw(16)<<apo
                                     <<std::setw(16)<<apo*apo/nb_info_[i].r_min2;
                        }
                        std::cout<<std::endl;
                        std::cout<<"I_binary_strong: "<<i_bin_strong<<std::endl;
#endif
                        PS::S32 insert_group=-1, insert_index=-1;
                        if(used_mask[i]>=0) {
                            insert_group = used_mask[i];
                            insert_index = j;
                        }
                        else if(used_mask[j]>=0) {
                            insert_group = used_mask[j];
                            insert_index = i;
                        }
                        // the case of merging group
                        if(insert_group>=0) {
                            // shift the first index in last group to last 
                            PS::S32 last_group_offset = _new_group_offset[n_new_group-1];
                            _new_group_member_index[offset] = _new_group_member_index[last_group_offset];
                            _new_group_member_adr  [offset] = _new_group_member_adr  [last_group_offset];
                            offset++;

                            // in the case insert_group is not the last group
                            if(insert_group<n_new_group-1) {
                                // shift first index in each group to the end to allow to insert j in i_group
                                for (PS::S32 k=n_new_group-1; k>insert_group; k--) {
                                    PS::S32 k_group_offset = _new_group_offset[k];
                                    PS::S32 k0_group_offset = _new_group_offset[k-1];
                                    _new_group_member_index[k_group_offset] = _new_group_member_index[k0_group_offset];
                                    _new_group_member_adr  [k_group_offset] = _new_group_member_adr  [k0_group_offset];
                                    _new_group_offset[k]++;
                                }
                            }
                            // replace the first position of insert_group with insert ptcl
                            PS::S32 insert_group_offset=_new_group_offset[insert_group];
                            _new_group_member_index[insert_group_offset] = insert_index;
                            _new_group_member_adr[insert_group_offset] = ptcl_ptr_[insert_index];
                            used_mask[insert_index] = insert_group;
                        }
                        else {   // new group case
                            _new_group_offset[n_new_group] = offset;
                            _new_group_member_index[offset] = i;
                            _new_group_member_adr[offset] = ptcl_ptr_[i];
                            used_mask[i] = n_new_group;
                            offset++;

                            _new_group_member_index[offset] = j;
                            _new_group_member_adr[offset] = ptcl_ptr_[j];
                            used_mask[j] = n_new_group;
                            offset++;

                            n_new_group++;
                        }
                    }
                }
            }
        }
        // for total number of members
        _new_group_offset[n_new_group] = offset;
#ifdef HERMITE_DEBUG
        assert(offset<=ptcl_.size());
#endif
        return n_new_group;
    }
    
    void moveCM(const PS::F64 dt) {
        pcm_.pos += pcm_.vel * dt;
        //std::cerr<<"pcm.pos "<<pcm_.pos<<" pcm.vel "<<pcm_.vel<<" dt "<<dt<<std::endl;
    }

    //!shift ptcl back to original frame
    void shiftToOrigin() {
        shiftToOrigin(pcm_,ptcl_.getPointer(),ptcl_.size());
    }

    //!shift ptcl to c.m. frame and save c.m.
    void shiftToCM() {
        calcAndShiftToCM(pcm_, ptcl_.getPointer(), ptcl_.size());
    }

    //! Reserve memory
    /*! To reserve memory for use
      @param[in] _n: maximum particle number to reserve
     */
    void reserveMem(const PS::S32 _n) {
        ptcl_.reserve(_n);
        ptcl_ptr_.reserve(_n);
        nb_info_.reserve(_n);
        pred_.reserve(_n);
        force_.reserve(_n);
        adr_dt_sorted_.reserve(_n);
        time_next_.reserve(_n);
    }

    ////! Copy particle data to local array
    ///* @param[in] _ptcl: particle data array
    //   @param[in] _n_ptcl: number of particle need to be added
    //   @param[in] _ptcl_list: particle address in _ptcl
    // */
    //template <class Tptcl>
    //void addInitPtcl(Tptcl * _ptcl, 
    //             const PS::S32 _n_ptcl, 
    //             const PS::S32* _ptcl_list) {
    //    for (PS::S32 i=0; i<_n_ptcl; i++) {
    //        ptcl_.increaseSize(1);
    //        ptcl_.back().DataCopy(_ptcl[_ptcl_list[i]]);
    //        //ptcl_.pushBackNoCheck(ptcl_org[ptcl_list[i]]);
    //    }
    //}
    // 
    ////! Copy particle data to local array
    ///* @param[in] _ptcl: particle data array
    //   @param[in] _n_ptcl: number of particles need to be added
    // */
    //template <class Tptcl>
    //void addInitPtcl(Tptcl * _ptcl, 
    //             const PS::S32 _n_ptcl) {
    //    for (PS::S32 i=0; i<_n_ptcl; i++) {
    //        ptcl_.increaseSize(1);
    //        ptcl_.back().DataCopy(_ptcl[i]);
    //    }
    //}

    //! Insert one particle
    /*! insert one particle to local array at position i, notice i is added in the front of adr_dt_sorted_.
      @param[in] _i: insert position
      @param[in] _ptcl: particle to be inserted
      @param[in] _n_list_correct: number of ptcl (count from first) need to correct the neighbor list
      @param[in] _time_sys: new ptcl time
     */
    template <class Tptcl>
    void insertOnePtcl(const PS::S32 _i,
                       Tptcl &_ptcl,
                       const PS::S32 _n_list_correct,
                       const PS::F64 _time_sys) {
#ifdef HERMITE_DEBUG
        assert(_i<=ptcl_.size());
        assert(ptcl_.size()+1<=n_nb_off_);
        assert(ptcl_.capacity()>=ptcl_.size()+1);
#endif 
        // increase all data size by 1
        ptcl_.increaseSize(1);
        ptcl_ptr_.increaseSize(1);
        pred_.increaseSize(1);
        force_.increaseSize(1);
        adr_dt_sorted_.increaseSize(1);
        time_next_.increaseSize(1);
        nb_info_.increaseSize(1);
        nb_info_.back().reserveNBFlag(ptcl_.capacity());
        nb_info_.back().reserveNBList(ptcl_.capacity());

#ifdef HERMITE_DEBUG
        assert(ptcl_.size()==pred_.size());
        assert(ptcl_.size()==force_.size());
        assert(ptcl_.size()==time_next_.size());
        assert(ptcl_.size()==nb_info_.size());
#endif

        // last just add
        if (_i==ptcl_.size()-1) {
            ptcl_.back().DataCopy(_ptcl);
            ptcl_.back().acc0 = ptcl_.back().acc1 = force_.back().acc0 = force_.back().acc1 = 0.0;
            ptcl_.back().time = _time_sys;
            ptcl_.back().dt = 0.0;
            ptcl_ptr_.back()=(Tptcl*)&_ptcl;

            // shift adr_dt_sorted
            for (PS::S32 i=ptcl_.size()-1; i>0; i--) {
                adr_dt_sorted_[i] = adr_dt_sorted_[i-1];
            }
            adr_dt_sorted_[0] = _i;
        }
        // shift _i to last and add
        else {
            PS::S32 inew = ptcl_.size()-1;
#ifdef HERMITE_DEBUG
            assert(&(ptcl_.back())==&(ptcl_[inew]));
#endif
            ptcl_.back() = ptcl_[_i];
            ptcl_ptr_.back() = ptcl_ptr_[_i];
            ptcl_[_i].DataCopy(_ptcl);
            ptcl_[_i].acc0 = ptcl_[_i].acc1 = force_[_i].acc0 = force_[_i].acc1 = 0.0;
            ptcl_[_i].time = _time_sys;
            ptcl_[_i].dt = 0.0;
            ptcl_ptr_[_i]=(Tptcl*)&_ptcl;
            
            pred_.back()       = pred_[_i];
            force_.back()      = force_[_i];
            time_next_.back()  = time_next_[_i];
            nb_info_.back()    = nb_info_[_i];

            // find sorted index and change
            for (PS::S32 i=inew; i>0; i--) {
                adr_dt_sorted_[i] = adr_dt_sorted_[i-1];
                if(adr_dt_sorted_[i]==_i) {
                    adr_dt_sorted_[i] = inew;
                }
            }
            adr_dt_sorted_[0] = _i;
            n_act_++;

            // correct neighbor information
            for (PS::S32 i=0; i<_n_list_correct; i++) 
                nb_info_[i].replaceIndex(_i, inew, true);
            for (PS::S32 i=_n_list_correct; i<inew; i++) 
                nb_info_[i].replaceIndex(_i, inew, false);
        }
    }

    //! Replace one particle
    /*! Replace one particle of _i and shift the index in adr_dt_softed_ to front
      @param[in] _i: modified index in Hint.ptcl_
      @param[in] _ptcl: new ptcl information
      @param[in] _time_sys: new ptcl time
    */
    template <class Tptcl>
    void modOnePtcl(const PS::S32 _i,
                    const Tptcl &_ptcl,
                    const PS::F64 _time_sys ) {
#ifdef HERMITE_DEBUG
        assert(_i<=ptcl_.size());
#endif 
        ptcl_[_i].DataCopy(_ptcl);
        ptcl_[_i].acc0 = ptcl_[_i].acc1 = force_[_i].acc0 = force_[_i].acc1 = 0;
        ptcl_[_i].time = _time_sys;
        ptcl_[_i].dt = 0.0;
        ptcl_ptr_[_i] = (Tptcl*)&_ptcl;

        // find index and shift to front
        const PS::S32 adr_size=adr_dt_sorted_.size();

        PS::S32 i=0;
        while(adr_dt_sorted_[i]!=_i&&i<adr_size) i++;
        if(i>=adr_size) adr_dt_sorted_.increaseSize(1);

        // shift
        for (PS::S32 k=i; k>0; k--) 
            adr_dt_sorted_[k] = adr_dt_sorted_[k-1];
        adr_dt_sorted_[0] = _i;
        if(i>=n_act_) n_act_++;
    }


    //! Add a list of particles at the end
    /*! Add a list of particles at the end of Hint.ptcl_, notice all new particles have index in front of adr_dt_sorted_.
      @param[in] _ptcl_origin: particle array to be added
      @param[in] _list_origin: adding particle index in _ptcl_origin
      @param[in] _n_list: number of adding particles
      @param[in] _n_nb_correct: number of ptcl (count from first): need to correct the neighbor list
      @param[in] _time_sys: new ptcl time
      @param[in] _adr_dt_front_flag: add new ptcl index in front of adr_dt_sorted_ (true) or end (false)
     */
    template <class Tptcl>
    void addPtclList(Tptcl* _ptcl_origin,
                     const PS::S32* _list_origin,
                     const PS::S32 _n_list,
                     const PS::S32 _n_nb_correct,
                     const PS::F64 _time_sys,
                     const bool _adr_dt_front_flag) {

        // quite if no new
        if(_n_list==0) return;
#ifdef HERMITE_DEBUG
        assert(_n_list>0);
#endif

        // original size
        const PS::S32 n_org = ptcl_.size();
        const PS::S32 n_adr_dt_org = adr_dt_sorted_.size();

#ifdef HARD_DEBUG
        assert(ptcl_.capacity()>=n_org+_n_list);
#endif        
        // increase all data size by _n_list
        ptcl_.increaseSize(_n_list);
        ptcl_ptr_.increaseSize(_n_list);
        pred_.increaseSize(_n_list);
        force_.increaseSize(_n_list);
        adr_dt_sorted_.increaseSize(_n_list);
        time_next_.increaseSize(_n_list);
        nb_info_.increaseSize(_n_list);
        for (PS::S32 i=0; i<_n_list; i++) {
            nb_info_[i+n_org].reserveNBFlag(ptcl_.capacity());
            nb_info_[i+n_org].reserveNBList(ptcl_.capacity());
        }

#ifdef HERMITE_DEBUG
        assert(ptcl_.size()<ARRAY_ALLOW_LIMIT);
        assert(ptcl_.size()==ptcl_ptr_.size());
        assert(ptcl_.size()==pred_.size());
        assert(ptcl_.size()==force_.size());
        assert(ptcl_.size()==time_next_.size());
        assert(ptcl_.size()==nb_info_.size());
#endif

        // add ptcl in order
        if (_list_origin==NULL) {
            for (PS::S32 i=0; i<_n_list; i++) {
                const PS::S32 inew= n_org+i;
                ptcl_[inew].DataCopy(_ptcl_origin[i]);
                ptcl_[inew].time = _time_sys;
                ptcl_[inew].acc0 = ptcl_[inew].acc1 = force_[inew].acc0 = force_[inew].acc1 = 0.0;
                ptcl_[inew].dt = 0.0;
                ptcl_ptr_[inew] = (Tptcl*)&_ptcl_origin[i];
                //_single_index_origin[_n_single+i] = i;
            }
        }
        // add ptcl from list
        else {
            for (PS::S32 i=0; i<_n_list; i++) {
                const PS::S32 iadr= _list_origin[i];
                const PS::S32 inew= n_org+i;
                ptcl_[inew].DataCopy(_ptcl_origin[iadr]);
                ptcl_[inew].time = _time_sys;
                ptcl_[inew].acc0 = ptcl_[inew].acc1 = force_[inew].acc0 = force_[inew].acc1 = 0.0;
                ptcl_[inew].dt = 0.0;
                ptcl_ptr_[inew] = (Tptcl*)&_ptcl_origin[iadr];
                //_single_index_origin[_n_single+i] = iadr;
            }
        }
        // nb_list disp

        //_n_single += _n_list;
        
        // add new ptcl in front of adr_dt_sorted
        if(_adr_dt_front_flag) {
            // shift adr_dt_sorted
            for (PS::S32 i=n_adr_dt_org-1; i>=0; i--) 
                adr_dt_sorted_[i+_n_list] = adr_dt_sorted_[i];
            // add new ptcl to adr_dt_sorted
            for (PS::S32 i=0; i<_n_list; i++) 
                adr_dt_sorted_[i] = n_org+i;
        }
        // add at then end of adr_dt_sorted
        else {
            for (PS::S32 i=n_adr_dt_org; i<n_adr_dt_org+_n_list; i++) 
                adr_dt_sorted_[i] = i;
        }
        n_act_ += _n_list;

        // add new ptcl to neighbor list of c.m.
        for (PS::S32 i=0; i<_n_nb_correct; i++) { 
            // escape the suppressed case
            if(ptcl_[i].mass==0&&nb_info_[i].getNBN()==0) continue;
            for (PS::S32 j=n_org; j<n_org+_n_list; j++) 
                nb_info_[i].addNB(j);
        }
    }
                    
    //! Remove a list of particle
    /*! Remove particles from _list, keep first _n_correct in ptcl_ undeleted (but mass to zero) and update their neighbor lists
      @param[in] _list: removing particle index for ptcl_
      @param[in] _n_list: number of removing particles
      @param[in] _n_group: number of c.m. ptcl (count from first): ptcl will not be delected but only remove from adr_dt_sorted_, also the offset of the index between Hint.ptcl_ _single_index_origin
      @param[in] _n_nb_correct: number of ptcl (count from first) to correct neighbor list
     */
    void removePtclList(const PS::S32* _list,
                        const PS::S32 _n_list,
                        const PS::S32 _n_group,
                        const PS::S32 _n_nb_correct) {
        // quit if no deleted ones
        if(_n_list==0) return;
#ifdef HERMITE_DEBUG
        assert(_n_list>0);
#endif

        const PS::S32 n_org=ptcl_.size();
        const PS::S32 n_new = n_org - _n_list; // new ptcl number
        // moving trace (initial -1)
        PS::S32 trace[n_org];
        for (PS::S32 i=0; i<n_org; i++) trace[i] = -1;
        PS::S32 ilast = n_org-1;

        // create moving table (trace)
        // record the new position if the ptcl is moved, if deleted, set to ptcl_.size()
        PS::S32 n_decrease=0;
        for (PS::S32 i=0; i<_n_list; i++) {
            PS::S32 k=_list[i];
            // check no delete case
            if(k<_n_group) {
                trace[k]=n_org;
                // set suppressed group mass to zero
                ptcl_[k].mass = 0.0;
                // set status to -20 to identify the suppressed ptcl
                ptcl_[k].status = -20;
                // set neighbor to zero to avoid issue
                nb_list_n_[k] = 0;
                // set r_min_index to -1
                nb_info_[k].r_min_index = -1;
                continue;
            }

            PS::S32 idel = k;
            // set ilast >= idel for safety.
            ilast = std::max(ilast, idel);
            
            // check whether position is already moved, if so, check the moved new position and set current to delete
            if(trace[k]>=0) {
                idel=trace[k];
                trace[k]=n_org;
            }
#ifdef HERMITE_DEBUG
            assert(k<n_org);
            assert(idel<n_org);
#endif

            // check the last avaiable particle that can be moved to idel
            // If the ilast is already moved, check whether the new moved position trace[ilast] is before the current idel.
            // If trace[ilast] is after the current idel, update the trace[ilast] to idel and set trace[ilast] as new idel to check, until trace[ilast]=-1 or idel >= ilast
            while (idel>=0) {
                PS::S32 itrlast = -1;
                //cond: last is moved && ( moved pos before new n_ptcl || last is del ) &&  ilast > idel 
                while(trace[ilast]>=0 && (trace[ilast]<n_new || trace[ilast]==n_org) && ilast>idel) ilast--;

                // if ilast is not yet moved or (new pos of ilast > idel and ilast is not del, move ilast to idel)
                if(trace[ilast]<0||(idel<trace[ilast]&&trace[ilast]<n_org)) {
                    itrlast=trace[ilast];
                    trace[ilast]=idel;
                }
                // idel is already at last, remove it
                trace[idel]=n_org; 
                idel = itrlast;
            }

            n_decrease++;
        }

        // move data
        for (PS::S32 i=0; i<n_org; i++) {
            // no change case
            if(trace[i]<0) continue;
            // remove case
            else if(trace[i]<n_org) {
                const PS::S32 inew=trace[i];
                ptcl_[inew] = ptcl_[i];
                ptcl_ptr_[inew] = ptcl_ptr_[i];
                pred_[inew] = pred_[i];
                force_[inew] = force_[i];
                time_next_[inew] = time_next_[i];
                nb_info_[inew] = nb_info_[i];
                //_single_index_origin[inew-_n_group] = _single_index_origin[i-_n_group];

                // shift nb_list
                nb_list_n_[inew] = nb_list_n_[i];
                nb_list_n_[i] = 0;
                PS::S32 disp_i = nb_list_disp_[i];
                PS::S32 disp_inew = nb_list_disp_[inew];
                for (PS::S32 i=0; i<nb_list_n_[inew]; i++) 
                    nb_list_[disp_inew+i] = nb_list_[disp_i+i];
            }
        }
        ptcl_.decreaseSize(n_decrease);
        ptcl_ptr_.decreaseSize(n_decrease);
        pred_.decreaseSize(n_decrease);
        force_.decreaseSize(n_decrease);
        time_next_.decreaseSize(n_decrease);
        nb_info_.decreaseSize(n_decrease);
        nb_list_n_.decreaseSize(n_decrease);
        nb_list_disp_.decreaseSize(n_decrease);
        nb_list_.decreaseSize(n_decrease*n_nb_off_);
        //_n_single -= n_decrease;

        // check nb_info
        for (PS::S32 i=0; i<ptcl_.size(); i++) {
            const PS::S32 i_min=nb_info_[i].r_min_index;
            if(i_min>=0) {
                const PS::S32 jtr=trace[i_min];
                if(jtr>=0) {
                    if(jtr<n_org) nb_info_[i].r_min_index=jtr;
                    else nb_info_[i].r_min_index=-1;
                }
            }
        }

        // update adr_dt_sorted
        PS::S32 adr_dt_sorted_new_size=modifyList(adr_dt_sorted_.getPointer(), adr_dt_sorted_.size(), trace, n_org);
        adr_dt_sorted_.decreaseSize(_n_list);

        assert(adr_dt_sorted_new_size==adr_dt_sorted_.size());

        
        // correct neighbor list
        for (PS::S32 i=0; i<_n_nb_correct; i++) {
            // escape suppressed one
            if (trace[i]==n_org) continue;
            
            PS::S32 nb_list_i_new_size=modifyList(nb_list_.getPointer(nb_list_disp_[i]), nb_list_n_[i], trace, n_org);
            nb_list_n_[i] = nb_list_i_new_size;
#ifdef HERMITE_DEBUG
            assert(nb_list_n_[i]>=0);
#endif
        }
            
    }                       

    //! Write back particle data to ptcl
    /* @param[out] _ptcl: ptcl array to write
       @param[in] _ptcl_list: particle address in _ptcl
       @param[in] _n_ptcl: number of particles need for copy
       @param[in] _i_start: start index in local particle array for copy
     */
    template <class Tptcl>
    void writeBackPtcl(Tptcl * _ptcl, 
                       const PS::S32* _ptcl_list,
                       const PS::S32 _n_ptcl, 
                       const PS::S32 _i_start) {
#ifdef HERMITE_DEBUG
        assert(_i_start+_n_ptcl<=ptcl_.size());
#endif
        for (PS::S32 i=0; i<_n_ptcl; i++) {
#ifdef HERMITE_DEBUG
            assert(_ptcl[_ptcl_list[i]].id==ptcl_[i+_i_start].id);
            assert(!std::isnan(ptcl_[i+_i_start].pos[0]));
            assert(!std::isnan(ptcl_[i+_i_start].pos[1]));
            assert(!std::isnan(ptcl_[i+_i_start].pos[2]));
            assert(!std::isnan(ptcl_[i+_i_start].vel[0]));
            assert(!std::isnan(ptcl_[i+_i_start].vel[1]));
            assert(!std::isnan(ptcl_[i+_i_start].vel[2]));
#endif
            _ptcl[_ptcl_list[i]].DataCopy(ptcl_[i+_i_start]);
        }
    }

    //! Write back particle data to original array
    /*!
      @param[in] _i_start: starting index in ptcl_ to copy
    */
    void writeBackPtcl(const PS::S32 _i_start) {
        for (PS::S32 i=_i_start; i<ptcl_.size(); i++) {
#ifdef HERMITE_DEBUG
            assert(!std::isnan(ptcl_[i].pos[0]));
            assert(!std::isnan(ptcl_[i].pos[1]));
            assert(!std::isnan(ptcl_[i].pos[2]));
            assert(!std::isnan(ptcl_[i].vel[0]));
            assert(!std::isnan(ptcl_[i].vel[1]));
            assert(!std::isnan(ptcl_[i].vel[2]));
//            assert(ptcl_[i].pos[0]==ptcl_[i].pos[0]);
//            assert(ptcl_[i].pos[1]==ptcl_[i].pos[1]);
//            assert(ptcl_[i].pos[2]==ptcl_[i].pos[2]);
//            assert(ptcl_[i].vel[0]==ptcl_[i].vel[0]);
//            assert(ptcl_[i].vel[1]==ptcl_[i].vel[1]);
//            assert(ptcl_[i].vel[2]==ptcl_[i].vel[2]);
            assert(ptcl_ptr_[i]->id==ptcl_[i].id);
#endif
            ptcl_ptr_[i]->DataCopy(ptcl_[i]);
        }
    }

    //! Write back a list of particle data to original array 
    /*!
      @param[in] _ptcl_list: particle index in ptcl_
      @param[in] _n_ptcl: number of particles need for copy
      @param[in] _n_avoid: if particle index is below _n_avoid, no write (to avoid write group c.m.)
    */
    void writeBackPtcl(const PS::S32* _ptcl_list,
                       const PS::S32 _n_ptcl,
                       const PS::S32 _n_avoid) {
        for (PS::S32 i=0; i<_n_ptcl; i++) {
            const PS::S32 k = _ptcl_list[i];
            if(k<_n_avoid) continue;
#ifdef HERMITE_DEBUG
            assert(k<ptcl_.size()&&k>=0);
            assert(!std::isnan(ptcl_[k].pos[0]));
            assert(!std::isnan(ptcl_[k].pos[1]));
            assert(!std::isnan(ptcl_[k].pos[2]));
            assert(!std::isnan(ptcl_[k].vel[0]));
            assert(!std::isnan(ptcl_[k].vel[1]));
            assert(!std::isnan(ptcl_[k].vel[2]));
//            assert(ptcl_[k].pos[0]==ptcl_[k].pos[0]);
//            assert(ptcl_[k].pos[1]==ptcl_[k].pos[1]);
//            assert(ptcl_[k].pos[2]==ptcl_[k].pos[2]);
//            assert(ptcl_[k].vel[0]==ptcl_[k].vel[0]);
//            assert(ptcl_[k].vel[1]==ptcl_[k].vel[1]);
//            assert(ptcl_[k].vel[2]==ptcl_[k].vel[2]);
            assert(ptcl_ptr_[k]->id==ptcl_[k].id);
#endif            
            ptcl_ptr_[k]->DataCopy(ptcl_[k]);
        }
    }

    //! Write back one particle data to original array 
    /*!
      @param[in] _i: particle index in ptcl_
    */
    void writeBackOnePtcl(const PS::S32 _i) {
#ifdef HERMITE_DEBUG
        assert(_i<ptcl_.size()&&_i>=0);
        assert(!std::isnan(ptcl_[_i].pos[0]));
        assert(!std::isnan(ptcl_[_i].pos[1]));
        assert(!std::isnan(ptcl_[_i].pos[2]));
        assert(!std::isnan(ptcl_[_i].vel[0]));
        assert(!std::isnan(ptcl_[_i].vel[1]));
        assert(!std::isnan(ptcl_[_i].vel[2]));
//        assert(ptcl_[_i].pos[0]==ptcl_[_i].pos[0]);
//        assert(ptcl_[_i].pos[1]==ptcl_[_i].pos[1]);
//        assert(ptcl_[_i].pos[2]==ptcl_[_i].pos[2]);
//        assert(ptcl_[_i].vel[0]==ptcl_[_i].vel[0]);
//        assert(ptcl_[_i].vel[1]==ptcl_[_i].vel[1]);
//        assert(ptcl_[_i].vel[2]==ptcl_[_i].vel[2]);
        assert(ptcl_ptr_[_i]->id==ptcl_[_i].id);
#endif            
        ptcl_ptr_[_i]->DataCopy(ptcl_[_i]);
    }

    //! Search perturber and neighbor
    /*! Assume the suppressed ptcl has zero mass
      @param[in] _apo_bin: apo-center distance of binary 
      @param[in] _n_bin: number of binaries (locate at begining of ptcl_)
      @param[in] _n_search: number of particles to search perturber
     */
    void searchPerturberBin(PS::F64 _apo_bin[], const PS::S32 _n_bin, const PS::S32 _n_search) {
        PS::S32 n = ptcl_.size();
        
        // find perturber
        for(PS::S32 i=0; i<_n_search; i++) {
            PS::S32 disp=nb_list_disp_[i];
            PS::S32 n_pert=0;
            nb_info_[i].r_min2 = PS::LARGE_FLOAT;
            nb_info_[i].min_mass = PS::LARGE_FLOAT;
            for(PS::S32 j=0; j<_n_bin; j++) {
                PS::F64vec dr = ptcl_[i].pos-ptcl_[j].pos;
                PS::F64 r2 = dr*dr;
                PS::F64 r_search = std::max(ptcl_[i].r_search,ptcl_[j].r_search) + _apo_bin[j];
                PS::F64 r_search2 = r_search*r_search;
                if (r2<r_search2&&i!=j&&ptcl_[j].mass>0.0) {
                    nb_list_[disp+n_pert]=j;
                    n_pert++;
                    if(r2<nb_info_[i].r_min2) {
                        nb_info_[i].r_min2 = r2;
                        nb_info_[i].r_min_index = j;
                    }
                    nb_info_[i].min_mass = std::min(nb_info_[i].min_mass, ptcl_[j].mass);
                }
            }
            for(PS::S32 j=_n_bin; j<n; j++) {
                PS::F64vec dr = ptcl_[i].pos-ptcl_[j].pos;
                PS::F64 r2 = dr*dr;
                PS::F64 r_search = std::max(ptcl_[i].r_search,ptcl_[j].r_search);
                PS::F64 r_search2 = r_search*r_search;
                if (r2<r_search2&&i!=j&&ptcl_[j].mass>0.0) {
                    nb_list_[disp+n_pert]=j;
                    n_pert++;
                    if(r2<nb_info_[i].r_min2) {
                        nb_info_[i].r_min2 = r2;
                        nb_info_[i].r_min_index = j;
                    }
                    nb_info_[i].min_mass = std::min(nb_info_[i].min_mass, ptcl_[j].mass);
                }
            }
//            nb_list_disp_[i] = disp;
            nb_list_n_[i] = n_pert;
#ifdef HERMITE_DEBUG
            assert(n_pert<=n_nb_off_);
#endif
        }
    }

    //! Search perturber and find nearest particles for the first N particles in time step sorted list
    /*!
      @param[in] _n_search: number of particles to search perturber
     */
    void searchNeighborWithInfoN(const PS::S32 _n_search) {
        for(PS::S32 i=0; i<_n_search; i++) searchNeighborWithInfoOne(i);
    }

    //! Search neighbor and find nearest particles for one particle
    /*! Search perturber and neighbor for one particle from adr_dt_sorted list.
        In this case, the suppressed particle can be avoid
      @param[in] _i: index of particle index in adr_dt_sorted to search perturber
     */
    inline void searchNeighborWithInfoOne(const PS::S32 _i) {
        PS::S32 n = adr_dt_sorted_.size();
        // find perturber
        PS::S32 disp=nb_list_disp_[_i];
        PS::S32 n_pert=0;
        nb_info_[_i].r_min2 = PS::LARGE_FLOAT;
        nb_info_[_i].min_mass = PS::LARGE_FLOAT;
        for(PS::S32 k=0; k<n; k++) {
            PS::S32 j = adr_dt_sorted_[k];
            PS::F64vec dr = ptcl_[_i].pos-ptcl_[j].pos;
            PS::F64 r2 = dr*dr;
            PS::F64 r_search = std::max(ptcl_[_i].r_search,ptcl_[j].r_search);
            PS::F64 r_search2 = r_search*r_search;
            if (r2<r_search2&&_i!=j) {
#ifdef HERMITE_DEBUG
                assert(ptcl_[j].mass>0);
#endif
                nb_list_[disp+n_pert]=j;
                n_pert++;
                if(r2<nb_info_[_i].r_min2) {
                    nb_info_[_i].r_min2 = r2;
                    nb_info_[_i].r_min_index = j;
                }
                nb_info_[_i].min_mass = std::min(nb_info_[_i].min_mass, ptcl_[j].mass);
            }
        }
//            nb_list_disp_[i] = disp;
        nb_list_n_[_i] = n_pert;
#ifdef HERMITE_DEBUG
        assert(n_pert<=n_nb_off_);
#endif
    }

    //! update rsearch of ptcl
    /*! 
      @param[in] _i_start: start index to calculate
      @param[in] _dt_tree: tree time step
      @param[in] _v_max: maximum velocity used to calcualte r_search
      \return the maximum research
     */
    PS::F64 updateRSearch(const PS::S32 _i_start, const PS::F64 _dt_tree, const PS::F64 _v_max) {
        PS::F64 dt_reduce_factor=1.0;
        for(PS::S32 i=_i_start; i<ptcl_.size(); i++) {
            PS::F64 dt_reduce_fi = ptcl_[i].calcRSearch(_dt_tree, _v_max);
            dt_reduce_factor = std::max(dt_reduce_fi, dt_reduce_factor);
        }
        return dt_reduce_factor;
    }

    const PS::S32* getPertList(const PS::S32 i) {
        return &nb_list_[nb_list_disp_[i]];
    }

    PS::S32 getPertN(const PS::S32 i) const {
        return nb_list_n_[i];
    }

    PS::S32 getPertListSize() const {
        return nb_list_.size();
    }

    const PtclH4* getPtcl() const {
        return ptcl_.getPointer();
    }
    
    const PtclForce* getForce() const {
        return force_.getPointer();
    }

    PS::S32 getPtclN() const {
        return ptcl_.size();
    }

    const Tgptcl** getPtclAdr() const {
        return ptcl_ptr_.getPointer();
    }

    const Tneighbor& getNbInfo(const PS::S32 i) const {
        return nb_info_[i];
    }

#ifdef HERMITE_DEBUG_PRINT
    void writePtcl(FILE* _fout, const PS::S32 _i_start) const{
        for (PS::S32 i=_i_start; i<ptcl_.size(); i++) {
            ptcl_[i].ParticleBase::writeAscii(_fout);
        }
    }
#endif

    //! Set neighbor list offset (maximum number)
    /*!
      @param[in] _n_nb_off: neighbor list offset (maximum number)
     */
    void setNBOff(const PS::S32 _n_nb_off) {
        n_nb_off_ = _n_nb_off;
    }

    //! Get Minimum mass of ptcl
    /*!
      \return minimum mass
     */
    const PS::F64 getMassMin(){
        PS::F64 mass_min = PS::LARGE_FLOAT;
        for(PS::S32 i=0; i<ptcl_.size(); i++){
            if(mass_min > ptcl_[i].mass)  mass_min = ptcl_[i].mass;
        }
        return mass_min;
    }

    //! Initial Hermite ptcl force and step 
    /*! Initial f, fdot and step
      @param[in] _list: ptcl index to initialize
      @param[in] _n_list: number of ptcl to initialize
      @param[in] _time_sys: current set time
      @param[in] _int_pars: integration parameters (include changeover functions)
      @param[in,out] _Aint: ARC integrator class (resolve and shift is used)
      @param[in] _start_flag: indicate whether it is the initial step; true: set Nact to zero; false: predict particles for force calculation
     */
    template <class ARCint, class Tpars>
    bool initial(const PS::S32* _list,
                 const PS::S32 _n_list,
                 const PS::F64 _time_sys,
                 const Tpars * _int_pars,
                 ARCint* _Aint,
                 const bool _start_flag = true) {

        // if no particle need initial, quit
        if(_n_list==0) return false;
#ifdef HERMITE_DEBUG
        assert(_n_list>0);
#endif

        const PS::S32* ptcl_list=_list;
        if(ptcl_list==NULL) ptcl_list = adr_dt_sorted_.getPointer();

        PS::F64 mass_min = PS::LARGE_FLOAT;
        for(PS::S32 i=0; i<_n_list; i++){
            PS::S32 iadr = ptcl_list[i];
            pred_[iadr].r_search = ptcl_[iadr].r_search;
            ptcl_[iadr].time = _time_sys;
            ptcl_[iadr].dt   = _time_sys;
            ptcl_[iadr].acc0 = ptcl_[iadr].acc1 = 0.0;
            mass_min = std::min(mass_min, ptcl_[iadr].mass);
        }

        time_step_.setTimeAndCalcNextDtMax(_time_sys);

        if(_start_flag) {
            n_act_ = 0;
            time_step_.calcAcc0OffsetSq(mass_min, _int_pars->changeover.getRout());
        }
        else PredictAll(pred_.getPointer(), ptcl_.getPointer(), ptcl_.size(), _time_sys);

        if(_Aint!=NULL) {
            _Aint->updateCM(ptcl_.getPointer());
            _Aint->resolve();
        }

        // force::acc0,acc1, neighbor list updated
//        if(_calc_full_flag) 
        bool fail_flag=CalcActAcc0Acc1(force_.getPointer(), 
                                       nb_info_.getPointer(),
                                       nb_list_.getPointer(), nb_list_n_.getPointer(), 
                                       nb_list_disp_.getPointer(), 
                                       ptcl_.getPointer(), ptcl_.size(), 
                                       ptcl_list, _n_list, _int_pars,
                                       _Aint);
#ifdef HERMITE_DEBUG_DUMP
        if(fail_flag) return true;
#endif
//        // only neighbor force calculation
//        else CalcAcc0Acc1ActNb(force_.getPointer(), 
//                               ptcl_.getPointer(), 
//                               ptcl_list, _n_list,
//                               nb_list_.getPointer(), nb_list_disp_.getPointer(), nb_list_n_.getPointer(), 
//                               r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, _Aint, _n_group);
    
        // store predicted force
        for(PS::S32 i=0; i<_n_list; i++){
            PS::S32 iadr = ptcl_list[i];
            ptcl_[iadr].acc0 = force_[iadr].acc0;
            ptcl_[iadr].acc1 = force_[iadr].acc1;
        }

        if(_Aint!=NULL) _Aint->shift2CM();

        fail_flag=calcDt2ndAct(ptcl_.getPointer(), force_.getPointer(), adr_dt_sorted_.getPointer(), _n_list, time_step_);

#ifdef HERMITE_DEBUG_DUMP
        if(fail_flag) return true;
#endif
        for(PS::S32 i=0; i<_n_list; i++){
            PS::S32 iadr = ptcl_list[i];
            time_next_[iadr] = ptcl_[iadr].time + ptcl_[iadr].dt;
        }

        return fail_flag;
    }
    
//    template<class Energy>
//    void CalcEnergy(Energy & eng) {
//        CalcEnergy(ptcl_.getPointer(), ptcl_.size(), eng, r_in_, r_out_, eps_sq_);
//    }

    PS::F64 getNextTime() const {
        return time_next_[adr_dt_sorted_[0]];
    }

    PS::F64 getTime() const {
        return time_step_.getTime();
    }

    //PS::F64* getNextTimeList() const {
    //    return time_next_.getPointer();
    //}
    
    PS::F64 getOneTime(const std::size_t i) const {
        return ptcl_[i].time;
    }

    PS::F64 getOneDt(const std::size_t i) const {
        return ptcl_[i].dt;
    }

    //! Integration active particles
    /*! Integrated to time_sys
      @param[in] _int_pars: integration parameters (include changeover functions)
      @param[in,out] _Aint: ARC integrator class (resolve and shift is used)
      \return fail_flag: If step size < dt_min, return true
     */
    template <class ARCint, class Tpars>
    bool integrateOneStepAct(const Tpars* _int_pars,
                             ARCint* _Aint) {
        PS::S32 n_group = 0;

        // get next time
        PS::F64 time_sys = getNextTime();
        time_step_.setTimeAndCalcNextDtMax(time_sys);

        // pred::mass,pos,vel updated
        PredictAll(pred_.getPointer(), ptcl_.getPointer(), ptcl_.size(), time_sys);
        if(_Aint!=NULL) {
            n_group = _Aint->getNGroups();
            _Aint->updateCM(pred_.getPointer());
            _Aint->resolve();
        }

        // force::acc0,acc1, neighbor list updated
//        if(_calc_full_flag) 
        bool fail_flag=CalcAcc0Acc1Act(force_.getPointer(), 
                                       nb_info_.getPointer(),
                                       nb_list_.getPointer(), nb_list_n_.getPointer(), 
                                       nb_list_disp_.getPointer(), 
                                       pred_.getPointer(), ptcl_.size(), 
                                       adr_dt_sorted_.getPointer(), n_act_, 
                                       _int_pars,
                                       _Aint);
#ifdef HERMITE_DEBUG_DUMP
        if (fail_flag) return true;
#endif
//        // only neighbor force calculation
//        else CalcAcc0Acc1ActNb(force_.getPointer(), 
//                               pred_.getPointer(), 
//                               adr_dt_sorted_.getPointer(), n_act_, 
//                               nb_list_.getPointer(), nb_list_disp_.getPointer(), nb_list_n_.getPointer(), 
//                               r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, _Aint, _n_group);

        // ptcl_org::pos,vel; pred::time,dt,acc0,acc1,acc2,acc3 updated
        fail_flag=CorrectAndCalcDt4thAct(ptcl_.getPointer(), force_.getPointer(), adr_dt_sorted_.getPointer(), n_act_, time_step_);
#ifdef HERMITE_DEBUG_DUMP
        if (fail_flag) return true;
#endif

        for(PS::S32 i=0; i<n_act_; i++){
            PS::S32 adr = adr_dt_sorted_[i];
            time_next_[adr] = ptcl_[adr].time + ptcl_[adr].dt;
            // update new perturber list
            if(adr<n_group) _Aint->updatePertOneGroup(adr, ptcl_.getPointer(), force_.getPointer(), getPertList(adr), getPertN(adr));
        }

        if(_Aint!=NULL) {
            // shift member to c.m. frame
            _Aint->shift2CM();
            // go back to original time
            _Aint->updateCM(ptcl_.getPointer());
            //Aint.updateCM(Hint.getPtcl(), group_act_list.getPointer(), group_act_n);
        }
        
        return fail_flag;
    }

    //! Integration a list of particle to current time (ingore dt)
    /*! Integrate a list of particle to current time.
      @param[in] _ptcl_list: particle index list to integrate
      @param[in] _n_list: number of particles
      @param[in] _int_pars: integrator paramters
      @param[in,out] _Aint: ARC integrator class (resolve and shift is used)
      \return fail_flag: If step size < dt_min, return true
     */
    template <class ARCint, class Tpars>
    bool integrateOneStepListNoPred(const PS::S32* _ptcl_list,
                                    const PS::S32 _n_list,
                                    const Tpars* _int_pars,
                                    ARCint* _Aint) {

        if(_n_list==0) return false;
#ifdef HERMITE_DEBUG
#ifdef HERMITE_DEBUG_DUMP
        if(_n_list<0) {
            std::cerr<<"Error: number of particle in list is negative!\n";
            return true;
        }
#else
        assert(_n_list>0);
#endif
#endif
        PS::S32 n_group = 0;

        // get next time
        PS::F64 time_sys = time_step_.getTime();

        if(_Aint!=NULL) n_group = _Aint->getNGroups();

        // adjust dt
        PS::S32 int_list[_n_list];
        PS::S32 group_int_list[n_group+1];
        PS::S32 n_group_int=0;
        PS::S32 n_int=0;
        for (PS::S32 i=0; i<_n_list; i++) {
            PS::S32 iadr = _ptcl_list[i];
            if(ptcl_[iadr].time>=time_sys) continue;
            ptcl_[iadr].dt = time_sys - ptcl_[iadr].time;
            int_list[n_int++] = iadr;
            if(iadr<n_group) group_int_list[n_group_int++] = iadr;
        }

        if (n_int==0) return false;

        // update next time step limit

        if(_Aint!=NULL) _Aint->resolve();

        // force::acc0,acc1, neighbor list updated
//        if(_calc_full_flag) 
        bool fail_flag=CalcActAcc0Acc1(force_.getPointer(), 
                                       nb_info_.getPointer(),
                                       nb_list_.getPointer(), nb_list_n_.getPointer(), 
                                       nb_list_disp_.getPointer(), 
                                       pred_.getPointer(), ptcl_.size(), 
                                       int_list, n_int,
                                       _int_pars,
                                       _Aint);
//        // only neighbor force calculation
//        else CalcAcc0Acc1ActNb(force_.getPointer(), 
//                               pred_.getPointer(), 
//                               int_list, n_int,
//                               nb_list_.getPointer(), nb_list_disp_.getPointer(), nb_list_n_.getPointer(), 
//                               r_in_, r_out_, r_oi_inv_, r_A_, eps_sq_, _Aint, _n_group);

#ifdef HERMITE_DEBUG_DUMP
        if (fail_flag) return true;
#endif

        fail_flag=CorrectAndCalcDt4thAct(ptcl_.getPointer(), 
                                         force_.getPointer(), 
                                         int_list, n_int, time_step_);


#ifdef HERMITE_DEBUG_DUMP
        if (fail_flag) return true;
#endif

        // shift member to c.m. frame
        if(_Aint!=NULL) {
            _Aint->shift2CM();
            // update c.m. of _i
            _Aint->updateCM(ptcl_.getPointer(), group_int_list, n_group_int);
        }

        for(PS::S32 i=0; i<n_int; i++){
            PS::S32 iadr=int_list[i];
            // update time table
            time_next_[iadr] = ptcl_[iadr].time + ptcl_[iadr].dt;
        }

        // update new perturber list
        for(PS::S32 i=0; i<n_group_int; i++) {
            PS::S32 iadr=group_int_list[i];
            _Aint->updatePertOneGroup(iadr, ptcl_.getPointer(), force_.getPointer(), getPertList(iadr), getPertN(iadr));
        }

        return false;
    }

    void SortAndSelectIp(PS::S32 group_act_list[],
                         PS::S32 &group_act_n,
                         const PS::S32 n_groups) {
        SortAndSelectIp(adr_dt_sorted_.getPointer(), time_next_.getPointer(), n_act_, time_next_.size(), group_act_list, group_act_n, n_groups);
        
    }

    void SortAndSelectIp() {
        SortAndSelectIp(adr_dt_sorted_.getPointer(), time_next_.getPointer(), n_act_, adr_dt_sorted_.size());
    }
    
    PS::S32 getNact() const{
        return n_act_;
    }

    const PS::S32* getActList() const{
        return adr_dt_sorted_.getPointer();
    }

#ifdef HERMITE_DEBUG
    template <class ARCint>
    bool checkAdrList(const ARCint& _Aint) {
        bool fail_flag = false;
        const PS::S32 n_ptcl = ptcl_.size();
        const PS::S32 n_adr = adr_dt_sorted_.size();
        PS::S32 ncheck[n_ptcl];
        for (PS::S32 i=0; i<n_ptcl; i++) ncheck[i]=0;
        for (PS::S32 i=0; i<n_adr; i++) {
            PS::S32 k =adr_dt_sorted_[i];
#ifdef HERMITE_DEBUG_DUMP
            if (time_next_[k] != ptcl_[k].time + ptcl_[k].dt) fail_flag = true;
#else
            assert(time_next_[k]==ptcl_[k].time + ptcl_[k].dt);
#endif
            PS::F64 tr=PS::S32(ptcl_[k].time/ptcl_[k].dt);
#ifdef HERMITE_DEBUG_DUMP
            if (tr*ptcl_[k].dt != ptcl_[k].time) fail_flag = true;
#else
            assert(tr*ptcl_[k].dt==ptcl_[k].time);
#endif
            ncheck[adr_dt_sorted_[i]]++;
        }
        const PS::S32 n_group = _Aint.getNGroups();
        for (PS::S32 i=0; i<n_group; i++) {
#ifdef HERMITE_DEBUG_DUMP
            if(_Aint.getMask(i)) { 
                if (ncheck[i]!=0) fail_flag = true;
            }
            else if (ncheck[i]!=1) 
                fail_flag = true;
#else
            if(_Aint.getMask(i)) assert(ncheck[i]==0);
            else assert(ncheck[i]==1);
#endif
        }
        for (PS::S32 i=n_group; i<n_ptcl; i++) {
#ifdef HERMITE_DEBUG_DUMP
            if (ncheck[i]!=1) fail_flag = true;
#else
            assert(ncheck[i]==1);
#endif
        }
        for (PS::S32 i=0; i<n_adr-1; i++) {
#ifdef HERMITE_DEBUG_DUMP
            if (ptcl_[adr_dt_sorted_[i]].dt>ptcl_[adr_dt_sorted_[i+1]].dt) fail_flag = true;
            if (time_next_[adr_dt_sorted_[i]]>time_next_[adr_dt_sorted_[i+1]]) fail_flag = true;
#else
            assert(ptcl_[adr_dt_sorted_[i]].dt<=ptcl_[adr_dt_sorted_[i+1]].dt);
            assert(time_next_[adr_dt_sorted_[i]]<=time_next_[adr_dt_sorted_[i+1]]);
#endif

        }
        return fail_flag;
    }

    void printStepHist(){
        std::map<PS::F64, PS::S32> stephist;
        for(int i=0; i<adr_dt_sorted_.size(); i++) {
            PS::S32 k = adr_dt_sorted_[i];
            std::map<PS::F64, PS::S32>::iterator p = stephist.find(ptcl_[k].dt);
            if (p==stephist.end()) stephist[ptcl_[k].dt]=1;
            else stephist[ptcl_[k].dt]++;
        }
        std::cerr<<"Step hist:\n";
        for(auto i=stephist.begin(); i!=stephist.end(); i++) {
            std::cerr<<std::setw(24)<<i->first;
        }
        std::cerr<<std::endl;
        for(auto i=stephist.begin(); i!=stephist.end(); i++) {
            std::cerr<<std::setw(24)<<i->second;
        }
        std::cerr<<std::endl;
    }
#endif

};

