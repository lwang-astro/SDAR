#pragma once

#include "Common/Float.h"
#include "particle.h"
#include "perturber.h"

namespace H4{

    //! hermite interaction class 
    class HermiteInteraction{
    public:
        template<class Tptcl>
        inline bool CalcAccAndNeighbor(PtclForce & _fi,
                                       Float & _r2,
                                       const Tptcl& _pi,
                                       const Tptcl& _pj,
                                       const HermiteParams& _pars) {
            const Floatvec dR = _pi.pos - _pj.pos;
            r2 = dR*dR;
            const Float dr2_eps = r2 + pars.eps_sq;
            const Float r_out = _pars.changeover.getRout();
            if(dr2_eps <= r_out*r_out){
                const Floatvec dV = _pi.vel - _pj.vel;
                const Float drdv = dR * dV;
                const Float dr_eps = sqrt(dr2_eps);
                const Float dri = 1.0/dr_eps;
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
                //const Float R = 1.0 / sqrt(r2_eps);
                const Float dri2 = dri*dri;
                const Float dri3 = dri2*dri;
                const Float k = pars.changeover.calcAcc0W(r_eps);
                const Float dk = pars.changeover.calcAcc1W(r_eps) * pars.changeover.getNorm();
                const Float mri3 = _pj.mass*dri3;
                const Floatvec F0 = -mri3*k*dR;
                const Floatvec F1 = -mri3*k*dV - 3.0*drdv*(dri2*F0 + mri3*dk*dri*dR);
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
                                    const int _n_tot,
                                    const HermiteParams& _int_pars,
                                    const ARCint* _Aint = NULL) {
            if(_Aint!=NULL) n_group = _Aint->getNGroups();
            for(int j=0; j<n_group; j++) {
                if(pi.id==_ptcl[j].id) continue;
                if(_Aint->getMask(j)) continue;
#ifdef HERMITE_DEBUG
                Float mcmcheck =0.0;
#endif
                const auto* pj = _Aint->getGroupPtcl(j);
                Float sdj = 1.0/_Aint->getSlowDown(j);
                Floatvec vcmsdj = (1.0-sdj)*_ptcl[j].vel;

                for(int k=0; k<_Aint->getGroupN(j); k++) {
                    Float r2 = 0.0;

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
                    Float rs=std::max(_pi.r_search, _pj[k].r_search);
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
            for(int j=n_group; j<_n_tot; j++){
                if(_pi.id==_ptcl[j].id) continue;

                Float r2 = 0.0;
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
                Float rs = std::max(_pi.r_search, _ptcl[j].r_search);
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
                                    const int _n_tot,
                                    const int _list_act[],
                                    const int _n_act,
                                    const HermiteParams& _int_pars,
                                    const ARCint* _Aint) {
            bool fail_flag = false;
            int n_group = 0;
            if(_Aint!=NULL) n_group = _Aint->getNGroups();

            const Floatvec vzero = Floatvec(0.0);
            // PS::ReallocatableArray< std::pair<int, int> > & merge_pair ){

            // active iparticles loop
            for(int i=0; i<_n_act; i++){
                const int iadr = _list_act[i];
                _force[iadr].clear();
            
                _nb_info[iadr].resetRmin();
                _nb_info[iadr].resetMassMin();
                _nb_info[iadr].resetNBFlag();
            
                // for group particle
                if (iadr<n_group) {
#ifdef HERMITE_DEBUG
                    Float mcmcheck =0.0;
#endif
                    const int ni = _Aint->getGroupN(iadr);             // number of members
#ifdef HERMITE_DEBUG
                    assert(ni>0);
#endif
                    const auto* pi = _Aint->getGroupPtcl(iadr);
                    const Float sdi = 1.0/_Aint->getSlowDown(iadr);      // slowdown factor
                    const Floatvec vcmsdi = (1.0-sdi)*_ptcl[iadr].vel;   // slowdown velocity
//                Float r2min[_n_tot]={PS::LARGE_FLOAT};

                    for (int j=0; j<ni; j++) {
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

    };
}
