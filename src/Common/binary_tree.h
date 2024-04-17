#pragma once

#include<cmath>
#include<algorithm>
#include"Common/Float.h"
#include"Common/matrix.h"

namespace COMM{
    const Float PI = 4.0*atan(1.0);

    //! Binary parameter class
    class Binary{
    public:
        Float semi;         // semi: semi-major axis
        Float ecc;          // ecc: eccentricity
        Float incline;      // incline: inclination angle, Eular angle beta
        Float rot_horizon;  // rotational angle in horizon, Eular angle alpha
        Float rot_self;     // rotational angle in orbital plane, Eular angle gamma
        Float t_peri;       // t_peri: time to peri-center
        Float period;       // period: period
        Float ecca;         // ecca: eccentricty anomaly (-pi, pi)
        Float m1;           // m1: mass 1
        Float m2;           // m2: mass 2
        Float r;            // distance between too members
        Vector3<Float> am;  // specific angular momentum
        Float stab;         // stability factor 

        //! Constructor
        /*! Set all variables to numeric maximum of Float, convenient for debugging
         */
        Binary(): semi(NUMERIC_FLOAT_MAX), ecc(NUMERIC_FLOAT_MAX), 
                  incline(NUMERIC_FLOAT_MAX), rot_horizon(NUMERIC_FLOAT_MAX),
                  rot_self(NUMERIC_FLOAT_MAX), t_peri(NUMERIC_FLOAT_MAX),
                  period(NUMERIC_FLOAT_MAX), ecca(NUMERIC_FLOAT_MAX),
                  m1(NUMERIC_FLOAT_MAX), m2(NUMERIC_FLOAT_MAX), r(NUMERIC_FLOAT_MAX),
                  am(NUMERIC_FLOAT_MAX,NUMERIC_FLOAT_MAX,NUMERIC_FLOAT_MAX), stab(NUMERIC_FLOAT_MAX) {}
                  

        //! Orbit to position and velocity
        /*! refer to the P3T code developed by Iwasawa M.
          @param[out] _p1: particle 1
          @param[out] _p2: particle 2
          @param[in] _bin: binary parameter
          @param[in] _ecca: eccentric anomaly (-pi, pi)
          @param[in] _G: gravitational constant
        */
        template <class Tptcl>
        static void orbitToParticle(Tptcl& _p1, Tptcl& _p2, const Binary& _bin, const Float& _ecca, const Float _G) {
            Float m_tot = _bin.m1 + _bin.m2;
            Vector3<Float> pos_star, vel_star;
            // hyper bolic orbit
            if (_bin.semi<0) {
                Float coshu = cosh(_ecca);
                Float sinhu = sinh(_ecca);
                Float eccfac = sqrt(_bin.ecc*_bin.ecc - 1.0);
                pos_star.x = -_bin.semi*(_bin.ecc - coshu);   // -a (e - cosh(E))
                pos_star.y = -_bin.semi*eccfac*sinhu;         // -a sqrt(e^2-1) sinh(E)
                pos_star.z = 0.0;
                Float afac = sqrt(_G*m_tot / -_bin.semi);     // mean_motion * a
                Float fac  = afac/(_bin.ecc*coshu-1.0);  // a/r * sqrt(GM/a)
                vel_star.x = -fac*sinhu;                 // - sqrt(GM/-a) sinh(E)/(e cosh(E)-1)
                vel_star.y =  fac*eccfac*coshu;          //   sqrt(GM(e^2-1)/-a) cosh(E)/(e cosh(E)-1)
                vel_star.z = 0.0;
            }
            else{ // ellipse orbit
                Float cosu = cos(_ecca);
                Float sinu = sin(_ecca);
                Float eccfac = sqrt(1.0 - _bin.ecc*_bin.ecc);     
                pos_star.x = _bin.semi*(cosu - _bin.ecc);         // a (cos(E) - e)
                pos_star.y = _bin.semi*eccfac*sinu;               // a sqrt(1-e^2) sin(E)
                pos_star.z = 0.0;
                Float afac = sqrt(_G*m_tot / _bin.semi);
                Float fac = afac/(1.0-_bin.ecc*cosu);   // a/r * sqrt(GM/a)
                vel_star.x = -fac*sinu;                 // - sqrt(GM/a) sin(E)/(1 - e cos(E))
                vel_star.y =  fac*eccfac*cosu;          //   sqrt(GM(1-e^2)/a) cos(E)/(1 - e cos(E))
                vel_star.z = 0.0;
            }
            Matrix3<Float> rot;
            rot.rotation(_bin.incline, _bin.rot_horizon, _bin.rot_self);
            Vector3<Float> pos_red = rot*pos_star;
            Vector3<Float> vel_red = rot*vel_star;
            _p1.mass = _bin.m1;
            _p2.mass = _bin.m2;
            
            Float m2_mt = _p2.mass / m_tot;
            Float m1_mt = _p1.mass / m_tot;
            _p1.pos[0] = -m2_mt * pos_red.x;
            _p1.pos[1] = -m2_mt * pos_red.y;
            _p1.pos[2] = -m2_mt * pos_red.z;

            _p2.pos[0] =  m1_mt * pos_red.x;
            _p2.pos[1] =  m1_mt * pos_red.y;
            _p2.pos[2] =  m1_mt * pos_red.z;

            _p1.vel[0] = -m2_mt * vel_red.x;
            _p1.vel[1] = -m2_mt * vel_red.y;
            _p1.vel[2] = -m2_mt * vel_red.z;

            _p2.vel[0] =  m1_mt * vel_red.x;
            _p2.vel[1] =  m1_mt * vel_red.y;
            _p2.vel[2] =  m1_mt * vel_red.z;
        }

        //! position velocity to orbit
        /* refer to the P3T code developed by Iwasawa M.
           @param[out] _bin: binary parameter
           @param[in]  _p1: particle 1
           @param[in]  _p2: particle 2
           @param[in] _G: gravitational constant
        */
        template <class Tptcl>
        static void particleToOrbit(Binary& _bin, const Tptcl& _p1, const Tptcl& _p2, const Float _G){
            _bin.m1 = _p1.mass;
            _bin.m2 = _p2.mass;
            Float m_tot = _p1.mass + _p2.mass;
            Float Gm_tot = _G*m_tot;
            Vector3<Float> pos_red(_p2.pos[0] - _p1.pos[0], _p2.pos[1] - _p1.pos[1], _p2.pos[2] - _p1.pos[2]);
            Vector3<Float> vel_red(_p2.vel[0] - _p1.vel[0], _p2.vel[1] - _p1.vel[1], _p2.vel[2] - _p1.vel[2]);
            Float r_sq = pos_red * pos_red;
            _bin.r = sqrt(r_sq);
            Float inv_dr = 1.0 / _bin.r;
            Float v_sq = vel_red * vel_red;
            _bin.semi = 1.0 / (2.0*inv_dr - v_sq /  Gm_tot);
            _bin.am = pos_red ^ vel_red;
            _bin.incline = atan2( sqrt(_bin.am.x*_bin.am.x+_bin.am.y*_bin.am.y), _bin.am.z);
            _bin.rot_horizon = _bin.am.y!=0.0? atan2(_bin.am.x, -_bin.am.y) : 0.0;

            // position and velocity in rest frame
            Vector3<Float> pos_bar, vel_bar;
            // OMG: rot_horizon; omg: rot_self; inc: incline
            Float cosOMG = cos(_bin.rot_horizon);
            Float sinOMG = sin(_bin.rot_horizon);
            Float cosinc = cos(_bin.incline);
            Float sininc = sin(_bin.incline);
            pos_bar.x =   pos_red.x*cosOMG + pos_red.y*sinOMG;
            pos_bar.y = (-pos_red.x*sinOMG + pos_red.y*cosOMG)*cosinc + pos_red.z*sininc;
            pos_bar.z = 0.0;
            vel_bar.x =   vel_red.x*cosOMG + vel_red.y*sinOMG;
            vel_bar.y = (-vel_red.x*sinOMG + vel_red.y*cosOMG)*cosinc + vel_red.z*sininc;
            vel_bar.z = 0.0;
            // h: specific angular momentum value (|r x v|)
            Float h = sqrt(_bin.am*_bin.am);
            // eccentricity vector: (v x h)/GM - r/|r|
            // angle between AM and vel is 90 degree
            Float ecccosomg =  h/Gm_tot*vel_bar.y - pos_bar.x*inv_dr;
            Float eccsinomg = -h/Gm_tot*vel_bar.x - pos_bar.y*inv_dr;
            _bin.ecc = sqrt( ecccosomg*ecccosomg + eccsinomg*eccsinomg );
            _bin.rot_self = atan2(eccsinomg, ecccosomg);
            // angle of position referring to the rest frame: f + omg (f: true anomaly)
            Float phi = atan2(pos_bar.y, pos_bar.x); 
            Float true_anomaly = phi - _bin.rot_self;
            // eccentric anomaly
            //Float sin_ecca = _bin.r*sin(true_anomaly) / (_bin.semi*sqrt(1.0 - _bin.ecc*_bin.ecc));
            //Float cos_ecca = (_bin.r*cos(true_anomaly) / _bin.semi) + _bin.ecc;
            if (_bin.semi>0) {
                _bin.ecca = atan2(sin(true_anomaly)*sqrt(1.0 - _bin.ecc*_bin.ecc), _bin.ecc+ cos(true_anomaly));  
                Float mean_motion = sqrt(Gm_tot/(_bin.semi*_bin.semi*_bin.semi)); 
                _bin.period = 2*PI/mean_motion;
                Float mean_anomaly = _bin.ecca - _bin.ecc*sin(_bin.ecca); 
                _bin.t_peri = mean_anomaly / mean_motion; 
            }
            else {
                _bin.ecca = atanh(sin(true_anomaly)*sqrt(_bin.ecc*_bin.ecc-1)/(_bin.ecc+ cos(true_anomaly)));  
                Float mean_motion = sqrt(Gm_tot/(-_bin.semi*_bin.semi*_bin.semi)); 
                _bin.period = 2*PI/mean_motion;
                Float mean_anomaly = _bin.ecc*sinh(_bin.ecca) - _bin.ecca;
                _bin.t_peri = mean_anomaly / mean_motion; 
            }
        }

        //! from period calculate semi-major axis
        /*!
          @param[in] _period: period
          @param[in] _mtot: total mass of binary
          @param[in] _G: gravitational constant
         */
        static Float periodToSemi(const Float& _period, const Float& _mtot, const Float& _G) {
            return pow(_period*_period*_G*_mtot/(4*PI*PI),1.0/3.0);
        }

        //! from semi-major axis calculate period;
        /*!
          @param[in] _semi: semi-major axis
          @param[in] _mtot: total mass of binary
          @param[in] _G: gravitational constant
         */
        static Float semiToPeriod(const Float& _semi, const Float& _mtot, const Float& _G) {
            Float mean_motion = sqrt(_G*_mtot/(fabs(_semi*_semi*_semi))); 
            return 2*PI/mean_motion;
        }
        

        //! position velocity to orbit semi-major axis and eccentricity
        /*! @param[out] _semi: semi-major axis
           @param[out] _ecc:  eccentricity
           @param[out] _r: distance between two particles
           @param[out] _rv: relative position dot velocity
           @param[in]  _p1: particle 1
           @param[in]  _p2: particle 2
           @param[in]  _G: gravitational constant
        */
        template <class Tpi, class Tpj>
        static void particleToSemiEcc(Float& _semi, Float& _ecc, Float& _r, Float& _rv, const Tpi& _p1, const Tpj& _p2, const Float _G){
            Float m_tot = _p1.mass + _p2.mass;
            Float Gm_tot = _G*m_tot;
            Vector3<Float> pos_red(_p2.pos[0] - _p1.pos[0], _p2.pos[1] - _p1.pos[1], _p2.pos[2] - _p1.pos[2]);
            Vector3<Float> vel_red(_p2.vel[0] - _p1.vel[0], _p2.vel[1] - _p1.vel[1], _p2.vel[2] - _p1.vel[2]);
            Float r_sq = pos_red * pos_red;
            _r = sqrt(r_sq);
            Float v_sq = vel_red * vel_red;
            _rv = pos_red * vel_red;
            _semi = 1.0 / (2.0 / _r - v_sq / Gm_tot);
            Float p = 1.0 - _r/_semi;
            _ecc = sqrt(p*p + _rv*_rv/_semi/Gm_tot);
        }

        //! calcualte semi, ecc and period
        /*! 
           @param[in,out]  _p1: particle 1
           @param[in,out]  _p2: particle 2
           @param[in]  _G: gravitational constant
        */
        template <class Tpi, class Tpj>
        void particleToSemiEccPeriod(const Tpi& _p1, const Tpj& _p2, const Float _G) {
            Float r,rv;
            particleToSemiEcc(semi, ecc, r, rv, _p1, _p2, _G);
            Float mtot=m1+m2;
            period = semiToPeriod(semi, mtot, _G);
        }

        //! calculate eccentric anomaly from separation (0-pi)
        /*!
          @param[in] _r: seperation
          \return eccentric anomaly
         */
        Float calcEccAnomaly(const Float _r) {
            if (semi>0) {
                Float cos_ecca = (1.0-_r/semi)/ecc;
                return acos(cos_ecca);
            }
            else {
                Float cosh_ecca = (1.0-_r/semi)/ecc;
                return acosh(cosh_ecca);
            }
        }

        //! calculate mean anomaly from eccentric anomaly (0-pi)
        /*!
           @param[in] _ecca:  eccentric anomaly
           @param[in] _ecc:  eccentricity
           \return mean anomaly
         */
        static Float calcMeanAnomaly(const Float _ecca, const Float _ecc) {
            if (_ecc<1.0) return _ecca - _ecc*sin(_ecca); 
            else return _ecc*sinh(_ecca) - _ecca;
        }

        //! calculate eccentric anomaly from mean anomaly
        /* refer to the P3T code developed by Iwasawa M.
          @param[in] _mean_anomaly: mean_anomaly
          @param[in] _ecc: eccentricity
          \return eccentric anomaly
        */
        static Float calcEccAnomaly(const Float _mean_anomaly,
                                    const Float _ecc){
            // a: semi-major axis
            // l: mean anomaly
            // e: eccentricity
            // u: eccentric anomaly
            // n: mean mortion
            Float u0 = _mean_anomaly;
            Float u1;
            int loop = 0;
            while(1){
                loop++;
                Float su0 = sin(u0);
                Float cu0 = sqrt(1.0 - su0*su0);
                //u1 = u0 - keplereq(l, e, u0)/keplereq_dot(e, u0);
                u1 = u0 - ((u0- _ecc*su0-_mean_anomaly)/(1.0 - _ecc*cu0));
                //if( fabs(u1-u0) < 1e-13 ){ return u1; }
                if( fabs(u1-u0) < 1e-15 ){ return u1; }
                else{ u0 = u1; }
                if (loop>1e5) {
                    std::cerr<<"Error: kepler solver cannot converge to find correct eccentricity anomaly!\n";
                    abort();
                }
            }
        }

        //! solve kepler orbit after dt
        /*! refer to the P3T code developed by Iwasawa M.
          @param[in,out] _bin: binary orbit
          @param[in] _dt: evolution time
        */
        static void solveKepler(Binary& _bin,
                                const Float _dt) {
            Float freq = sqrt( (_bin.m1+_bin.m2) / (_bin.semi*_bin.semi*_bin.semi) );
            Float mean_anomaly_old = _bin.ecca - _bin.ecc * sin(_bin.ecca);
            Float dt_tmp = _dt - to_int(_dt/_bin.period)*_bin.period;
            Float mean_anomaly_new = freq * dt_tmp + mean_anomaly_old; // mean anomaly
            _bin.ecca = calcEccAnomaly(mean_anomaly_new, _bin.ecc); // eccentric anomaly
        }
    
        //! calculate kepler Orbit from particles
        /*! 
          @param[out] _p1: particle 1
          @param[out] _p2: particle 2
          @param[in] _G: gravitational constant
        */    
        template <class Tptcl>
        void calcOrbit(const Tptcl& _p1, const Tptcl& _p2, const Float _G) {
            particleToOrbit(*this, _p1, _p2, _G);
        }

        //! calculate two components from kepler Orbit 
        /*! 
          @param[out] _p1: particle 1
          @param[out] _p2: particle 2
          @param[in] _G: gravitational constant
        */    
        template <class Tptcl>
        void calcParticles(Tptcl& _p1, Tptcl& _p2, const Float& _G) {
            orbitToParticle(_p1, _p2, *this, this->ecca, _G);
        }

        //! from period calculate semi-major axis
        /*!
          @param[in] _G: gravitational constant
         */
        void calcSemiFromPeriod(const Float& _G) {
            Float mtot = m1+m2;
            semi = periodToSemi(period, mtot, _G);
        }

        //! calculate two components from kepler Orbit with input eccentricity anomaly
        /*! 
          @param[out] _p1: particle 1
          @param[out] _p2: particle 2
          @param[in] _ecca: eccentricity anomaly
          @param[in] _G: gravitational constant
        */    
        template <class Tptcl>
        void calcParticlesEcca(Tptcl& _p1, Tptcl& _p2, const Float _ecca, const Float _G) const {
            orbitToParticle(_p1, _p2, *this, _ecca, _G);
        }

        //! rotate position vector from binary rest-frame to original frame based on three angles
        /*!
          @param[out] _vec: vector to rotate
         */
        void rotateToOriginalFrame(Float* _vec ) {
            Matrix3<Float> rot;
            rot.rotation(incline, rot_horizon, rot_self);
            Vector3<Float> vec (_vec[0], _vec[1], _vec[2]);
            vec = rot*vec;
            _vec[0] = vec.x;
            _vec[1] = vec.y;
            _vec[2] = vec.z;
        }

        //! Solve kepler motion for dt 
        /*! Evolve current orbit for dt 
          @param[in] _dt: evolution time
        */
        void evolve(const Float _dt) {
            solveKepler(*this, _dt);
        }

        // IO functions
        void print(std::ostream& _os) const{
            _os<<"semi:        semi-major axis: "<<semi<<std::endl
               <<"ecc:         eccentricity:    "<<ecc<<std::endl
               <<"incline:     inclination angle:                 "<<incline<<std::endl
               <<"rot_horizon: rotational angle in horizon:       "<<rot_horizon<<std::endl
               <<"rot_self:    rotational angle in orbital plane: "<<rot_self<<std::endl
               <<"t_peri:      phase to peri-center:  "<<t_peri<<std::endl
               <<"period:      period:               "<<period<<std::endl
               <<"ecca:        eccentricty anomaly: "<<ecca<<std::endl
               <<"m1:          mass 1: "<<m1<<std::endl
               <<"m2:          mass 2: "<<m2<<std::endl
               <<"r:           separation: "<<r<<std::endl
               <<"L:           Angular momentum: "<<am.x<<" "<<am.y<<" "<<am.z<<std::endl
               <<"stab:        stability factor: "<<stab<<std::endl;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        static void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"semi"
                 <<std::setw(_width)<<"ecc"
                 <<std::setw(_width)<<"incline"
                 <<std::setw(_width)<<"rot_horizon"
                 <<std::setw(_width)<<"rot_self"
                 <<std::setw(_width)<<"t_peri"
                 <<std::setw(_width)<<"period"
                 <<std::setw(_width)<<"ecca"
                 <<std::setw(_width)<<"m1"
                 <<std::setw(_width)<<"m2"
                 <<std::setw(_width)<<"r"
                 <<std::setw(_width)<<"L.x"
                 <<std::setw(_width)<<"L.y"
                 <<std::setw(_width)<<"L.z"
                 <<std::setw(_width)<<"stab";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<semi
                 <<std::setw(_width)<<ecc
                 <<std::setw(_width)<<incline
                 <<std::setw(_width)<<rot_horizon
                 <<std::setw(_width)<<rot_self
                 <<std::setw(_width)<<t_peri
                 <<std::setw(_width)<<period
                 <<std::setw(_width)<<ecca
                 <<std::setw(_width)<<m1
                 <<std::setw(_width)<<m2
                 <<std::setw(_width)<<r
                 <<std::setw(_width)<<am.x
                 <<std::setw(_width)<<am.y
                 <<std::setw(_width)<<am.z
                 <<std::setw(_width)<<stab;
        }


        //! write class data to file with binary format
        /*! @param[in] _fp: FILE type file for output
         */
        void writeBinary(FILE *_fp) const {
            fwrite(this, sizeof(*this),1,_fp);
        }

        //! read class data to file with binary format
        /*! @param[in] _fin: FILE type file for reading
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(this, sizeof(*this),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }

        //! write class data to file with ASCII format
        /*! @param[in] _fout: std:osteram file for output
         */
        void writeAscii(std::ostream& _fout) const {
            _fout<<semi<<" "
                 <<ecc<<" "  
                 <<incline<<" "   
                 <<rot_horizon<<" "
                 <<rot_self<<" "
                 <<t_peri<<" "
                 <<period<<" "
                 <<ecca<<" "
                 <<m1<<" "  
                 <<m2<<" "    
                 <<r<<" "    
                 <<am.x<<" "
                 <<am.y<<" "  
                 <<am.z<<" "
                 <<stab<<" ";
        }

        //! read class data to file with ASCII format
        /*! @param[in] _fin: std::istream file for input
         */
        void readAscii(std::istream&  _fin) {
            _fin>>semi>>ecc>>incline>>rot_horizon>>rot_self>>t_peri>>period>>ecca>>m1>>m2>>r>>am.x>>am.y>>am.z>>stab;
        }

        //void PosVel2SemiEcc(Float & semi,
        //                    Float & ecc,
        //                    const Vector3<Float> & pos0, const Vector3<Float> & pos1,
        //                    const Vector3<Float> & vel0, const Vector3<Float> & vel1,
        //                    const Float mass0, const Float mass1) {
        //    Float m_tot = mass0 + mass1;
        //    Vector3<Float> pos_red = pos1 - pos0;
        //    Vector3<Float> vel_red = vel1 - vel0;
        //    Float r_sq = pos_red * pos_red;
        //    Float r = sqrt(r_sq);
        //    Float inv_dr = 1.0 / r;
        //    Float v_sq = vel_red * vel_red;
        //    semi = 1.0 / (2.0*inv_dr - v_sq / m_tot);
        // 
        //    Float rs = 1.0 - r/semi;
        //    Float rv = pos_red * vel_red;
        //    ecc = std::sqrt(rs*rs + rv*rv/(m_tot*semi));
        //}
    };

//! Binary tree cell 
/*! 
*/
    template <class Tptcl, class Tbinary>
    class BinaryTree: public Tptcl, public Tbinary {
    private: 
        typedef std::pair<Float, int> R2Index;  // First is distance square, second is index
        typedef BinaryTree<Tptcl,Tbinary> BinaryTreeLocal; // local defined binarytree 
        /*                          Example for level & branch
                          -------------------------------------------------
                               level               branch 
                                0                      0                   
                                                      / \                 
                                1                    0   1                 
                                                    / \ / \
                                2                  0  1 2  3               
                          -------------------------------------------------*/
        int n_members; ///> Total member number belong to the tree (can be more than 2)
        int member_index[2];        ///> particle index of original list, -1 indicate member is binarytree, -2 indicates member is particle and its memory is locally allocated; -3 initial values
        Tptcl* member[2]; ///> member pointer

    public:
        int level;     ///> binary tree level, start from zero for root tree. Then each subtree increase level by one.
        int branch;    ///> the branch in each level, count from the left member to right. 

    private:
        static bool pairLess(const R2Index& a, const R2Index & b)  {
            return a.first < b.first;
        }

        //! Find the minimum distant particle and same the sqrt distance at first and index at second of r2_list
        /* @param[in,out] _member_list: group particle member index, will be reordered by the minimum distance chain.
           @param[in,out] _r2_list: a pair. first is the square distance of closes neighbor (next particle i+1 in _member_list); second is equal to particle current index i
           @param[in] _n_members: number of members
           @param[in] _ptcl_org: original particle data
        */
        static void calcMinDisList(int _ptcl_list[], 
                                   R2Index _r2_list[],
                                   const int _n_members,
                                   const Tptcl* _ptcl) {
            for (int i=0; i<_n_members-1; i++) {
                int k = _ptcl_list[i];
                int jc=-1;
                Float r2min = NUMERIC_FLOAT_MAX;
                for(int j=i+1; j<_n_members; j++) {
                    const int kj =_ptcl_list[j];
                    Vector3<Float> dr = {_ptcl[k].pos[0] - _ptcl[kj].pos[0],
                                         _ptcl[k].pos[1] - _ptcl[kj].pos[1],
                                         _ptcl[k].pos[2] - _ptcl[kj].pos[2]};
                    Float r2 = dr*dr;
                    if(r2<r2min) {
                        r2min = r2;
                        jc = j;
                    }
                }
#ifdef BINARY_DEBUG
                ASSERT(jc>=0);
#endif
                if (jc!=i+1) {
                    int jtmp = _ptcl_list[i+1];
                    _ptcl_list[i+1] = _ptcl_list[jc];
                    _ptcl_list[jc]  = jtmp;
                }
                _r2_list[i].first = r2min;
                _r2_list[i].second = i;
            }
        }
    

    public:

        //! constructor
        BinaryTree(): Tptcl(), Tbinary(), n_members(-1), member_index{-3,-3}, member{NULL,NULL}, level(-1), branch(-1) {}

        //! allocate memory for two members as particles
        /*! allocate memory as particle types for two members
            set member index to -2
         */
        void allocateParticleMember() {
            ASSERT(member_index[0]!=-2);
            ASSERT(member_index[1]!=-2);
            member[0] = new Tptcl();
            member[1] = new Tptcl();
            member_index[0] = -2;
            member_index[1] = -2;
        }

        //! clear memory allocation
        void clear() {
            for (int i=0; i<2; i++) {
                if (member_index[i] == -2) { 
                    delete member[i];
                    member[i] = NULL;
                    member_index[i] = -3;
                }
            }
        }

        //! destructor
        ~BinaryTree() {
            clear();
        }

        // copy function
        //BinaryTreeLocal & operator = (const BinaryTreeLocal& _bin) {
        //    std::memcpy(this, &_bin, sizeof(Tbinary)+sizeof(Tptcl)+sizeof(int));
        //    return *this;
        //}

        //! Three-body stability function for hierarchical binary tree 
        /* Use Myllaeri et al. (2018, MNRAS, 476, 830) 3-body stability criterion to check whether any B-S/B-B sub-system is stable, return the maximim stabliity factor
           @param[in] _bin: binary tree to check
           @param[in] _t_crit: time interval for stable check
           \return maximum stability factor <1 stable; >1 unstable; =0 binary system
        */
        static Float& stableCheckIter(BinaryTreeLocal& _bin, const Float _t_crit) {
            if (_bin.semi<0) _bin.stab=NUMERIC_FLOAT_MAX;
            else _bin.stab = 0.0;
            Float mout = _bin.mass;

            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    BinaryTreeLocal* bin_in = _bin.getMemberAsTree(k);
                    _bin.stab = std::max(_bin.stab,stableCheckIter(*bin_in, _t_crit));
                    if (_bin.stab<1.0) {
                        // if zero mass exist, set unstable
                        if (_bin.getMember(1-k)->mass ==0.0) {
                            _bin.stab = std::max(_bin.stab, Float(1000.0));
                        }
                        else {

                            Float incline=acos(std::min(Float(1.0), _bin.am*bin_in->am/sqrt((_bin.am*_bin.am)*(bin_in->am*bin_in->am))));
                            
                            Float cosinc = cos(incline);
                            Float fac = 1.0 - 2.0*bin_in->ecc/3.0 * (1.0 - 0.5*bin_in->ecc*bin_in->ecc) 
                                - 0.3*cosinc*(1.0 - 0.5*bin_in->ecc + 2.0*cosinc*(1.0 - 2.5*pow(bin_in->ecc,1.5) - cosinc));

                            Float min = bin_in->mass;
                            Float g = sqrt(std::max(bin_in->m1,bin_in->m2) /min)*(1.0 + mout/min);
    
                            Float qst = 1.52*pow(sqrt(_t_crit/_bin.period)/(1.0 - _bin.ecc),1.0/6.0)*pow(fac*g,1.0/3.0);

                            Float peri_out = _bin.semi * (_bin.ecc + 1.0);
                            Float q = peri_out/bin_in->semi;
                    
                            bin_in->stab = std::max(bin_in->stab, qst/q);
                            _bin.stab = std::max(_bin.stab, qst/q);
                        }
                    }
                }
            }
            return _bin.stab;
        }

        //! generate level and branch
        /*! @param[in] _bin: binary tree
            @param[in] _reset: if true, set _bin as root tree (level=0, branch=0), otherwise use the stored values.
         */
        static void generateLevelBranch(BinaryTreeLocal& _bin, const bool _reset) {
            if (_reset) {
                _bin.level=0;
                _bin.branch=0;
            }
            else {
                ASSERT(_bin.level>=0&&_bin.branch>=0);
            }
            for (int i=0; i<2; i++) {
                if(_bin.isMemberTree(i)) {
                    auto* bini = _bin.getMemberAsTree(i);
                    bini->level = _bin.level+1;
                    bini->branch = 2*_bin.branch + i;
                    generateLevelBranch(*bini, false);
                }
            }
        };

        //! check whether a given binarytree data is in the same branch 
        /*! @param[in] _bin: binary tree
          \return 0: no relation; 1: same binary; -2: the given _bin is in the upper root tree; 2: the given _bin is in the leaf tree
         */        
        int isSameBranch(const BinaryTreeLocal& _bin) const {
            ASSERT(level>=0&&branch>=0&&_bin.level>=0&&_bin.branch>=0);
            int dlevel = _bin.level - level;
            if (dlevel>0) { // check leaf
                int left_edge = (branch<<dlevel);
                int right_edge = ((branch+1)<<dlevel);
                if (_bin.branch>=left_edge&&_bin.branch<right_edge) return 2;
                else return 0;
            }
            else if (dlevel<0) { // check mother
                int left_edge = (_bin.branch<<-dlevel);
                int right_edge = ((_bin.branch+1)<<-dlevel);
                if (branch>=left_edge&&branch<right_edge) return -2;
                else return 0;
            }
            else if(_bin.branch==branch) return 1;
            else return 0;
        }

        //! Generate kepler binary tree for a group of particles
        /* 
           @param[in,out] _bins: binary tree, size of n_members
           @param[in,out] _ptcl_list: group particle member index in _ptcl, will be reordered by the minimum distance chain.
           @param[in] _n: number of particles in _ptcl_list
           @param[in] _ptcl: particle data array
           @param[in] _G: gravtational constant
           \return binary tree number
        */
        static void generateBinaryTree(BinaryTreeLocal _bins[],  // make sure bins.size = n_members-1!
                                       int _ptcl_list[],   // make sure list.size = n_members!
                                       const int _n,
                                       Tptcl* _ptcl,
                                       const Float _G) {
            ASSERT(_n>=2);

            R2Index r2_list[_n];
            // reorder _ptcl_list by minimum distance of each particles, and save square minimum distance r2min and index i in _ptcl_list (not particle index in _ptcl) in r2_list 
            calcMinDisList(_ptcl_list, r2_list, _n, _ptcl);
            // sort r2_list by r2min 
            if(_n>2) std::sort(r2_list, r2_list+_n-1, pairLess);

            // tree root for each binary pair
            BinaryTreeLocal* bin_host[_n]; 
            for(auto &p : bin_host) p=NULL;
    
            // check binary from the closest pair to longest pair
            for(int i=0; i<_n-1; i++) {
                // get the pair index k in _ptcl_list with sorted r2_list, i represent the sorted r2min order.
                int k = r2_list[i].second;
                Tptcl* p[2];
                int pindex[2]={-1,-1};
                // if no tree root assign, set member 1 to particle and their host to current bins i
                if(bin_host[k]==NULL) {
#ifdef BINARY_DEBUG
                    ASSERT(_ptcl_list[k]>=0);
#endif
                    p[0] = &_ptcl[_ptcl_list[k]];
                    pindex[0] = _ptcl_list[k];
                    bin_host[k] = &_bins[i];
                }

                // if tree root already exist, member 1 assigned to tree root bin_host[k], tree members' bin_host -> current bins i
                else {
                    p[0] = bin_host[k];
                    int ki = k;
                    while(bin_host[ki]==p[0]&&ki>=0) bin_host[ki--] = &_bins[i];
                }

                // if no tree root assign, set member 2 to particle and their host to current bins i
                if(bin_host[k+1]==NULL) {
#ifdef BINARY_DEBUG
                    ASSERT(_ptcl_list[k+1]>0);
#endif
                    p[1] = &_ptcl[_ptcl_list[k+1]];
                    pindex[1] = _ptcl_list[k+1];
                    bin_host[k+1] = &_bins[i];
                }

                // if tree root already exist, member 2 assigned to tree root bin_host[k], tree members' bin_host -> current bins i
                else {
                    p[1] = bin_host[k+1];
                    int ki = k+1;
                    while(bin_host[ki]==p[1]&&ki>=0) bin_host[ki++] = &_bins[i];
                }

                // calculate binary parameter
                _bins[i].setMembers(p[0], p[1], pindex[0], pindex[1]);
                // calculate kepler orbit
                _bins[i].calcOrbit(_G);
                // calculate center-of-mass 
                _bins[i].calcCenterOfMass();
            }

            // generate level branch
            generateLevelBranch(_bins[_n-2],true);

#ifdef BINARY_DEBUG
            for(int i=0; i<_n; i++) ASSERT(bin_host[i]==&_bins[_n-2]); // check whether all bin_host point to the last of bins
            ASSERT(_bins[_n-2].getMemberN() == _n); // makesure the particle number match the binarytree recored member
#endif
        }

        //! collect BinaryTree data iteratively and save to a BinaryTree array
        /*! The member number should be correct 
          @param[out] _bin: binarytree array to store the data (should be size of n_members-1), the root is put last (n_members-1)
        */
        void getherBinaryTreeIter(BinaryTreeLocal _bin[]) {
            copyDataIter(_bin, n_members-2);
        }

        //! calculate center-of-mass from members
        void calcCenterOfMass() {
#ifdef BINARY_DEBUG
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
#endif
            Float m1 = member[0]->mass;
            Float m2 = member[1]->mass;
            this->pos[0] = m1*member[0]->pos[0] + m2*member[1]->pos[0];
            this->pos[1] = m1*member[0]->pos[1] + m2*member[1]->pos[1];
            this->pos[2] = m1*member[0]->pos[2] + m2*member[1]->pos[2];
            this->vel[0] = m1*member[0]->vel[0] + m2*member[1]->vel[0];
            this->vel[1] = m1*member[0]->vel[1] + m2*member[1]->vel[1];
            this->vel[2] = m1*member[0]->vel[2] + m2*member[1]->vel[2];
            this->mass = m1+m2;
            this->pos[0] /=this->mass;
            this->pos[1] /=this->mass;
            this->pos[2] /=this->mass;
            this->vel[0] /=this->mass;
            this->vel[1] /=this->mass;
            this->vel[2] /=this->mass;
        }

        //! calculate Kepler orbit semi, ecc and period only
        void calcSemiEccPeriod(const Float& _G) {
            Tbinary::particleToSemiEccPeriod(*member[0], *member[1], _G);
        }

        //! calculate Kepler orbit from members
        void calcOrbit(const Float _G) {
            Tbinary::calcOrbit(*member[0], *member[1], _G);
        }

        //! calculate particles from orbit
        void calcParticles(const Float _G) {
            Tbinary::calcParticles(*member[0], *member[1], _G);
        }

        //! calc total number of members
        int calcMemberN() {
#ifdef BINARY_DEBUG
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
#endif
            n_members = 0;
            if (member_index[0]==-1) n_members += ((BinaryTreeLocal*)member[0])->getMemberN();
            else n_members++;
            if (member_index[1]==-1) n_members += ((BinaryTreeLocal*)member[1])->getMemberN();
            else n_members++;
            return n_members;
        }

        //! get total number of members
        int getMemberN() const {
            return n_members;
        }

        //! calc total number of members iteratively
        int calcMemberNIter() {
#ifdef BINARY_DEBUG
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
#endif
            n_members = 0;
            if (member_index[0]==-1) n_members += (BinaryTreeLocal*)member[0]->calcMemberNIter();
            else n_members++;
            if (member_index[1]==-1) n_members += (BinaryTreeLocal*)member[1]->calcMemberNIter();
            else n_members++;
            return n_members;
        }

        //! iteratively leaf processing function template
        template <class T>
        using ProcessFunctionLeaf = void (*) (T&, Tptcl* &);
    
        //! Process (bottom) leaf data with extra dat iteratively (from left to right)
        /*!
          @param[in] _dat: data for processing loop from left to right
          @param[in] _f: function to process leaf information
          \return: new data generated by processing function
        */
        template <class T>
        void processLeafIter(T& _dat, ProcessFunctionLeaf<T> _f) {
#ifdef BINARY_DEBUG
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
#endif
            if (member_index[0]==-1) ((BinaryTreeLocal*)member[0])->processLeafIter(_dat, _f);
            else _f(_dat, member[0]);
            if (member_index[1]==-1) ((BinaryTreeLocal*)member[1])->processLeafIter(_dat, _f);
            else _f(_dat, member[1]);
        }

        //! iteratively root processing function template
        template <class T> 
        using ProcessFunctionRoot = T (*) (T&, BinaryTreeLocal& );

        //! Process root data and return result iteratively
        /*! The process go from top root to leafs from left to right
          @param[in] _dat: data return from upper level processing function 
          @param[in] _f: function to process root information
          \return: new data generated by processing function
        */
        template <class T>
        T processRootIter(T& _dat, ProcessFunctionRoot<T> _f){ 
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
            T dat_new = _f(_dat, *this);
            for (int k=0; k<2; k++) 
                if (member_index[k]==-1) dat_new = ((BinaryTreeLocal*)member[k])->processRootIter(dat_new, _f);
            return dat_new;
        }

        //! Process root and leaf data and return result iteratively
        /*! The process go from top root to leafs from left to right
          @param[in] _dat_root: data return from upper root processing function 
          @param[in] _f_root: function to process root information
          @param[in] _dat_leaf: data for processing loop from left to right
          @param[in] _f_leaf: function to process leaf information
          \return: new data generated by processing function
        */
        template <class Troot, class Tleaf>
        Troot processRootLeafIter(Troot& _dat_root, ProcessFunctionRoot<Troot> _f_root, Tleaf& _dat_leaf, ProcessFunctionLeaf<Tleaf> _f_leaf){ 
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
            Troot dat_new = _f_root(_dat_root, *this);
            for (int k=0; k<2; k++)  {
                if (member_index[k]==-1) dat_new = ((BinaryTreeLocal*)member[k])->processRootLeafIter(dat_new, _f_root, _dat_leaf, _f_leaf);
                else _f_leaf(_dat_leaf, member[k]);
            }
            return dat_new;
        }

        //! iteratively tree processing function template
        template <class T, class Tr>
        using ProcessFunctionTree = Tr (*) (T&, const Tr&, const Tr&, BinaryTreeLocal& );

        //! Process tree data with extra dat iteratively (from bottom to top)
        /*! 
          @param[in] _dat: general data type for process using
          @param[in] _res1: return from branch 1 
          @param[in] _res2: return from branch 2 
          @param[in] _f: function to process root information
          \return: new data generated by processing function
        */
        template <class T, class Tr>
        Tr processTreeIter(T& _dat, const Tr& _res1, const Tr& _res2, ProcessFunctionTree<T,Tr> _f){ 
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
            Tr res[2] = {_res1, _res2};
            for (int k=0; k<2; k++) 
                if (member_index[k]==-1) res[k] = ((BinaryTreeLocal*)member[k])->processTreeIter(_dat, _res1, _res2, _f);
            return _f(_dat, res[0], res[1], *this);
        }
    

        //! copy a binary tree to a binaryTree array iteratively
        /*! 
          @param [in] _bin_out: BinaryTree array to store the data
          @param [in] _i: index in _bin_out for copying the current tree cell
        */
        int copyDataIter(BinaryTreeLocal _bin_out[], const int _i) {
            int n_iter = 1;
            int inow = _i; // array index to fill (backwards moving)
            _bin_out[_i] = *this;
            for (int k=0; k<2; k++) {
                if (member_index[k]==-1) {
                    inow--;
                    ASSERT(inow>=0);
                    _bin_out[_i].member[k] = &_bin_out[inow];
                    n_iter += ((BinaryTreeLocal*)member[k])->copyDataIter(_bin_out, inow);
                }
            }
            return n_iter;
        }

        //! set members
        /*!
          @param[in] _p1: particle 1
          @param[in] _p2: particle 2
          @param[in] _i1: particle 1 index (should be >=0 if it is not binary tree)
          @param[in] _i2: particle 2 index (should be >=0 if it is not binary tree)
         */
        void setMembers(Tptcl* _p1, Tptcl* _p2, const int _i1, const int _i2) {
            ASSERT(member_index[0]!=-2);
            ASSERT(member_index[1]!=-2);
            member[0] = _p1;
            member[1] = _p2;
            member_index[0] = _i1;
            member_index[1] = _i2;
            calcMemberN();
        }

        //! get member
        Tptcl* getMember(const size_t i) const {
            ASSERT(i<2);
            return member[i];
        }

        //! get member as tree
        BinaryTreeLocal* getMemberAsTree(const size_t i) const {
            ASSERT(i<2);
            return (BinaryTreeLocal*)member[i];
        }

        //! get member index
        int getMemberIndex(const size_t i) const {
            ASSERT(i<2);
            return member_index[i];
        }

        bool isMemberTree(const size_t i) const {
            ASSERT(i<2);
            return (member_index[i]==-1);
        }

        //! get left member
        Tptcl* getLeftMember() const {
            return member[0];
        }

        //! get left member as tree
        BinaryTreeLocal* getLeftMemberAsTree() const {
            return (BinaryTreeLocal*)member[0];
        }

        //! get right member
        Tptcl* getRightMember() const {
            return member[1];
        }

        //! get right member as tree
        BinaryTreeLocal* getRightMemberAsTree() const {
            return (BinaryTreeLocal*)member[1];
        }

        ////! copy operator = 
        ///*! Copy the data and member address. If member address point to a particle (id>0), copy address; else if it is a binary tree, copy the address and also correct the address difference, this assume the whole BinaryTree is stored in a continuing array, thus correct address difference also make the new binary tree have consistent member address in new array.
        // */
        //BinaryTreeLocal & operator = (const BinaryTreeLocal& _bin) {
        //    // copy everything first
        //    std::memcpy(this, &_bin, sizeof(BinaryTreeLocal));
        //    // get Different of address
        //    const BinaryTreeLocal* adr_diff = (const BinaryTreeLocal*)this - &_bin; 
        //    // correct member address
        //    for (int k=0; k<2; k++) {
        //        if (_bin.member[k]!=NULL)  {
        //            // if a member is another binarytree
        //            if (_bin.member_index[k]==-1) {
        //                // Add the address difference to make the new member address consistent in the new binary tree array
        //                member[k] += adr_diff;
        //            }
        //        }
        //    }
        //    return *this;
        //}

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        static void printColumnTitle(std::ostream & _fout, const int _width=20) {
            Tptcl::printColumnTitle(_fout, _width);
            Tbinary::printColumnTitle(_fout, _width);
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            Tptcl::printColumn(_fout, _width);
            Tbinary::printColumn(_fout, _width);
        }

        //! write class data to file with ASCII format
        /*! @param[in] _fout: std:osteram file for output
         */
        void writeAscii(std::ostream& _fout) const {
            Tptcl::writeAscii(_fout);
            Tbinary::writeAscii(_fout);
        }

        //! read class data to file with ASCII format
        /*! @param[in] _fin: std::istream file for input
         */
        void readAscii(std::istream&  _fin) { 
            Tptcl::readAscii(_fin);
            Tbinary::readAscii(_fin);
        }       

        //! print binary and member information
        void printBinaryTreeIter(std::ostream & _fout, const int _width=20){
            Tbinary::printColumn(_fout, _width);
            for (int k=0; k<2; k++) member[k]->printColumn(_fout,_width);
            for (int k=0; k<2; k++) 
                if (member_index[k]==-1) ((BinaryTreeLocal*)member[k])->printBinaryTreeIter(_fout, _width);
        }

        //! print binary function
        void printBinaryIter(std::ostream & _fout, const int _width=20){
            printColumn(_fout, _width);
            for (int k=0; k<2; k++) {
                if (member_index[k]==-1) ((BinaryTreeLocal*)member[k])->printBinaryIter(_fout, _width);
            }
        }
    };
}
