#pragma once

#include<cmath>
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
        Vector3<Float> am;        // angular momentum


    private:
        //! Orbit to position and velocity
        /*! @param[out]: _p1: particle 1
          @param[out]: _p2: particle 2
          @param[in]: _bin: binary parameter
        */
        template <class Tptcl>
        void orbitToParticle(Tptcl& _p1, Tptcl& _p2, const Binary& _bin, const Float& _ecca) const {
            Float m_tot = _bin.m1 + _bin.m2;
            Float n = sqrt( m_tot / (_bin.semi*_bin.semi*_bin.semi) );
            Float cosu = cos(_ecca);
            Float sinu = sin(_ecca);
            Float c0 = sqrt(1.0 - _bin.ecc*_bin.ecc);
            Vector3<Float> pos_star(_bin.semi*(cosu - _bin.ecc), _bin.semi*c0*sinu, 0.0);
            Vector3<Float> vel_star(-_bin.semi*n*sinu/(1.0-_bin.ecc*cosu), _bin.semi*n*c0*cosu/(1.0-_bin.ecc*cosu), 0.0);
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
        /* @param[out]: _bin: binary parameter
           @param[in]:  _p1: particle 1
           @param[in]:  _p2: particle 2
        */
        template <class Tptcl>
        void particleToOrbit(Binary& _bin, const Tptcl& _p1, const Tptcl& _p2){
            _bin.m1 = _p1.mass;
            _bin.m2 = _p2.mass;
            Float m_tot = _p1.mass + _p2.mass;
            Vector3<Float> pos_red(_p2.pos[0] - _p1.pos[0], _p2.pos[1] - _p1.pos[1], _p2.pos[2] - _p1.pos[2]);
            Vector3<Float> vel_red(_p2.vel[0] - _p1.vel[0], _p2.vel[1] - _p1.vel[1], _p2.vel[2] - _p1.vel[2]);
            Float r_sq = pos_red * pos_red;
            r = sqrt(r_sq);
            Float inv_dr = 1.0 / r;
            Float v_sq = vel_red * vel_red;
            _bin.semi = 1.0 / (2.0*inv_dr - v_sq / m_tot);
            //    ASSERT(semi > 0.0);
            _bin.am = pos_red ^ vel_red;
            _bin.incline = atan2( sqrt(_bin.am.x*_bin.am.x+_bin.am.y*_bin.am.y), _bin.am.z);
            _bin.rot_horizon = atan2(_bin.am.x, -_bin.am.y);
     
            Vector3<Float> pos_bar, vel_bar;
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
            Float h = sqrt(_bin.am*_bin.am);
            Float ecccosomg =  h/m_tot*vel_bar.y - pos_bar.x*inv_dr;
            Float eccsinomg = -h/m_tot*vel_bar.x - pos_bar.y*inv_dr;
            _bin.ecc = sqrt( ecccosomg*ecccosomg + eccsinomg*eccsinomg );
            _bin.rot_self = atan2(eccsinomg, ecccosomg);
            Float phi = atan2(pos_bar.y, pos_bar.x); // f + omg (f: true anomaly)
            Float f = phi - _bin.rot_self;
            Float sinu = r*sin(f) / (_bin.semi*sqrt(1.0 - _bin.ecc*_bin.ecc));
            Float cosu = (r*cos(f) / _bin.semi) + _bin.ecc;
            _bin.ecca = atan2(sinu, cosu); // eccentric anomaly
            Float n = sqrt(m_tot/(_bin.semi*_bin.semi*_bin.semi)); // mean motion
            _bin.period = 8.0*std::atan(1.0)/n;
            Float l = _bin.ecca - _bin.ecc*sin(_bin.ecca);  // mean anomaly
            _bin.t_peri = l / n; 
        }

        //! calculate eccentricy anomaly from mean anomaly
        /*
          @param[in] _mean_anomaly: mean_anomaly
          @param[in] _ecc: eccentricity
          \return eccentric anomaly
        */
        Float calcEccAnomaly(const Float _mean_anomaly,
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
        /*!
          @param[in,out] _bin: binary orbit
          @param[in] _dt: evolution time
        */
        void solveKepler(Binary& _bin,
                         const Float _dt) {
            Float freq = sqrt( (_bin.m1+_bin.m2) / (_bin.semi*_bin.semi*_bin.semi) );
            Float mean_anomaly_old = _bin.ecca - _bin.ecc * sin(_bin.ecca);
            Float dt_tmp = _dt - to_int(_dt/_bin.period)*_bin.period;
            Float mean_anomaly_new = freq * dt_tmp + mean_anomaly_old; // mean anomaly
            _bin.ecca = calcEccAnomaly(mean_anomaly_new, _bin.ecc); // eccentric anomaly
        }
    
    public:

        //! calculate kepler Orbit from particles
        /*! 
          @param[out]: _p1: particle 1
          @param[out]: _p2: particle 2
        */    
        template <class Tptcl>
        void calcOrbit(const Tptcl& _p1, const Tptcl& _p2) {
            particleToOrbit(*this, _p1, _p2);
        }

        //! calculate two components from kepler Orbit
        /*! 
          @param[out]: _p1: particle 1
          @param[out]: _p2: particle 2
        */    
        template <class Tptcl>
        void calcParticles(Tptcl& _p1, Tptcl& _p2) {
            orbitToParticle(_p1, _p2, *this, this->ecca);
        }

        //! calculate two components from kepler Orbit with input eccentricity anomaly
        /*! 
          @param[out]: _p1: particle 1
          @param[out]: _p2: particle 2
          @param[in]: _ecca: eccentricity anomaly
        */    
        template <class Tptcl>
        void calcParticlesEcca(Tptcl& _p1, Tptcl& _p2, const Float _ecca) const {
            orbitToParticle(_p1, _p2, *this, _ecca);
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
               <<"r:           separation: "<<r<<std::endl;
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
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
                 <<std::setw(_width)<<"r";
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
                 <<std::setw(_width)<<r;
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
/*! Require id as a member in Tptcl.
  Binary id will be the negative value of first member id
*/
    template <class Tptcl>
    class BinaryTree: public Tptcl, public Binary {
    private: 
        typedef std::pair<Float, int> R2Index;  // First is distance square, second is index

        int n_members; ///> Total member number belong to the tree (can be more than 2)
        Tptcl* member[2]; ///> member pointer
    
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
#ifdef KEPLER_DEBUG
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
        BinaryTree(): Tptcl(), Binary(), n_members(-1), member{NULL,NULL} {}


        // copy function
        //BinaryTree<Tptcl> & operator = (const BinaryTree<Tptcl>& _bin) {
        //    std::memcpy(this, &_bin, sizeof(Binary)+sizeof(Tptcl)+sizeof(int));
        //    return *this;
        //}


        //! Generate kepler binary tree for a group of particles
        /* Notice binarytree id is (- first member id), thus particle id must be positive to avoid error
           @param[in,out] _bins: binary tree, size of n_members
           @param[in,out] _ptcl_list: group particle member index in _ptcl, will be reordered by the minimum distance chain.
           @param[in] _n: number of particles in _ptcl_list
           @param[in] _ptcl: particle data array
           \return binary tree number
        */
        static void generateBinaryTree(BinaryTree<Tptcl> _bins[],  // make sure bins.size = n_members-1!
                                       int _ptcl_list[],   // make sure list.size = n_members!
                                       const int _n,
                                       Tptcl* _ptcl) {

            R2Index r2_list[_n];
            // reorder _ptcl_list by minimum distance of each particles, and save square minimum distance r2min and index i in _ptcl_list (not particle index in _ptcl) in r2_list 
            calcMinDisList(_ptcl_list, r2_list, _n, _ptcl);
            // sort r2_list by r2min 
            if(_n>2) std::sort(r2_list, r2_list+_n-1, pairLess);

            // tree root for each binary pair
            BinaryTree<Tptcl>* bin_host[_n]; 
            for(auto &p : bin_host) p=NULL;
    
            // check binary from the closest pair to longest pair
            for(int i=0; i<_n-1; i++) {
                // get the pair index k in _ptcl_list with sorted r2_list, i represent the sorted r2min order.
                int k = r2_list[i].second;
                Tptcl* p[2];
                // if no tree root assign, set member 1 to particle and their host to current bins i
                if(bin_host[k]==NULL) {
                    p[0] = &_ptcl[_ptcl_list[k]];
#ifdef KEPLER_DEBUG
                    ASSERT(p[0]->id>0);
#endif
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
                    p[1] = &_ptcl[_ptcl_list[k+1]];
#ifdef KEPLER_DEBUG
                    ASSERT(p[1]->id>0);
#endif
                    bin_host[k+1] = &_bins[i];
                }

                // if tree root already exist, member 2 assigned to tree root bin_host[k], tree members' bin_host -> current bins i
                else {
                    p[1] = bin_host[k+1];
                    int ki = k+1;
                    while(bin_host[ki]==p[1]&&ki>=0) bin_host[ki++] = &_bins[i];
                }

                // calculate binary parameter
                _bins[i].setMembers(p[0], p[1]);
                // calculate kepler orbit
                _bins[i].calcOrbit();
                // calculate center-of-mass 
                _bins[i].calcCenterOfMass();
            }

#ifdef KEPLER_DEBUG
            for(int i=0; i<_n; i++) ASSERT(bin_host[i]==&_bins[_n-2]); // check whether all bin_host point to the last of bins
            ASSERT(_bins[_n-2].getMemberN() == _n); // makesure the particle number match the binarytree recored member
#endif
        }


        //! collect BinaryTree data iteratively and save to a BinaryTree array
        /*! The member number should be correct 
          @param[out] _bin: binarytree array to store the data (should be size of n_members-1), the root is put last (n_members-1)
        */
        void getherBinaryTreeIter(BinaryTree<Tptcl> _bin[]) {
            copyDataIter(_bin, n_members-2);
        }

        //! calculate center-of-mass from members
        void calcCenterOfMass() {
#ifdef KEPLER_DEBUG
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
            ASSERT(member[0]->id!=0);
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
            this->id = - std::abs(member[0]->id);
        }

        //! calculate Kepler orbit from members
        void calcOrbit() {
            Binary::calcOrbit(*member[0], *member[1]);
        }

        //! calc total number of members
        int calcMemberN() {
#ifdef KEPLER_DEBUG
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
#endif
            n_members = 0;
            if (member[0]->id<0) n_members += ((BinaryTree<Tptcl>*)member[0])->getMemberN();
            else n_members++;
            if (member[1]->id<0) n_members += ((BinaryTree<Tptcl>*)member[1])->getMemberN();
            else n_members++;
            return n_members;
        }

        //! get total number of members
        int getMemberN() const {
            return n_members;
        }

        //! calc total number of members iteratively
        int calcMemberNIter() {
#ifdef KEPLER_DEBUG
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
#endif
            n_members = 0;
            if (member[0]->id<0) n_members += (BinaryTree<Tptcl>*)member[0]->calcMemberNIter();
            else n_members++;
            if (member[1]->id<0) n_members += (BinaryTree<Tptcl>*)member[1]->calcMemberNIter();
            else n_members++;
            return n_members;
        }

        template <class T>
        using ProcessFunctionLeaf = void (*) (T&, Tptcl* &);
    
        template <class T>
        void processLeafIter(T& _c, ProcessFunctionLeaf<T> _f) {
#ifdef KEPLER_DEBUG
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
#endif
            if (member[0]->id<0) ((BinaryTree<Tptcl>*)member[0])->processLeafIter(_c, _f);
            else _f(_c, member[0]);
            if (member[1]->id<0) ((BinaryTree<Tptcl>*)member[1])->processLeafIter(_c, _f);
            else _f(_c, member[1]);
        }

        //! iteratively root processing function template
        template <class T> 
        using ProcessFunctionRoot = T (*) (const T&, BinaryTree<Tptcl>& );

        //! Process root data and return result iteratively
        /*! The process go from top root to leafs from left to right
          @param[in] _dat: data return from upper level processing function 
          @param[in] _f: function to process root information
          \return: new data generated by processing function
        */
        template <class T>
        T processRootIter(const T& _dat, ProcessFunctionRoot<T> _f){ 
            ASSERT(member[0]!=NULL);
            ASSERT(member[1]!=NULL);
            T dat_new = _f(_dat, *this);
            for (int k=0; k<2; k++) 
                if (member[k]->id<0) dat_new = ((BinaryTree<Tptcl>*)member[k])->processRootIter(dat_new, _f);
            return dat_new;
        }

    

        // copy a binary tree to a binaryTree array iteratively
        /*! 
          @param [in] _bin_out: BinaryTree array to store the data
          @param [in] _i: index in _bin_out for copying the current tree cell
        */
        int copyDataIter(BinaryTree<Tptcl> _bin_out[], const int _i) {
            int n_iter = 1;
            int inow = _i; // array index to fill (backwards moving)
            _bin_out[_i] = *this;
            for (int k=0; k<2; k++) {
                if (member[k]->id<0) {
                    inow--;
                    ASSERT(inow>=0);
                    _bin_out[_i].member[k] = &_bin_out[inow];
                    n_iter += ((BinaryTree<Tptcl>*)member[k])->copyDataIter(_bin_out, inow);
                }
            }
            return n_iter;
        }

        //! set members
        void setMembers(Tptcl* _p1, Tptcl* _p2) {
            member[0] = _p1;
            member[1] = _p2;
            calcMemberN();
        }

        //! get member
        Tptcl* getMember(const size_t i) const {
            ASSERT(i<2);
            return member[i];
        }

        //! get left member
        Tptcl* getLeftMember() const {
            return member[0];
        }

        //! get right member
        Tptcl* getRightMember() const {
            return member[1];
        }

        ////! copy operator = 
        ///*! Copy the data and member address. If member address point to a particle (id>0), copy address; else if it is a binary tree, copy the address and also correct the address difference, this assume the whole BinaryTree is stored in a continuing array, thus correct address difference also make the new binary tree have consistent member address in new array.
        // */
        //BinaryTree<Tptcl> & operator = (const BinaryTree<Tptcl>& _bin) {
        //    // copy everything first
        //    std::memcpy(this, &_bin, sizeof(BinaryTree<Tptcl>));
        //    // get Different of address
        //    const BinaryTree<Tptcl>* adr_diff = (const BinaryTree<Tptcl>*)this - &_bin; 
        //    // correct member address
        //    for (int k=0; k<2; k++) {
        //        if (_bin.member[k]!=NULL)  {
        //            // if a member is another binarytree
        //            if (_bin.member[k]->id<0) {
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
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            Binary::printColumnTitle(_fout, _width);
            Tptcl::printColumnTitle(_fout, _width);
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            Binary::printColumn(_fout, _width);
            Tptcl::printColumn(_fout, _width);
        }

        //! print function
        void printColumnIter(std::ostream & _fout, const int _width=20){
            printColumn(_fout, _width);
            _fout<<std::endl;
            for (int k=0; k<2; k++) {
                if (member[k]->id<0) ((BinaryTree<Tptcl>*)member[k])->printColumnIter(_fout, _width);
            }
        }
    };
}
