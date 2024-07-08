#pragma once
#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "AR/slow_down.h"

namespace AR {

    //! binary parameter with slowdown
    class BinarySlowDown: public COMM::Binary {
    public:
        SlowDown slowdown;
        Float stab_check_time;

        //! write class data to file with binary format
        /*! @param[in] _fp: FILE type file for output
         */
        void writeBinary(FILE *_fp) const {
            fwrite(this, sizeof(*this),1,_fp);
        }

        //! read class data to file with binary format
        /*! @param[in] _fin FILE type file for reading
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(this, sizeof(*this),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        static void printColumnTitle(std::ostream & _fout, const int _width=20) {
            Binary::printColumnTitle(_fout, _width);
            SlowDown::printColumnTitle(_fout,_width);
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            Binary::printColumn(_fout, _width);
            slowdown.printColumn(_fout,_width);
        }

        //! write class data to file with ASCII format
        /*! @param[in] _fout: std:osteram file for output
         */
        //void writeAscii(std::ostream& _fout) const {
        //    COMM::BinaryTree<Tparticle>::writeAscii(_fout);
        //    SlowDown::writeAscii(_fout);
        //}

        //! read class data to file with ASCII format
        /*! @param[in] _fin: std::istream file for input
         */
        //void readAscii(std::istream&  _fin) { 
        //    COMM::BinaryTree<Tparticle>::readAscii(_fin);
        //    SlowDown::readAscii(_fin);
        //}       

    };
    
    //! define ar binary tree
    template <class Tparticle>
    using BinaryTree=COMM::BinaryTree<Tparticle,BinarySlowDown>;

    //! Fix step options for integration with adjusted step (not for time sychronizatio phase)
    /*! always: use the given step without change \n
        later: fix step after a few adjustment of initial steps due to energy error
        none: don't fix step
     */
    enum class FixStepOption {always, later, none};

    //! A class contains information (e.g. parameters, binary tree, indices) about the particle group
    /*! The member of this class should not be the data that must be recored and should be possible calculated any time based on the main class (TimeTransformedSymplecticIntegrator) data. 
      This class must be inherited when a different Information class is applied in the main class, since the binary tree is used in the integration for slowdown factor
      The basic members are used in the integration are \\
      ds: integration step size \\
      binarytree: the Kepler orbital parameters of the hierarchical systems and slowdown factors \\
      fix_step_option: option to control whether the adjustment of step sizes are used \\
     */
    template <class Tparticle, class Tpcm>
    class Information{
    public:
        Float ds;  ///> initial step size for integration
        Float time_offset; ///> offset of time to obtain real physical time (real time = TimeTransformedSymplecticIntegrator:time_ + info.time_offset)
        Float r_break_crit;    // group break radius criterion
        FixStepOption fix_step_option; ///> fix step option for integration
        COMM::List<BinaryTree<Tparticle>> binarytree; ///> a list of binary tree that contain the hierarchical orbital parameters of the particle group.
#ifdef AR_DEBUG_DUMP
        bool dump_flag; ///> for debuging dump
#endif

        //! initializer, set ds to zero, fix_step_option to none
        Information(): ds(0.0), time_offset(0.0), r_break_crit(-1.0), fix_step_option(AR::FixStepOption::none), binarytree() {
#ifdef AR_DEBUG_DUMP
            dump_flag = false;
#endif
        }

        //! check whether parameters values are correct initialized
        /*! \return true: all correct
         */
        bool checkParams() {
            ASSERT(r_break_crit>=0.0);
            ASSERT(binarytree.getSize()>0);
            return true;
        }

        //! reserve memory of binarytree list
        void reserveMem(const int _nmax) {
            binarytree.setMode(COMM::ListMode::local);
            binarytree.reserveMem(_nmax);
        }

        //! get the root of binary tree
        BinaryTree<Tparticle>& getBinaryTreeRoot() const {
            int n = binarytree.getSize();
            ASSERT(n>0);
            return binarytree[n-1];
        }

        //! calculate ds for an elliptic orbit
        /*!
          For an ellpitic orbit, step ds=dt*G*m1*m2/r, the estimation of ds for one orbit is: 2*pi*sqrt(G*semi/(m1+m2))*m1*m2 
          @param[in] _bin: binary tree to check
          @param[in] _G: gravitational constant
          @param[in] _coff: coefficient for ds, default is 2*pi/32 (1/32 orbit)
         */
        inline Float calcDsElliptic(BinaryTree<Tparticle>& _bin, const Float& _G, const Float _coff = 0.19634954084) {
            return _coff*sqrt(_G*_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
        }

        //! calculate ds for a hyperbolic orbit
        /*!
          For an hyperbolic orbit, step ds=dt*G*m1*m2/r, the estimation of ds for one orbit is: 2*pi*sqrt(-G*semi/(m1+m2))*m1*m2 
          @param[in] _bin: binary tree to check
          @param[in] _G: gravitational constant
          @param[in] _coff: coefficient for ds, default is 2*pi/256 (1/256 orbit)
         */
        inline Float calcDsHyperbolic(BinaryTree<Tparticle>& _bin, const Float& _G, const Float _coff = 0.0245436926) {
            return _coff*sqrt(-_G*_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
        }

#ifdef AR_TIME_FUNCTION_MUL_POT
        //! iteration function to calculate summation and production of ds for all inner kepler orbits of a binary tree
        /*! 
          @param[out] ds_sum: summation of inner binaries' ds
          @param[out] ds_prod: production of inner binaries' ds
          @param[in] _bin: binary tree to check
          @param[in] _int_order: symplectic integrator accurate order
          @param[in] _G: gravitational constant
         */
        /*
        void calcSumProdDsKeplerIter(Float& ds_sum, Float& ds_prod, BinaryTree<Tparticle>& _bin, const int _int_order, const Float& _G) {
            if (_bin.getMemberN()>2) { // if not inner most kepler orbit, go inside members
                for (int k=0; k<2; k++) 
                    if (_bin.isMemberTree(k)) calcSumProdDsKeplerIter(ds_sum, ds_prod, *_bin.getMemberAsTree(k), _int_order, _G);
            }
            else {
                // zero mass cause ds=0
                if (_bin.m1>0 && _bin.m2>0) { 
                    // perturbation ratio
                    Float pert_ratio = (_bin.slowdown.pert_out>0&&_bin.slowdown.pert_in>0)? _bin.slowdown.pert_in/_bin.slowdown.pert_out: 1.0;
                    // scale step based on perturbation and sym method order
                    Float scale_factor = std::min(Float(1.0),pow(1e-1*pert_ratio,1.0/Float(_int_order)));

                    Float ds;
                    if (_bin.semi>0) ds = calcDsElliptic(_bin, _G, 1.0)*scale_factor;
                    else ds = calcDsHyperbolic(_bin, _G, 1.0);
                    ASSERT(ds<NUMERIC_FLOAT_MAX && ds>0);

                    ds_sum += ds;
                    ds_prod *= ds;
                }
            }
        }
        */

        //! calculate ds for a binary tree with multiple inner binaries/hyperbolics 
        /*! use calcSumProdDsKeplerIter with additional _ds_scale factor
          @param[in] _bin: binary tree to check
          @param[in] _int_order: symplectic integrator accurate order
          @param[in] _G: gravitational constant
          @param[in] _ds_scale: global scaling for final ds returned from calcProdDsKeplerIter (default: 1.0)
          
          \return ds: cumulative production of ds for all inner most kepler orbit of _bin
         */
        /*
        Float calcDsKeplerBinaryTree(BinaryTree<Tparticle>& _bin, const int _int_order, const Float& _G, const Float& _ds_scale = 1.0) {
            Float ds_sum = 0.0;
            Float ds_prod = 1.0;
            calcSumProdDsKeplerIter(ds_sum, ds_prod, _bin, _int_order, _G);
            // cofficient of 2*pi/32
            return 0.19634954084*_ds_scale*ds_prod/ds_sum;
        }
        */

        //! calculate ds for a binary tree with multiple inner binaries/hyperbolics 
        /*! calculate ds via formula: 1/ds = 1/ds1 + 1/ds2 + ..., where ds1, ds2 are ds estimated fro each inner most kepler orbit
          @param[in] _bin: binary tree to check
          @param[in] _int_order: symplectic integrator accurate order
          @param[in] _G: gravitational constant
          @param[in] _ds_scale: global scaling for final ds (default: 1.0)
          
          \return ds: ds estimated based on all inner most kepler orbit of _bin
         */
        Float calcDsKeplerBinaryTree(BinaryTree<Tparticle>& _bin, const int _int_order, const Float& _G, const Float& _ds_scale = 1.0) {
            Float ds = 0.0;
            if (_bin.getMemberN()>2) { // if not inner most kepler orbit, go inside members
                for (int k=0; k<2; k++) 
                    if (_bin.isMemberTree(k)) {
                        Float dsk = calcDsKeplerBinaryTree(*_bin.getMemberAsTree(k), _int_order, _G);
                        if (ds == 0.0) ds = dsk;
                        else if (dsk > 0) ds = ds*dsk/(ds+dsk);
                    }
            }
            else {
                // zero mass cause ds=0
                if (_bin.m1>0 && _bin.m2>0) { 
                    // perturbation ratio
                    Float pert_ratio = (_bin.slowdown.pert_out>0&&_bin.slowdown.pert_in>0)? _bin.slowdown.pert_in/_bin.slowdown.pert_out: 1.0;
                    // scale step based on perturbation and sym method order
                    Float scale_factor = std::min(Float(1.0),pow(1e-1*pert_ratio,1.0/Float(_int_order)));

                    if (_bin.semi>0) ds = calcDsElliptic(_bin, _G)*scale_factor;
                    else ds = calcDsHyperbolic(_bin, _G);
                    ASSERT(ds<NUMERIC_FLOAT_MAX && ds>0);
                }
            }

            // cofficient of 2*pi/32
            return _ds_scale*ds;
        }

#else        

        //! calculate ds for a binary tree with multiple inner binaries/hyperbolics 
        /*! calculate ds by select the minimum ds of binary/hyperbolic
          @param[in] _bin: binary tree to check
          @param[in] _int_order: symplectic integrator accurate order
          @param[in] _G: gravitational constant
          @param[in] _ds_scale: global scaling for final ds (default: 1.0), only multiplied once in the root level

          \return ds: minium ds * _ds_scale
         */
        Float calcDsKeplerBinaryTree(BinaryTree<Tparticle>& _bin, const int _int_order, const Float& _G, const Float& _ds_scale = 1.0) {
            Float ds_min = NUMERIC_FLOAT_MAX;
            if (_bin.getMemberN()>2) { // if not inner most kepler orbit, go inside members
                for (int k=0; k<2; k++) 
                    if (_bin.isMemberTree(k)) {
                        Float dsk = calcDsKeplerBinaryTree(*_bin.getMemberAsTree(k), _int_order, _G);
                        if (dsk < ds_min) ds_min = dsk;
                    }
            }
            else {
                // zero mass cause ds=0
                if (_bin.m1>0 && _bin.m2>0) { 
                    // perturbation ratio
                    Float pert_ratio = (_bin.slowdown.pert_out>0&&_bin.slowdown.pert_in>0)? _bin.slowdown.pert_in/_bin.slowdown.pert_out: 1.0;
                    // scale step based on perturbation and sym method order
                    Float scale_factor = std::min(Float(1.0),pow(1e-1*pert_ratio,1.0/Float(_int_order)));

                    if (_bin.semi>0) ds = calcDsElliptic(_bin, _G)*scale_factor;
                    else ds = calcDsHyperbolic(_bin, _G)*scale_factor;
                    ASSERT(ds<NUMERIC_FLOAT_MAX && ds>0);
                    if (ds < ds_min) ds_min = ds;
                }
            }
            return _ds_scale*ds_min;
        }

        //! iteration function to calculate average kepler ds for a binary tree
        /*
        void calcDsOverEbinMinKeplerIter(Float& _ds_over_ebin_min_bin, Float& _ds_min_bin, Float& _ds_min_hyp, Float& _etot_sd, const Float& _G, const Float& _nest_sd_up, BinaryTree<Tparticle>& _bin, const int _intergrator_order) {
            Float nest_sd = _nest_sd_up * _bin.slowdown.getSlowDownFactor();
            // perturbation ratio
            Float pert_ratio = (_bin.slowdown.pert_out>0&&_bin.slowdown.pert_in>0)? _bin.slowdown.pert_in/_bin.slowdown.pert_out: 1.0;
            // scale step based on perturbation and sym method order
            Float scale_factor = std::min(Float(1.0),pow(1e-1*pert_ratio,1.0/Float(_intergrator_order)));
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    calcDsOverEbinMinKeplerIter(_ds_over_ebin_min_bin, _ds_min_bin, _ds_min_hyp, _etot_sd, _G, nest_sd, *_bin.getMemberAsTree(k), _intergrator_order);
                }
            }
            // zero mass cause ds=0
            if (_bin.m1>0&&_bin.m2>0) {
                if (_bin.semi>0) {
                    Float dsi = calcDsElliptic(_bin, _G);

                    // scale by slowdown energy
                     Float ebin_sd = _G*(_bin.m1*_bin.m2)/(2*_bin.semi*nest_sd);
                    ASSERT(dsi>0&&ebin_sd>0);
                    Float ds_over_ebin = dsi*scale_factor/ebin_sd;
                    if (ds_over_ebin<_ds_over_ebin_min_bin) {
                        _ds_over_ebin_min_bin = ds_over_ebin;
                        _ds_min_bin = dsi*scale_factor;
                    }
                    _etot_sd += ebin_sd;
                }
                else {
                    Float dsi = calcDsHyperbolic(_bin, _G);
                    ASSERT(dsi>0);
                    //Float factor = std::min(Float(1.0), pow(nest_sd_org,Float(1.0/3.0)));
                    _ds_min_hyp = std::min(dsi, _ds_min_hyp);
                }
            }
        }
        */

        //! calculate average kepler ds iterately for a binary tree
        /*! use calcDsMinKeplerIter
         */
        /*
        Float calcDsKeplerBinaryTree(BinaryTree<Tparticle>& _bin, const int _int_order, const Float& _G, const Float& _ds_scale) {
            Float ds_over_ebin_min=NUMERIC_FLOAT_MAX;
            Float ds_min_hyp=NUMERIC_FLOAT_MAX;
            Float ds_min_bin=NUMERIC_FLOAT_MAX;
            Float etot_sd = 0.0;
            calcDsMinKeplerIter(ds_over_ebin_min, ds_min_bin, ds_min_hyp, etot_sd, _G, 1.0, _bin, _int_order);
            //Float ds_min_bin = ds_over_ebin_min*etot_sd/(bin_root.getMemberN()-1);
            ASSERT(ds_min_hyp<NUMERIC_FLOAT_MAX||ds_min_bin<NUMERIC_FLOAT_MAX);
            return _ds_scale*std::min(ds_min_bin,ds_min_hyp);
        }
        */

#endif

        //! calculate ds from the inner most binary with minimum period, determine the fix step option
        /*! Estimate ds first from the inner most binary orbit (eccentric anomaly), set fix_step_option to later
          @param[in] _int_order: accuracy order of the symplectic integrator.
          @param[in] _G: gravitational constant
          @param[in] _ds_scale: scaling factor to determine ds
         */
        void calcDsAndStepOption(const int _int_order, const Float& _G, const Float& _ds_scale) {
            auto& bin_root = getBinaryTreeRoot();
            ds = calcDsKeplerBinaryTree(bin_root, _int_order, _G, _ds_scale);
            ASSERT(ds>0);

            // Avoid too small step
            //if (_sd_org<1.0) ds *= std::max(1.0/8.0*pow(_sd_org, 1.0/Float(_int_order)),0.125);
            //auto& bin_root = getBinaryTreeRoot();
            const int n_particle = bin_root.getMemberN();

            // determine the fix step option
            //fix_step_option = FixStepOption::later;
            fix_step_option = FixStepOption::none;
            //// for two-body case, determine the step at begining then fix
            if (n_particle==2||bin_root.stab<1.0) fix_step_option = AR::FixStepOption::later;
            //// for multiple case, check whether outer peri-center is close to inner apo-center, if not, use fix step
            //if (n_particle>2) {
            //    Float apo_in_max = 0;
            //    for (int j=0; j<2; j++) {
            //        if (bin_root.isMemberTree(j)) {
            //            auto* bin_sub = (COMM::BinaryTree<Tparticle>*) bin_root.getMember(j);
            //            apo_in_max = std::max(apo_in_max,bin_sub->semi*(1+bin_sub->ecc));
            //        }
            //    }
            //    Float peri_out = bin_root.semi*(1-bin_root.ecc);
            //    if (peri_out>3*apo_in_max) fix_step_option = FixStepOption::later;
            //}
        }

        //! generate binary tree for the particle group
        /*! 
          - Construct the binary tree based on particle positions and velocities\n
          - If particle mass is 0.0, set to unused particles and put at the outermost orbits\n
          - If USE_CM_FRAME is used, switch each coordinate frame of the binary members to their center-of-the-mass frame\n

          @param[in] _particles: particle group 
          @param[in] _G: gravitational constant
         */
        void generateBinaryTree(COMM::ParticleGroup<Tparticle, Tpcm>& _particles, const Float _G) {
            // If the existing binarytree is not in original frame, first shift to original frame
#ifdef USE_CM_FRAME
            if (binarytree.getSize()>0) {
                auto& bin_root = getBinaryTreeRoot();
                if (!bin_root.isOriginFrame()) bin_root.shiftToOriginFrame();
            }
#endif

            const int n_particle = _particles.getSize();
            const int binary_tree_index = static_cast<int>(COMM::BinaryTreeMemberIndexTag::binarytree);
            ASSERT(n_particle>1);
            binarytree.resizeNoInitialize(n_particle-1);
            int particle_index_local[n_particle];
            int particle_index_unused[n_particle];
            int n_particle_real = 0;
            int n_particle_unused=0;
            for (int i=0; i<n_particle; i++) {
                if (_particles[i].mass>0.0) particle_index_local[n_particle_real++] = i;
                else particle_index_unused[n_particle_unused++] = i;
            }
            if (n_particle_real>1) 
                BinaryTree<Tparticle>::generateBinaryTree(binarytree.getDataAddress(), particle_index_local, n_particle_real, _particles.getDataAddress(), _G);

            // Add unused particles to the outmost orbit
            if (n_particle_real==0) {
                binarytree[0].setMembers(&(_particles[0]), &(_particles[1]), 0, 1);
                binarytree[0].mass = 0.0;
                binarytree[0].m1 = 0.0;
                binarytree[0].m2 = 0.0;  
                for (int i=2; i<n_particle; i++) {
                    binarytree[i-1].setMembers((Tparticle*)&(binarytree[i-2]), &( _particles[i]), binary_tree_index, i);
                    binarytree[i-1].mass = 0.0;
                    binarytree[i-1].m1 = 0.0;
                    binarytree[i-1].m2 = 0.0;
                }
            }
            else if (n_particle_real==1) {
                int i1 = particle_index_local[0];
                int i2 = particle_index_unused[0];
                binarytree[0].setMembers(&(_particles[i1]), &(_particles[i2]), i1 ,i2);
                binarytree[0].m1 = _particles[i1].mass;
                binarytree[0].m2 = 0.0;
                binarytree[0].mass = binarytree[0].m1;
                binarytree[0].pos[0] = _particles[i1].pos[0];
                binarytree[0].pos[1] = _particles[i1].pos[1];
                binarytree[0].pos[2] = _particles[i1].pos[2];
                binarytree[0].vel[0] = _particles[i1].vel[0];
                binarytree[0].vel[1] = _particles[i1].vel[1];
                binarytree[0].vel[2] = _particles[i1].vel[2];
                for (int i=1; i<n_particle_unused; i++) {
                    int k = particle_index_unused[i];
                    binarytree[i].setMembers((Tparticle*)&(binarytree[i-1]), &(_particles[k]), binary_tree_index, k);
                    binarytree[i].m1 = binarytree[i-1].mass;
                    binarytree[i].m2 = 0.0;
                    binarytree[i].mass = binarytree[i].m1;
                    binarytree[i].pos[0] = binarytree[i-1].pos[0];
                    binarytree[i].pos[1] = binarytree[i-1].pos[1];
                    binarytree[i].pos[2] = binarytree[i-1].pos[2];
                    binarytree[i].vel[0] = binarytree[i-1].vel[0];
                    binarytree[i].vel[1] = binarytree[i-1].vel[1];
                    binarytree[i].vel[2] = binarytree[i-1].vel[2];
                }
            }
            else {
                for (int i=0; i<n_particle_unused; i++) {
                    int ilast = n_particle_real-1+i;
                    ASSERT(ilast<n_particle-1);
                    int k = particle_index_unused[i];
                    binarytree[ilast].setMembers((Tparticle*)&(binarytree[ilast-1]), &(_particles[k]), binary_tree_index, k);
                    binarytree[ilast].m1 = binarytree[ilast-1].mass;
                    binarytree[ilast].m2 = 0.0;
                    binarytree[ilast].mass = binarytree[ilast].m1;
                    binarytree[ilast].pos[0] = binarytree[ilast-1].pos[0];
                    binarytree[ilast].pos[1] = binarytree[ilast-1].pos[1];
                    binarytree[ilast].pos[2] = binarytree[ilast-1].pos[2];
                    binarytree[ilast].vel[0] = binarytree[ilast-1].vel[0];
                    binarytree[ilast].vel[1] = binarytree[ilast-1].vel[1];
                    binarytree[ilast].vel[2] = binarytree[ilast-1].vel[2];
                }
            }

#ifdef USE_CM_FRAME
            getBinaryTreeRoot().shiftToCenterOfMassFrame();
#endif

        }

        //! check binary tree member pair id, if consisent, return ture. otherwise set the member pair id
        /*! 
          @param[in] _bin: binary tree to check
          @param[in] _reset_flag: if true, reset pair id to zero
        */
        bool checkAndSetBinaryPairIDIter(BinaryTree<Tparticle>& _bin, const bool _reset_flag) {
            bool return_flag=true;
            Tparticle* p[2] = {_bin.getLeftMember(), _bin.getRightMember()};

            for (int i=0; i<2; i++) {
                if (_bin.isMemberTree(i)) {
                    return_flag = return_flag & checkAndSetBinaryPairIDIter(*_bin.getMemberAsTree(i),_reset_flag);
                }
            }
            for (int i=0; i<2; i++) {
                if (!_bin.isMemberTree(i)) {
                    auto pair_id = p[1-i]->id;
                    return_flag = return_flag & (p[i]->getBinaryPairID()==pair_id);
                    if (_reset_flag) p[i]->setBinaryPairID(0);
                    else p[i]->setBinaryPairID(pair_id);
                }
            }
            if (p[0]->id<p[1]->id) _bin.id = p[0]->id;
            else _bin.id = p[1]->id;
            return return_flag;
        }

        //! clear function
        void clear() {
            ds=0.0;
            time_offset = 0.0;
            r_break_crit=-1.0;
            fix_step_option = FixStepOption::none;
            binarytree.clear();
#ifdef AR_DEBUG_DUMP
            dump_flag=false;
#endif
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"ds";
            _fout<<std::setw(_width)<<"Time_offset";
            _fout<<std::setw(_width)<<"r_break_crit";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<ds;
            _fout<<std::setw(_width)<<time_offset;
            _fout<<std::setw(_width)<<r_break_crit;
        }

        //! write class data to file with binary format
        /*! @param[in] _fout: FILE type file for output
         */
        void writeBinary(FILE *_fout) const {
            fwrite(&ds, sizeof(int),1,_fout);
            fwrite(&time_offset, sizeof(Float),1,_fout);
            fwrite(&r_break_crit, sizeof(Float),1,_fout);
            fwrite(&fix_step_option, sizeof(FixStepOption),1,_fout);
        }

        //! read class data to file with binary format
        /*! @param[in] _fin: FILE type file for reading
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(&ds, sizeof(int),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
            rcount = fread(&time_offset, sizeof(Float),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
            rcount = fread(&r_break_crit, sizeof(Float),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
            rcount = fread(&fix_step_option, sizeof(FixStepOption),1,_fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }    
    };

}
