#pragma once
#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "AR/slow_down.h"

namespace AR {

    //! binary parameter with slowdown
    class BinarySlowDown: public COMM::Binary {
    public:
        SlowDown slowdown;

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

    //! A class contains the Kepler orbital parameters and initial step size
    /*! The class has members to generate the Kepler orbital binary tree of the particle group and store the information in binarytree. \\
      It also provides the initial step size estimator and the option of fix step.
     */
    template <class Tparticle, class Tpcm>
    class Information{
    private:
        //! calculate average kepler ds for ARC
        void calcDsMinKeplerIter(Float& _ds_over_ebin_min_bin, Float& _ds_min_hyp, Float& _etot_sd, const Float& _G, const Float& _nest_sd_up, BinaryTree<Tparticle>& _bin) {
            Float nest_sd = _nest_sd_up * _bin.slowdown.getSlowDownFactor();
            for (int k=0; k<2; k++) {
                if (_bin.isMemberTree(k)) {
                    calcDsMinKeplerIter(_ds_over_ebin_min_bin, _ds_min_hyp, _etot_sd, _G, nest_sd, *_bin.getMemberAsTree(k));
                }
            }
            if (_bin.semi>0) {
                //kepler orbit, step ds=dt*G*m1*m2/r estimation (1/32 orbit): 2*pi/32*sqrt(semi/(m1+m2))*m1*m2 
                Float dsi = 0.19634954084*sqrt(_G*_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
                // scale by /Ebin_sd
                Float ebin_sd = _G*(_bin.m1*_bin.m2)/(2*_bin.semi*nest_sd);
                ASSERT(dsi>0&&ebin_sd>0);
                _ds_over_ebin_min_bin = std::min(dsi/ebin_sd, _ds_over_ebin_min_bin);
                _etot_sd += ebin_sd;
            }
            else {
                //hyperbolic orbit, step ds=dt*m1*m2/r estimation (1/256 orbit): pi/128*sqrt(semi/(m1+m2))*m1*m2
                Float dsi = 0.0245436926*sqrt(-_G*_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
                _ds_min_hyp = std::min(dsi, _ds_min_hyp);
            }
        }

    public:
        Float ds;  ///> initial step size for integration
        FixStepOption fix_step_option; ///> fix step option for integration
        COMM::List<BinaryTree<Tparticle>> binarytree; ///> a list of binary tree that contain the hierarchical orbital parameters of the particle group.

        //! initializer, set ds to zero, fix_step_option to none
        Information(): ds(Float(0.0)),  fix_step_option(AR::FixStepOption::none), binarytree() {}

        //! check whether parameters values are correct initialized
        /*! \return true: all correct
         */
        bool checkParams() {
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

        //! calculate ds from the inner most binary with minimum period, determine the fix step option
        /*! Estimate ds first from the inner most binary orbit (eccentric anomaly), set fix_step_option to later
          @param[in] _int_order: accuracy order of the symplectic integrator.
          @param[in] _G: gravitational constant
         */
        void calcDsAndStepOption(const int _int_order, const Float& _G) {
            auto& bin_root = getBinaryTreeRoot();
            Float ds_over_ebin_min=NUMERIC_FLOAT_MAX, ds_min_hyp=NUMERIC_FLOAT_MAX;
            Float etot_sd = 0.0;
            calcDsMinKeplerIter(ds_over_ebin_min, ds_min_hyp, etot_sd, _G, 1.0, bin_root);
            Float ds_min_bin = ds_over_ebin_min*etot_sd/(bin_root.getMemberN()-1);
            if (ds_min_bin!=0) ds = ds_min_bin;
            else {
                ASSERT(bin_root.semi<0);
                ds = ds_min_hyp;
            }
            ASSERT(ds>0);

            // Avoid too small step
            //if (_sd_org<1.0) ds *= std::max(1.0/8.0*pow(_sd_org, 1.0/Float(_int_order)),0.125);
            //auto& bin_root = getBinaryTreeRoot();
            //const int n_particle = bin_root.getMemberN();

            // determine the fix step option
            fix_step_option = FixStepOption::later;
            //fix_step_option = FixStepOption::none;
            //// for two-body case, determine the step at begining then fix
            //if (n_particle==2) fix_step_option = AR::FixStepOption::later;
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
        /*! @param[in] _particles: particle group 
         */
        void generateBinaryTree(COMM::ParticleGroup<Tparticle, Tpcm>& _particles, const Float _G) {
            const int n_particle = _particles.getSize();
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
                for (int i=2; i<n_particle; i++) {
                    binarytree[i-1].setMembers((Tparticle*)&(binarytree[i-2]), &( _particles[i]), -1, i);
                }
            }
            else if (n_particle_real==1) {
                int i1 = particle_index_local[0];
                int i2 = particle_index_unused[0];
                binarytree[0].setMembers(&(_particles[i1]), &(_particles[i2]), i1 ,i2);
                for (int i=1; i<n_particle_unused; i++) {
                    int k = particle_index_unused[i];
                    binarytree[i].setMembers((Tparticle*)&(binarytree[i-1]), &(_particles[k]), -1, k);
                }
            }
            else {
                for (int i=0; i<n_particle_unused; i++) {
                    int ilast = n_particle_real-1+i;
                    ASSERT(ilast<n_particle-1);
                    int k = particle_index_unused[i];
                    binarytree[ilast].setMembers((Tparticle*)&(binarytree[ilast-1]), &(_particles[k]), -1, k);
                }
            }
        }

        //! clear function
        void clear() {
            ds=0.0;
            fix_step_option = FixStepOption::none;
            binarytree.clear();
        }

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"ds";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<ds;
        }

        //! write class data to file with binary format
        /*! @param[in] _fp: FILE type file for output
         */
        void writeBinary(FILE *_fout) const {
            fwrite(&ds, sizeof(int),1,_fout);
            fwrite(&fix_step_option, sizeof(FixStepOption),1,_fout);
        }

        //! read class data to file with binary format
        /*! @param[in] _fp: FILE type file for reading
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(&ds, sizeof(int),1,_fin);
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
