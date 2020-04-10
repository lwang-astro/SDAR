#pragma once
#include "Common/Float.h"
#include "Common/binary_tree.h"

namespace AR {
    //! Fix step options for integration with adjusted step (not for time sychronizatio phase)
    /*! always: use the given step without change \n
        later: fix step after a few adjustment of initial steps due to energy error
        none: don't fix step
     */
    enum class FixStepOption {always, later, none};

    template <class Tparticle, class Tpcm>
    //! A class contains the Kepler orbital parameters and initial step size
    /*! The class has members to generate the Kepler orbital binary tree of the particle group and store the information in binarytree. \\
      It also provides the initial step size estimator and the option of fix step.
     */
    class Information{
    private:
        //! calculate average kepler ds for ARC
        template <class Tds>
        static Tds calcDsKeplerIter(Tds& _ds, COMM::BinaryTree<Tparticle>& _bin) {
            Tds ds;
            ds.G   = _ds.G;
#ifdef AR_TTL_GT_MULTI
            // ecca [v/r]
            Float tov2 = 0.03855314219*_bin.semi*_bin.semi*_bin.semi/(_ds.G*_bin.m1+_bin.m2);
            // semi cum
            ds.tov = std::min(_ds.tov,tov2);
            ds.r   = 2.0*_ds.r*_bin.semi;
#else
            //kepler orbit, step ds=dt*m1*m2/r estimation (1/32 orbit): 2*pi/32*sqrt(semi/(m1+m2))*m1*m2 
            if (_bin.semi>0) ds.min = 0.19634954084*sqrt(_ds.G*_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
            //hyperbolic orbit, step ds=dt*m1*m2/r estimation (1/256 orbit): pi/128*sqrt(semi/(m1+m2))*m1*m2
            else ds.min = 0.0245436926*sqrt(-_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
            ds.min = std::min(ds.min, _ds.min);
#endif
            return ds;
        }

    public:
        Float ds;  ///> initial step size for integration
        FixStepOption fix_step_option; ///> fix step option for integration
        COMM::List<COMM::BinaryTree<Tparticle>> binarytree; ///> a list of binary tree that contain the hierarchical orbital parameters of the particle group.

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
        COMM::BinaryTree<Tparticle>& getBinaryTreeRoot() const {
            int n = binarytree.getSize();
            ASSERT(n>0);
            return binarytree[n-1];
        }

        //! calculate ds from the inner most binary with minimum period, determine the fix step option
        /*! Estimate ds first from the inner most binary orbit (eccentric anomaly), set fix_step_option to later
          @param[in] _sd_org: original slowdown factor of the group
          @param[in] _int_order: accuracy order of the symplectic integrator.
          @param[in] _G: gravitational constant
         */
        void calcDsAndStepOption(const Float _sd_org, const int _int_order, const Float _G) {
            auto& bin_root = getBinaryTreeRoot();
#ifdef AR_TTL_GT_MULTI
            struct {Float G, tov, r; } ds_iter = {_G, NUMERIC_FLOAT_MAX, 1.0};
            ds_iter = bin_root.processRootIter(ds_iter, calcDsKeplerIter);
            ds = sqrt(ds_dat.tov)/ds_dat.r;
#else
            struct {Float G, min;} ds_iter = {_G, NUMERIC_FLOAT_MAX};
            ds_iter = bin_root.processRootIter(ds_iter, calcDsKeplerIter);
            ds = ds_iter.min;
#endif
            // Avoid too small step
            //if (_sd_org<1.0) ds *= std::max(1.0/8.0*pow(_sd_org, 1.0/Float(_int_order)),0.125);
            //auto& bin_root = getBinaryTreeRoot();
            //const int n_particle = bin_root.getMemberN();

            // determine the fix step option
            //fix_step_option = FixStepOption::later;
            fix_step_option = FixStepOption::none;
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
                COMM::BinaryTree<Tparticle>::generateBinaryTree(binarytree.getDataAddress(), particle_index_local, n_particle_real, _particles.getDataAddress(), _G);

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
