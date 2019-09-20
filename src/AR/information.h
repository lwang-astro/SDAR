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
    class Information{
    private:
        //! calculate minimum kepler ds for ARC
        static Float calcDsKeplerIter(Float& _ds, COMM::BinaryTree<Tparticle>& _bin) {
            Float ds;
            //kepler orbit, step ds=dt*m1*m2/r estimation (1/32 orbit): 2*pi/32*sqrt(semi/(m1+m2))*m1*m2 
            if (_bin.semi>0) ds = 0.19634954084*sqrt(_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
            //hyperbolic orbit, step ds=dt*m1*m2/r estimation (1/256 orbit): pi/128*sqrt(semi/(m1+m2))*m1*m2
            else ds = 0.0245436926*sqrt(-_bin.semi/(_bin.m1+_bin.m2))*(_bin.m1*_bin.m2);
            return std::min(ds, _ds);
        }

    public:
        Float ds;  ///> estimated step size for AR integration
        FixStepOption fix_step_option; ///> fixt step option for integration
        COMM::List<COMM::BinaryTree<Tparticle>> binarytree;

        Information(): ds(Float(0.0)),  fix_step_option(AR::FixStepOption::none), binarytree() {}

        //! check whether parameters values are correct initialized
        /*! \return true: all correct
         */
        bool checkParams() {
            ASSERT(binarytree.getSize()>0);
            return true;
        }

        //! reserve memory
        void reserveMem(const int _nmax) {
            binarytree.setMode(COMM::ListMode::local);
            binarytree.reserveMem(_nmax);
        }

        //! get binary tree root 
        COMM::BinaryTree<Tparticle>& getBinaryTreeRoot() const {
            int n = binarytree.getSize();
            ASSERT(n>0);
            return binarytree[n-1];
        }

        //! calculate ds from binary tree information and adjust by slowdown, determine the fix step option
        /*! Estimate ds first from binary orbit (eccentricity anomaly), then adjust based on original slowdown factor and order of integrator
          @param[in] _sd_org: original slowdown factor
          @param[in] _int_order: integrator order 
         */
        void calcDsAndStepOption(const Float _sd_org, const int _int_order) {
            ds = NUMERIC_FLOAT_MAX;
            ds = getBinaryTreeRoot().processRootIter(ds, calcDsKeplerIter);
            // Avoid too small step
            //if (_sd_org<1.0) ds *= std::max(1.0/8.0*pow(_sd_org, 1.0/Float(_int_order)),0.125);
            //auto& bin_root = getBinaryTreeRoot();
            //const int n_particle = bin_root.getMemberN();

            // determine the fix step option
            fix_step_option = FixStepOption::later;
            //// for two-body case, determine the step at begining then fix
            //if (n_particle==2) fix_step_option = AR::FixStepOption::later;
            //// for multiple case, check whether outer peri-center is close to inner apo-center, if not, use fix step
            //if (n_particle>2) {
            //    Float apo_in_max = 0;
            //    for (int j=0; j<2; j++) {
            //        if (bin_root.getMember(j)->id<0) {
            //            auto* bin_sub = (COMM::BinaryTree<Tparticle>*) bin_root.getMember(j);
            //            apo_in_max = std::max(apo_in_max,bin_sub->semi*(1+bin_sub->ecc));
            //        }
            //    }
            //    Float peri_out = bin_root.semi*(1-bin_root.ecc);
            //    if (peri_out>3*apo_in_max) fix_step_option = FixStepOption::later;
            //}
        }

        //! generate binary tree for a particle group
        /*! @param[in] _particles: particle group 
         */
        void generateBinaryTree(COMM::ParticleGroup<Tparticle, Tpcm>& _particles) {
            const int n_particle = _particles.getSize();
            ASSERT(n_particle>1);
            binarytree.resizeNoInitialize(n_particle-1);
            int particle_index_local[n_particle];
            for (int i=0; i<n_particle; i++) particle_index_local[i] = i;
            COMM::BinaryTree<Tparticle>::generateBinaryTree(binarytree.getDataAddress(), particle_index_local, n_particle, _particles.getDataAddress());
        }

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
