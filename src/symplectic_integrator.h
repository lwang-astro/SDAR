#pragma once

#include "symplectic_step.h"
#include "force.h"
#include "slow_down.h"
#include "particle_group.h"
#include "list.h"
#include "binary_tree.h"

//! Algorithmic regularization chain (ARC) namespace
/*!
  All major ARC classes and related acceleration functions (typedef) are defined
*/
namespace AR {
    //! Symplectic integrator manager
    /*! Tmethod is the class contain the interaction function, see sample of interaction.h:\n
     */
    template <class Tmethod>
    class SymplecticManager {
    public:
        Tmethod interaction; ///> class contain interaction function
        SymplecticStep step;  ///> class to manager kick drift step
    };
    
    //! Symplectic integrator class for a group of particles
    /*! The basic steps to integrate \n
      1. Add particles (particles.addParticle/particles.linkParticleList)  \n
      2. Initial system (initial) \n
      3. Integration (integrateOneStep/integrateToTime) \n
      Requirement for Tparticle class, public memebers: pos[3], vel[3], mass\n
      Template dependence: Tparticle: particle type;  Tpert: perturber class type, Tmethod: interaction class;
    */
    template <class Tparticle, class Tpert, class Tmethod>
    class SymplecticIntegrator {
    private:
        // intergrated variables
        Float time_;   ///< integrated time (not real physical time if slowdown is on)
        Float etot_; ///< integrated system energy

        // calculated varaiables
        Float ekin_;   ///< kinetic energy
        Float epot_;   ///< potential

#ifdef AR_TTL
        // transformation factors
        Float gt_inv_;     ///< integrated time transformation factor for drift 
#endif

        // force array
        List<Force> force_; ///< acceleration array 

    public:
        SymplecticManager<Tmethod>* manager; ///< integration manager
        Tpert   perturber; ///< perturber class 
        SlowDown slowdown; ///< slowdown controller
        ParticleGroup<Tparticle> particles; ///< particle group manager
        BinaryTree binarytree;   ///< particle chain
        
        //! Constructor
#ifdef AR_TTL
        SymplecticIntegrator(): time_(0), etot_(0), ekin_(0), epot_(0), gt_inv_(0), force_(), manager(NULL), perturber(), slowdown(), particles(), binarytree() {}
#else
        SymplecticIntegrator(): time_(0), etot_(0), ekin_(0), epot_(0), force_(), manager(NULL), perturber(), slowdown(), particles(), binarytree() {}
#endif

        //! reserve memory for force
        /*! The size of force depends on the particle data size.Thus particles should be added first before call this function
          @param[in] _nmax: maximum number of particles for memory allocation, if not given (0), use particles reserved size
        */
        void reserveMem() {
            // force array always allocated local memory
            force_.setMode(ListMode::local);
            int nmax = particles.getSizeMax();
            assert(nmax>0);
            force_.reserveMem(nmax);
        }

        //! Clear function
        /*! Free dynamical memory space allocated
         */
        void clear() {
            force_.clear();
            perturber.clear();
            slowdown.clear();
            particles.clear();
            binarytree.clear();
        }

        //! destructor
        ~SymplecticIntegrator() {
            clear();
        }

        //! operator = 
        /*! Copy function will remove the local data and also copy the particle data or the link
         */
        SymplecticIntegrator& operator = (const SymplecticIntegrator& _sym) {
            clear();
            time_   = _sym.time_;
            etot_   = _sym.etot_;

            ekin_   = _sym.ekin_;
            epot_   = _sym.epot_;
#ifdef AR_TTL
            gt_inv_ = _sym.gt_inv_;
#endif
            force_  = _sym.force_;
            manager = _sym.manager;
            slowdown = _sym.slowdown;
            particles = _sym.particles;
            binarytree = _sym.binarytree;

            return *this;
        }

    private:
        //! Calculate kinetic energy
        inline void calcEKin(){
            ekin_ = Float(0.0);
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            for (int i=0; i<num; i++) {
                const Float *vi=pdat[i].vel;
                ekin_ += 0.5 * pdat[i].mass * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
            }
        }

        //! kick velocity
        /*! First time step will be calculated, the velocities are kicked
          @param[in] _dt: time size
        */
        inline void kickVel(const Float _dt) {
            const int num = particles.getSize();            
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            const Float kappa = slowdown.getSlowDownFactor();
            for (int i=0; i<num; i++) {
                // kick velocity
                Float* vel = pdat[i].vel;
                Float* acc = force[i].acc_in;
                Float* pert= force[i].acc_pert;
                // half dv 
                vel[0] += _dt * (acc[0] + kappa*pert[0]);
                vel[1] += _dt * (acc[1] + kappa*pert[1]);
                vel[2] += _dt * (acc[2] + kappa*pert[2]);
            }
        }

        //! drift time and position
        /*! First (real) time is drifted, then positions are drifted
          @param[in] _dt: time step
        */
        inline void driftTimeAndPos(const Float _dt) {
            // drift time 
            time_ += _dt;

            // update real time
            slowdown.driftRealTime(_dt);

            // drift position
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            for (int i=0; i<num; i++) {
                Float* pos = pdat[i].pos;
                Float* vel = pdat[i].vel;
                pos[0] += _dt * vel[0];
                pos[1] += _dt * vel[1];
                pos[2] += _dt * vel[2];
            }
        }


#ifdef AR_TTL
        //! kick energy and time transformation function for drift
        /*!
          @param[in] _dt: time step
        */

        inline void kickEtotAndGTDrift(const Float _dt) {
            Float de = Float(0.0);
            Float dgt = Float(0.0);
            const Float kappa = slowdown.getSlowDownFactor();
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].vel;
                Float* pert= force[i].acc_pert;
                Float* gtgrad=force[i].gtgrad;
                de += mass * (vel[0] * kappa*pert[0] +
                              vel[1] * kappa*pert[1] +
                              vel[2] * kappa*pert[2]);
                dgt +=  (vel[0] * gtgrad[0] +
                         vel[1] * gtgrad[1] +
                         vel[2] * gtgrad[2]);
            }
            etot_ += _dt * de;
            gt_inv_ += _dt * dgt;
        }

#else
        //! kick energy
        /*!
          @param[in] _dt: time step
        */
        inline void kickEtot(const Float _dt) {
            Float de = 0.0;
            const Float kappa = slowdown.getSlowDownFactor();
            const int num = particles.getSize();
            Tparticle* pdat = particles.getDataAddress();
            Force* force = force_.getDataAddress();
            for (int i=0;i<num;i++) {
                Float  mass= pdat[i].mass;
                Float* vel = pdat[i].vel;
                Float* pert= force[i].acc_pert;
                de += mass * (vel[0] * kappa*pert[0] +
                              vel[1] * kappa*pert[1] +
                              vel[2] * kappa*pert[2]);
            }
            etot_ += _dt * de;
        }
#endif
            

        //! Backup integration data 
        /*! Backup #time_, #etot_, #ekin_, $epot_, #gt_drift_, $gt_kick_, #particles, $slowdown to one Float data array
          \return backup array size
        */
        int backupIntData(Float* _bk) {
            _bk[0] = time_;
            _bk[1] = etot_;
            _bk[2] = ekin_;
            _bk[3] = epot_;
#ifdef AR_TTL
            _bk[4] = gt_inv_;
            int bk_size = 5; 
#else
            int bk_size = 4;
#endif
            bk_size += particles.backupParticlePosVel(&_bk[bk_size]); 
            bk_size +=  slowdown.backupSlowDownFactorAndTimeReal(&_bk[bk_size]); // time_real, slowdownfactor
            return bk_size;
        }

        //! Restore integration data
        /*! restore #time_, #etot_, #ekin_, $epot_, #gt_drift_, $gt_kick_, #particles, $slowdown from one Float data array
          \return backup array size
        */
        int restoreIntData(Float* _bk) {
            time_    = _bk[0];
            etot_    = _bk[1];
            ekin_    = _bk[2];
            epot_    = _bk[3];
#ifdef AR_TTL
            gt_inv_= _bk[4];
            int bk_size = 5; 
#else
            int bk_size = 4; 
#endif
            bk_size += particles.restoreParticlePosVel(&_bk[bk_size]);
            bk_size +=  slowdown.restoreSlowDownFactorAndTimeReal(&_bk[bk_size]);
            return bk_size;
        }

    public:

        //! initialization for integration
        /*! initialize the system. Acceleration, energy and time transformation factors are updated. If the center-of-mass is not yet calculated, the system will be shifted to center-of-mass frame.
          @param[in] _time_real: real physical time to initialize
        */
        void initial(const Float _time_real) {
            // Initial intgrt value t (avoid confusion of real time when slowdown is used)
            time_ = Float(0.0);

            // slowdown initial time
            slowdown.setRealTime(_time_real);

            // reset particle modification flag
            particles.setModifiedFalse();

            // check the center-of-mass initialization
            if(particles.isOriginFrame()) {
                particles.calcCenterOfMass();
                particles.shiftToCenterOfMassFrame();
            }

            Tparticle* particle_data = particles.getDataAddress();
            Force* force_data = force_.getDataAddress();
            const int n_particle = particles.getSize();

            if (n_particle==2) {

#ifdef AR_TTL
                // calculate acceleration, potential, time transformation function gradient and time transformation factor 
                Float gt_kick = manager->interaction.calcAccPotAndGTKickTwo(force_data, epot_, particle_data, n_particle); 

                // initially gt_drift 
                gt_inv_ = 1.0/gt_kick;
#else 
                // calculate acceleration, potential and time transformation factor 
                manager->interaction.calcAccPotAndGTKickTwo(force_data, epot_, particle_data, n_particle); 
#endif

                // two-body case, also initialslowdown
                // calcualte perturbation force (cumulative to acc) and slowdown perturbation estimation
                Float sd_pert = manager->interaction.calcAccAndSlowDownPert(force_data, particle_data, n_particle, particles.cm, perturber);
                // slowdown factor
                slowdown.calcSlowDownFactorBinary(etot_, sd_pert);

            }
            else {
#ifdef AR_TTL
                // calculate acceleration, potential, time transformation function gradient and time transformation factor 
                Float gt_kick = manager->interaction.calcAccPotAndGTKick(force_data, epot_, particle_data, n_particle); 

                // initially gt_drift 
                gt_inv_ = 1.0/gt_kick;
#else 
                // calculate acceleration, potential and time transformation factor 
                manager->interaction.calcAccPotAndGTKick(force_data, epot_, particle_data, n_particle); 
#endif
            }

            // calculate kinetic energy
            calcEKin();

            // initial total energy
            etot_ = ekin_ + epot_;

        }

        //! integration for one step
        /*!
          @param[in] _ds: step size
        */
        void integrateOneStep(const Float _ds) {
            assert(!particles.isModified());
            assert(_ds>0);

            // symplectic step coefficent group n_particleber
            const int nloop = manager->step.getCDPairSize();

            // particle pointer and size
            const int n_particle = particles.getSize();
            Tparticle* particle_data = particles.getDataAddress();
            Force* force_data = force_.getDataAddress();
            
            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds_drift = manager->step.getCK(i)*_ds;

                // time transformation factor for drift
#ifdef AR_TTL
                Float gt_drift = 1.0/gt_inv_;
#else 
                Float gt_drift = manager->interaction.calcGTDrift(ekin_-etot_); // pt = -etot
#endif

                // drift
                Float dt_drift = ds_drift*gt_drift;

                // drift time and postion
                driftTimeAndPos(dt_drift);

                // step for kick
                Float ds_kick = manager->step.getDK(i)*_ds;

                // calculate acceleration, potential and time transfromation factor
                // Notice in TTL method, time transformation function gradient is also updated 
                Float gt_kick = manager->interaction.calcAccPotAndGTKick(force_data, epot_, particle_data, n_particle); 

                // calcualte perturbation force (cumulative to acc)
                manager->interaction.calcAccAndSlowDownPert(force_data, particle_data, n_particle, particles.cm, perturber);

                // time step for kick
                Float dt_kick = ds_kick* gt_kick;

                // kick half step for velocity
                kickVel(0.5*dt_kick);

#ifdef AR_TTL                
                // kick total energy and time transformation factor for drift
                kickEtotAndGTDrift(dt_kick);
#else
                // kick total energy 
                kickEtot(dt_kick);
#endif
                // kick half step for velocity
                kickVel(0.5*dt_kick);
                
                // calculate kinetic energy
                calcEKin();
            }
        }

        //! integration for two body one step
        /*! For two-body problem the calculation can be much symplified to improve performance. 
          Besides, the slow-down factor calculation is embedded in the Drift (for time) and Kick (for perturbation)
         */
        void integrateTwoOneStep(const Float _ds) {
            assert(!particles.isModified());
            assert(_ds>0);

            // symplectic step coefficent group number
            const int nloop = manager->step.getCDPairSize();
            
            const int n_particle = particles.getSize();
            assert(n_particle==2);

            Tparticle* particle_data = particles.getDataAddress();
            Float mass1 = particle_data[0].mass;
            Float* pos1 = particle_data[0].pos;
            Float* vel1 = particle_data[0].vel;

            Float mass2 = particle_data[1].mass;
            Float* pos2 = particle_data[1].pos;
            Float* vel2 = particle_data[1].vel;

            Force* force_data = force_.getDataAddress();
            Float* acc1 = force_data[0].acc_in;
            Float* pert1= force_data[0].acc_pert;

            Float* acc2 = force_data[1].acc_in;
            Float* pert2= force_data[1].acc_pert;
#ifdef AR_TTL
            Float* gtgrad1 = force_data[0].gtgrad;
            Float* gtgrad2 = force_data[1].gtgrad;
#endif
            for (int i=0; i<nloop; i++) {
                // step for drift
                Float ds = manager->step.getCK(i)*_ds;
                // time transformation factor for drift
#ifdef AR_TTL
                Float gt = 1.0/gt_inv_;
#else 
                Float gt = manager->interaction.calcGTDrift(ekin_-etot_); // pt = -etot
#endif
                // drift
                Float dt = ds*gt;
                
                // drift time 
                time_ += dt;

                // update real time
                slowdown.driftRealTime(dt);

                // drift position
                pos1[0] += dt * vel1[0];
                pos1[1] += dt * vel1[1];
                pos1[2] += dt * vel1[2];

                pos2[0] += dt * vel2[0];
                pos2[1] += dt * vel2[1];
                pos2[2] += dt * vel2[2];

                // step for kick
                ds = manager->step.getDK(i)*_ds;

                // calculate acceleration, potential, and time transformation factor for kick
                // Notice in TTL method, time transformation function gradient is also updated 
                gt = manager->interaction.calcAccPotAndGTKickTwo(force_data, epot_, particle_data, n_particle); 

                // time step for kick
                dt = 0.5*ds*gt;

                // calcualte perturbation force (cumulative to acc) and slowdown perturbation estimation
                Float sd_pert = manager->interaction.calcAccAndSlowDownPert(force_data, particle_data, n_particle, particles.cm, perturber);

                // slowdown factor
                const Float kappa = slowdown.calcSlowDownFactorBinary(etot_, sd_pert);

                // kick half step for velocity
                Float dvel1[3], dvel2[3];
                Float kpert1[3], kpert2[3];
                kpert1[0] = kappa*pert1[0];
                kpert1[1] = kappa*pert1[1];
                kpert1[2] = kappa*pert1[2];

                dvel1[0] = dt * (acc1[0] + kpert1[0]);
                dvel1[1] = dt * (acc1[1] + kpert1[1]);
                dvel1[2] = dt * (acc1[2] + kpert1[2]);

                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];

                kpert2[0] = kappa*pert2[0];
                kpert2[1] = kappa*pert2[1];
                kpert2[2] = kappa*pert2[2];

                dvel2[0] = dt * (acc2[0] + kappa*pert2[0]);
                dvel2[1] = dt * (acc2[1] + kappa*pert2[1]);
                dvel2[2] = dt * (acc2[2] + kappa*pert2[2]);

                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];

                // kick total energy and time transformation factor for drift
                etot_ += 2.0*dt * (mass1* (vel1[0] * kpert1[0] + 
                                           vel1[1] * kpert1[1] + 
                                           vel1[2] * kpert1[2]) +
                                   mass2* (vel2[0] * kpert2[0] + 
                                           vel2[1] * kpert2[1] + 
                                           vel2[2] * kpert2[2]));

#ifdef AR_TTL                
                gt_inv_ +=  2.0*dt* (vel1[0] * gtgrad1[0] +
                                     vel1[1] * gtgrad1[1] +
                                     vel1[2] * gtgrad1[2] +
                                     vel2[0] * gtgrad2[0] +
                                     vel2[1] * gtgrad2[1] +
                                     vel2[2] * gtgrad2[2]);
#endif

                // kick half step for velocity
                vel1[0] += dvel1[0];
                vel1[1] += dvel1[1];
                vel1[2] += dvel1[2];
                
                vel2[0] += dvel2[0];
                vel2[1] += dvel2[1];
                vel2[2] += dvel2[2];
                
                // calculate kinetic energy
                ekin_ = 0.5 * (mass1 * (vel1[0]*vel1[0]+vel1[1]*vel1[1]+vel1[2]*vel1[2]) +
                               mass2 * (vel2[0]*vel2[0]+vel2[1]*vel2[1]+vel2[2]*vel2[2]));
            }
        }
        
        //! Get current physical time
        /*! \return current physical time
         */
        Float getTime() const {
            return slowdown.getRealTime();
        }

        //! Get current kinetic energy
        /*! \return current kinetic energy
         */
        Float getEkin() const {
            return ekin_;
        }

        //! Get current potential energy
        /*! \return current potetnial energy (negative value for bounded systems)
         */
        Float getEpot() const {
            return epot_;
        }

        //! Get current total integrated energy 
        /*! \return total integrated energy 
         */
        Float getEtot() const {
            return etot_;
        }

#ifdef AR_TTL
        //! Get integrated inverse time transformation factor
        /*! In TTF case, it is calculated by integrating \f$ \frac{dg}{dt} = \sum_k \frac{\partial g}{\partial \vec{r_k}} \bullet \vec{v_k} \f$.
          Notice last step is the sub-step in one symplectic loop
          \return time transformation factor for drift
        */
        Float getGTIntInv() const {
            return gt_inv_;
        }
#endif

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
                 <<std::setw(_width)<<"Epot";
#ifdef AR_TTL
            _fout<<std::setw(_width)<<"Gt_inv";
#endif
            slowdown.printColumnTitle(_fout, _width);
            particles.printColumnTitle(_fout, _width);
        }
        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<time_
                 <<std::setw(_width)<<ekin_+epot_-etot_
                 <<std::setw(_width)<<etot_
                 <<std::setw(_width)<<ekin_
                 <<std::setw(_width)<<epot_;
#ifdef AR_TTL
            _fout<<std::setw(_width)<<gt_inv_;
#endif
            slowdown.printColumn(_fout, _width);
            particles.printColumn(_fout, _width);
        }

    };
}
