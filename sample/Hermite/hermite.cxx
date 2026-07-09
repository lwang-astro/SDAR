#include <iostream>
#include <fstream>
//#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <cassert>

#define ASSERT(expr) assert(expr)
#define DATADUMP(x) abort()

#include "Common/io.h"
#include "AR/symplectic_integrator.h"
#include "Hermite/hermite_integrator.h"
#include "particle.h"
#include "hermite_perturber.h"
#include "ar_interaction.h"
#include "hermite_interaction.h"
#include "hermite_information.h"


using namespace H4;

typedef HermiteIntegrator<Particle, Particle, HermitePerturber, Neighbor<Particle>, HermiteInteraction, ARInteraction, HermiteInformation> H4Int;

//! IO parameters for Hermite integration
class IOParamsH4{
public:
    COMM::IOParamsContainer input_par_store;

    COMM::IOParams<int>     print_width;
    COMM::IOParams<int>     print_precision;
    COMM::IOParams<int>     nstep_max;
    COMM::IOParams<int>     sym_order;
    COMM::IOParams<int>     dt_min_power_index;
    COMM::IOParams<int>     dt_max_power_index;
    COMM::IOParams<int>     dt_out_power_index;
    COMM::IOParams<int>     n_neighbor_max;
    COMM::IOParams<double>  ds_scale;
    COMM::IOParams<int>     interrupt_detection_option;
    COMM::IOParams<double>  energy_error;
    COMM::IOParams<double>  time_error;
    COMM::IOParams<double>  time_zero;
    COMM::IOParams<double>  time_end;
    COMM::IOParams<double>  r_group;
    COMM::IOParams<double>  r_neighbor_over_group;
    COMM::IOParams<double>  eta_4th;
    COMM::IOParams<double>  eta_2nd;
    COMM::IOParams<double>  eps_sq;
    COMM::IOParams<double>  grav_const;
    COMM::IOParams<double>  slowdown_ref;
#ifdef SLOWDOWN_MASSRATIO
    COMM::IOParams<double>  slowdown_mass_ref;
#endif
    COMM::IOParams<double>  slowdown_timescale_max;
#ifdef USE_MPFRC
    COMM::IOParams<int>     mpfr_digits;
#endif
    COMM::IOParams<std::string> filename_par;

    IOParamsH4()
        : input_par_store()
        , print_width         (input_par_store, WRITE_WIDTH,        "print-width",          "print width of value")
        , print_precision     (input_par_store, WRITE_PRECISION,    "print-precision",      "print digital precision")
        , nstep_max           (input_par_store, 1000000,            "n-step-max",           "number of maximum step for AR integration")
        , sym_order           (input_par_store, -6,                 "k",                    "Symplectic integrator order, should be even number")
        , dt_min_power_index  (input_par_store, 40,                 "dt-min-power",         "power index to calculate mimimum hermite time step: dt_max*0.5^n")
        , dt_max_power_index  (input_par_store, 2,                  "dt-max-power",         "power index of 0.5 for maximum hermite time step")
        , dt_out_power_index  (input_par_store, 2,                  "o",                    "power index of 0.5 for output time interval")
        , n_neighbor_max      (input_par_store, -1,                 "n-neighbor-max",       "maximum number of neighbors for group","same as N")
        , ds_scale            (input_par_store, 1.0,                "ds-scale",             "step size scaling factor for Ar integration")
        , interrupt_detection_option(input_par_store, 0,            "i",                    "modify orbits and check interruption; 0: turn off; 1: modify the binary orbits based on detection criterion; 2. only record the binary information when interruption criterion is triggered")
        , energy_error        (input_par_store, 1e-10,              "e",                    "relative energy error limit for AR")
        , time_error          (input_par_store, 0.0,                "time-error",           "time synchronization absolute error limit for AR","default is 0.25*dt-min")
        , time_zero           (input_par_store, 0.0,                "time-start",           "initial physical time")
        , time_end            (input_par_store, 1.0,                "t",                    "ending physical time ")
        , r_group             (input_par_store, 1e-3,               "r-group",                    "distance criterion (group radius) reference for switching AR and Hermite;  the final radius is scaled by max(1,(mass/<mass>)^(1/3))")
        , r_neighbor_over_group(input_par_store, 2.0,                "r-neighbor-over-group", "coefficient to compute neighbor radius from group radius")
        , eta_4th             (input_par_store, 0.1,                "eta-4th",              "time step coefficient for 4th order")
        , eta_2nd             (input_par_store, 0.001,              "eta-2nd",              "time step coefficient for 2nd order")
        , eps_sq              (input_par_store, 0.0,                "eps",                  "softerning parameter")
        , grav_const          (input_par_store, 1.0,                "G",                    "gravitational constant")
        , slowdown_ref        (input_par_store, 1e-6,               "slowdown-ref",         "slowdown perturbation ratio reference")
#ifdef SLOWDOWN_MASSRATIO
        , slowdown_mass_ref   (input_par_store, 0.0,                "slowdown-mass-ref",    "slowdowm mass reference","averaged mass")
#endif
        , slowdown_timescale_max(input_par_store, 0.0,              "slowdown-timescale-max", "maximum timescale for maximum slowdown factor","time-end")
#ifdef USE_MPFRC
        , mpfr_digits         (input_par_store, 30,                 "mpfr-dights",          "dights for MPFR precison")
#endif
        , filename_par        (input_par_store, "",                 "p",                    "filename to load manager parameters","input name")
    {}

    int read(int argc, char* argv[], const char* bin_name) {
        static int h4_flag = -1;
        static struct option long_options[] = {
            {print_width.key,              required_argument, &h4_flag, 1},
            {print_precision.key,          required_argument, &h4_flag, 2},
            {nstep_max.key,                required_argument, &h4_flag, 3},
            {dt_min_power_index.key,       required_argument, &h4_flag, 4},
            {dt_max_power_index.key,       required_argument, &h4_flag, 5},
            {n_neighbor_max.key,           required_argument, &h4_flag, 6},
            {ds_scale.key,                 required_argument, &h4_flag, 7},
            {time_error.key,               required_argument, &h4_flag, 8},
            {time_zero.key,                required_argument, &h4_flag, 9},
            {eta_4th.key,                  required_argument, &h4_flag, 10},
            {eta_2nd.key,                  required_argument, &h4_flag, 11},
            {eps_sq.key,                   required_argument, &h4_flag, 12},
            {slowdown_ref.key,             required_argument, &h4_flag, 13},
            {slowdown_timescale_max.key,   required_argument, &h4_flag, 14},
#ifdef SLOWDOWN_MASSRATIO
            {slowdown_mass_ref.key,        required_argument, &h4_flag, 15},
#endif
#ifdef USE_MPFRC
            {mpfr_digits.key,              required_argument, &h4_flag, 16},
#endif
            {"r-neighbor-over-group",      required_argument, &h4_flag, 17},
            {"r-group",                     required_argument, &h4_flag, 18},
            {"help",                       no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };

        int opt_used = 0;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "t:k:G:e:o:i:p:h", long_options, &option_index)) != -1)
            switch (copt) {
            case 0:
                switch (h4_flag) {
                case 1:
                    print_width.value = atoi(optarg);
                    opt_used += 2;
                    break;
                case 2:
                    print_precision.value = atoi(optarg);
                    opt_used += 2;
                    break;
                case 3:
                    nstep_max.value = atoi(optarg);
                    opt_used += 2;
                    break;
                case 4:
                    dt_min_power_index.value = atoi(optarg);
                    opt_used += 2;
                    break;
                case 5:
                    dt_max_power_index.value = atoi(optarg);
                    opt_used += 2;
                    break;
                case 6:
                    n_neighbor_max.value = atoi(optarg);
                    opt_used += 2;
                    break;
                case 7:
                    ds_scale.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 8:
                    time_error.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 9:
                    time_zero.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 10:
                    eta_4th.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 11:
                    eta_2nd.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 12:
                    eps_sq.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 13:
                    slowdown_ref.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 14:
                    slowdown_timescale_max.value = atof(optarg);
                    opt_used += 2;
                    break;
#ifdef SLOWDOWN_MASSRATIO
                case 15:
                    slowdown_mass_ref.value = atof(optarg);
                    opt_used += 2;
                    break;
#endif
#ifdef USE_MPFRC
                case 16:
                    mpfr_digits.value = atoi(optarg);
                    opt_used += 2;
                    break;
#endif
                case 17:
                    r_neighbor_over_group.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 18:
                    r_group.value = atof(optarg);
                    opt_used += 2;
                    break;
                }
                break;
            case 't':
                time_end.value = atof(optarg);
                opt_used++;
                break;
            case 'k':
                sym_order.value = atoi(optarg);
                opt_used++;
                break;
            case 'G':
                grav_const.value = atof(optarg);
                opt_used++;
                break;
            case 'e':
                energy_error.value = atof(optarg);
                opt_used++;
                break;
            case 'o':
                dt_out_power_index.value = atoi(optarg);
                opt_used++;
                break;
            case 'i':
                interrupt_detection_option.value = atoi(optarg);
                opt_used++;
                break;
            case 'p':
                filename_par.value = optarg;
                {
                    FILE* fpar_in;
                    if( (fpar_in = fopen(filename_par.value.c_str(),"r")) == NULL) {
                        fprintf(stderr,"Error: Cannot open file %s.\n", filename_par.value.c_str());
                        abort();
                    }
                    input_par_store.readAscii(fpar_in);
                    fclose(fpar_in);
                }
                opt_used++;
                break;
            case 'h':
                std::cout<<bin_name<<" [option] data_filename\n"
                         <<"Input data file format: \n"
                         <<"  First   line:  number of particles(N)\n"
                         <<"  2-(N+1) line:  mass, x, y, z, vx, vy, vz, radius\n"
                         <<"  last    line:  N_group, group_offset_index_lst[N_group], group_member_particle_index[N_member_total]\n";
                input_par_store.printHelp(std::cout);
                std::cout<<"Size of integrator: (bytes)"<<sizeof(H4Int)<<std::endl;
                return -1;
            default:
                std::cerr<<"Unknown argument. check '-h' for help.\n";
                abort();
            }
        return opt_used;
    }
};

int main(int argc, char **argv){

    //unsigned int oldcw;
    //fpu_fix_start(&oldcw);

    IOParamsH4 iop;
    
    // Check whether all options are defined
    std::vector<COMM::IOParamsContainer*> all_pars;
    all_pars.push_back(&iop.input_par_store);
    std::vector<std::string> known_options;
    known_options.push_back("help");
    known_options.push_back("h");
    FindUndefinedOptions(all_pars, argc, argv, &known_options);

    if (argc==1) {
        std::cerr<<"Please provide particle data filename\n";
        abort();
    }

    int opt_used = iop.read(argc, argv, "hermite");
    if (opt_used < 0) return 0;

    // data file name
    char* filename = argv[argc-1];

#ifdef USE_MPFRC
    setMPFRPrec(iop.mpfr_digits.value);
#endif

    // manager
    HermiteManager<HermiteInteraction> manager;
    AR::TimeTransformedSymplecticManager<ARInteraction> ar_manager;

    manager.step.eta_4th = iop.eta_4th.value;
    manager.step.eta_2nd = iop.eta_2nd.value;
    Float dt_max = pow(Float(0.5), Float(iop.dt_max_power_index.value));
    manager.step.setDtRange(dt_max, iop.dt_min_power_index.value - iop.dt_max_power_index.value);
    manager.interaction.eps_sq = iop.eps_sq.value;
    manager.interaction.gravitational_constant = iop.grav_const.value;
    ar_manager.interaction.eps_sq = iop.eps_sq.value;
    ar_manager.interaction.gravitational_constant = iop.grav_const.value;
    ar_manager.time_step_min = manager.step.getDtMin();
    ar_manager.ds_scale = iop.ds_scale.value;
    if (iop.time_error.value == 0.0) ar_manager.time_error_max = 0.25*ar_manager.time_step_min;
    else ar_manager.time_error_max = iop.time_error.value;

    ASSERT(ar_manager.time_error_max>1e-14);
    // time error cannot be smaller than round-off error
    ar_manager.energy_error_relative_max = iop.energy_error.value; 
    ar_manager.slowdown_pert_ratio_ref = iop.slowdown_ref.value;
    if (iop.slowdown_timescale_max.value>0.0) ar_manager.slowdown_timescale_max = iop.slowdown_timescale_max.value;
    else ar_manager.slowdown_timescale_max = iop.time_end.value;
    ar_manager.step_count_max = iop.nstep_max.value;
    // set symplectic order
    ar_manager.step.initialSymplecticCofficients(iop.sym_order.value);
    ar_manager.interaction.interrupt_detection_option = iop.interrupt_detection_option.value;

    // store input parameters
    std::string fpar_out = std::string(filename) + ".par";
    std::FILE* fout = std::fopen(fpar_out.c_str(),"w");
    if (fout==NULL) {
        std::cerr<<"Error: data file "<<fpar_out<<" cannot be open!\n";
        abort();
    }
    iop.input_par_store.writeAscii(fout);
    fclose(fout);

    // interrupt file output
    std::ofstream finterrupt;
    if (iop.interrupt_detection_option.value>0) {
        std::string finterrupt_name = std::string(filename) + ".interrupt";
        finterrupt.open(finterrupt_name.c_str(),std::ofstream::out);
        AR::InterruptBinary<Particle>::printColumnTitleAscii(finterrupt,20,true);
        finterrupt<<std::endl;
    }

    // integrator
    H4Int h4_int;
    h4_int.manager = &manager;
    h4_int.ar_manager = &ar_manager;

    std::fstream fin;
    fin.open(filename,std::fstream::in);
    if(!fin.is_open()) {
        std::cerr<<"Error: data file "<<filename<<" cannot be open!\n";
        abort();
    }
    h4_int.particles.setMode(COMM::ListMode::local);
    h4_int.particles.readMemberAscii(fin);
    for (int i=0; i<h4_int.particles.getSize(); i++) h4_int.particles[i].id = i+1;
    h4_int.particles.calcCenterOfMass();
    h4_int.particles.shiftToCenterOfMassFrame();
    h4_int.particles.calcCenterOfMass();

    if (iop.n_neighbor_max.value <=0) manager.n_neighbor_max = h4_int.particles.getSize();
    else manager.n_neighbor_max = iop.n_neighbor_max.value;
        
    Float m_ave = h4_int.particles.cm.mass/h4_int.particles.getSize();
    // initialize per-particle group and neighbor radii with mass-dependent weighting
    Float r_neighbor_sum = 0.0;
    for (int i=0; i<h4_int.particles.getSize(); i++) {
        h4_int.particles[i].setRGroupAndNeighbor(iop.r_group.value, iop.r_neighbor_over_group.value, m_ave);
        r_neighbor_sum += h4_int.particles[i].getRNeighbor();
    }
    Float r_neighbor_ave = r_neighbor_sum / h4_int.particles.getSize();
    manager.step.calcAcc0OffsetSq(m_ave, r_neighbor_ave, iop.grav_const.value);
    h4_int.step = manager.step;

#ifdef SLOWDOWN_MASSRATIO
    if (iop.slowdown_mass_ref.value<=0.0) ar_manager.slowdown_mass_ref = m_ave;
    else ar_manager.slowdown_mass_ref = iop.slowdown_mass_ref.value;
#endif
    // print parameters
    manager.print(std::cerr);
    ar_manager.print(std::cerr);


    std::cerr<<"CM: after shift ";
    h4_int.particles.cm.printColumnAscii(std::cerr, 22);
    std::cerr<<std::endl;

    h4_int.groups.setMode(COMM::ListMode::local);
    h4_int.groups.reserveMem(h4_int.particles.getSize());
    h4_int.reserveIntegratorMem();
    // initial system 
    h4_int.initialSystemSingle(iop.time_zero.value);
    h4_int.readGroupConfigureAscii(fin);

    // initialization 
    h4_int.initialIntegration(); // get neighbors and min particles
    const int n_group_init = h4_int.getNGroup();
    // AR inner slowdown number
    int n_group_sub_init[n_group_init], n_group_sub_tot_init=0;
    for (int i=0; i<n_group_init; i++) {
#ifdef AR_SLOWDOWN_ARRAY
        n_group_sub_init[i] = h4_int.groups[i].binary_slowdown.getSize();
#elif AR_SLOWDOWN_TREE
        n_group_sub_init[i] = h4_int.groups[i].info.binarytree.getSize();
#endif
        n_group_sub_tot_init += n_group_sub_init[i];
    }
    h4_int.adjustGroups(true);
    h4_int.initialIntegration();
    h4_int.sortDtAndSelectActParticle();

    // precision
    std::cout<<std::setprecision(iop.print_precision.value);

    // get initial energy
    h4_int.calcEnergySlowDown(true);
    // cm
    h4_int.particles.calcCenterOfMass();
    std::cerr<<"CM:";
    h4_int.particles.cm.printColumnAscii(std::cerr, 22);
    std::cerr<<std::endl;

    //print column title
    h4_int.printColumnTitleAscii(std::cout, iop.print_width.value, n_group_sub_init, n_group_init, n_group_sub_tot_init);
    std::cout<<std::endl;

    //print initial data
    h4_int.printColumnAscii(std::cout, iop.print_width.value, n_group_sub_init, n_group_init, n_group_sub_tot_init);
    std::cout<<std::endl;
    
    // dt_out
    Float dt_out = pow(Float(0.5),Float(iop.dt_out_power_index.value));
    Float time_out = iop.time_zero.value + dt_out;

    // integration loop
    while (h4_int.getTime()<iop.time_end.value) {
        h4_int.integrateGroupsOneStep();
        int n_interrupt = h4_int.getNInterrupt();
        for (int i=0; i<n_interrupt; i++) {
            auto& interrupt_info = h4_int.getInterruptInfo(i);
            std::cerr<<"Interrupt "<<i<<" : ";
            switch (interrupt_info.status) {
            case AR::InterruptStatus::change:
                std::cerr<<" Change";
                break;
            case AR::InterruptStatus::merge:
                std::cerr<<" Merge";
                break;
            case AR::InterruptStatus::destroy:
                std::cerr<<" Destroy";
                break;
            case AR::InterruptStatus::none:
                break;
            }
            std::cerr<<std::endl;
            interrupt_info.printColumnTitleAscii(std::cerr);
            std::cerr<<std::endl;
            interrupt_info.printColumnAscii(std::cerr);
            std::cerr<<std::endl;
            if (iop.interrupt_detection_option.value>0) {
                interrupt_info.printColumnAscii(finterrupt, 20, true);
                finterrupt<<std::endl;
            }
        }
        h4_int.integrateSingleOneStepAct();
        h4_int.adjustGroups(false);
        h4_int.initialIntegration();
        h4_int.modifySingleParticles();
        h4_int.sortDtAndSelectActParticle();

        if (h4_int.getTime()>=time_out) {
            h4_int.calcEnergySlowDown(false);
            
            h4_int.particles.calcCenterOfMass();
            std::cerr<<"CM:";
            h4_int.particles.cm.printColumnAscii(std::cerr, 22);
            std::cerr<<std::endl;

            // Notice in energy calculation, writeBackGroupMembers() is already done;
            h4_int.printColumnAscii(std::cout, iop.print_width.value, n_group_sub_init, n_group_init, n_group_sub_tot_init);
            std::cout<<std::endl;
            h4_int.printStepHist();

            time_out += dt_out;
        }
    }

    //fpu_fix_end(&oldcw);

    return 0;
}

