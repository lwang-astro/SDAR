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

#define ASSERT(x) assert(x)
#define DATADUMP(expr) 

#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "Common/io.h"
#include "AR/symplectic_integrator.h"
#include "AR/information.h"
#include "particle.h"
#include "perturber.h"
#include "interaction.h"

using namespace AR;

typedef TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, Interaction, Information<Particle,Particle>> ARInt;

//! IO parameters for AR integration
class IOParamsAR{
public:
    COMM::IOParamsContainer input_par_store;

    COMM::IOParams<int>     print_width;
    COMM::IOParams<int>     print_precision;
    COMM::IOParams<int>     nstep_max;
    COMM::IOParams<int>     sym_order;
    COMM::IOParams<double>  energy_error;
    COMM::IOParams<double>  time_error;
    COMM::IOParams<double>  time_zero;
    COMM::IOParams<double>  time_end;
    COMM::IOParams<double>  r_break;
    COMM::IOParams<int>     nstep;
    COMM::IOParams<double>  s;
    COMM::IOParams<double>  ds_scale;
    COMM::IOParams<double>  gravitational_constant;
    COMM::IOParams<double>  dt_min;
    COMM::IOParams<double>  dt_out;
    COMM::IOParams<double>  slowdown_ref;
#ifdef AR_SLOWDOWN_MASSRATIO
    COMM::IOParams<double>  slowdown_mass_ref;
#endif
    COMM::IOParams<double>  slowdown_timescale_max;
    COMM::IOParams<int>     interrupt_detection_option;
    COMM::IOParams<int>     fix_step_option;
#ifdef USE_MPFRC
    COMM::IOParams<int>     mpfr_digits;
#endif
#ifdef AR_G_FUNC
    COMM::IOParams<int>     g_func_option;
    COMM::IOParams<std::string> g_func_switch_option;
#endif
    COMM::IOParams<std::string> filename_par;
    COMM::IOParams<std::string> filename_out;
    COMM::IOParams<std::string> integration_mode;
    COMM::IOParams<int> load_flag;

    IOParamsAR()
        : input_par_store()
        , print_width         (input_par_store, WRITE_WIDTH,        "print-width",     "print width of value")
        , print_precision     (input_par_store, WRITE_PRECISION,    "print-precision", "print digital precision")
        , nstep_max           (input_par_store, 1000000,            "n-step-max",      "number of maximum (integrate/output) step for AR integration")
        , sym_order           (input_par_store, -6,                 "k",               "Symplectic integrator order, should be even number, positive value for Yoshida 1st method (can be arbitrary precision);  negative value for Yoshida 2nd method (only limited to double precision)")
        , energy_error        (input_par_store, 1e-10,              "e",               "relative energy error limit for AR")
        , time_error          (input_par_store, 0.0,                "time-error",      "time synchronization absolute error limit for AR","default is 0.25*dt-min")
        , time_zero           (input_par_store, 0.0,                "time-start",      "initial physical time")
        , time_end            (input_par_store, 0.0,                "t",               "ending physical time")
        , r_break             (input_par_store, 1e-3,               "r",               "distance criterion for checking stability")
        , nstep               (input_par_store, 0,                  "n",               "number of integration steps (higher priority than time_end)")
        , s                   (input_par_store, 0.0,                "s",               "step size, not physical time step;  <=0: auto, try to achieve min binary period/32;   >0: fixed")
        , ds_scale            (input_par_store, 1.0,                "ds-scale",        "step size scaling factor")
        , gravitational_constant(input_par_store, 1.0,              "G",               "gravitational constant")
        , dt_min              (input_par_store, 1e-13,              "dt-min",          "minimum physical time step")
        , dt_out              (input_par_store, 0.0,                "o",               "output time interval")
        , slowdown_ref        (input_par_store, 1e-6,               "slowdown-ref",    "slowdown perturbation ratio reference")
#ifdef AR_SLOWDOWN_MASSRATIO
        , slowdown_mass_ref   (input_par_store, 0.0,                "slowdown-mass-ref", "slowdowm mass reference","averaged mass")
#endif
        , slowdown_timescale_max(input_par_store, 0.0,              "slowdown-timescale-max", "maximum timescale for maximum slowdown factor","time-end")
        , interrupt_detection_option(input_par_store, 0,            "i",               "modify orbits and check interruption;  0: turn off;  1: modify the binary orbits based on interruption criterion;  2. recored binary parameters based on interruption criterion")
        , fix_step_option     (input_par_store, -1,                 "fix-step-option", "fix step options: always, later, none","auto")
#ifdef USE_MPFRC
        , mpfr_digits         (input_par_store, 30,                 "mpfr-dights",     "dights for MPFR precison")
#endif
#ifdef AR_G_FUNC_MUL_POT
        , g_func_option      (input_par_store, 0,                  "g-func",          "time transformation (g) function mode;  0=standard LogH;  1=BLogH (innermost binaries);  2=normalized BLogH;  3=all pairs;  4=BTLogH (tree-level product)")
#elif AR_G_FUNC_MAX_POT
        , g_func_option      (input_par_store, 0,                  "g-func",          "g-function mode: 0=standard LogH, 1=use maximum pair potential","0")
#elif AR_G_FUNC_ADD_POT
        , g_func_option      (input_par_store, 0,                  "g-func",          "g-function mode: 0=standard LogH, 1=use summation of innermost pair potentials","0")
#endif
#ifdef AR_G_FUNC
        , g_func_switch_option(input_par_store, "fixed",           "g-func-switch",   "how to apply --g-func: fixed (always use it), auto (switch between --g-func and 0 based on perturbation)","fixed")
#endif
        , filename_par        (input_par_store, "",                 "p",               "filename to load manager parameters","input name")
        , filename_out        (input_par_store, "",                 "f",               "filename to output snapshots in BINARY format;  if not given, print directly in standard output","input name")
        , integration_mode    (input_par_store, "base",             "m",               "Integration mode;  base: no time synchronization, orbital parameters and step size update;  orbit: no time synchronization, update orbital parameters and step size every step;  full: use time synchronization (integrateToTime) with full orbital parameters and step size update")
        , load_flag           (input_par_store, 0,                  "l",               "Load dumped data for restart (if used, the input file is dumped data)")
    {}

    int read(int argc, char* argv[], const char* bin_name) {
        static int ar_flag = -1;
        static struct option long_options[] = {
            {print_width.key,              required_argument, &ar_flag, 1},
            {print_precision.key,          required_argument, &ar_flag, 2},
            {nstep_max.key,                required_argument, &ar_flag, 3},
            {sym_order.key,                required_argument, &ar_flag, 4},
            {energy_error.key,             required_argument, &ar_flag, 5},
            {time_error.key,               required_argument, &ar_flag, 6},
            {time_zero.key,                required_argument, &ar_flag, 7},
            {time_end.key,                 required_argument, &ar_flag, 8},
            {r_break.key,                  required_argument, &ar_flag, 9},
            {nstep.key,                    required_argument, &ar_flag, 10},
            {s.key,                        required_argument, &ar_flag, 11},
            {ds_scale.key,                 required_argument, &ar_flag, 12},
            {gravitational_constant.key,   required_argument, &ar_flag, 13},
            {dt_min.key,                   required_argument, &ar_flag, 14},
            {slowdown_ref.key,             required_argument, &ar_flag, 16},
#ifdef AR_SLOWDOWN_MASSRATIO
            {slowdown_mass_ref.key,        required_argument, &ar_flag, 17},
#endif
            {slowdown_timescale_max.key,   required_argument, &ar_flag, 18},
            {interrupt_detection_option.key, required_argument, &ar_flag, 24},
            {fix_step_option.key,          required_argument, &ar_flag, 19},
#ifdef USE_MPFRC
            {mpfr_digits.key,              required_argument, &ar_flag, 20},
#endif
#ifdef AR_G_FUNC
            {g_func_option.key,            required_argument, &ar_flag, 21},
            {g_func_switch_option.key,     required_argument, &ar_flag, 27},
#endif
            {filename_par.key,             required_argument, &ar_flag, 22},
            {filename_out.key,             required_argument, &ar_flag, 23},
            {load_flag.key,                no_argument,       &ar_flag, 26},
            {"help",                       no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };

        int opt_used = 0;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-n:t:r:s:m:k:G:e:p:f:i:o:lh", long_options, &option_index)) != -1)
            switch (copt) {
            case 0:
                switch (ar_flag) {
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
                    sym_order.value = atoi(optarg);
                    opt_used += 2;
                    break;
                case 5:
                    energy_error.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 6:
                    time_error.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 7:
                    time_zero.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 8:
                    time_end.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 9:
                    r_break.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 10:
                    nstep.value = atoi(optarg);
                    opt_used += 2;
                    break;
                case 11:
                    s.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 12:
                    ds_scale.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 13:
                    gravitational_constant.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 14:
                    dt_min.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 16:
                    slowdown_ref.value = atof(optarg);
                    opt_used += 2;
                    break;
#ifdef AR_SLOWDOWN_MASSRATIO
                case 17:
                    slowdown_mass_ref.value = atof(optarg);
                    opt_used += 2;
                    break;
#endif
                case 18:
                    slowdown_timescale_max.value = atof(optarg);
                    opt_used += 2;
                    break;
                case 19:
                    if (!strcmp(optarg,"none")) fix_step_option.value = 2;
                    else if (!strcmp(optarg,"always")) fix_step_option.value = 0;
                    else if (!strcmp(optarg,"later")) fix_step_option.value = 1;
                    else {
                        std::cerr<<"Error: fix step option unknown ("<<optarg<<"), should be always, later, none\n";
                        abort();
                    }
                    opt_used += 2;
                    break;
#ifdef USE_MPFRC
                case 20:
                    mpfr_digits.value = atoi(optarg);
                    opt_used += 2;
                    break;
#endif
#ifdef AR_G_FUNC
                case 21:
                    g_func_option.value = atoi(optarg);
                    opt_used += 2;
                    break;
                case 27:
                    g_func_switch_option.value = optarg;
                    opt_used += 2;
                    break;
#endif
                case 22:
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
                    opt_used += 2;
                    break;
                case 23:
                    filename_out.value = optarg;
                    opt_used += 2;
                    break;
                case 24:
                    interrupt_detection_option.value = atoi(optarg);
                    opt_used += 2;
                    break;
                case 26:
                    load_flag.value = 1;
                    opt_used++;
                    break;
                }
                break;
            case 'n':
                nstep.value = atoi(optarg);
                opt_used++;
                break;
            case 'o':
                dt_out.value = atof(optarg);
                opt_used++;
                break;
            case 't':
                time_end.value = atof(optarg);
                opt_used++;
                break;
            case 'r':
                r_break.value = atof(optarg);
                opt_used++;
                break;
            case 's':
                s.value = atof(optarg);
                opt_used++;
                break;
            case 'k':
                sym_order.value = atoi(optarg);
                opt_used++;
                break;
            case 'G':
                gravitational_constant.value = atof(optarg);
                opt_used++;
                break;
            case 'e':
                energy_error.value = atof(optarg);
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
            case 'f':
                filename_out.value = optarg;
                opt_used++;
                break;
            case 'i':
                interrupt_detection_option.value = atoi(optarg);
                opt_used++;
                break;
            case 'm':
                integration_mode.value = optarg;
                opt_used += 2;
                break;
            case 'l':
                load_flag.value = 1;
                opt_used++;
                break;
            case 'h':
                std::cout<<bin_name<<" [option] data_filename\n"
                         <<"Input data file format: \n"
                         <<"    header line: number_of_particle\n"
                         <<"    following lines: mass, x, y, z, vx, vy, vz, radius\n";
                input_par_store.printHelp(std::cout);
                std::cout<<"Size of integrator class: (bytes) "<<sizeof(ARInt)<<std::endl;
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }
        return opt_used;
    }
};

int main(int argc, char **argv){

    //unsigned int oldcw;
    //fpu_fix_start(&oldcw);

    IOParamsAR iop;

    // Check whether all options are defined
    std::vector<COMM::IOParamsContainer*> all_pars;
    all_pars.push_back(&iop.input_par_store);
    std::vector<std::string> known_options;
    known_options.push_back("help");
    known_options.push_back("h");
    FindUndefinedOptions(all_pars, argc, argv, &known_options);

    FILE* fsnap = NULL;

#ifdef AR_TTL
    std::string bin_name("ar.ttl");
#else
    std::string bin_name("ar.logh");
#endif
#ifdef AR_SLOWDOWN_ARRAY
    bin_name += ".sd.a";
#elif AR_SLOWDOWN_TREE
    bin_name += ".sd.t";
#endif

    int opt_used = iop.read(argc, argv, bin_name.c_str());
    if (opt_used < 0) return 0;

    // Open output file if filename_out was specified
    if (!iop.filename_out.value.empty()) {
        if( (fsnap = fopen(iop.filename_out.value.c_str(),"r")) == NULL) {
            fprintf(stderr,"Error: Cannot open file %s.\n", iop.filename_out.value.c_str());
            abort();
        }
    }

    if (argc==1) {
        std::cerr<<"Please provide particle data filename\n";
        abort();
    }

    // data file name
    char* filename = argv[argc-1];

#ifdef USE_MPFRC
    setMPFRPrec(iop.mpfr_digits.value);
#endif

    // manager
    TimeTransformedSymplecticManager<Interaction> manager;
    manager.interaction.gravitational_constant = iop.gravitational_constant.value;
    manager.time_step_min = iop.dt_min.value;
    manager.ds_scale = iop.ds_scale.value;
    if (iop.time_error.value>0.0)  manager.time_error_max = iop.time_error.value;
    else manager.time_error_max = 0.25*iop.dt_min.value;
    manager.energy_error_relative_max = iop.energy_error.value; 
    if (iop.slowdown_timescale_max.value>0.0) manager.slowdown_timescale_max = iop.slowdown_timescale_max.value;
    else if (iop.time_end.value>0.0) manager.slowdown_timescale_max = iop.time_end.value;
    else manager.slowdown_timescale_max = NUMERIC_FLOAT_MAX;
    manager.slowdown_pert_ratio_ref = iop.slowdown_ref.value;
    manager.step_count_max = iop.nstep_max.value;
    // set symplectic order
    manager.step.initialSymplecticCofficients(iop.sym_order.value);

    manager.interaction.interrupt_detection_option = iop.interrupt_detection_option.value;

    // store input parameters
    std::string fpar_out = std::string(filename) + ".par";
    std::FILE* fout = std::fopen(fpar_out.c_str(),"w");
    if (fout==NULL) {
        std::cerr<<"Error: data file "<<fpar_out<<" cannot be open!\n";
        abort();
    }
    iop.input_par_store.writeAscii(fout);
    fclose(fout);
    
    // integrator
    ARInt sym_int;
    sym_int.manager = &manager;

    if(iop.load_flag.value) {
        std::FILE* fin = std::fopen(filename,"r");
        if (fin==NULL) {
            std::cerr<<"Error: data file "<<filename<<" cannot be open!\n";
            abort();
        }
        sym_int.readBinary(fin);
        fclose(fin);
    }
    else {
        std::fstream fin;
        fin.open(filename,std::fstream::in);
        if(!fin.is_open()) {
            std::cerr<<"Error: data file "<<filename<<" cannot be open!\n";
            abort();
        }
        sym_int.particles.setMode(COMM::ListMode::local);
        sym_int.particles.readMemberAscii(fin);
        sym_int.reserveIntegratorMem();
        fin.close();

    }
    sym_int.particles.calcCenterOfMass();
#ifdef AR_SLOWDOWN_MASSRATIO
    Float m_ave = sym_int.particles.cm.mass/sym_int.particles.getSize();
    if (iop.slowdown_mass_ref.value<=0.0) manager.slowdown_mass_ref = m_ave;
    else manager.slowdown_mass_ref = iop.slowdown_mass_ref.value;
#endif
    manager.print(std::cerr);

    for (int i=0; i<sym_int.particles.getSize(); i++) sym_int.particles[i].id = i+1;

    sym_int.info.reserveMem(sym_int.particles.getSize());
    sym_int.info.generateBinaryTree(sym_int.particles,manager.interaction.gravitational_constant);

#ifdef AR_G_FUNC
    sym_int.g_func_user = iop.g_func_option.value;
    if (iop.g_func_switch_option.value == "auto")
        sym_int.g_func_switch = AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, Interaction, AR::Information<Particle, Particle>>::GFUNC_AUTO;
    else
        sym_int.g_func_switch = AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, Interaction, AR::Information<Particle, Particle>>::GFUNC_FIXED;
#endif

    // r_break
    sym_int.info.r_break_crit = iop.r_break.value;

    // no initial when both parameters and data are load
    if(!iop.load_flag.value) {
        // initialization 
        sym_int.initialIntegration(iop.time_zero.value);
#ifdef AR_G_FUNC
        sym_int.info.calcDsAndStepOption(manager.step.getOrder(), manager.interaction.gravitational_constant, manager.ds_scale, sym_int.g_func);
#else
        sym_int.info.calcDsAndStepOption(manager.step.getOrder(), manager.interaction.gravitational_constant, manager.ds_scale);
#endif
    }

    // use input fix step option
    if (iop.fix_step_option.value>=0) {
        switch (iop.fix_step_option.value) {
        case 2:
            sym_int.info.fix_step_option = FixStepOption::none;
            break;
        case 0:
            sym_int.info.fix_step_option = FixStepOption::always;
            break;
        case 1:
            sym_int.info.fix_step_option = FixStepOption::later;
            break;
        }
    }

    // use input ds
    if (iop.s.value>0.0) sym_int.info.ds = iop.s.value;

    // precision
    std::cout<<std::setprecision(iop.print_precision.value);

#ifdef AR_SLOWDOWN_ARRAY
    int n_sd = sym_int.binary_slowdown.getSize();
#elif AR_SLOWDOWN_TREE
    int n_sd = sym_int.info.binarytree.getSize();
#else
    int n_sd = 0;
#endif
    //print column title
    sym_int.printColumnTitleAscii(std::cout, iop.print_width.value, n_sd);
    std::cout<<std::endl;

    //print initial data
    sym_int.printColumnAscii(std::cout, iop.print_width.value, n_sd);
    std::cout<<std::endl;

    
    // integration loop
    const int n_particle = sym_int.particles.getSize();
    if (iop.integration_mode.value != "full") {
        bool do_orbit_update = (iop.integration_mode.value == "orbit");
        Float time_out = iop.time_zero.value + iop.dt_out.value;
        Float time_table[manager.step.getCDPairSize()];
        sym_int.profile.step_count = 1;
        auto IntegrateOneStep = [&] (){
            if (do_orbit_update) {
                bool update_flag = sym_int.updateBinarySemiEccPeriodIter(
                    sym_int.info.getBinaryTreeRoot(),
                    manager.interaction.gravitational_constant,
                    sym_int.getTime());
                if (update_flag) {
                    auto& root = sym_int.info.getBinaryTreeRoot();
                    root.stableCheckIter(root, 10000 * root.period);
                }
#ifdef AR_SLOWDOWN_TREE
                sym_int.syncTreeSlowDownAndDs(true, false);
#else
                if (update_flag) {
                    sym_int.info.calcDsAndStepOption(
                        manager.step.getOrder(),
                        manager.interaction.gravitational_constant,
                        manager.ds_scale
    #ifdef AR_G_FUNC
                        , sym_int.g_func
    #endif
                    );
                }
#endif
            }
#ifdef SDAR_TIME_MEASURE
            sym_int.profile.prof_tot.start();
#endif
            if(n_particle==2) sym_int.integrateTwoOneStep(sym_int.info.ds, time_table);
            else sym_int.integrateOneStep(sym_int.info.ds, time_table);
#ifdef SDAR_TIME_MEASURE
            sym_int.profile.prof_tot.end();
#endif
            if (sym_int.getTime()>=time_out) {
                if (fsnap != NULL) 
                    sym_int.writeBinary(fsnap);
                else {
                    sym_int.printColumnAscii(std::cout, iop.print_width.value, n_sd);
                    std::cout<<std::endl;
                }
                time_out += iop.dt_out.value;
            }
            sym_int.profile.step_count_sum++;
        };
        if (iop.nstep.value>0) for (int i=0; i<iop.nstep.value; i++) IntegrateOneStep();
        else while (sym_int.getTime()<iop.time_end.value) IntegrateOneStep();
    }
    else {
        if (iop.dt_out.value>0.0) iop.nstep.value = int(iop.time_end.value/iop.dt_out.value+0.5);
        else if (iop.nstep.value>0) iop.dt_out.value = iop.time_end.value/iop.nstep.value;
        for (int i=1; i<=iop.nstep.value; i++) {
            auto bin_interrupt = sym_int.integrateToTime(iop.dt_out.value*i);
            if (bin_interrupt.status!=InterruptStatus::none) {
                std::cerr<<"Interrupt condition triggered! ";
                switch (bin_interrupt.status) {
                case InterruptStatus::change:
                    std::cerr<<" Change";
                    break;
                case InterruptStatus::merge:
                    std::cerr<<" merge";
                    break;
                case InterruptStatus::destroy:
                    std::cerr<<" Destroy";
                    break;
                case InterruptStatus::none:
                    break;
                }
                std::cerr<<std::endl;
                bin_interrupt.printColumnTitleAscii(std::cerr);
                std::cerr<<std::endl;
                bin_interrupt.printColumnAscii(std::cerr);
                std::cerr<<std::endl;

                Particle* p1 = bin_interrupt.getBinaryTreeAddress()->getLeftMember();
                Particle* p2 = bin_interrupt.getBinaryTreeAddress()->getRightMember();
                // merger case, quit integration
                if (n_particle==2&&(p1->mass==0||p2->mass==0)) {
                    sym_int.printColumnAscii(std::cout, iop.print_width.value, n_sd);
                    std::cout<<std::endl;
                    break;
                }
            }
#ifndef USE_CM_FRAME
            sym_int.info.generateBinaryTree(sym_int.particles, manager.interaction.gravitational_constant);
#endif
            if (fsnap != NULL) 
                sym_int.writeBinary(fsnap);
            else {
                sym_int.printColumnAscii(std::cout, iop.print_width.value, n_sd);
                std::cout<<std::endl;
            }
        }
    }


    // dump final data
    std::string fdata_out = std::string(filename) + ".last";
    fout = std::fopen(fdata_out.c_str(),"w");
    if (fout==NULL) {
        std::cerr<<"Error: data file "<<fdata_out<<" cannot be open!\n";
        abort();
    }
    sym_int.writeBinary(fout);
    fclose(fout);

#ifdef USE_CM_FRAME
    sym_int.info.getBinaryTreeRoot().shiftToOriginFrame();
#endif

    //fpu_fix_end(&oldcw);

    return 0;
}

