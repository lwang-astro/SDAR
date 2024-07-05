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

int main(int argc, char **argv){

    //unsigned int oldcw;
    //fpu_fix_start(&oldcw);

    COMM::IOParamsContainer input_par_store;

    COMM::IOParams<int> print_width    (input_par_store, WRITE_WIDTH,     "print width of value"); //print width
    COMM::IOParams<int> print_precision(input_par_store, WRITE_PRECISION, "print digital precision"); //print digital precision
    COMM::IOParams<int> nstep_max      (input_par_store, 1000000, "number of maximum (integrate/output) step for AR integration"); // maximum time step allown for tsyn integration
    COMM::IOParams<int> sym_order      (input_par_store, -6,   "Symplectic integrator order, should be even number, positive value for Yoshida 1st method (can be arbitrary precision); negative value for Yoshida 2nd method (only limited to double precision)"); // symplectic integrator order
    COMM::IOParams<double> energy_error (input_par_store, 1e-10,"relative energy error limit for AR"); // phase error requirement
    COMM::IOParams<double> time_error   (input_par_store, 0.0,  "time synchronization absolute error limit for AR","default is 0.25*dt-min"); // time synchronization error
    COMM::IOParams<double> time_zero    (input_par_store, 0.0,  "initial physical time");    // initial physical time
    COMM::IOParams<double> time_end     (input_par_store, 0.0,  "ending physical time"); // ending physical time
    COMM::IOParams<double> r_break      (input_par_store, 1e-3, "distance criterion for checking stability"); // binary break criterion
    COMM::IOParams<int>   nstep        (input_par_store,  0, "number of integration steps (higher priority than time_end)"); // total step size
    COMM::IOParams<double> s            (input_par_store, 0.0,  "step size, not physical time step","auto");    // step size
    COMM::IOParams<double> ds_scale     (input_par_store, 1.0,  "step size scaling factor");    // step size scaling factor
    COMM::IOParams<double> gravitational_constant   (input_par_store, 1.0, "gravitational constant"); // gravitational constant
    COMM::IOParams<double> dt_min       (input_par_store, 1e-13,"minimum physical time step"); // minimum physical time step
    COMM::IOParams<double> dt_out       (input_par_store, 0.0,"output time interval"); // output time interval
    COMM::IOParams<double> slowdown_ref (input_par_store, 1e-6, "slowdown perturbation ratio reference"); // slowdown reference factor
#ifdef AR_SLOWDOWN_MASSRATIO
    COMM::IOParams<double> slowdown_mass_ref (input_par_store, 0.0, "slowdowm mass reference","averaged mass"); // slowdown mass reference
#endif
    COMM::IOParams<double> slowdown_timescale_max (input_par_store, 0.0, "maximum timescale for maximum slowdown factor","time-end"); // slowdown timescale
    COMM::IOParams<int>   interrupt_detection_option(input_par_store, 0, "modify orbits and check interruption: 0: turn off; 1: modify the binary orbits based on interruption criterion; 2. recored binary parameters based on interruption criterion");  // modify orbit or check interruption using modifyAndInterruptIter function
    COMM::IOParams<int>   fix_step_option (input_par_store, -1, "fix step options: always, later, none","auto"); // if true; use input fix step option
#ifdef USE_MPFRC
    COMM::IOParams<int>   mpfr_digits     (input_par_store, 30, "dights for MPFR precison");
#endif
    COMM::IOParams<std::string> filename_par (input_par_store, "", "filename to load manager parameters","input name"); // par dumped filename
    bool load_flag=false;  // if true; load dumped data
    bool synch_flag=false; // if true, switch on time synchronization

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

    int copt;
    static struct option long_options[] = {
        {"time-start", required_argument, 0, 0},
        {"time-end", required_argument, 0, 't'},
        {"r-break", required_argument, 0, 'r'},
        {"fix-step-option", required_argument, 0, 2},
        {"energy-error",required_argument, 0, 'e'},
        {"time-error",required_argument, 0, 4},
        {"dt-min",required_argument, 0, 5},
        {"n-step-max",required_argument, 0, 6},
        {"slowdown-ref",required_argument, 0, 7},
#ifdef AR_SLOWDOWN_MASSRATIO
        {"slowdown-mass-ref",required_argument, 0, 8},
#endif
        {"slowdown-timescale-max",required_argument, 0, 9},
        {"print-width",required_argument, 0, 10},
        {"print-precision",required_argument, 0, 11},
        {"ds-scale",required_argument, 0, 12},
#ifdef USE_MPFRC
        {"mpfr-dights", required_argument, 0, 13},
#endif
        {"load-data",no_argument, 0, 'l'},
        {"help",no_argument, 0, 'h'},
        {0,0,0,0}
    };
  
    int option_index;
    while ((copt = getopt_long(argc, argv, "N:n:t:r:s:Sk:G:e:p:o:i:lh", long_options, &option_index)) != -1)
        switch (copt) {
        case 0:
            time_zero.value = atof(optarg);
            break;
        case 2:
            if (!strcmp(optarg,"none")) fix_step_option.value=2;
            else if (!strcmp(optarg,"always")) fix_step_option.value=0;
            else if (!strcmp(optarg,"later")) fix_step_option.value=1;
            else {
                std::cerr<<"Error: fix step option unknown ("<<optarg<<"), should be always, later, none\n";
                abort();
            }
            break;
        case 4:
            time_error.value = atof(optarg);
            break;
        case 5:
            dt_min.value = atof(optarg);
                break;
        case 6:
            nstep_max.value = atoi(optarg);
            break;
        case 7:
            slowdown_ref.value = atof(optarg);
            break;
#ifdef AR_SLOWDOWN_MASSRATIO
        case 8:
            slowdown_mass_ref.value = atof(optarg);
            break;
#endif
        case 9:
            slowdown_timescale_max.value = atof(optarg);
            break;
        case 10:
            print_width.value = atof(optarg);
            break;
        case 11:
            print_precision.value = atoi(optarg);
            break;
        case 12:
            ds_scale.value = atof(optarg);
            break;
#ifdef USE_MPFRC
        case 13:
            mpfr_digits.value = atoi(optarg);
#endif
        case 'G':
            gravitational_constant.value = atof(optarg);
            break;
        case 'n':
            nstep.value = atoi(optarg);
            break;
        case 't':
            time_end.value = atof(optarg);
            break;
        case 'r':
            r_break.value = atof(optarg);
            break;
        case 's':
            s.value = atof(optarg);
            break;
        case 'S':
            synch_flag = true;
            break;
        case 'k':
            sym_order.value = atoi(optarg);
            break;
        case 'e':
            energy_error.value = atof(optarg);
            break;
        case 'l':
            load_flag = true;
            break;
        case 'p':
            filename_par.value = optarg;
            FILE* fpar_in;
            if( (fpar_in = fopen(filename_par.value.c_str(),"r")) == NULL) {
                fprintf(stderr,"Error: Cannot open file %s.\n", filename_par.value.c_str());
                abort();
            }
            input_par_store.readAscii(fpar_in);
            fclose(fpar_in);
            break;
        case 'o':
            dt_out.value = atof(optarg);
            break;
        case 'i':
            interrupt_detection_option.value = atoi(optarg);
            break;
        case 'h':
            std::cout<<bin_name<<" [option] data_filename\n"
                     <<"Input data file format: \n"
                     <<"    header line: number_of_particle\n"
                     <<"    following lines: mass, x, y, z, vx, vy, vz, radius\n"
                     <<"Options: (*) show defaulted values\n"
                     <<"          --dt-min          [int]  :  "<<dt_min<<"\n"
                     <<"          --ds-scale        [Float]:  "<<ds_scale<<"\n"
                     <<"    -e [Float]:  "<<energy_error<<"\n"
                     <<"          --energy-error    [Float]:  same as -e\n"
                     <<"          --fix-step-option [char] :  "<<fix_step_option<<"\n"
                     <<"    -G [Float]:  "<<gravitational_constant<<"\n"
                     <<"    -k [int]:    "<<sym_order<<"\n"
                     <<"    -i [string]: "<<interrupt_detection_option<<"\n"
                     <<"    -l :          load dumped data for restart (if used, the input file is dumped data)\n"
                     <<"          --load-data (same as -l)\n"
#ifdef USE_MPFRC
                     <<"          --mpfr-digits     [int]  :  "<<mpfr_digits<<"\n"
#endif
                     <<"    -n [int]:    "<<nstep<<"\n"
                     <<"    -o [float]:  "<<dt_out<<"\n"
                     <<"    -p [string]: "<<filename_par<<"\n"
                     <<"          --print-width     [int]  : "<<print_width<<"\n"
                     <<"          --print-precision [int]  : "<<print_precision<<"\n"
                     <<"    -r [Float]:  "<<r_break<<"\n"
                     <<"          --r-break      [Float]: same as -r\n"
                     <<"    -s [Float]:  "<<s<<"\n"
                     <<"          --slowdown-ref:           [Float]: "<<slowdown_ref<<"\n"
#ifdef AR_SLOWDOWN_MASSRATIO
                     <<"          --slowdown-mass-ref       [Float]: "<<slowdown_mass_ref<<"\n"
#endif
                     <<"          --slowdown-timescale-max: [Float]: "<<slowdown_timescale_max<<"\n"
                     <<"    -S :         Switch on time synchronization (use with -o or -n)\n"
                     <<"    -t [Float]:  "<<time_end<<"\n"
                     <<"          --time-start      [Float]:  "<<time_zero<<"\n"
                     <<"          --time-end        [Float]:  same as -t\n"
                     <<"          --time-error      [Float]:  "<<time_error<<"\n"
                     <<"    -h :          print option information\n"
                     <<"          --help (same as -h)\n";
            std::cout<<"Size of integrator class: (bytes) "<<sizeof(ARInt)<<std::endl;
            return 0;
        default:
            std::cerr<<"Unknown argument. check '-h' for help.\n";
            abort();
        }

    if (argc==1) {
        std::cerr<<"Please provide particle data filename\n";
        abort();
    }

    // data file name
    char* filename = argv[argc-1];

#ifdef USE_MPFRC
    setMPFRPrec(mpfr_digits.value);
#endif

    // manager
    TimeTransformedSymplecticManager<Interaction> manager;
    manager.interaction.gravitational_constant = gravitational_constant.value;
    manager.time_step_min = dt_min.value;
    manager.ds_scale = ds_scale.value;
    if (time_error.value>0.0)  manager.time_error_max = time_error.value;
    else manager.time_error_max = 0.25*dt_min.value;
    manager.energy_error_relative_max = energy_error.value; 
    if (slowdown_timescale_max.value>0.0) manager.slowdown_timescale_max = slowdown_timescale_max.value;
    else if (time_end.value>0.0) manager.slowdown_timescale_max = time_end.value;
    else manager.slowdown_timescale_max = NUMERIC_FLOAT_MAX;
    manager.slowdown_pert_ratio_ref = slowdown_ref.value;
    manager.step_count_max = nstep_max.value;
    // set symplectic order
    manager.step.initialSymplecticCofficients(sym_order.value);

    manager.interaction.interrupt_detection_option = interrupt_detection_option.value;


    // store input parameters
    std::string fpar_out = std::string(filename) + ".par";
    std::FILE* fout = std::fopen(fpar_out.c_str(),"w");
    if (fout==NULL) {
        std::cerr<<"Error: data file "<<fpar_out<<" cannot be open!\n";
        abort();
    }
    input_par_store.writeAscii(fout);
    fclose(fout);
    
    // integrator
    ARInt sym_int;
    sym_int.manager = &manager;

    if(load_flag) {
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
    if (slowdown_mass_ref.value<=0.0) manager.slowdown_mass_ref = m_ave;
    else manager.slowdown_mass_ref = slowdown_mass_ref.value;
#endif
    manager.print(std::cerr);

    for (int i=0; i<sym_int.particles.getSize(); i++) sym_int.particles[i].id = i+1;

    sym_int.info.reserveMem(sym_int.particles.getSize());
    sym_int.info.generateBinaryTree(sym_int.particles,manager.interaction.gravitational_constant);

    // r_break
    sym_int.info.r_break_crit = r_break.value;

    // no initial when both parameters and data are load
    if(!load_flag) {
        // initialization 
        sym_int.initialIntegration(time_zero.value);
        sym_int.info.calcDsAndStepOption(manager.step.getOrder(), manager.interaction.gravitational_constant, manager.ds_scale);
    }

    // use input fix step option
    if (fix_step_option.value>=0) {
        switch (fix_step_option.value) {
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
    if (s.value>0.0) sym_int.info.ds = s.value;

    // precision
    std::cout<<std::setprecision(print_precision.value);

#ifdef AR_SLOWDOWN_ARRAY
    int n_sd = sym_int.binary_slowdown.getSize();
#elif AR_SLOWDOWN_TREE
    int n_sd = sym_int.info.binarytree.getSize();
#else
    int n_sd = 0;
#endif
    //print column title
    sym_int.printColumnTitle(std::cout, print_width.value, n_sd);
    std::cout<<std::endl;

    //print initial data
    sym_int.printColumn(std::cout, print_width.value, n_sd);
    std::cout<<std::endl;

    
    // integration loop
    const int n_particle = sym_int.particles.getSize();
    if (!synch_flag) {
        float time_out = time_zero.value + dt_out.value;
        Float time_table[manager.step.getCDPairSize()];
        sym_int.profile.step_count = 1;
        auto IntegrateOneStep = [&] (){
#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
            sym_int.updateSlowDownAndCorrectEnergy(true, false);
#endif
            if(n_particle==2) sym_int.integrateTwoOneStep(sym_int.info.ds, time_table);
            else sym_int.integrateOneStep(sym_int.info.ds, time_table);
            if (sym_int.getTime()>=time_out) {
                sym_int.printColumn(std::cout, print_width.value, n_sd);
                std::cout<<std::endl;
                time_out += dt_out.value;
            }
            sym_int.profile.step_count_sum++;
        };

        if (nstep.value>0) for (int i=0; i<nstep.value; i++) IntegrateOneStep();
        else while (sym_int.getTime()<time_end.value) IntegrateOneStep();
    }
    else {
        if (dt_out.value>0.0) nstep.value = int(time_end.value/dt_out.value+0.5);
        else if (nstep.value>0) dt_out.value = time_end.value/nstep.value;
        for (int i=1; i<=nstep.value; i++) {
            auto bin_interrupt = sym_int.integrateToTime(dt_out.value*i);
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
                bin_interrupt.printColumnTitle(std::cerr);
                std::cerr<<std::endl;
                bin_interrupt.printColumn(std::cerr);
                std::cerr<<std::endl;

                Particle* p1 = bin_interrupt.getBinaryTreeAddress()->getLeftMember();
                Particle* p2 = bin_interrupt.getBinaryTreeAddress()->getRightMember();
                // merger case, quit integration
                if (n_particle==2&&(p1->mass==0||p2->mass==0)) {
                    sym_int.printColumn(std::cout, print_width.value, n_sd);
                    std::cout<<std::endl;
                    break;
                }
            }
#ifndef USE_CM_FRAME
            sym_int.info.generateBinaryTree(sym_int.particles, manager.interaction.gravitational_constant);
#endif
            sym_int.printColumn(std::cout, print_width.value, n_sd);
            std::cout<<std::endl;
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

