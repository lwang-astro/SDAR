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

int main(int argc, char **argv){

    //unsigned int oldcw;
    //fpu_fix_start(&oldcw);

    COMM::IOParamsContainer input_par_store;

    COMM::IOParams<int> print_width    (input_par_store, WRITE_WIDTH,     "print width of value"); //print width
    COMM::IOParams<int> print_precision(input_par_store, WRITE_PRECISION, "print digital precision"); //print digital precision
    COMM::IOParams<int> nstep_max      (input_par_store, 1000000, "number of maximum (integrate/output) step for AR integration"); // maximum time step allown for tsyn integration
    COMM::IOParams<int> sym_order      (input_par_store, -6,   "Symplectic integrator order, should be even number"); // symplectic integrator order
    COMM::IOParams<double> energy_error (input_par_store, 1e-10,"relative energy error limit for AR"); // phase error requirement
    COMM::IOParams<double> time_error   (input_par_store, 0.0,  "time synchronization absolute error limit for AR","default is 0.25*dt-min"); // time synchronization error
    COMM::IOParams<double> time_zero    (input_par_store, 0.0,  "initial physical time");    // initial physical time
    COMM::IOParams<double> time_end     (input_par_store, 0.0,  "ending physical time"); // ending physical time
    COMM::IOParams<int>   nstep        (input_par_store,  0, "number of integration steps (higher priority than time_end)"); // total step size
    COMM::IOParams<double> s            (input_par_store, 0.0,  "step size, not physical time step","auto");    // step size
    COMM::IOParams<double> gravitational_constant   (input_par_store, 1.0, "gravitational constant"); // gravitational constant
    COMM::IOParams<double> dt_min       (input_par_store, 1e-13,"minimum physical time step"); // minimum physical time step
    COMM::IOParams<double> dt_out       (input_par_store, 0.0,"output time interval"); // output time interval
    COMM::IOParams<double> slowdown_ref (input_par_store, 1e-6, "slowdown perturbation ratio reference"); // slowdown reference factor
#ifdef AR_SLOWDOWN_MASSRATIO
    COMM::IOParams<double> slowdown_mass_ref (input_par_store, 0.0, "slowdowm mass reference","averaged mass"); // slowdown mass reference
#endif
    COMM::IOParams<double> slowdown_timescale_max (input_par_store, 0.0, "maximum timescale for maximum slowdown factor","time-end"); // slowdown timescale
    COMM::IOParams<int>   fix_step_option (input_par_store, -1, "always, later, none","auto"); // if true; use input fix step option
    COMM::IOParams<std::string> filename_par (input_par_store, "", "filename to load manager parameters","input name"); // par dumped filename
    bool load_flag=false;  // if true; load dumped data
    bool synch_flag=false; // if true, switch on time synchronization

#ifdef AR_TTL
    std::string bin_name("ar.ttl");
#else
    std::string bin_name("ar.logh");
#endif
#ifdef AR_SLOWDOWN_INNER
    bin_name += ".sd";
#endif

    int copt;
    static struct option long_options[] = {
        {"time-start", required_argument, 0, 0},
        {"time-end", required_argument, 0, 't'},
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
        {"load-data",no_argument, 0, 'l'},
        {"help",no_argument, 0, 'h'},
        {0,0,0,0}
    };
  
    int option_index;
    while ((copt = getopt_long(argc, argv, "N:n:t:s:Sk:G:e:p:o:lh", long_options, &option_index)) != -1)
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
        case 'G':
            gravitational_constant.value = atof(optarg);
            break;
        case 'n':
            nstep.value = atoi(optarg);
            break;
        case 't':
            time_end.value = atof(optarg);
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
        case 'h':
            std::cout<<bin_name<<" [option] data_filename\n"
                     <<"Input data file format: each line: mass, x, y, z, vx, vy, vz\n"
                     <<"Options: (*) show defaulted values\n"
                     <<"          --dt-min          [int]  :  "<<dt_min<<"\n"
                     <<"    -e [Float]:  "<<energy_error<<"\n"
                     <<"          --energy-error    [Float]:  same as -e\n"
                     <<"          --fix-step-option [char] :  "<<fix_step_option<<"\n"
                     <<"    -G [Float]:  "<<gravitational_constant<<"\n"
                     <<"    -l :          load dumped data for restart (if used, the input file is dumped data)\n"
                     <<"          --load-data (same as -l)\n"
                     <<"    -k [int]:    "<<sym_order<<"\n"
                     <<"    -n [int]:    "<<nstep<<"\n"
                     <<"    -o [float]:  "<<dt_out<<"\n"
                     <<"    -p [string]: "<<filename_par<<"\n"
                     <<"          --print-width     [int]  : "<<print_width<<"\n"
                     <<"          --print-precision [int]  : "<<print_precision<<"\n"
                     <<"    -s [Float]:  "<<s<<"\n"
                     <<"          --slowdown-ref:           [Float]: "<<slowdown_ref<<"\n"
#ifdef AR_SLOWDOWN_MASSRATIO
                     <<"          --slowdown-mass-ref       [Float]: "<<slowdown_mass_ref<<"\n"
#endif
                     <<"          --slowdown-timescale-max: [Float]: "<<slowdown_timescale_max<<"\n"
                     <<"    -S :         Switch on time synchronization\n"
                     <<"    -t [Float]:  "<<time_end<<"\n"
                     <<"          --time-start      [Float]:  "<<time_zero<<"\n"
                     <<"          --time-end        [Float]:  same as -t\n"
                     <<"          --time-error      [Float]:  "<<time_error<<"\n"
                     <<"    -h :          print option information\n"
                     <<"          --help (same as -h)\n";
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

    // manager
    SymplecticManager<Interaction> manager;
    manager.interaction.gravitational_constant = gravitational_constant.value;
    manager.time_step_real_min = dt_min.value;
    if (time_error.value>0.0)  manager.time_error_max_real = time_error.value;
    else manager.time_error_max_real = 0.25*dt_min.value;
    manager.energy_error_relative_max = energy_error.value; 
    if (slowdown_timescale_max.value>0.0) manager.slowdown_timescale_max = slowdown_timescale_max.value;
    else if (time_end.value>0.0) manager.slowdown_timescale_max = time_end.value;
    else manager.slowdown_timescale_max = NUMERIC_FLOAT_MAX;
    manager.slowdown_pert_ratio_ref = slowdown_ref.value;
    manager.step_count_max = nstep_max.value;
    // set symplectic order
    manager.step.initialSymplecticCofficients(sym_order.value);


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
    SymplecticIntegrator<Particle, Particle, Perturber, Interaction, Information<Particle,Particle>> sym_int;
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

    // no initial when both parameters and data are load
    if(!load_flag) {
        // initialization 
        sym_int.initialIntegration(time_zero.value);
        sym_int.info.calcDsAndStepOption(sym_int.slowdown.getSlowDownFactorOrigin(), manager.step.getOrder(), manager.interaction.gravitational_constant);
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

#ifdef AR_SLOWDOWN_INNER
    int n_sd = sym_int.slowdown_inner.getSize();
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
            sym_int.updateSlowDownAndCorrectEnergy();
            if(n_particle==2) sym_int.integrateTwoOneStep(sym_int.info.ds, time_table);
            else sym_int.integrateOneStep(sym_int.info.ds, time_table);
            if (sym_int.slowdown.getRealTime()>=time_out) {
                sym_int.printColumn(std::cout, print_width.value, n_sd);
                std::cout<<std::endl;
                time_out += dt_out.value;
            }
            sym_int.profile.step_count_sum++;
        };

        if (nstep.value>0) for (int i=0; i<nstep.value; i++) IntegrateOneStep();
        else while (sym_int.slowdown.getRealTime()<time_end.value) IntegrateOneStep();
    }
    else {
        int nstep_per_out = int(time_end.value/dt_out.value+0.5);
        for (int i=1; i<=nstep_per_out; i++) {
            sym_int.integrateToTime(dt_out.value*i);
            sym_int.info.generateBinaryTree(sym_int.particles, manager.interaction.gravitational_constant);
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

    //fpu_fix_end(&oldcw);

    return 0;
}

