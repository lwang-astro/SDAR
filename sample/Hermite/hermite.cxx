#include <iostream>
#include <fstream>
//#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <iomanip>
#include <cmath>

#define ASSERT(expr) assert(expr)
#define DATADUMP(x) abort()

#include "AR/symplectic_integrator.h"
#include "Hermite/hermite_integrator.h"
#include "particle.h"
#include "hermite_perturber.h"
#include "ar_interaction.h"
#include "hermite_interaction.h"
#include "hermite_information.h"


using namespace H4;

int main(int argc, char **argv){

    //unsigned int oldcw;
    //fpu_fix_start(&oldcw);

    int print_width=22; //print width
    int print_precision=14; //print digital precision
    int nstep_max=1000000; // maximum time step allown for tsyn integration
    int sym_order=-6; // symplectic integrator order
    Float energy_error=1e-10; // phase error requirement
    Float time_error=1e-6; // time synchronization error
    Float time_zero=0.0;    // initial physical time
    Float time_end=1.0; // ending physical time
    Float dt_output = 0.25; // output time interval
    Float dt_min = 1e-24; // minimum physical time step
    char* filename_par=NULL; // par dumped filename
    Float r_break = 1e-3; // binary break criterion
    Float eta_4th = 0.02; // time step coefficient 
    Float eta_2nd = 0.005; // time step coefficient for 2nd order
    Float eps_sq = 0.0;    // softening parameter
    Float G = 1.0;      // gravitational constant
    bool load_flag=false; // if true; load dumped data

    int copt;
    static struct option long_options[] = {
        {"time-start", required_argument, 0, 0},
        {"time-end", required_argument, 0, 't'},
        {"r-break", required_argument, 0, 0},
        {"energy-error",required_argument, 0, 'e'},
        {"time-error",required_argument, 0, 0},
        {"dt-min",required_argument, 0, 0},
        {"n-step-max",required_argument, 0, 0},
        {"eta-4th",required_argument, 0, 0},
        {"eta-2nd",required_argument, 0, 0},
        {"eps",required_argument, 0, 0},
        {"print-width",required_argument, 0, 0},
        {"print-precision",required_argument, 0, 0},
        {"load-data",no_argument, 0, 'l'},
        {"load-par",required_argument, 0, 0},    
        {"help",no_argument, 0, 'h'},
        {0,0,0,0}
    };
  
    int option_index;
    while ((copt = getopt_long(argc, argv, "t:k:G:e:o:lh", long_options, &option_index)) != -1)
        switch (copt) {
        case 0:
#ifdef DEBUG
            std::cerr<<"option "<<option_index<<" "<<long_options[option_index].name<<" "<<optarg<<std::endl;
#endif
            switch (option_index) {
            case 0:
                time_zero = atof(optarg);
                break;
            case 2:
                r_break = atof(optarg);
                break;
            case 4:
                time_error = atof(optarg);
                break;
            case 5:
                dt_min = atof(optarg);
                break;
            case 6:
                nstep_max = atoi(optarg);
                break;
            case 7:
                eta_4th = atof(optarg);
                break;
            case 8:
                eta_2nd = atof(optarg);
                break;
            case 9:
                eps_sq = atof(optarg);
                break;
            case 10:
                print_width = atof(optarg);
                break;
            case 11:
                print_precision = atoi(optarg);
                break;
            case 13:
                filename_par = optarg;
                break;
            default:
                std::cerr<<"Unknown option. check '-h' for help.\n";
                abort();
            }
            break;
        case 't':
            time_end = atof(optarg);
            break;
        case 'k':
            sym_order = atoi(optarg);
            break;
        case 'G':
            G = atof(optarg);
            break;
        case 'e':
            energy_error = atof(optarg);
            break;
        case 'o':
            dt_output = atof(optarg);
            break;
        case 'l':
            load_flag = true;
            break;
        case 'h':
            std::cout<<"chain [option] data_filename\n"
                     <<"Input data file format: each line: mass, x, y, z, vx, vy, vz\n"
                     <<"Options: (*) show defaulted values\n"
                     <<"          --time-start [Float]:  initial physical time ("<<time_zero<<")\n"
                     <<"    -t [Float]:  ending physical time; if set, -n will be invalid (unset)\n"
                     <<"          --time-end (same as -t)\n"
                     <<"    -o [Float]:  output time interval ("<<dt_output<<")\n"
                     <<"    -k [int]:  Symplectic integrator order, should be even number ("<<sym_order<<")\n"
                     <<"    -e [Float]:  relative energy error limit for AR ("<<energy_error<<")\n"
                     <<"          --energy-error (same as -e)\n"
                     <<"          --time-error [Float]:    time synchronization relative error limit for AR ("<<time_error<<")\n"
                     <<"    -G [Float]: gravitational constant ("<<G<<")\n"
                     <<"          --eta-4th:   time step coefficient for 4th order ("<<eta_4th<<")\n"
                     <<"          --eta-2nd:   time step coefficient for 2nd order ("<<eta_2nd<<")\n"
                     <<"          --eps:       softerning parameter ("<<eps_sq<<")\n"
                     <<"          --print-width [int]:     print width of value ("<<print_width<<")\n"
                     <<"          --print-precision [int]: print digital precision ("<<print_precision<<")\n"
                     <<"    -l :          load dumped data for restart (if used, the input file is dumped data)\n"
                     <<"          --load-data (same as -l)\n"
                     <<"          --load-par    [char]:    filename to load manager parameters\n"
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

#ifdef DEBUG
    // parameter list
    std::cerr<<"Options:\n"
             <<"steps: "<<nstep<<std::endl
             <<"time-start: "<<time_zero<<std::endl
             <<"time-end: "<<time_end<<std::endl
             <<"step size: "<<s<<std::endl
             <<"energy-error: "<<energy_error<<std::endl
             <<"time-error: "<<time_error<<std::endl
             <<"print width: "<<print_width<<std::endl
             <<"print precision: "<<print_precision<<std::endl;
#endif

    // data file name
    char* filename = argv[argc-1];

    // open data file
    std::fstream fin;
    fin.open(filename,std::fstream::in);
    if(!fin.is_open()) {
        std::cerr<<"Error: Filename "<<filename<<" not found\n";
        abort();
    }

    // manager
    HermiteManager<HermiteInteraction> manager;
    manager.r_break_crit = r_break;
    manager.step.eta_4th = eta_4th;
    manager.step.eta_2nd = eta_2nd;
    manager.interaction.eps_sq = eps_sq;
    manager.interaction.G = G;
    AR::SymplecticManager<ARInteraction> ar_manager;
    ar_manager.interaction.eps_sq = eps_sq;
    ar_manager.interaction.G = G;
    ar_manager.time_step_real_min = dt_min;
    ar_manager.time_error_relative_max_real = time_error;
    ar_manager.energy_error_relative_max = energy_error; 
    ar_manager.step_count_max = nstep_max;
    // set symplectic order
    ar_manager.step.initialSymplecticCofficients(sym_order);

    // integrator
    HermiteIntegrator<Particle, Particle, HermitePerturber, HermiteInteraction, ARInteraction, HermiteInformation> h4_int;
    h4_int.manager = &manager;
    h4_int.ar_manager = &ar_manager;

    if(load_flag) {
        
    }
    else {
        h4_int.particles.setMode(COMM::ListMode::local);
        h4_int.particles.readAscii(fin);
        for (int i=0; i<h4_int.particles.getSize(); i++) h4_int.particles[i].id = i+1;
        h4_int.reserveMem(h4_int.particles.getSize());
    }

    // initialization 
    h4_int.initialSystem(time_zero);
    h4_int.info.time = h4_int.getTime();

    // add primordial groups
    h4_int.addGroupsAscii(fin);

    // initial integrator
    h4_int.adjustSystemAfterModify();

    // precision
    std::cout<<std::setprecision(print_precision);

    // get initial energy
    h4_int.info.calcEnergy(h4_int.particles, manager.interaction, true);

    //print column title
    h4_int.info.printColumnTitle(std::cout, print_width);
    h4_int.particles.printColumnTitle(std::cout, print_width);
    std::cout<<std::endl;

    //print initial data
    h4_int.info.printColumn(std::cout, print_width);
    h4_int.particles.printColumn(std::cout, print_width);
    std::cout<<std::endl;
    
    // integration loop
    while (h4_int.info.time<time_end) {
        h4_int.integrateOneStepAct();
        h4_int.adjustGroups();
        h4_int.sortDtAndSelectActParticle();
        h4_int.info.time = h4_int.getTime();

        if (fmod(h4_int.info.time, dt_output)==0.0) {
            h4_int.writeBackGroupMembers();
            h4_int.info.calcEnergy(h4_int.particles, manager.interaction, false);
            
            h4_int.info.printColumn(std::cout, print_width);
            h4_int.particles.printColumn(std::cout, print_width);
            std::cout<<std::endl;
        }
    }

    //fpu_fix_end(&oldcw);

    return 0;
}
