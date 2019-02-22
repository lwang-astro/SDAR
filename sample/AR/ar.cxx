#include <iostream>
#include <fstream>
//#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <cassert>

#define ASSERT(x) assert(x)

#include "AR/symplectic_integrator.h"
#include "particle.h"
#include "perturber.h"
#include "interaction.h"
#include "information.h"

using namespace AR;

int main(int argc, char **argv){

    //unsigned int oldcw;
    //fpu_fix_start(&oldcw);

    int print_width=22; //print width
    int print_precision=14; //print digital precision
    int nstep=1000; // total step size
    int nstep_max=1000000; // maximum time step allown for tsyn integration
    int sym_order=-6; // symplectic integrator order
    Float energy_error=1e-10; // phase error requirement
    Float time_error=5e-14; // time synchronization error
    Float s=0.5;    // step size
    Float time_zero=0.0;    // initial physical time
    Float time_end=-1.0; // ending physical time
    Float dt_min = 1e-13; // minimum physical time step
    char* filename_par=NULL; // par dumped filename
    FixStepOption fix_step_option=FixStepOption::none; // fix step option
    bool load_flag=false; // if true; load dumped data

    int copt;
    static struct option long_options[] = {
        {"time-start", required_argument, 0, 0},
        {"time-end", required_argument, 0, 't'},
        {"fix-step-option", required_argument, 0, 0},
        {"energy-error",required_argument, 0, 'e'},
        {"time-error",required_argument, 0, 0},
        {"dt-min",required_argument, 0, 0},
        {"n-step-max",required_argument, 0, 0},
        {"print-width",required_argument, 0, 0},
        {"print-precision",required_argument, 0, 0},
        {"load-data",no_argument, 0, 'l'},
        {"load-par",required_argument, 0, 0},    
        {"help",no_argument, 0, 'h'},
        {0,0,0,0}
    };
  
    int option_index;
    while ((copt = getopt_long(argc, argv, "N:n:t:s:k:e:lh", long_options, &option_index)) != -1)
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
                if (strcmp(optarg,"none")) fix_step_option=FixStepOption::none;
                else if (strcmp(optarg,"always")) fix_step_option=FixStepOption::always;
                else if (strcmp(optarg,"later")) fix_step_option=FixStepOption::later;
                else {
                    std::cerr<<"Error: fix step option unknown ("<<optarg<<"), should be always, later, none\n";
                    abort();
                }
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
                print_width = atof(optarg);
                break;
            case 8:
                print_precision = atoi(optarg);
                break;
            case 10:
                filename_par = optarg;
                break;
            default:
                std::cerr<<"Unknown option. check '-h' for help.\n";
                abort();
            }
            break;
        case 'n':
            nstep = atoi(optarg);
            break;
        case 't':
            time_end = atof(optarg);
            break;
        case 's':
            s = atof(optarg);
            break;
        case 'k':
            sym_order = atoi(optarg);
            break;
        case 'e':
            energy_error = atof(optarg);
            break;
        case 'l':
            load_flag = true;
            break;
        case 'h':
            std::cout<<"chain [option] data_filename\n"
                     <<"Input data file format: each line: mass, x, y, z, vx, vy, vz\n"
                     <<"Options: (*) show defaulted values\n"
                     <<"    -n [int]:     number of integration steps ("<<nstep<<")\n"
                     <<"          --time-start [Float]:  initial physical time ("<<time_zero<<")\n"
                     <<"    -t [Float]:  ending physical time; if set, -n will be invalid (unset)\n"
                     <<"          --time-end (same as -t)\n"
                     <<"    -s [Float]:  step size, not physical time step ("<<s<<")\n"
                     <<"          --fix-step-option [char]:  always, later, none (none) \n"
                     <<"    -k [int]:  Symplectic integrator order,should be even number ("<<sym_order<<")\n"
                     <<"    -e [Float]:  relative energy error limit ("<<energy_error<<")\n"
                     <<"          --energy-error (same as -e)\n"
                     <<"          --time-error [Float]:    time synchronization absolute error limit ("<<time_error<<")\n"
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
    SymplecticManager<Interaction> manager;
    manager.time_step_real_min = dt_min;
    manager.time_error_max_real = time_error;
    manager.energy_error_relative_max = energy_error; 
    manager.step_count_max = nstep_max;
    
    // set symplectic order
    manager.step.initialSymplecticCofficients(sym_order);

    // integrator
    SymplecticIntegrator<Particle, Particle, Perturber, Interaction, Information> sym_int;
    sym_int.manager = &manager;

    if(load_flag) {
        
    }
    else {
        sym_int.particles.setMode(COMM::ListMode::local);
        sym_int.particles.readAscii(fin);
        sym_int.reserveIntegratorMem();
    }

    // initialization 
    sym_int.initialIntegration(time_zero);

    // precision
    std::cout<<std::setprecision(print_precision);

    //print column title
    sym_int.printColumnTitle(std::cout, print_width);
    std::cout<<std::endl;

    //print initial data
    sym_int.printColumn(std::cout, print_width);
    std::cout<<std::endl;

    
    // integration loop
    const int n_particle = sym_int.particles.getSize();
    if (time_end<0.0) {
        Float time_table[manager.step.getCDPairSize()];
        for (int i=0; i<nstep; i++) {
            if(n_particle==2) sym_int.integrateTwoOneStep(s, time_table);
            else sym_int.integrateOneStep(s, time_table);
            sym_int.printColumn(std::cout, print_width);
            std::cout<<std::endl;
        }
    }
    else {
        sym_int.integrateToTime(s, time_end, fix_step_option);
        sym_int.printColumn(std::cout, print_width);
        std::cout<<std::endl;
    }

    //fpu_fix_end(&oldcw);

    return 0;
}

