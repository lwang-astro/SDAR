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

    COMM::IOParams<int> print_width    (input_par_store, 22,   "print width of value"); //print width
    COMM::IOParams<int> print_precision(input_par_store, 14,   "print digital precision"); //print digital precision
    COMM::IOParams<int> nstep_max      (input_par_store, 1000000, "number of maximum (integrate/output) step for AR integration"); // maximum time step allown for tsyn integration
    COMM::IOParams<int> sym_order      (input_par_store, -6,   "Symplectic integrator order, should be even number"); // symplectic integrator order
    COMM::IOParams<Float> energy_error (input_par_store, 1e-10,"relative energy error limit for AR"); // phase error requirement
    COMM::IOParams<Float> time_error   (input_par_store, 0.0,  "time synchronization absolute error limit for AR","default is 0.25*dt-min"); // time synchronization error
    COMM::IOParams<Float> time_zero    (input_par_store, 0.0,  "initial physical time");    // initial physical time
    COMM::IOParams<Float> time_end     (input_par_store, 1.0,  "ending physical time"); // ending physical time
    COMM::IOParams<int>   nstep        (input_par_store, 1000, "number of integration steps"); // total step size
    COMM::IOParams<Float> s            (input_par_store, 0.0,  "step size, not physical time step","auto");    // step size
    COMM::IOParams<Float> dt_min       (input_par_store, 1e-13,"minimum physical time step"); // minimum physical time step
    COMM::IOParams<int>   fix_step_option (input_par_store, -1, "always, later, none","auto"); // if true; use input fix step option
    COMM::IOParams<std::string> filename_par (input_par_store, "", "filename to load manager parameters","input name"); // par dumped filename
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
                time_zero.value = atof(optarg);
                break;
            case 2:
                if (strcmp(optarg,"none")) fix_step_option.value=0;
                else if (strcmp(optarg,"always")) fix_step_option.value=1;
                else if (strcmp(optarg,"later")) fix_step_option.value=2;
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
                print_width.value = atof(optarg);
                break;
            case 8:
                print_precision.value = atoi(optarg);
                break;
            default:
                std::cerr<<"Unknown option. check '-h' for help.\n";
                abort();
            }
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
        case 'h':
            std::cout<<"chain [option] data_filename\n"
                     <<"Input data file format: each line: mass, x, y, z, vx, vy, vz\n"
                     <<"Options: (*) show defaulted values\n"
                     <<"          --dt-min          [int]  :  "<<dt_min<<"\n"
                     <<"    -e [Float]:  "<<energy_error<<"\n"
                     <<"          --energy-error    [Float]:  same as -e\n"
                     <<"          --fix-step-option [char] :  "<<fix_step_option<<"\n"
                     <<"    -l :          load dumped data for restart (if used, the input file is dumped data)\n"
                     <<"          --load-data (same as -l)\n"
                     <<"    -k [int]:    "<<sym_order<<"\n"
                     <<"    -n [int]:    "<<nstep<<"\n"
                     <<"    -p [string]: "<<filename_par<<"\n"
                     <<"          --print-width     [int]  : "<<print_width<<"\n"
                     <<"          --print-precision [int]  : "<<print_precision<<"\n"
                     <<"    -s [Float]:  "<<s<<"\n"
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
    manager.time_step_real_min = dt_min.value;
    if (time_error.value>0.0)  manager.time_error_max_real = time_error.value;
    else manager.time_error_max_real = 0.25*dt_min.value;
    manager.energy_error_relative_max = energy_error.value; 
    manager.slowdown_timescale_max = 1.0;
    manager.slowdown_mass_ref = 1.0;
    manager.slowdown_pert_ratio_ref = 1e-6;
    manager.step_count_max = nstep_max.value;
    // set symplectic order
    manager.step.initialSymplecticCofficients(sym_order.value);

    manager.print(std::cerr);

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

    for (int i=0; i<sym_int.particles.getSize(); i++) sym_int.particles[i].id = i+1;

    sym_int.info.reserveMem(sym_int.particles.getSize());
    sym_int.info.generateBinaryTree(sym_int.particles);

    // no initial when both parameters and data are load
    if(!load_flag) {
        // initialization 
        sym_int.initialIntegration(time_zero.value);
        sym_int.info.calcDsAndStepOption(sym_int.slowdown.getSlowDownFactorOrigin(), manager.step.getOrder());
    }

    // use input fix step option
    if (fix_step_option.value>=0) {
        switch (fix_step_option.value) {
        case 0:
            sym_int.info.fix_step_option = FixStepOption::none;
            break;
        case 1:
            sym_int.info.fix_step_option = FixStepOption::always;
            break;
        case 2:
            sym_int.info.fix_step_option = FixStepOption::later;
            break;
        }
    }

    // use input ds
    if (s.value>0.0) sym_int.info.ds = s.value;
    
    // precision
    std::cout<<std::setprecision(print_precision.value);

    //print column title
    sym_int.printColumnTitle(std::cout, print_width.value);
    std::cout<<std::endl;

    //print initial data
    sym_int.printColumn(std::cout, print_width.value);
    std::cout<<std::endl;

    
    // integration loop
    const int n_particle = sym_int.particles.getSize();
    if (time_end.value<0.0) {
        Float time_table[manager.step.getCDPairSize()];
        for (int i=0; i<nstep.value; i++) {
            if(n_particle==2) sym_int.integrateTwoOneStep(sym_int.info.ds, time_table);
            else sym_int.integrateOneStep(sym_int.info.ds, time_table);
            sym_int.printColumn(std::cout, print_width.value);
            std::cout<<std::endl;
        }
    }
    else {
        Float time_step = time_end.value/nstep.value;
        for (int i=1; i<=nstep.value; i++) {
            sym_int.integrateToTime(time_step*i);
            sym_int.printColumn(std::cout, print_width.value);
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

