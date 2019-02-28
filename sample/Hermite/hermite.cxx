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
    int dt_min_power_index = 40; // power index to calculate minimum physical time step
    Float energy_error=1e-10; // phase error requirement
    Float time_error=0;; // time synchronization error
    Float time_zero=0.0;    // initial physical time
    Float time_end=1.0; // ending physical time
    Float dt_output = 0.25; // output time interval
    Float dt_max = 0.25; // maximum physical time step
    Float r_break = 1e-3; // binary break criterion
    Float r_search = 5.0; // neighbor search radius for AR
    Float eta_4th = 0.1; // time step coefficient 
    Float eta_2nd = 0.001; // time step coefficient for 2nd order
    Float eps_sq = 0.0;    // softening parameter
    Float slowdown_ref =1e-6; // slowdown reference factor
    Float slowdown_max =1e4; // slowdown reference factor
    Float G = 1.0;      // gravitational constant
    char* filename_par=NULL; // par dumped filename

    int copt;
    static struct option long_options[] = {
        {"time-start", required_argument, 0, 0},
        {"time-end", required_argument, 0, 't'},
        {"r-break", required_argument, 0, 'r'},
        {"energy-error",required_argument, 0, 'e'},
        {"time-error",required_argument, 0, 0},
        {"dt-max",required_argument, 0, 0},
        {"dt-min-power",required_argument, 0, 0},
        {"n-step-max",required_argument, 0, 0},
        {"eta-4th",required_argument, 0, 0},
        {"eta-2nd",required_argument, 0, 0},
        {"eps",required_argument, 0, 0},
        {"slowdown-ref",required_argument, 0, 0},
        {"slowdown-max",required_argument, 0, 0},
        {"print-width",required_argument, 0, 0},
        {"print-precision",required_argument, 0, 0},
        {"load-par",required_argument, 0, 0},    
        {"help",no_argument, 0, 'h'},
        {0,0,0,0}
    };
  
    int option_index;
    while ((copt = getopt_long(argc, argv, "t:R:r:k:G:e:o:lh", long_options, &option_index)) != -1)
        switch (copt) {
        case 0:
#ifdef DEBUG
            std::cerr<<"option "<<option_index<<" "<<long_options[option_index].name<<" "<<optarg<<std::endl;
#endif
            switch (option_index) {
            case 0:
                time_zero = atof(optarg);
                break;
            case 4:
                time_error = atof(optarg);
                break;
            case 5:
                dt_max = atof(optarg);
                break;
            case 6:
                dt_min_power_index = atoi(optarg);
                break;
            case 7:
                nstep_max = atoi(optarg);
                break;
            case 8:
                eta_4th = atof(optarg);
                break;
            case 9:
                eta_2nd = atof(optarg);
                break;
            case 10:
                eps_sq = atof(optarg);
                break;
            case 11:
                slowdown_ref = atof(optarg);
                break;
            case 12:
                slowdown_max = atof(optarg);
                break;
            case 13:
                print_width = atof(optarg);
                break;
            case 14:
                print_precision = atoi(optarg);
                break;
            case 15:
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
        case 'r':
            r_break = atof(optarg);
            break;
        case 'R':
            r_search = atof(optarg);
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
        case 'h':
            std::cout<<"chain [option] data_filename\n"
                     <<"Input data file format: each line: mass, x, y, z, vx, vy, vz\n"
                     <<"Options: (*) show defaulted values\n"
                     <<"          --time-start [Float]:  initial physical time ("<<time_zero<<")\n"
                     <<"    -t [Float]:  ending physical time ("<<time_end<<")\n"
                     <<"          --time-end (same as -t)\n"
                     <<"          --dt-max       [Float]: maximum hermite time step ("<<dt_max<<")\n"
                     <<"          --dt-min-power [int]  : power index to calculate mimimum hermite time step ("<<dt_min_power_index<<")\n"
                     <<"    -r [Float]:  distance criterion for switching AR and Hermite ("<<r_break<<")\n"
                     <<"          --r-break             : same as -r\n"
                     <<"    -R [Float]:  neighbor search radius ("<<r_search<<")\n"
                     <<"    -o [Float]:  output time interval ("<<dt_output<<")\n"
                     <<"    -k [int]:  Symplectic integrator order, should be even number ("<<sym_order<<")\n"
                     <<"    -e [Float]:  relative energy error limit for AR ("<<energy_error<<")\n"
                     <<"          --energy-error (same as -e)\n"
                     <<"          --time-error [Float]:    time synchronization absolute error limit for AR, default is 0.25*dt-min ("<<time_error<<")\n"
                     <<"          --n-step-max [int]  :    number of maximum step for AR integration ("<<nstep_max<<")\n"
                     <<"    -G [Float]: gravitational constant ("<<G<<")\n"
                     <<"          --eta-4th:   time step coefficient for 4th order ("<<eta_4th<<")\n"
                     <<"          --eta-2nd:   time step coefficient for 2nd order ("<<eta_2nd<<")\n"
                     <<"          --eps:       softerning parameter ("<<eps_sq<<")\n"
                     <<"          --slowdown-ref:  slowdown perturbation ratio reference ("<<slowdown_ref<<")\n"
                     <<"          --slowdown-max:  slowdown maximum factor ("<<slowdown_max<<")\n"
                     <<"          --print-width [int]:     print width of value ("<<print_width<<")\n"
                     <<"          --print-precision [int]: print digital precision ("<<print_precision<<")\n"
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


    // manager
    HermiteManager<HermiteInteraction> manager;
    AR::SymplecticManager<ARInteraction> ar_manager;

    if (filename_par!=NULL) {
        std::FILE* fp = std::fopen(filename_par,"r");
        if (fp==NULL) {
            std::cerr<<"Error: parameter file "<<filename_par<<" cannot be open!\n";
            abort();
        }
        manager.readBinary(fp);
        ar_manager.readBinary(fp);
        fclose(fp);
    }
    else {
        manager.r_break_crit = r_break;
        manager.r_neighbor_crit = r_search;
        manager.step.eta_4th = eta_4th;
        manager.step.eta_2nd = eta_2nd;
        manager.step.setDtRange(dt_max, dt_min_power_index);
        manager.interaction.eps_sq = eps_sq;
        manager.interaction.G = G;
        ar_manager.interaction.eps_sq = eps_sq;
        ar_manager.interaction.G = G;
        ar_manager.time_step_real_min = manager.step.getDtMin();
        if (time_error == 0.0) time_error = 0.25*ar_manager.time_step_real_min;
        ASSERT(time_error>1e-14);
        ar_manager.time_error_max_real = time_error;
        // time error cannot be smaller than round-off error
        ar_manager.energy_error_relative_max = energy_error; 
        ar_manager.slowdown_pert_ratio_ref = slowdown_ref;
        ar_manager.slowdown_factor_max = slowdown_max;
        ar_manager.step_count_max = nstep_max;
        // set symplectic order
        ar_manager.step.initialSymplecticCofficients(sym_order);
    }

    // integrator
    HermiteIntegrator<Particle, Particle, HermitePerturber, Neighbor<Particle>, HermiteInteraction, ARInteraction, HermiteInformation> h4_int;
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
        
    Float m_ave = h4_int.particles.cm.mass/h4_int.particles.getSize();
    manager.step.calcAcc0OffsetSq(m_ave, r_search);
    ar_manager.slowdown_mass_ref = m_ave;

    std::cerr<<"CM: after shift ";
    h4_int.particles.cm.printColumn(std::cerr, print_width);
    std::cerr<<std::endl;

    h4_int.groups.setMode(COMM::ListMode::local);
    h4_int.groups.reserveMem(h4_int.particles.getSize());
    h4_int.reserveIntegratorMem();
    // initial system 
    h4_int.initialSystemSingle(time_zero);
    h4_int.readGroupConfigureAscii(fin);

    // no initial when both parameters and data are load
    // initialization 
    h4_int.initialIntegration(); // get neighbors and min particles
    const int n_group_init = h4_int.getNGroup();
    h4_int.adjustGroups(true);
    h4_int.initialIntegration();
    h4_int.sortDtAndSelectActParticle();
    h4_int.info.time = h4_int.getTime();

    // precision
    std::cout<<std::setprecision(print_precision);

    // get initial energy
    h4_int.writeBackGroupMembers();
    h4_int.info.calcEnergy(h4_int.particles, manager.interaction, true);
    // cm
    h4_int.particles.calcCenterOfMass();
    std::cerr<<"CM:";
    h4_int.particles.cm.printColumn(std::cerr, print_width);
    std::cerr<<std::endl;

    //print column title
    h4_int.info.printColumnTitle(std::cout, print_width);
    std::cout<<std::setw(print_width)<<"Ngroup";
    for (int i=0; i<n_group_init; i++) h4_int.groups[i].slowdown.printColumnTitle(std::cout, print_width);
    h4_int.particles.printColumnTitle(std::cout, print_width);
    std::cout<<std::endl;

    //print initial data
    h4_int.info.printColumn(std::cout, print_width);
    h4_int.particles.printColumn(std::cout, print_width);
    std::cout<<std::endl;
    
    // integration loop
    while (h4_int.info.time<time_end) {
        h4_int.integrateOneStepAct();
        h4_int.adjustGroups(false);
        h4_int.initialIntegration();
        h4_int.sortDtAndSelectActParticle();
        h4_int.info.time = h4_int.getTime();

        if (dt_output==0.0||fmod(h4_int.info.time, dt_output)==0.0) {
            h4_int.writeBackGroupMembers();
            h4_int.info.calcEnergy(h4_int.particles, manager.interaction, false);
            
            h4_int.particles.calcCenterOfMass();
            std::cerr<<"CM:";
            h4_int.particles.cm.printColumn(std::cerr, print_width);
            std::cerr<<std::endl;

            h4_int.info.printColumn(std::cout, print_width);
            std::cout<<std::setw(print_width)<<n_group_init;
            for (int i=0; i<n_group_init; i++) h4_int.groups[i].slowdown.printColumn(std::cout, print_width);
            h4_int.particles.printColumn(std::cout, print_width);
            std::cout<<std::endl;
            h4_int.printStepHist();
        }
    }

    std::string fpar_out = std::string(filename) + ".par";
    std::FILE* fout = std::fopen(fpar_out.c_str(),"w");
    if (fout==NULL) {
        std::cerr<<"Error: data file "<<fpar_out<<" cannot be open!\n";
        abort();
    }
    manager.writeBinary(fout);
    ar_manager.writeBinary(fout);
    fclose(fout);

    //fpu_fix_end(&oldcw);

    return 0;
}

