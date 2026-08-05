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
    COMM::IOParams<int>         load_flag;
    COMM::IOParams<std::string> filename_chkpt;

    IOParamsH4()
        : input_par_store()
        , print_width         (input_par_store, WRITE_WIDTH,        "print-width",          "print width of value")
        , print_precision     (input_par_store, WRITE_PRECISION,    "print-precision",      "print digital precision")
        , nstep_max           (input_par_store, 1000000,            "n-step-max",           "number of maximum step for AR integration")
        , sym_order           (input_par_store, -6,                 "k",                    "Symplectic integrator order, should be even number")
        , dt_min_power_index  (input_par_store, 40,                 "dt-min-power",         "power index to calculate mimimum hermite time step: 0.5^n")
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
        , load_flag           (input_par_store, 0,                  "l",                    "Load dumped data for restart (if used, the input file is dumped data)")
        , filename_chkpt      (input_par_store, "",                 "c",                    "filename for the checkpoint at the last output time for restart (binary format)", "<data_filename>.last")
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
            {filename_chkpt.key,           required_argument, &h4_flag, 19},
            {load_flag.key,                no_argument,       &h4_flag, 20},
            {"help",                       no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };

        int opt_used = 0;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "t:k:G:e:o:i:p:c:lh", long_options, &option_index)) != -1)
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
                case 19:
                    filename_chkpt.value = optarg;
                    opt_used += 2;
                    break;
                case 20:
                    load_flag.value = 1;
                    opt_used++;
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
            case 'c':
                filename_chkpt.value = optarg;
                opt_used++;
                break;
            case 'l':
                load_flag.value = 1;
                opt_used++;
                break;
            case 'h':
                std::cout<<bin_name<<" [option] data_filename\n"
                         <<"Input data file format: \n"
                         <<"  First   line:  number of particles(N)\n"
                         <<"  2-(N+1) line:  mass, x, y, z, vx, vy, vz, radius\n"
                         <<"  last    line:  N_group, group_offset_index_lst[N_group], group_member_particle_index[N_member_total]\n"
                         <<"Restart: a binary checkpoint is written at every output time.\n"
                         <<"  Default name: <data_filename>.last;  -c <name> overrides it.\n"
                         <<"  To resume, run with -l <checkpoint> (e.g. -l <data_filename>.last)\n";
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

    // Checkpoint file name: -c overrides; default is <input>.last
    std::string chkpt_filename = iop.filename_chkpt.value.empty()
        ? std::string(filename) + ".last"
        : iop.filename_chkpt.value;
    iop.filename_chkpt.value = chkpt_filename;

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

    // defensive checks
    // 1) keep the AR time synchronization tolerance a small fraction of the minimum block
    //    step (hard limit time_error_max<=0.5*dt_min is enforced in
    //    HermiteIntegrator::checkParams; the default is 0.25*dt_min)
    if (ar_manager.time_error_max > 0.25*manager.step.getDtMin()) {
        std::cerr<<"Warning: time-error ("<<ar_manager.time_error_max
                 <<") > 0.25*dt_min ("<<0.25*manager.step.getDtMin()<<"). "
                 <<"The AR time synchronization tolerance is unusually large relative to the "
                 <<"minimum block step, group times may be inconsistent with the block grid.\n";
    }
    // 2) dt_min must stay resolvable in floating point up to the end of the run. Once it
    //    approaches ~eps*t_end, the block time grid (correctTimeRoundOff) and the AR time
    //    resolution near the end degrade.
    {
        const Float t_end = iop.time_end.value > 1.0 ? iop.time_end.value : 1.0;
        const Float dt_min_resolvable = 8.0*std::numeric_limits<Float>::epsilon()*t_end;
        if (manager.step.getDtMin() < dt_min_resolvable) {
            std::cerr<<"Warning: dt_min ("<<manager.step.getDtMin()
                     <<") < 8*eps*t_end ("<<dt_min_resolvable<<"). "
                     <<"The block time-step grid is no longer resolvable in the current Float "
                     <<"precision by the end of the run (t_end="<<iop.time_end.value<<"). "
                     <<"Consider increasing dt-min, shortening the run, or enabling "
                     <<"-mpfr-digits for higher precision.\n";
        }
    }

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
    // Set effective default filenames BEFORE writing the .par file, so the saved
    // file contains real (non-empty) values for p and c (empty values would
    // misalign the ASCII parameter parsing on restart with -p).
    // filename_chkpt's effective value is already computed above (chkpt_filename).
    if (iop.filename_par.value.empty())
        iop.filename_par.value = std::string(filename) + ".par";
    std::string fpar_out = iop.filename_par.value;
    std::FILE* fout = std::fopen(fpar_out.c_str(),"w");
    if (fout==NULL) {
        std::cerr<<"Error: data file "<<fpar_out<<" cannot be open!\n";
        abort();
    }
    iop.input_par_store.writeAscii(fout);
    fclose(fout);

    // interrupt file output (append mode for restart)
    std::ofstream finterrupt;
    if (iop.interrupt_detection_option.value>0) {
        std::string finterrupt_name = std::string(filename) + ".interrupt";
        auto ios_mode = iop.load_flag.value ? std::ofstream::app : std::ofstream::out;
        finterrupt.open(finterrupt_name.c_str(), ios_mode);
        if (!iop.load_flag.value) {
            AR::InterruptBinary<Particle>::printColumnTitleAscii(finterrupt,20,true);
            finterrupt<<std::endl;
        }
    }

    // integrator
    H4Int h4_int;
    h4_int.manager = &manager;
    h4_int.ar_manager = &ar_manager;

    typedef H4::ParticleH4<Particle> H4Particle;
    COMM::List<H4Particle> restart_particles;
    COMM::List<int> restart_group_offsets;
    COMM::List<int> restart_group_members;
    Float restart_time_int = iop.time_zero.value;  // internal clock to restore
    Float restart_time_offset = 0.0;               // time offset to restore
    int n_ptcl = 0, n_group = 0;

    if (iop.load_flag.value) {
        // === RESTART: read binary checkpoint ===
        std::FILE* fin_bin = std::fopen(filename, "r");
        if (fin_bin == NULL) {
            std::cerr << "Error: data file " << filename << " cannot be open!\n";
            abort();
        }
        uint32_t magic;
        size_t rcount = fread(&magic, sizeof(uint32_t), 1, fin_bin);
        if (rcount < 1 || magic != 0x48434B50U) {
            std::cerr << "Error: Not a valid Hermite checkpoint file (bad magic).\n";
            abort();
        }
        rcount = 0;
        rcount += fread(&n_ptcl, sizeof(int), 1, fin_bin);
        restart_particles.setMode(COMM::ListMode::local);
        restart_particles.reserveMem(n_ptcl);
        restart_particles.resizeNoInitialize(n_ptcl);
        for (int i = 0; i < n_ptcl; i++) {
            size_t rp = fread(&restart_particles[i], sizeof(H4Particle), 1, fin_bin);
            if (rp < 1) {
                std::cerr << "Error: checkpoint particle reading fails!\n";
                abort();
            }
        }

        rcount += fread(&n_group, sizeof(int), 1, fin_bin);
        size_t n_expect = 4;  // n_ptcl, n_group, time_int, time_off
        if (n_group > 0) {
            restart_group_offsets.setMode(COMM::ListMode::local);
            restart_group_offsets.reserveMem(n_group + 1);
            restart_group_offsets.resizeNoInitialize(n_group + 1);
            for (int i = 0; i < n_group + 1; i++)
                rcount += fread(&restart_group_offsets[i], sizeof(int), 1, fin_bin);
            int n_members = restart_group_offsets[n_group];
            restart_group_members.setMode(COMM::ListMode::local);
            restart_group_members.reserveMem(n_members);
            restart_group_members.resizeNoInitialize(n_members);
            for (int i = 0; i < n_members; i++)
                rcount += fread(&restart_group_members[i], sizeof(int), 1, fin_bin);
            n_expect += (size_t)(n_group + 1) + (size_t)n_members;
        }
        // internal clock: real time = restart_time_int + restart_time_offset
        rcount += fread(&restart_time_int, sizeof(Float), 1, fin_bin);
        rcount += fread(&restart_time_offset, sizeof(Float), 1, fin_bin);
        std::fclose(fin_bin);
        if (rcount != n_expect) {
            std::cerr << "Error: Hermite checkpoint file is truncated.\n";
            abort();
        }

        // --- Reconstruct integrator state from checkpoint ---
        h4_int.particles.setMode(COMM::ListMode::local);
        h4_int.particles.reserveMem(n_ptcl);
        h4_int.particles.resizeNoInitialize(n_ptcl);
        for (int i = 0; i < n_ptcl; i++) h4_int.particles[i] = restart_particles[i];

        h4_int.particles.calcCenterOfMass();

        if (iop.n_neighbor_max.value <= 0) manager.n_neighbor_max = n_ptcl;
        else manager.n_neighbor_max = iop.n_neighbor_max.value;

        Float m_ave = h4_int.particles.cm.mass / n_ptcl;
        Float r_neighbor_sum = 0.0;
        for (int i = 0; i < n_ptcl; i++)
            r_neighbor_sum += h4_int.particles[i].getRNeighbor();
        Float r_neighbor_ave = r_neighbor_sum / n_ptcl;
        manager.step.calcAcc0OffsetSq(m_ave, r_neighbor_ave, iop.grav_const.value);
        h4_int.step = manager.step;

#ifdef SLOWDOWN_MASSRATIO
        if (iop.slowdown_mass_ref.value <= 0.0) ar_manager.slowdown_mass_ref = m_ave;
        else ar_manager.slowdown_mass_ref = iop.slowdown_mass_ref.value;
#endif

        h4_int.groups.setMode(COMM::ListMode::local);
        h4_int.groups.reserveMem(n_ptcl);
        h4_int.reserveIntegratorMem();

        // Standard init flow: needed to build groups, neighbors, forces
        // Restore the internal clock and its offset so that real time = time_ + time_offset_
        // matches the checkpoint and the (small, shifted) internal particle times are
        // consistent with time_. AR groups pick up time_offset_ inside initialIntegration.
        h4_int.setTimeOffset(restart_time_offset);
        h4_int.initialSystemSingle(restart_time_int);
        if (n_group > 0)
            h4_int.addGroups(restart_group_members.getDataAddress(),
                             restart_group_offsets.getDataAddress(), n_group);
        h4_int.initialIntegration();

        // initialIntegration zeroes acc0/acc1/dt for "init" particles -
        // restore the saved values from the checkpoint so the predictor
        // step in the next integration uses the correct derivatives.
        // NOTE: AR integrator state inside groups is rebuilt from scratch;
        // this means restart is NOT bit-exact, but results are within
        // the integration accuracy (~1e-7 relative error).
        for (int i = 0; i < n_ptcl; i++) {
            h4_int.particles[i].acc0[0] = restart_particles[i].acc0[0];
            h4_int.particles[i].acc0[1] = restart_particles[i].acc0[1];
            h4_int.particles[i].acc0[2] = restart_particles[i].acc0[2];
            h4_int.particles[i].acc1[0] = restart_particles[i].acc1[0];
            h4_int.particles[i].acc1[1] = restart_particles[i].acc1[1];
            h4_int.particles[i].acc1[2] = restart_particles[i].acc1[2];
            h4_int.particles[i].pot    = restart_particles[i].pot;
            h4_int.particles[i].time   = restart_particles[i].time;
            h4_int.particles[i].dt     = restart_particles[i].dt;
        }
    } else {
        // === NORMAL START: read ASCII input ===
        std::fstream fin;
        fin.open(filename, std::fstream::in);
        if (!fin.is_open()) {
            std::cerr << "Error: data file " << filename << " cannot be open!\n";
            abort();
        }
        h4_int.particles.setMode(COMM::ListMode::local);
        h4_int.particles.readMemberAscii(fin);
        for (int i = 0; i < h4_int.particles.getSize(); i++) h4_int.particles[i].id = i + 1;
        h4_int.particles.calcCenterOfMass();
        h4_int.particles.shiftToCenterOfMassFrame();
        h4_int.particles.calcCenterOfMass();

        if (iop.n_neighbor_max.value <= 0) manager.n_neighbor_max = h4_int.particles.getSize();
        else manager.n_neighbor_max = iop.n_neighbor_max.value;

        Float m_ave = h4_int.particles.cm.mass / h4_int.particles.getSize();
        Float r_neighbor_sum = 0.0;
        for (int i = 0; i < h4_int.particles.getSize(); i++) {
            h4_int.particles[i].setRGroupAndNeighbor(iop.r_group.value, iop.r_neighbor_over_group.value, m_ave);
            r_neighbor_sum += h4_int.particles[i].getRNeighbor();
        }
        Float r_neighbor_ave = r_neighbor_sum / h4_int.particles.getSize();
        manager.step.calcAcc0OffsetSq(m_ave, r_neighbor_ave, iop.grav_const.value);
        h4_int.step = manager.step;

#ifdef SLOWDOWN_MASSRATIO
        if (iop.slowdown_mass_ref.value <= 0.0) ar_manager.slowdown_mass_ref = m_ave;
        else ar_manager.slowdown_mass_ref = iop.slowdown_mass_ref.value;
#endif

        h4_int.groups.setMode(COMM::ListMode::local);
        h4_int.groups.reserveMem(h4_int.particles.getSize());
        h4_int.reserveIntegratorMem();
        h4_int.initialSystemSingle(iop.time_zero.value);
        h4_int.readGroupConfigureAscii(fin);
        h4_int.initialIntegration();
    }

    // --- Common post-init ---
    manager.print(std::cerr);
    ar_manager.print(std::cerr);
    std::cerr << "CM: after shift ";
    h4_int.particles.cm.printColumnAscii(std::cerr, 22);
    std::cerr << std::endl;

    // Group slowdown tracking (for output column formatting)
    const int n_group_init = h4_int.getNGroup();
    int n_group_sub_tot_init = 0;
    COMM::List<int> n_group_sub_init_lst;
    n_group_sub_init_lst.setMode(COMM::ListMode::local);
    n_group_sub_init_lst.reserveMem(n_group_init);
    n_group_sub_init_lst.resizeNoInitialize(n_group_init);
    for (int i = 0; i < n_group_init; i++) {
#ifdef AR_SLOWDOWN_ARRAY
        n_group_sub_init_lst[i] = h4_int.groups[i].binary_slowdown.getSize();
#elif AR_SLOWDOWN_TREE
        n_group_sub_init_lst[i] = h4_int.groups[i].info.binarytree.getSize();
#else
        n_group_sub_init_lst[i] = 0;
#endif
        n_group_sub_tot_init += n_group_sub_init_lst[i];
    }

    bool use_adjust_true = !iop.load_flag.value;
    h4_int.adjustGroups(use_adjust_true);
    h4_int.initialIntegration();
    h4_int.sortDtAndSelectActParticle();

    // precision
    std::cout << std::setprecision(iop.print_precision.value);

    if (!iop.load_flag.value) {
        // get initial energy (fresh start only)
        h4_int.calcEnergySlowDown(true);
        h4_int.particles.calcCenterOfMass();
        std::cerr << "CM:";
        h4_int.particles.cm.printColumnAscii(std::cerr, 22);
        std::cerr << std::endl;

        const int* sd_arr = n_group_sub_init_lst.getDataAddress();
        h4_int.printColumnTitleAscii(std::cout, iop.print_width.value, sd_arr, n_group_init, n_group_sub_tot_init);
        std::cout << std::endl;
        h4_int.printColumnAscii(std::cout, iop.print_width.value, sd_arr, n_group_init, n_group_sub_tot_init);
        std::cout << std::endl;
    }

    // dt_out and first output time
    Float dt_out = pow(Float(0.5), Float(iop.dt_out_power_index.value));
    Float time_out;
    if (iop.load_flag.value) {
        time_out = iop.time_zero.value + dt_out;
        while (time_out <= h4_int.getTime() + Float(1e-15)) time_out += dt_out;
    } else {
        time_out = iop.time_zero.value + dt_out;
    }

    const int* sd_arr = n_group_sub_init_lst.getDataAddress();

    // BINARY CHECKPOINT WRITER
    auto writeCheckpoint = [&](const char* chkpt_path) {
        std::FILE* fchkpt = std::fopen(chkpt_path, "w");
        if (fchkpt == NULL) return;
        uint32_t magic = 0x48434B50U;
        fwrite(&magic, sizeof(uint32_t), 1, fchkpt);
        int n_p = h4_int.particles.getSize();
        fwrite(&n_p, sizeof(int), 1, fchkpt);
        h4_int.writeBackGroupMembers();
        for (int i = 0; i < n_p; i++) fwrite(&h4_int.particles[i], sizeof(H4Particle), 1, fchkpt);
        int n_g = h4_int.getNGroup();
        fwrite(&n_g, sizeof(int), 1, fchkpt);
        if (n_g > 0) {
            COMM::List<int> offsets, members;
            offsets.setMode(COMM::ListMode::local);
            offsets.reserveMem(n_g + 1);
            offsets.resizeNoInitialize(n_g + 1);
            offsets[0] = 0;
            int tot = 0;
            for (int i = 0; i < n_g; i++) {
                int nm = h4_int.groups[i].particles.getSize();
                tot += nm;
                offsets[i + 1] = tot;
            }
            members.setMode(COMM::ListMode::local);
            members.reserveMem(tot);
            members.resizeNoInitialize(tot);
            int idx = 0;
            for (int g = 0; g < n_g; g++) {
                int nm = h4_int.groups[g].particles.getSize();
                for (int j = 0; j < nm; j++)
                    members[idx++] = h4_int.groups[g].info.particle_index[j];
            }
            for (int i = 0; i < n_g + 1; i++) fwrite(&offsets[i], sizeof(int), 1, fchkpt);
            for (int i = 0; i < tot; i++) fwrite(&members[i], sizeof(int), 1, fchkpt);
        }
        // Save the internal clock state so restart can reconstruct both time_ and
        // time_offset_ (real time = time_ + time_offset_). This is required after the
        // periodic time-origin shift keeps the internal clock small while time_offset_
        // grows; otherwise a restart would set time_ to the (large) real time while the
        // particle .time fields are the (small) internal times.
        Float time_int = h4_int.getTimeInt();
        Float time_off = h4_int.getTime() - time_int;
        fwrite(&time_int, sizeof(Float), 1, fchkpt);
        fwrite(&time_off, sizeof(Float), 1, fchkpt);
        std::fclose(fchkpt);
    };

    // integration loop
    while (h4_int.getTime() < iop.time_end.value) {
        // periodically re-anchor the internal clock to keep |time_| small (reset to
        // [0, dt_max) once it reaches a dt_max multiple). The real physical time
        // (h4_int.getTime()) is unchanged, but the internal clock stays far from the
        // machine round-off limit even for very long runs. The shift amount is a
        // multiple of dt_max (hence of dt_min), so the block-time-step grid stays aligned.
        {
            const Float dt_max_shift = manager.step.getDtMax();
            const Float time_int = h4_int.getTimeInt();
            if (time_int >= dt_max_shift) {
                const Float shift = dt_max_shift * floor(time_int/dt_max_shift);
                h4_int.shiftTimeOrigin(shift);
            }
        }
        h4_int.integrateGroupsOneStep();
        int n_interrupt = h4_int.getNInterrupt();
        for (int i = 0; i < n_interrupt; i++) {
            auto& interrupt_info = h4_int.getInterruptInfo(i);
            std::cerr << "Interrupt " << i << " : ";
            switch (interrupt_info.status) {
            case AR::InterruptStatus::change:
                std::cerr << " Change";
                break;
            case AR::InterruptStatus::merge:
                std::cerr << " Merge";
                break;
            case AR::InterruptStatus::destroy:
                std::cerr << " Destroy";
                break;
            case AR::InterruptStatus::none:
                break;
            }
            std::cerr << std::endl;
            interrupt_info.printColumnTitleAscii(std::cerr);
            std::cerr << std::endl;
            interrupt_info.printColumnAscii(std::cerr);
            std::cerr << std::endl;
            if (iop.interrupt_detection_option.value > 0) {
                interrupt_info.printColumnAscii(finterrupt, 20, true);
                finterrupt << std::endl;
            }
        }
        h4_int.integrateSingleOneStepAct();
        h4_int.adjustGroups(false);
        h4_int.initialIntegration();
        h4_int.modifySingleParticles();
        h4_int.sortDtAndSelectActParticle();

        if (h4_int.getTime() >= time_out) {
            h4_int.calcEnergySlowDown(false);
            h4_int.particles.calcCenterOfMass();
            std::cerr << "CM:";
            h4_int.particles.cm.printColumnAscii(std::cerr, 22);
            std::cerr << std::endl;

            h4_int.printColumnAscii(std::cout, iop.print_width.value, sd_arr, n_group_init, n_group_sub_tot_init);
            std::cout << std::endl;
            h4_int.printStepHist();

            // Write restart checkpoint at every output time (always enabled)
            writeCheckpoint(chkpt_filename.c_str());

            time_out += dt_out;
        }
    }

    //fpu_fix_end(&oldcw);

    return 0;
}

