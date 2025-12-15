// Test the search neighbor function in ParticleMesh class
#include <cassert>
#define ASSERT(expr) assert(expr)
#include <vector>
#include <iostream>
#include <cmath>

#include "Common/particle_mesh.h"
#include "Common/particle_group.h"
#include "Common/io.h"
#include <getopt.h>

class Particle{
public:
    int id;
    Float pos[3];
    Float r_search;

    Float getRSearch() const {
        return r_search;
    }
};

class Group{
public:
    Particle cm;
};

int main(int argc, char* argv[]) {

    // input parameters:
    COMM::IOParamsContainer input_par_store;

    COMM::IOParams<int> n_particles(input_par_store, 1000, "number of particles");
    COMM::IOParams<int> n_groups(input_par_store, 100, "number of groups");
    COMM::IOParams<Float> r_cluster(input_par_store, 1.0, "cluster radius");
    COMM::IOParams<Float> r_search_min(input_par_store, 0.01, "minimum search radius");
    COMM::IOParams<Float> r_search_max(input_par_store, 0.1, "maximum search radius");
    COMM::IOParams<int> mesh_n_cells_min(input_par_store, 10, "minimum number of mesh cells");
    COMM::IOParams<int> mesh_n_particles_per_cell_min(input_par_store, 4, "minimum number of particles per cell in mesh");
    COMM::IOParams<Float> mesh_max_particles_large_r_search_fraction(input_par_store, 0.1, "maximum fraction of particles with large search radius in mesh");
    COMM::IOParams<std::string> filename_par (input_par_store, "", "filename to load manager parameters","input name"); // par dumped filename

    static int opt_flag = -1;
    static struct option long_options[] = {
        {"n-particles", required_argument, &opt_flag, 'n'},
        {"n-groups", required_argument, &opt_flag, 'g'},
        {"r-cluster", required_argument, &opt_flag, 'r'},
        {"r-search-min", required_argument, &opt_flag, 1},
        {"r-search-max", required_argument, &opt_flag, 2},
        {"mesh-n-cells-min", required_argument, &opt_flag, 3},
        {"mesh-np-cell-min", required_argument, &opt_flag, 4},
        {"mesh-max-rs-frac", required_argument, &opt_flag, 5},
        {"help",no_argument, &opt_flag, 'h'},
        {0,0,0,0}
    };    
    
    int copt;
    int option_index;
    while ((copt = getopt_long(argc, argv, "n:g:r:h", long_options, &option_index)) != -1) {
        switch (copt) {
        case 'n':
            n_particles.value = atoi(optarg);
            break;
        case 'g':
            n_groups.value = atoi(optarg);
            break;
        case 'r':
            r_cluster.value = atof(optarg);
            break;
        case 'h':
            std::cout<<"particlemesstest [option]\n"
                     <<"Options: (*) show defaulted values\n"
                     <<"          --n-particles(-n) [int]:  "<<n_particles<<"\n"
                     <<"          --n-groups(-g)    [int]:  "<<n_groups<<"\n"
                     <<"          --r-cluster(-r)   [Float]: "<<r_cluster<<"\n"
                     <<"          --r-search-min    [Float]: "<<r_search_min<<"\n"
                     <<"          --r-search-max    [Float]: "<<r_search_max<<"\n"
                     <<"          --mesh-n-cells-min [int]: "<<mesh_n_cells_min<<"\n"
                     <<"          --mesh-np-cell-min [int]: "<<mesh_n_particles_per_cell_min<<"\n"
                     <<"          --mesh-max-rs-frac [Float]: "<<mesh_max_particles_large_r_search_fraction<<"\n"
                     <<"          --help(-h):                 help information\n";
            return 0;
        case 0:
            switch (opt_flag) {
            case 1:
                r_search_min.value = atof(optarg);
                break;
            case 2:
                r_search_max.value = atof(optarg);
                break;
            case 3:
                mesh_n_cells_min.value = atoi(optarg);
                break;
            case 4:
                mesh_n_particles_per_cell_min.value = atoi(optarg);
                break;
            case 5:
                mesh_max_particles_large_r_search_fraction.value = atof(optarg);
                break;
            default:
                std::cerr<<"Unknown argument. check '-h' for help.\n";
                abort();
            }
            break;
        default:
            std::cerr<<"Unknown argument. check '-h' for help.\n";
            abort();
        }
    }

    // print initial parameters
    std::cout<<"Initial parameters:\n"
             <<"  number of particles: "<<n_particles<<"\n"
             <<"  number of groups:    "<<n_groups<<"\n"
             <<"  cluster radius:      "<<r_cluster<<"\n"
             <<"  minimum search radius: "<<r_search_min<<"\n"
             <<"  maximum search radius: "<<r_search_max<<"\n"
             <<"  mesh minimum number of cells: "<<mesh_n_cells_min<<"\n"
             <<"  mesh minimum particles per cell: "<<mesh_n_particles_per_cell_min<<"\n"
             <<"  mesh maximum fraction of particles with large search radius: "<<mesh_max_particles_large_r_search_fraction<<"\n";
    
    // Create a particle group
    COMM::ParticleGroup<Particle,Particle> particles;
    particles.setMode(COMM::ListMode::local);
    particles.reserveMem(n_particles.value);
    for (int i = 0; i < n_particles.value; ++i) {
        Particle p;
        p.id = i;
        p.pos[0] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        p.pos[1] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        p.pos[2] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        p.r_search = r_search_min.value + static_cast<Float>(rand()) / RAND_MAX * (r_search_max.value - r_search_min.value); // random search radius between 5 and 15
        particles.addMember(p);
    }
    // create particle indices list, randomly select the particles with step size from 0 to 2 and keep sorted
    COMM::List<int> particle_indices;
    particle_indices.setMode(COMM::ListMode::local);
    particle_indices.reserveMem(n_particles.value);
    for (int i = 0; i < n_particles.value; ++i) {
        if (rand() % 3 == 0) {
            particle_indices.addMember(i);
        }
    }

    COMM::List<Group> groups;
    groups.setMode(COMM::ListMode::local);
    groups.reserveMem(n_groups.value);
    for (int i = 0; i < n_groups.value; ++i) {
        Group g;
        auto& gcm = g.cm;
        gcm.id = i + n_particles.value;
        gcm.pos[0] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        gcm.pos[1] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        gcm.pos[2] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        gcm.r_search = 2*(r_search_min.value + static_cast<Float>(rand()) / RAND_MAX * (r_search_max.value - r_search_min.value)); // random search radius between 10 and 30
        groups.addMember(g);
    }
    COMM::List<int> group_indices;
    group_indices.setMode(COMM::ListMode::local);
    group_indices.reserveMem(n_groups.value);
    for (int i = 0; i < n_groups.value; ++i) {
        if (rand() % 2 == 0) {
            group_indices.addMember(i);
        }
    }

    // print selected particle and group numbers
    std::cout << "Selected " << particle_indices.getSize() << " particles." << std::endl;
    std::cout << "Selected " << group_indices.getSize() << " groups." << std::endl;
    std::cout << "Total selected: " 
              << particle_indices.getSize() + group_indices.getSize() << std::endl;
    
    // Create a particle mesh
    COMM::ParticleMeshForSearchNeighbor mesh;

    // Set up the mesh with the particle group
    bool use_mesh = mesh.findOptimizedDivision(&particles, &particle_indices, &groups, &group_indices, 
                                               mesh_n_cells_min.value, 
                                               mesh_n_particles_per_cell_min.value, 
                                               mesh_max_particles_large_r_search_fraction.value);
    if (!use_mesh) {
        std::cout << "Mesh not used due to insufficient particles or large search radii." << std::endl;
        int n_div[3];
        mesh.getNDiv(n_div);
        Float box_size[3];
        mesh.getBoxSize(box_size);
        Float box_center[3];
        mesh.getBoxCenter(box_center);
        std::cout << "Box center: (" << box_center[0] << ", " << box_center[1] << ", " << box_center[2] << ")" << std::endl;
        std::cout << "Box size: (" << box_size[0] << ", " << box_size[1] << ", " << box_size[2] << ")" << std::endl;
        std::cout << "Mesh n_div: (" << n_div[0] << ", " << n_div[1] << ", " << n_div[2] << ")" << std::endl;
        return 0;
    }

    mesh.buildCells();    
    mesh.addParticleAndGroups(&particles, &particle_indices, &groups, &group_indices, particles.getSize());

    // Perform the search neighbor operation
    mesh.checkSearchNeighborForAllParticles(&particles, &particle_indices, &groups, &group_indices);

    return 0;
}
