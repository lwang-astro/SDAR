// Test the search neighbor function in ParticleKDTree class
#include <cassert>
#define ASSERT(expr) assert(expr)
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "Common/particle_kdtree.h"
#include "Common/particle_group.h"
#include "Common/io.h"
#include <getopt.h>

class Particle{
public:
    int id;
    Float pos[3];
    Float r_search;

    Float getRNeighbor() const {
        return r_search;
    }
};

class Group{
public:
    COMM::ParticleGroup<Particle,Particle> particles;
    
    Group(): particles() {
        particles.setMode(COMM::ListMode::local);
        particles.reserveMem(1);
    }
};

// Helper for brute force check
Float get_dist_sq(const Float p1[3], const Float p2[3]) {
    Float d2 = 0.0;
    for(int k=0; k<3; k++) d2 += (p1[k]-p2[k])*(p1[k]-p2[k]);
    return d2;
}

int main(int argc, char* argv[]) {

    // input parameters:
    COMM::IOParamsContainer input_par_store;

    COMM::IOParams<int> n_particles(input_par_store, 1000, "number of particles");
    COMM::IOParams<int> n_groups(input_par_store, 100, "number of groups");
    COMM::IOParams<Float> r_cluster(input_par_store, 1.0, "cluster radius");
    COMM::IOParams<Float> r_search_min(input_par_store, 0.01, "minimum search radius");
    COMM::IOParams<Float> r_search_max(input_par_store, 0.1, "maximum search radius");
    
    static int opt_flag = -1;
    static struct option long_options[] = {
        {"n-particles", required_argument, &opt_flag, 'n'},
        {"n-groups", required_argument, &opt_flag, 'g'},
        {"r-cluster", required_argument, &opt_flag, 'r'},
        {"r-search-min", required_argument, &opt_flag, 1},
        {"r-search-max", required_argument, &opt_flag, 2},
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
            std::cout<<"particlekdtreetest [option]\n"
                     <<"Options: (*) show defaulted values\n"
                     <<"          --n-particles(-n) [int]:  "<<n_particles<<"\n"
                     <<"          --n-groups(-g)    [int]:  "<<n_groups<<"\n"
                     <<"          --r-cluster(-r)   [Float: "<<r_cluster<<"\n"
                     <<"          --r-search-min    [Float]: "<<r_search_min<<"\n"
                     <<"          --r-search-max    [Float]: "<<r_search_max<<"\n"
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
             <<"  maximum search radius: "<<r_search_max<<"\n";
    
    // Create a particle group
    COMM::ParticleGroup<Particle,Particle> particles;
    particles.setMode(COMM::ListMode::local);
    particles.reserveMem(n_particles.value+100);
    for (int i = 0; i < n_particles.value; ++i) {
        Particle p;
        p.id = i;
        p.pos[0] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        p.pos[1] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        p.pos[2] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        p.r_search = r_search_min.value + static_cast<Float>(rand()) / RAND_MAX * (r_search_max.value - r_search_min.value); 
        particles.addMember(p);
    }
    
    // create particle indices list (subset)
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
    groups.reserveMem(n_groups.value+100);
    for (int i = 0; i < n_groups.value; ++i) {
        Group g;
        auto& gcm = g.particles.cm;
        gcm.id = i + n_particles.value;
        gcm.pos[0] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        gcm.pos[1] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        gcm.pos[2] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
        gcm.r_search = 2*(r_search_min.value + static_cast<Float>(rand()) / RAND_MAX * (r_search_max.value - r_search_min.value)); 
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
    
    std::cout << "Selected " << particle_indices.getSize() << " particles." << std::endl;
    std::cout << "Selected " << group_indices.getSize() << " groups." << std::endl;

    // ---------------------------------------------------------
    // KDTree Test
    // ---------------------------------------------------------
    std::cout << "Building KDTree..." << std::endl;
    
    COMM::ParticleKDTree kdtree;
    
    // Add particles and groups using the subset indices
    // CHANGED: Pass address of indices
    kdtree.addParticles(particles, &particle_indices);
    kdtree.addGroups(groups, &group_indices);

    std::cout << "KDTree built successfully." << std::endl;

    // ---------------------------------------------------------
    // Verification: 
    // ---------------------------------------------------------
    std::cout << "Verifying all particles and groups in subset..." << std::endl;

    std::vector<int> nb_p_list;
    std::vector<int> nb_g_list;
    int errors = 0;
    int total_checks = 0;

    // 1. Verify Particles in subset
    for (int k = 0; k < particle_indices.getSize(); ++k) {
        int target_idx = particle_indices[k];
        const Particle& target = particles[target_idx];
        total_checks++;

        // KDTree Search
        // CHANGED: Removed target_idx and is_group parameters
        kdtree.searchNeighbor(target, nb_p_list, nb_g_list);

        // Brute Force Search
        std::vector<int> bf_p_list;
        std::vector<int> bf_g_list;

        // Check against selected particles
        for (int j = 0; j < particle_indices.getSize(); ++j) {
            int idx = particle_indices[j];
            const Particle& p = particles[idx];
            // CHANGED: Removed self-check continue, include self in brute force
            
            Float d2 = get_dist_sq(target.pos, p.pos);
            Float r_sum = std::max(target.r_search, p.r_search);
            if (d2 < r_sum * r_sum) {
                bf_p_list.push_back(idx);
            }
        }

        // Check against selected groups
        for (int j = 0; j < group_indices.getSize(); ++j) {
            int idx = group_indices[j];
            const Group& g = groups[idx];
            Float d2 = get_dist_sq(target.pos, g.particles.cm.pos);
            Float r_sum = std::max(target.r_search, g.particles.cm.r_search);
            if (d2 < r_sum * r_sum) {
                bf_g_list.push_back(idx);
            }
        }

        // Compare
        std::sort(nb_p_list.begin(), nb_p_list.end());
        std::sort(bf_p_list.begin(), bf_p_list.end());
        std::sort(nb_g_list.begin(), nb_g_list.end());
        std::sort(bf_g_list.begin(), bf_g_list.end());

        if (nb_p_list != bf_p_list || nb_g_list != bf_g_list) {
            std::cerr << "Error at Particle target ID " << target.id << "\n";
            abort();
        }
    }

    // 2. Verify Groups in subset
    for (int k = 0; k < group_indices.getSize(); ++k) {
        int target_idx = group_indices[k];
        const Group& target = groups[target_idx];
        total_checks++;

        // KDTree Search
        // CHANGED: Removed target_idx and is_group parameters
        kdtree.searchNeighbor(target, nb_p_list, nb_g_list);

        // Brute Force Search
        std::vector<int> bf_p_list;
        std::vector<int> bf_g_list;

        // Check against selected particles
        for (int j = 0; j < particle_indices.getSize(); ++j) {
            int idx = particle_indices[j];
            const Particle& p = particles[idx];
            
            Float d2 = get_dist_sq(target.particles.cm.pos, p.pos);
            Float r_sum = std::max(target.particles.cm.r_search, p.r_search);
            if (d2 < r_sum * r_sum) {
                bf_p_list.push_back(idx);
            }
        }

        // Check against selected groups
        for (int j = 0; j < group_indices.getSize(); ++j) {
            int idx = group_indices[j];
            const Group& g = groups[idx];
            // CHANGED: Removed self-check continue

            Float d2 = get_dist_sq(target.particles.cm.pos, g.particles.cm.pos);
            Float r_sum = std::max(target.particles.cm.r_search, g.particles.cm.r_search);
            if (d2 < r_sum * r_sum) {
                bf_g_list.push_back(idx);
            }
        }

        // Compare
        std::sort(nb_p_list.begin(), nb_p_list.end());
        std::sort(bf_p_list.begin(), bf_p_list.end());
        std::sort(nb_g_list.begin(), nb_g_list.end());
        std::sort(bf_g_list.begin(), bf_g_list.end());

        if (nb_p_list != bf_p_list || nb_g_list != bf_g_list) {
            std::cerr << "Error at Group target ID " << target.particles.cm.id << "\n";
            abort();
        }
    }

    std::cout << "Checked " << total_checks << " targets." << std::endl;

    if (errors == 0) {
        std::cout << "Verification PASSED! All checks matched brute force results." << std::endl;
    } else {
        std::cout << "Verification FAILED with " << errors << " errors." << std::endl;
        return 1;
    }

    // ---------------------------------------------------------
    // Dynamic Update Test
    // ---------------------------------------------------------
    std::cout << "\n---------------------------------------------------------\n";
    std::cout << "Starting Dynamic Update Test (100 steps)..." << std::endl;
    
    int n_steps = 100;
    int update_errors = 0;

    for (int step = 0; step < n_steps; ++step) {
        if (step % 10 == 0) std::cout << "Step " << step << "..." << std::endl;

        // 1. Randomly move some particles
        int n_moved_p = 0;
        for (int k = 0; k < particle_indices.getSize(); ++k) {
            if (rand() % 10 == 0) { // 10% chance to move
                int idx = particle_indices[k];
                Particle& p = particles[idx];
                
                // Random small displacement
                Float dr = p.r_search * 0.2 * (static_cast<Float>(rand()) / RAND_MAX); // 0 to 0.2 * r_search
                Float theta = static_cast<Float>(rand()) / RAND_MAX * 2 * M_PI;
                Float phi = static_cast<Float>(rand()) / RAND_MAX * M_PI;
                
                p.pos[0] += dr * sin(phi) * cos(theta);
                p.pos[1] += dr * sin(phi) * sin(theta);
                p.pos[2] += dr * cos(phi);

                // Update KDTree
                kdtree.updateParticle(idx, p);
                n_moved_p++;
            }
        }

        // 2. Randomly move some groups
        int n_moved_g = 0;
        for (int k = 0; k < group_indices.getSize(); ++k) {
            if (rand() % 10 == 0) { // 10% chance to move
                int idx = group_indices[k];
                Group& g = groups[idx];
                
                Float dr = g.particles.cm.r_search * 0.2 * (static_cast<Float>(rand()) / RAND_MAX);
                Float theta = static_cast<Float>(rand()) / RAND_MAX * 2 * M_PI;
                Float phi = static_cast<Float>(rand()) / RAND_MAX * M_PI;
                
                g.particles.cm.pos[0] += dr * sin(phi) * cos(theta);
                g.particles.cm.pos[1] += dr * sin(phi) * sin(theta);
                g.particles.cm.pos[2] += dr * cos(phi);

                // Update KDTree
                kdtree.updateGroup(idx, g);
                n_moved_g++;
            }
        }

        // 3. Verify Correctness (Full Check)
        // Check ALL particles in subset
        for (int k = 0; k < particle_indices.getSize(); ++k) {
            int target_idx = particle_indices[k];
            const Particle& target = particles[target_idx];

            // CHANGED: Removed target_idx and is_group parameters
            kdtree.searchNeighbor(target, nb_p_list, nb_g_list);

            // Brute Force
            std::vector<int> bf_p_list;
            std::vector<int> bf_g_list;
            for (int j = 0; j < particle_indices.getSize(); ++j) {
                int idx = particle_indices[j];
                const Particle& p = particles[idx];
                // CHANGED: Removed self-check continue
                Float d2 = get_dist_sq(target.pos, p.pos);
                Float r_crit = std::max(target.r_search, p.r_search);
                if (d2 < r_crit * r_crit) bf_p_list.push_back(idx);
            }
            for (int j = 0; j < group_indices.getSize(); ++j) {
                int idx = group_indices[j];
                const Group& g = groups[idx];
                Float d2 = get_dist_sq(target.pos, g.particles.cm.pos);
                Float r_crit = std::max(target.r_search, g.particles.cm.r_search);
                if (d2 < r_crit * r_crit) bf_g_list.push_back(idx);
            }

            std::sort(nb_p_list.begin(), nb_p_list.end());
            std::sort(bf_p_list.begin(), bf_p_list.end());
            std::sort(nb_g_list.begin(), nb_g_list.end());
            std::sort(bf_g_list.begin(), bf_g_list.end());

            if (nb_p_list != bf_p_list || nb_g_list != bf_g_list) {
                std::cerr << "Update Error at Step " << step << ", Particle Target ID " << target.id << "\n";
                update_errors++;
                abort();
            }
        }

        // Check ALL groups in subset
        for (int k = 0; k < group_indices.getSize(); ++k) {
            int target_idx = group_indices[k];
            const Group& target = groups[target_idx];

            // CHANGED: Removed target_idx and is_group parameters
            kdtree.searchNeighbor(target, nb_p_list, nb_g_list);

            // Brute Force
            std::vector<int> bf_p_list;
            std::vector<int> bf_g_list;
            for (int j = 0; j < particle_indices.getSize(); ++j) {
                int idx = particle_indices[j];
                const Particle& p = particles[idx];
                Float d2 = get_dist_sq(target.particles.cm.pos, p.pos);
                Float r_crit = std::max(target.particles.cm.r_search, p.r_search);
                if (d2 < r_crit * r_crit) bf_p_list.push_back(idx);
            }
            for (int j = 0; j < group_indices.getSize(); ++j) {
                int idx = group_indices[j];
                const Group& g = groups[idx];
                // CHANGED: Removed self-check continue
                Float d2 = get_dist_sq(target.particles.cm.pos, g.particles.cm.pos);
                Float r_crit = std::max(target.particles.cm.r_search, g.particles.cm.r_search);
                if (d2 < r_crit * r_crit) bf_g_list.push_back(idx);
            }

            std::sort(nb_p_list.begin(), nb_p_list.end());
            std::sort(bf_p_list.begin(), bf_p_list.end());
            std::sort(nb_g_list.begin(), nb_g_list.end());
            std::sort(bf_g_list.begin(), bf_g_list.end());

            if (nb_p_list != bf_p_list || nb_g_list != bf_g_list) {
                std::cerr << "Update Error at Step " << step << ", Group Target ID " << target.particles.cm.id << "\n";
                update_errors++;
                abort();
            }
        }
    }

    if (update_errors == 0) {
        std::cout << "Dynamic Update Test PASSED! Completed " << n_steps << " steps." << std::endl;
    } else {
        std::cout << "Dynamic Update Test FAILED with " << update_errors << " errors." << std::endl;
    }

    // ---------------------------------------------------------
    // Remove and Insert Test
    // ---------------------------------------------------------
    std::cout << "\n---------------------------------------------------------\n";
    std::cout << "Starting Remove and Insert Test (10 steps)..." << std::endl;

    // Convert Lists to std::vector for easier manipulation in test
    std::vector<int> active_p;
    for(int i=0; i<particle_indices.getSize(); ++i) active_p.push_back(particle_indices[i]);
    
    std::vector<int> active_g;
    for(int i=0; i<group_indices.getSize(); ++i) active_g.push_back(group_indices[i]);

    int ri_errors = 0;
    int n_ri_steps = 10;

    for (int step = 0; step < n_ri_steps; ++step) {
        std::cout << "Step " << step << "..." << std::endl;

        // 1. Remove random particles
        for (int k = 0; k < 5; ++k) {
            if (active_p.empty()) break;
            int rand_idx = rand() % active_p.size();
            int p_idx = active_p[rand_idx];
            
            kdtree.removeParticle(p_idx);
            
            // Remove from active list (swap with back)
            active_p[rand_idx] = active_p.back();
            active_p.pop_back();
        }

        // 2. Remove random groups
        for (int k = 0; k < 5; ++k) {
            if (active_g.empty()) break;
            int rand_idx = rand() % active_g.size();
            int g_idx = active_g[rand_idx];
            
            kdtree.removeGroup(g_idx);
            
            active_g[rand_idx] = active_g.back();
            active_g.pop_back();
        }

        // 3. Insert new particles
        for (int k = 0; k < 5; ++k) {
            Particle p;
            p.id = 100000 + step * 100 + k; // Unique ID
            p.pos[0] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
            p.pos[1] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
            p.pos[2] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
            p.r_search = r_search_min.value + static_cast<Float>(rand()) / RAND_MAX * (r_search_max.value - r_search_min.value);
            
            particles.addMember(p);
            int new_idx = particles.getSize() - 1;
            
            kdtree.InsertParticle(new_idx, p);
            active_p.push_back(new_idx);
        }

        // 4. Insert new groups
        for (int k = 0; k < 5; ++k) {
            Group g;
            g.particles.cm.id = 200000 + step * 100 + k;
            g.particles.cm.pos[0] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
            g.particles.cm.pos[1] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
            g.particles.cm.pos[2] = static_cast<Float>(rand()) / RAND_MAX * r_cluster.value;
            g.particles.cm.r_search = 2*(r_search_min.value + static_cast<Float>(rand()) / RAND_MAX * (r_search_max.value - r_search_min.value));
            
            groups.addMember(g);
            int new_idx = groups.getSize() - 1;
            
            kdtree.InsertGroup(new_idx, g);
            active_g.push_back(new_idx);
        }

        // 5. Verify (Full Check)
        // Check ALL active particles
        for (size_t k = 0; k < active_p.size(); ++k) {
            int target_idx = active_p[k];
            const Particle& target = particles[target_idx];

            // CHANGED: Removed target_idx and is_group parameters
            kdtree.searchNeighbor(target, nb_p_list, nb_g_list);

            // Brute Force
            std::vector<int> bf_p_list;
            std::vector<int> bf_g_list;
            
            for (int idx : active_p) {
                const Particle& p = particles[idx];
                // CHANGED: Removed self-check continue
                Float d2 = get_dist_sq(target.pos, p.pos);
                Float r_crit = std::max(target.r_search, p.r_search);
                if (d2 < r_crit * r_crit) bf_p_list.push_back(idx);
            }
            for (int idx : active_g) {
                const Group& g = groups[idx];
                Float d2 = get_dist_sq(target.pos, g.particles.cm.pos);
                Float r_crit = std::max(target.r_search, g.particles.cm.r_search);
                if (d2 < r_crit * r_crit) bf_g_list.push_back(idx);
            }

            std::sort(nb_p_list.begin(), nb_p_list.end());
            std::sort(bf_p_list.begin(), bf_p_list.end());
            std::sort(nb_g_list.begin(), nb_g_list.end());
            std::sort(bf_g_list.begin(), bf_g_list.end());

            if (nb_p_list != bf_p_list || nb_g_list != bf_g_list) {
                std::cerr << "R/I Error at Step " << step << ", Particle Target ID " << target.id << "\n";
                ri_errors++;
                abort();
            }
        }

        // Check ALL active groups
        for (size_t k = 0; k < active_g.size(); ++k) {
            int target_idx = active_g[k];
            const Group& target = groups[target_idx];

            // CHANGED: Removed target_idx and is_group parameters
            kdtree.searchNeighbor(target, nb_p_list, nb_g_list);

            // Brute Force
            std::vector<int> bf_p_list;
            std::vector<int> bf_g_list;
            
            for (int idx : active_p) {
                const Particle& p = particles[idx];
                Float d2 = get_dist_sq(target.particles.cm.pos, p.pos);
                Float r_crit = std::max(target.particles.cm.r_search, p.r_search);
                if (d2 < r_crit * r_crit) bf_p_list.push_back(idx);
            }
            for (int idx : active_g) {
                const Group& g = groups[idx];
                // CHANGED: Removed self-check continue
                Float d2 = get_dist_sq(target.particles.cm.pos, g.particles.cm.pos);
                Float r_crit = std::max(target.particles.cm.r_search, g.particles.cm.r_search);
                if (d2 < r_crit * r_crit) bf_g_list.push_back(idx);
            }

            std::sort(nb_p_list.begin(), nb_p_list.end());
            std::sort(bf_p_list.begin(), bf_p_list.end());
            std::sort(nb_g_list.begin(), nb_g_list.end());
            std::sort(bf_g_list.begin(), bf_g_list.end());

            if (nb_p_list != bf_p_list || nb_g_list != bf_g_list) {
                std::cerr << "R/I Error at Step " << step << ", Group Target ID " << target.particles.cm.id << "\n";
                ri_errors++;
                abort();
            }
        }
    }

    if (ri_errors == 0) {
        std::cout << "Remove/Insert Test PASSED!" << std::endl;
    } else {
        std::cout << "Remove/Insert Test FAILED with " << ri_errors << " errors." << std::endl;
    }

    return 0;
}