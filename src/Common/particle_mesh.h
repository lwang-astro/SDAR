#pragma once

#include "Common/Float.h"
#include "particle_group.h"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

namespace COMM{
    class ParticleMesh{
    private:
        struct Cell {
            std::vector<size_t> indices; // particle indexes in this cell
            int group_offset = 0; // offset for group search
        };

        int n_div_[3]; // number of divisions in x,y,z
        Float box_size_[3]; // box sizes in x,y,z
        Float box_center_[3]; // box center in x,y,z
        Float r_cell_max; // maximum cell size
        bool is_built = false; // whether the mesh has been built
        std::vector<std::vector<std::vector<Cell>>> cells_; // 3D grid of cells storing particle indices
        std::vector<std::array<int,3>> cell_indices_; // list of all cell indices
        std::vector<size_t> indices_large_r_search_; // list of all particle with r_search larger than r_cell_max
        int indices_large_r_search_group_offset_; // offset for group search in large r_search list

    public:
        ParticleMesh(): n_div_{0,0,0},
                        box_size_{0.0, 0.0, 0.0}, 
                        box_center_{0.0, 0.0, 0.0}, 
                        r_cell_max(0.0), is_built(false),
                        cells_(), cell_indices_(), 
                        indices_large_r_search_(), indices_large_r_search_group_offset_(0) {}

        /*! Set box size and center from particle group
        */                        
        template <class Tparticle, class Tpcm>
        void setBoxSizeAndCenter(const ParticleGroup<Tparticle, Tpcm>& particles) {
            // determine box size and center from particles
            Float pos_min[3] = {-NUMERIC_FLOAT_MAX, -NUMERIC_FLOAT_MAX, -NUMERIC_FLOAT_MAX};
            Float pos_max[3] = {NUMERIC_FLOAT_MAX, NUMERIC_FLOAT_MAX, NUMERIC_FLOAT_MAX};
            for (int i=0; i<particles.getSize(); i++){
                const auto& pos = particles[i].pos;
                for (int j=0; j<3; j++){
                    if (pos[j] < pos_min[j]) pos_min[j] = pos[j];
                    if (pos[j] > pos_max[j]) pos_max[j] = pos[j];
                }
            }
            for (int i=0; i<3; i++){
                box_size_[i] = pos_max[i] - pos_min[i];
                box_center_[i] = 0.5 * (pos_max[i] + pos_min[i]);
            }
        }

        /*! Set number of divisions in each dimension based on (maximum) cell radius r_cell. 
            When box size in a dimension is smaller than r_cell, set number of division to 1 in that dimension
            @param[in] r_cell: desired (maximum) cell size
            @return total number of cells
        */        
        int setDivision(const Float r_cell) {
            for (int i=0; i<3; i++){
                // check whether one dimension size < r_cell, if so, set box size to r_cell
                if (box_size_[i] < r_cell) {
                    box_size_[i] = r_cell;
                    n_div_[i] = 1;
                }
                else {
                    n_div_[i] = static_cast<int>(std::floor(box_size_[i] / r_cell));
                }
            }
            r_cell_max = r_cell;
            return n_div_[0] * n_div_[1] * n_div_[2];
        }

        /*! find optimized n_div_ for a particle group
            Assume particle has functin getRSearch(), generate a sorted rsearch list of all particles.
            Then test from largest to smallest, find the optimized r_search that can lead to a suitable total cell numbers, 
            where cell numbers > 2 at least in one dimension and cell numbers in all dimensions is much less than total particle numbers
            All particles with rsearch larger than the optimized r_search are stored in indices_large_r_search_ for later processing.
            These number should be small enough to not affect performance much.
            @param[in] particles: particle group to be inserted into the mesh
            @param[in] min_total_cells: minimum total cell numbers required (defaulted 2)
            @param[in] min_particles_per_cell: minimum average particle numbers per cell (defaulted 10)
            @param[in] max_particles_large_r_search_fraction: maximum fraction of particles allowed to have r_search larger than r_cell (defaulted 0.2)
            @return true if optimized division is found, false otherwise
        */
        template <class Tparticle, class Tpcm>
        bool findOptimizedDivision(const ParticleGroup<Tparticle, Tpcm>& particles, 
                                   const int min_total_cells=2, 
                                   const int min_particles_per_cell = 10,
                                   const int max_particles_large_r_search_fraction = 0.2) {
            if (particles.getSize() < min_particles_per_cell) {
                // too less particles, set division to 1 in all dimensions
                n_div_[0] = n_div_[1] = n_div_[2] = 1;
                return;
            }
                                        
            std::vector<Float> rsearch_list;
            rsearch_list.reserve(particles.getSize());
            for (int i=0; i<particles.getSize(); i++){
                search_list.push_back(particles[i].getRSearch());
            }
            // sort by r_search in descending order
            std::sort(rsearch_list.begin(), rsearch_list.end(), std::greater<Float>());

            // find optimized r_search
            size_t step_size = std::max(rsearch_list.size() / min_particles_per_cell, size_t(1));
            int n_steps = static_cast<int>(particles.getSize()*max_particles_large_r_search_fraction);
            for (size_t i=0; i<n_steps; i+=step_size){
                const Float test_r_search = rsearch_list[i];
                const int total_cells = setDivision(test_r_search);
                if ( ( (n_div_[0] > 1) || (n_div_[1] > 1) || (n_div_[2] > 1) ) 
                     && (total_cells <= static_cast<int>(particles.getSize()/min_particles_per_cell)) 
                     && (total_cells >= min_total_cells) ) {
                        // found suitable division
                        return true;
                     }
            }

            return false;
        }

        /*! Build empty cells based on current division settings
        */
        void buildCells() {
            ASSERT(n_div_[0] > 0 && n_div_[1] > 0 && n_div_[2] > 0);
            // initialize 3D vector of cells
            cells_.resize(n_div_[0]);
            for (int i=0; i<n_div_[0]; i++){
                cells_[i].resize(n_div_[1]);
                for (int j=0; j<n_div_[1]; j++){
                    cells_[i][j].resize(n_div_[2]);
                }
            }
        }

        /*! Clear cells and free memory
        */
        void clearCells(){
            // clear 3D vector storage
            cells_.clear();
            cell_indices_.clear();
        }

        /*! Clear particle data in each cell
        */        
        void clearCellData(){
            for (int i=0; i<n_div_[0]; i++){
                for (int j=0; j<n_div_[1]; j++){
                    for (int k=0; k<n_div_[2]; k++){
                        cells_[i][j][k].indices.clear();
                        cells_[i][j][k].group_offset = 0;
                    }
                }
            }
            cell_indices_.clear();
            indices_large_r_search_.clear();
            indices_large_r_search_group_offset_ = 0;
            is_built = false;
        }

        ~ParticleMesh(){
            clearCells();
        }

        /*! Get cell index for a given position
            @param[in] pos: position array
            @param[out] cell_index: cell index array
        */
        void getCellIndex(const Float* pos, std::array<int,3>& cell_index) const {
            for (int i=0; i<3; i++){
                Float rel_pos = pos[i] - (box_center_[i] - 0.5*box_size_[i]);
                int idx = static_cast<int>(std::floor(rel_pos / box_size_[i] * n_div_[i]));
                if (idx < 0) idx = 0;
                if (idx >= n_div_[i]) idx = n_div_[i] - 1;
                cell_index[i] = idx;
            }
        }

        /*! Insert particle into the mesh
            @param[in] particle: particle to be inserted
            @param[in] particle_index: index of the particle to be inserted
            @param[in] is_group: flag to indicate if the particle is a group particle
        */
       template <class Tparticle>
        void insertParticleInOrder(const Tparticle& particle, const size_t particle_index, const bool is_group){
            if (particle.getRSearch() > r_cell_max){
                cell_indices_[particle_index] = {indices_large_r_search_.size(), -1, -1};
                auto it = std::lower_bound(indices_large_r_search_.begin(), indices_large_r_search_.end(), particle_index);
                ASSERT(it == indices_large_r_search_.end() || *it != particle_index); // index must not exist
                indices_large_r_search_.insert(it, particle_index);
                return;
            }
            std::array<int,3> cell_index;
            getCellIndex(particle.pos, cell_index);
            auto& cell = cells_[cell_index[0]][cell_index[1]][cell_index[2]];
            // insert particle index in order
            auto it = std::lower_bound(cell.indices.begin(), cell.indices.end(), particle_index);
            ASSERT(it == cell.indices.end() || *it != particle_index); // index must not exist
            cell.indices.insert(it, particle_index);
            if (!is_group) cell.group_offset ++;
            ASSERT((!is_group && int(it - cell.indices.begin()) < cell.group_offset) || (is_group && int(it - cell.indices.begin()) >= cell.group_offset));
            if (cell_indices_.size() <= particle_index)
                cell_indices_.resize(particle_index + 1);
            cell_indices_[particle_index] = cell_index;
        }
        
        /*! add particles into cells
            @param[in] particles: particle group to be inserted into the mesh
            @param[in] index_group_offset: index offset for group particles
        */
        template <class Tparticle, class Tpcm, TGroup>        
        void addParticles(const ParticleGroup<Tparticle, Tpcm>& particles, const COMM::List<TGroup>& groups, const int index_group_offset){
            ASSERT(!is_built);
            indices_large_r_search_.resize(0);
            indices_large_r_search_group_offset_ = 0;
            cell_indices_.resize(particles.getSize()+groups.getSize());
            for (int i=0; i<particles.getSize(); i++){
                const Float r_search = particles[i].getRSearch();
                if (r_search > r_cell_max){
                    cell_indices_[i] = {indices_large_r_search_.size(), -1, -1};
                    indices_large_r_search_.push_back(i);
                    if (i < index_group_offset) indices_large_r_search_group_offset_ ++;
                    continue;
                }
                const auto& pos = particles[i].pos;
                std::array<int,3> cell_index;
                getCellIndex(pos, cell_index);
                auto& cell = cells_[cell_index[0]][cell_index[1]][cell_index[2]];
                cell.indices.push_back(i);
                if (i < index_group_offset) cell.group_offset ++;
                cell_indices_[i] = cell_index;
            }
            is_built = true;
        }

        /*! Update particle position in the mesh
            @param[in] pos: new position array
            @param[in] particle_index: index of the particle to be updated
            @param[in] is_group: flag to indicate if the particle is a group particle
        */
        void updateParticle(const Float* pos, const size_t particle_index, const bool is_group){
            const auto& old_cell_index = cell_indices_[particle_index];
            if (old_cell_index[1] == -1) {
                // particle was in large r_search list, do nothing
                return;
            }
            std::array<int,3> new_cell_index;
            getCellIndex(pos, new_cell_index);
            if (old_cell_index != new_cell_index){
                // remove from old cell
                auto& old_cell = cells_[old_cell_index[0]][old_cell_index[1]][old_cell_index[2]].indices;
                auto it = std::lower_bound(old_cell.begin(), old_cell.end(), particle_index);
                ASSERT(it != old_cell.end()); // index must exist
                old_cell.erase(it);
                // insert into new cell
                auto& cell = cells_[new_cell_index[0]][new_cell_index[1]][new_cell_index[2]].indices;
                it = std::lower_bound(cell.begin(), cell.end(), particle_index);
                ASSERT(it == cell.end() || *it != particle_index); // index must not exist
                cell.insert(it, particle_index);
                cell_indices_[particle_index] = new_cell_index;
            }
        }

        /*! Search for particles within a given position range
            @param[in] pos_min: minimum position array
            @param[in] pos_max: maximum position array
            @param[out] result: vector to store found particle indices
            @param[out] result_group: vector to store found group particle indices
        */
        void searchPosRange(const Float pos_min[3], const Float pos_max[3], std::vector<size_t>& result, std::vector<size_t>& result_group) const {
            result.clear();
            result_group.clear();
            std::array<int,3> cell_min, cell_max;
            getCellIndex(pos_min, cell_min);
            getCellIndex(pos_max, cell_max);
            for (int i=cell_min[0]; i<=cell_max[0]; i++){ 
                for (int j=cell_min[1]; j<=cell_max[1]; j++){
                    for (int k=cell_min[2]; k<=cell_max[2]; k++){
                        const auto& cell_particles = cells_[i][j][k].indices;
                        result.insert(result.end(), cell_particles.begin(), cell_particles.begin() + cells_[i][j][k].group_offset);
                        result_group.insert(result_group.end(), cell_particles.begin() + cells_[i][j][k].group_offset, cell_particles.end());
                    }
                }
            }
        }

        /*! Search for particles within a given radius from a position
            @param[in] pos: center position array
            @param[in] radius: search radius
            @param[out] result: vector to store found particle indices
            @param[out] result_group: vector to store found group particle indices
        */
        void searchNeighbor(const Float pos[3], const Float radius, std::vector<size_t>& result, std::vector<size_t>& result_group) const {
            Float pos_min[3], pos_max[3];
            for (int i=0; i<3; i++){
                pos_min[i] = pos[i] - radius;
                pos_max[i] = pos[i] + radius;
            }
            searchPosRange(pos_min, pos_max, result, result_group);
        }
    };

} // namespace COMM
