#pragma once

#include "Common/Float.h"
#include "particle_group.h"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

namespace COMM{
    class ParticleMeshForSearchNeighbor{
    private:
        struct Cell {
            std::vector<size_t> indices; // particle indexes in this cell
            int group_offset = 0; // offset for group search
        };

        int n_div_[3]; // number of divisions in x,y,z
        Float box_size_[3]; // box sizes in x,y,z
        Float box_center_[3]; // box center in x,y,z
        Float r_cell_max; // maximum cell size
        std::vector<std::vector<std::vector<Cell>>> cells_; // 3D grid of cells storing particle indices
        std::vector<std::array<int,3>> cell_indices_; // list of all cell indices
        std::vector<size_t> particle_indices_large_r_search_; // list of all particle indices with r_search larger than r_cell_max
        int particle_indices_large_r_search_group_offset_; // offset for groups with large r_search 
        bool is_box_set = false; // whether the box size and center have been set
        bool is_n_div_set = false; // whether the number of divisions have been set
        bool is_cells_initialized = false; // whether the cells_ structure has been initialized
        bool is_cells_built = false; // whether the particles have been added to the cells

    public:
        ParticleMeshForSearchNeighbor(): n_div_{0,0,0},
                                         box_size_{0.0, 0.0, 0.0}, 
                                         box_center_{0.0, 0.0, 0.0}, 
                                         r_cell_max(0.0), 
                                         cells_(), cell_indices_(), 
                                         particle_indices_large_r_search_(), 
                                         particle_indices_large_r_search_group_offset_(0),
                                         is_box_set(false), is_n_div_set(false), is_cells_initialized(false),
                                         is_cells_built(false) {}

        /*! Set box size and center from particle group
            @param[in] particles: particle group to be inserted into the mesh
            @param[in] particle_indices: optional particle index list to be considered for box size and center calculation
            @param[in] groups: optional group list to be considered for box size and center calculation
            @param[in] group_indices: optional group index list to be considered for box size and center calculation
        */                        
        template <class Tparticle, class Tpcm, class TGroup>
        void setBoxSizeAndCenter(const ParticleGroup<Tparticle, Tpcm>* particles, 
                                 const COMM::List<int>* particle_indices = nullptr,
                                 const COMM::List<TGroup>* groups = nullptr,
                                 const COMM::List<int>* group_indices = nullptr) {
            // determine box size and center from particles
            Float pos_min[3] = {-NUMERIC_FLOAT_MAX, -NUMERIC_FLOAT_MAX, -NUMERIC_FLOAT_MAX};
            Float pos_max[3] = {NUMERIC_FLOAT_MAX, NUMERIC_FLOAT_MAX, NUMERIC_FLOAT_MAX};
            if (particle_indices != nullptr) {
                for (int i=0; i<particle_indices->getSize(); i++){
                    const int idx = (*particle_indices)[i];
                    const auto& pos = (*particles)[idx].pos;
                    for (int j=0; j<3; j++){
                        if (pos[j] < pos_min[j]) pos_min[j] = pos[j];
                        if (pos[j] > pos_max[j]) pos_max[j] = pos[j];
                    }
                }
            }
            else {
                for (int i=0; i<particles->getSize(); i++){
                    const auto& pos = (*particles)[i].pos;
                    for (int j=0; j<3; j++){
                        if (pos[j] < pos_min[j]) pos_min[j] = pos[j];
                        if (pos[j] > pos_max[j]) pos_max[j] = pos[j];
                    }
                }
            }

            if (groups != nullptr) {
                if (group_indices != nullptr) {
                    for (int i=0; i<group_indices->getSize(); i++){
                        const int idx = (*group_indices)[i];
                        const auto& pos = (*groups)[idx].cm.pos;
                        for (int j=0; j<3; j++){
                            if (pos[j] < pos_min[j]) pos_min[j] = pos[j];
                            if (pos[j] > pos_max[j]) pos_max[j] = pos[j];
                        }
                    }
                }
                else {
                    for (int i=0; i<groups->getSize(); i++){
                        const auto& pos = (*groups)[i].cm.pos;
                        for (int j=0; j<3; j++){
                            if (pos[j] < pos_min[j]) pos_min[j] = pos[j];
                            if (pos[j] > pos_max[j]) pos_max[j] = pos[j];
                        }
                    }
                }
            }

            for (int i=0; i<3; i++){
                box_size_[i] = pos_max[i] - pos_min[i];
                box_center_[i] = 0.5 * (pos_max[i] + pos_min[i]);
            }

            is_box_set = true;
        }

        /*! Set number of divisions in each dimension based on (maximum) cell radius r_cell. 
            When box size in a dimension is smaller than r_cell, set number of division to 1 in that dimension
            @param[in] r_cell: desired (maximum) cell size
            @return total number of cells
        */        
        int setDivision(const Float r_cell) {
            ASSERT(is_box_set);

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

            is_n_div_set = true;
        }

        /*! find optimized n_div_ for a particle group
            Assume particle has functin getRSearch(), generate a sorted rsearch list of all particles.
            Then test from largest to smallest, find the optimized r_search that can lead to a suitable total cell numbers, 
            where cell numbers > 2 at least in one dimension and cell numbers in all dimensions is much less than total particle numbers
            All particles with rsearch larger than the optimized r_search are stored in indices_large_r_search_ for later processing.
            These number should be small enough to not affect performance much.
            @param[in] particles: particle group to be inserted into the mesh
            @param[in] n_cells_min: minimum total cell numbers required (defaulted 2)
            @param[in] n_particles_per_cell_min: minimum average particle numbers per cell (defaulted 10)
            @param[in] max_particles_large_r_search_fraction: maximum fraction of particles allowed to have r_search larger than r_cell (defaulted 0.2)
            @return true if optimized division is found, false otherwise
        */
        template <class Tparticle, class Tpcm, class TGroup>
        bool findOptimizedDivision(const ParticleGroup<Tparticle, Tpcm>* particles, 
                                   const COMM::List<int>* particle_indices = nullptr,
                                   const COMM::List<TGroup>* groups = nullptr,
                                   const COMM::List<int>* group_indices = nullptr,
                                   const int n_cells_min=2, 
                                   const int n_particles_per_cell_min=10,
                                   const int max_particles_large_r_search_fraction=0.2) {

            int n_particles = particle_indices != nullptr ? particle_indices->getSize() : 0;
            int n_groups = group_indices != nullptr ? group_indices->getSize() : (groups != nullptr ? groups->getSize() : 0);
            int n_tot = n_particles + n_groups;
            if (n_tot < n_particles_per_cell_min) {
                // too less particles, set division to 1 in all dimensions
                n_div_[0] = n_div_[1] = n_div_[2] = 1;
                return false;
            }

            // first set box size and center            
            setBoxSizeAndCenter(particles, particle_indices, groups, group_indices);
                                        
            // check r_search of all particles
            std::vector<Float> rsearch_list;
            rsearch_list.reserve(n_tot);
            if (particle_indices != nullptr) {
                for (int i=0; i<particle_indices->getSize(); i++){
                    const int idx = (*particle_indices)[i];
                    rsearch_list.push_back((*particles)[idx].getRSearch());
                }
            }
            else {
                for (int i=0; i<particles->getSize(); i++){
                    rsearch_list.push_back((*particles)[i].getRSearch());
                }
            }
            if (groups != nullptr) {
                if (group_indices != nullptr) {
                    for (int i=0; i<group_indices->getSize(); i++){
                        const int idx = (*group_indices)[i];
                        rsearch_list.push_back((*groups)[idx].cm.getRSearch());
                    }
                }
                else {
                    for (int i=0; i<groups->getSize(); i++){
                        rsearch_list.push_back((*groups)[i].cm.getRSearch());
                    }
                }
            }
            // sort rsearch in descending order
            std::sort(rsearch_list.begin(), rsearch_list.end(), std::greater<Float>());

            // find optimized r_search
            size_t step_size = std::max(rsearch_list.size() / n_particles_per_cell_min, size_t(1));
            int n_steps = static_cast<int>(particles.getSize()*max_particles_large_r_search_fraction);
            for (size_t i=0; i<n_steps; i+=step_size){
                const Float test_r_search = rsearch_list[i];
                const int total_cells = setDivision(test_r_search);
                if ( ( (n_div_[0] > 1) || (n_div_[1] > 1) || (n_div_[2] > 1) ) 
                     && (total_cells <= static_cast<int>(particles.getSize()/n_particles_per_cell_min)) 
                     && (total_cells >= n_cells_min) ) {
                        // found suitable division
                        return true;
                     }
            }

            return false;
        }

        /*! Build empty cells based on current division settings
        */
        void buildCells() {
            ASSERT(is_n_div_set);

            // initialize 3D vector of cells
            cells_.resize(n_div_[0]);
            for (int i=0; i<n_div_[0]; i++){
                cells_[i].resize(n_div_[1]);
                for (int j=0; j<n_div_[1]; j++){
                    cells_[i][j].resize(n_div_[2]);
                }
            }

            is_cells_initialized = true;
        }

        /*! Clear cells and free memory
        */
        void clearCells(){
            // clear 3D vector storage
            cells_.clear();
            cell_indices_.clear();

            is_cells_initialized = false;
        }

        /*! Clear particle data in each cell
        */        
        void clearCellData(){
            ASSERT(is_cells_initialized);

            for (int i=0; i<n_div_[0]; i++){
                for (int j=0; j<n_div_[1]; j++){
                    for (int k=0; k<n_div_[2]; k++){
                        cells_[i][j][k].indices.clear();
                        cells_[i][j][k].group_offset = 0;
                    }
                }
            }
            cell_indices_.clear();
            particle_indices_large_r_search_.clear();
            particle_indices_large_r_search_group_offset_ = 0;
            is_cells_built = false;
        }

        ~ParticleMeshForSearchNeighbor(){
            clearCells();
        }

        /*! Get cell index for a given position
            @param[in] pos: position array
            @param[out] cell_index: cell index array
        */
        void getCellIndex(const Float* pos, std::array<int,3>& cell_index) const {
            ASSERT(is_box_set);
            for (int i=0; i<3; i++){
                Float rel_pos = pos[i] - (box_center_[i] - 0.5*box_size_[i]);
                int idx = static_cast<int>(std::floor(rel_pos / box_size_[i] * n_div_[i]));
                if (idx < 0) idx = 0;
                if (idx >= n_div_[i]) idx = n_div_[i] - 1;
                cell_index[i] = idx;
            }
        }

        /*! add particles into cells
            @param[in] particles: particle group to be inserted into the mesh
            @param[in] particle_indices: optional particle index list to be considered for insertion
            @param[in] groups: optional group list to be considered for insertion
            @param[in] group_indices: optional group index list to be considered for insertion
            @param[in] index_group_offset: index offset for group particles
        */
        template <class Tparticle, class Tpcm, class TGroup>
        void addParticleAndGroups(const ParticleGroup<Tparticle, Tpcm>& particles, 
                                  const COMM::List<int>* particle_indices = nullptr,
                                  const COMM::List<TGroup>* groups = nullptr,
                                  const COMM::List<int>* group_indices = nullptr,
                                  const int index_group_offset){
            ASSERT(!is_cells_built);
            particle_indices_large_r_search_.resize(0);
            particle_indices_large_r_search_group_offset_ = 0;
            ASSERT(index_group_offset >= particles.getSize());
            cell_indices_.resize(index_group_offset + groups.getSize());
            for (int i=0; i< cell_indices_.size(); i++){
                cell_indices_[i] = {-1, -1, -1};
            }

            if (particle_indices) {
                for (int i=0; i<particle_indices->getSize(); i++){
                    const int idx = (*particle_indices)[i];
                    const Float r_search = particles[idx].getRSearch();
                    if (r_search > r_cell_max){
                        cell_indices_[idx] = {particle_indices_large_r_search_.size(), -1, -1};
                        particle_indices_large_r_search_.push_back(idx);
                        if (idx < index_group_offset) particle_indices_large_r_search_group_offset_ ++;
                        continue;
                    }
                    const auto& pos = particles[idx].pos;
                    std::array<int,3> cell_index;
                    getCellIndex(pos, cell_index);
                    auto& cell = cells_[cell_index[0]][cell_index[1]][cell_index[2]];
                    cell.indices.push_back(idx);
                    if (idx < index_group_offset) cell.group_offset ++;
                    cell_indices_[idx] = cell_index;
                }
            }
            else {
                for (int i=0; i<particles.getSize(); i++){
                    const Float r_search = particles[i].getRSearch();
                    if (r_search > r_cell_max){
                        cell_indices_[i] = {particle_indices_large_r_search_.size(), -1, -1};
                        particle_indices_large_r_search_.push_back(i);
                        if (i < index_group_offset) particle_indices_large_r_search_group_offset_ ++;
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
            }

            if (groups) {
                if (group_indices) {
                    for (int i=0; i<group_indices->getSize(); i++){
                        const int idx = (*group_indices)[i];
                        const Float r_search = (*groups)[idx].cm.getRSearch();
                        if (r_search > r_cell_max){
                            cell_indices_[idx + index_group_offset] = {particle_indices_large_r_search_.size(), -1, -1};
                            particle_indices_large_r_search_.push_back(idx + index_group_offset);
                            continue;
                        }
                        const auto& pos = (*groups)[idx].cm.pos;
                        std::array<int,3> cell_index;
                        getCellIndex(pos, cell_index);
                        auto& cell = cells_[cell_index[0]][cell_index[1]][cell_index[2]];
                        cell.indices.push_back(idx + index_group_offset);
                        cell_indices_[idx + index_group_offset] = cell_index;
                    }
                }
                else {
                    for (int i=0; i<groups->getSize(); i++){
                        const Float r_search = (*groups)[i].cm.getRSearch();
                        if (r_search > r_cell_max){
                            cell_indices_[i + index_group_offset] = {particle_indices_large_r_search_.size(), -1, -1};
                            particle_indices_large_r_search_.push_back(i + index_group_offset);
                            continue;
                        }
                        const auto& pos = (*groups)[i].cm.pos;
                        std::array<int,3> cell_index;
                        getCellIndex(pos, cell_index);
                        auto& cell = cells_[cell_index[0]][cell_index[1]][cell_index[2]];
                        cell.indices.push_back(i + index_group_offset);
                        cell_indices_[i + index_group_offset] = cell_index;
                    }
                }
            }

            is_cells_built = true;
        }

        bool isCellsBuilt() const {
            return is_cells_built;
        }

        /*! Insert particle into the mesh
            @param[in] particle: particle to be inserted
            @param[in] particle_index: index of the particle to be inserted
            @param[in] is_group: flag to indicate if the particle is a group particle
        */
        template <class Tparticle>
        void insertParticleInOrder(const Tparticle& particle, const size_t particle_index, const bool is_group){
            if (particle.getRSearch() > r_cell_max){
                cell_indices_[particle_index] = {particle_indices_large_r_search_.size(), -1, -1};
                auto it = std::lower_bound(particle_indices_large_r_search_.begin(), particle_indices_large_r_search_.end(), particle_index);
                ASSERT(it == particle_indices_large_r_search_.end() || *it != particle_index); // index must not exist
                particle_indices_large_r_search_.insert(it, particle_index);
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

        /*! Remove particle
            @param[in] particle_index: index of the particle to be removed
            @param[in] is_group: flag to indicate if the particle is a group particle
        */
        void removeParticle(const size_t particle_index, const bool is_group){
            const auto& cell_index = cell_indices_[particle_index];
            ASSERT(cell_index[0] != -1); // particle must exist
            if (cell_index[1] == -1) {
                // particle was in large r_search list
                auto it = std::lower_bound(particle_indices_large_r_search_.begin(), particle_indices_large_r_search_.end(), particle_index);
                ASSERT(it != particle_indices_large_r_search_.end()); // index must exist
                particle_indices_large_r_search_.erase(it);
                return;
            }
            auto& cell = cells_[cell_index[0]][cell_index[1]][cell_index[2]].indices;
            auto it = std::lower_bound(cell.begin(), cell.end(), particle_index);
            ASSERT(it != cell.end()); // index must exist
            if (!is_group) cells_[cell_index[0]][cell_index[1]][cell_index[2]].group_offset --;
            ASSERT((is_group && (static_cast<int>(it - cell.begin()) < cells_[cell_index[0]][cell_index[1]][cell_index[2]].group_offset)) || (!is_group && (static_cast<int>(it - cell.begin()) >= cells_[cell_index[0]][cell_index[1]][cell_index[2]].group_offset)));
            cell.erase(it);
        }

        /*! Update particle position in the mesh, assume no r_search change
            @param[in] pos: new position array
            @param[in] particle_index: index of the particle to be updated
            @param[in] is_group: flag to indicate if the particle is a group particle
        */
        void updateParticle(const Float* pos, const size_t particle_index, const bool is_group){
            const auto& old_cell_index = cell_indices_[particle_index];
            ASSERT(old_cell_index[0] != -1); // particle must exist
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
