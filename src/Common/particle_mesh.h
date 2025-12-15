#pragma once

#include "Common/Float.h"
#include "particle_group.h"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace COMM{
    class ParticleMeshForSearchNeighbor{
    private:
        struct Cell {
            std::vector<int> indices; // particle indexes in this cell
            int group_offset = 0; // offset for group search
        };

        int n_div_[3]; // number of divisions in x,y,z
        Float box_size_[3]; // box sizes in x,y,z
        Float box_center_[3]; // box center in x,y,z
        Float r_cell_max; // maximum cell size
        std::vector<std::vector<std::vector<Cell>>> cells_; // 3D grid of cells storing particle indices
        std::vector<std::array<int,3>> cell_indices_; // list of all cell indices
        std::vector<int> particle_indices_large_r_search_; // list of all particle indices with r_search larger than r_cell_max
        int particle_indices_large_r_search_group_offset_; // offset for groups with large r_search 
        int index_group_offset_; // index offset for group particles
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
            Float pos_min[3] = {NUMERIC_FLOAT_MAX, NUMERIC_FLOAT_MAX, NUMERIC_FLOAT_MAX};
            Float pos_max[3] = {-NUMERIC_FLOAT_MAX, -NUMERIC_FLOAT_MAX, -NUMERIC_FLOAT_MAX};
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
            is_n_div_set = true;

            return n_div_[0] * n_div_[1] * n_div_[2];
        }

        /*! find optimized n_div_ for a particle group
            Assume particle has functin getRSearch(), generate a sorted rsearch list of all particles.
            Then test from largest to smallest, find the optimized r_search that can lead to a suitable total cell numbers, 
            where cell numbers > 2 at least in one dimension and cell numbers in all dimensions is much less than total particle numbers
            All particles with rsearch larger than the optimized r_search are stored in indices_large_r_search_ for later processing.
            These number should be small enough to not affect performance much.
            @param[in] particles: particle group to be inserted into the mesh
            @param[in] particle_indices: optional particle index list to be considered for finding optimized division
            @param[in] groups: optional group list to be considered for finding optimized division
            @param[in] group_indices: optional group index list to be considered for finding optimized division
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
                                   const int n_cells_min=10, 
                                   const int n_particles_per_cell_min=4,
                                   const Float max_particles_large_r_search_fraction=0.2) {

            int n_particles = particle_indices != nullptr ? particle_indices->getSize() : 0;
            int n_groups = group_indices != nullptr ? group_indices->getSize() : (groups != nullptr ? groups->getSize() : 0);
            int n_tot = n_particles + n_groups;
            if (n_tot < n_cells_min * n_particles_per_cell_min) {
                // too less particles, set division to 1 in all dimensions
#ifdef PARTICLE_MESH_DEBUG
                std::cout << "ParticleMeshForSearchNeighbor: total selected particles (" << n_tot 
                          << ") less than minimum required (" << n_cells_min * n_particles_per_cell_min 
                          << "), mesh not used." << std::endl; 
#endif                
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

            // first check the largest rsearch
            int total_cells = setDivision(rsearch_list[0]);
            int n_cells_mean = n_tot/n_particles_per_cell_min;
            // cells should be less than n_tot
            if (total_cells > n_tot) {
                Float test_r_search = rsearch_list[0]*std::cbrt(Float(total_cells)/Float(n_tot));
                total_cells = setDivision(test_r_search);
#ifdef PARTICLE_MESH_DEBUG
                std::cout << "ParticleMeshForSearchNeighbor: adjusted r_cell = " << test_r_search
                          << ", n_div = (" << n_div_[0] << ", " << n_div_[1] << ", " << n_div_[2] << ")"
                          << ", total_cells = " << total_cells << std::endl;
#endif
                return true;
            }
            else if (total_cells < n_cells_mean) {
                // too less cells, find optimized r_search
                const int n_steps = static_cast<int>(n_tot*max_particles_large_r_search_fraction);
                for (int i=1; i<n_steps; ++i){
                    const Float test_r_search = rsearch_list[i]; 
                    const int total_cells = setDivision(test_r_search);
                    if (total_cells >= n_cells_mean) {
                        // found suitable division
    #ifdef PARTICLE_MESH_DEBUG
                        std::cout << "ParticleMeshForSearchNeighbor: test r_search = " << test_r_search
                                << ", n_div = (" << n_div_[0] << ", " << n_div_[1] << ", " << n_div_[2] << ")"
                                << ", total_cells = " << total_cells << std::endl; 
    #endif                
                        return true;
                    }
                }
                if (total_cells < n_cells_min) {
                // still too less cells, do not use mesh    
#ifdef PARTICLE_MESH_DEBUG
                    std::cout << "ParticleMeshForSearchNeighbor: total cells (" << total_cells 
                              << ") less than minimum required (" << n_cells_min 
                              << "), mesh not used." << std::endl;  
#endif            
                    return false;
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
            index_group_offset_ = 0;
            is_cells_built = false;
        }

        // clear function
        void clear() {
            n_div_[0] = n_div_[1] = n_div_[2] = 0;
            box_size_[0] = box_size_[1] = box_size_[2] = 0.0;
            box_center_[0] = box_center_[1] = box_center_[2] = 0.0;
            r_cell_max = 0.0;
            clearCells();
            particle_indices_large_r_search_.clear();
            particle_indices_large_r_search_group_offset_ = 0;
            index_group_offset_ = 0;
            is_box_set = false;
            is_n_div_set = false;
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
        void addParticleAndGroups(const ParticleGroup<Tparticle, Tpcm>* particles, 
                                  const COMM::List<int>* particle_indices = nullptr,
                                  const COMM::List<TGroup>* groups = nullptr,
                                  const COMM::List<int>* group_indices = nullptr,
                                  const int index_group_offset = 0){
            ASSERT(!is_cells_built);
            particle_indices_large_r_search_.resize(0);
            particle_indices_large_r_search_group_offset_ = 0;
            ASSERT(index_group_offset >= particles->getSize());
            cell_indices_.resize(index_group_offset + groups->getSize());
            for (size_t i=0; i< cell_indices_.size(); i++){
                cell_indices_[i] = {-1, -1, -1};
            }

            if (particle_indices) {
#ifdef PARTICLE_MESH_DEBUG
                // check whether particle indices are sorted and unique
                for (int i=1; i<particle_indices->getSize(); i++){
                    ASSERT( (*particle_indices)[i] > (*particle_indices)[i-1] );
                }
#endif                
                for (int i=0; i<particle_indices->getSize(); i++){
                    const int idx = (*particle_indices)[i];
                    const Float r_search = (*particles)[idx].getRSearch();
                    // when r_search > r_cell_max, put in large r_search list
                    if (r_search > r_cell_max){
                        cell_indices_[idx] = {static_cast<int>(particle_indices_large_r_search_.size()), -1, -1};
                        particle_indices_large_r_search_.push_back(idx);
                        if (idx < index_group_offset) particle_indices_large_r_search_group_offset_ ++;
                        continue;
                    }
                    const auto& pos = (*particles)[idx].pos;
                    std::array<int,3> cell_index;
                    getCellIndex(pos, cell_index);
                    auto& cell = cells_[cell_index[0]][cell_index[1]][cell_index[2]];
                    cell.indices.push_back(idx);
                    if (idx < index_group_offset) cell.group_offset ++;
                    cell_indices_[idx] = cell_index;
                }
            }
            else {
                for (int i=0; i<particles->getSize(); i++){
                    const Float r_search = (*particles)[i].getRSearch();
                    if (r_search > r_cell_max){
                        cell_indices_[i] = {static_cast<int>(particle_indices_large_r_search_.size()), -1, -1};
                        particle_indices_large_r_search_.push_back(i);
                        if (i < index_group_offset) particle_indices_large_r_search_group_offset_ ++;
                        continue;
                    }
                    const auto& pos = (*particles)[i].pos;
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
#ifdef PARTICLE_MESH_DEBUG
                    // check whether group indices are sorted and unique
                    for (int i=1; i<group_indices->getSize(); i++){
                        ASSERT( (*group_indices)[i] > (*group_indices)[i-1] );
                    }                
#endif
                    for (int i=0; i<group_indices->getSize(); i++){
                        const int idx = (*group_indices)[i];
                        const Float r_search = (*groups)[idx].cm.getRSearch();
                        if (r_search > r_cell_max){
                            cell_indices_[idx + index_group_offset] = {static_cast<int>(particle_indices_large_r_search_.size()), -1, -1};
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
                            cell_indices_[i + index_group_offset] = {static_cast<int>(particle_indices_large_r_search_.size()), -1, -1};
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
#ifdef PARTICLE_MESH_DEBUG
            // print indices_large_r_search_ info
            std::cout << "ParticleMeshForSearchNeighbor: number of particles with r_search > r_cell (" 
                      << r_cell_max << ") = " << particle_indices_large_r_search_.size() 
                      << " out of total " << ( (particle_indices != nullptr ? particle_indices->getSize() : particles->getSize()) 
                                              + (group_indices != nullptr ? group_indices->getSize() : (groups != nullptr ? groups->getSize() : 0)) )
                      << ", fraction = " << static_cast<Float>(particle_indices_large_r_search_.size()) 
                                         / static_cast<Float>( (particle_indices != nullptr ? particle_indices->getSize() : particles->getSize()) 
                                                              + (group_indices != nullptr ? group_indices->getSize() : (groups != nullptr ? groups->getSize() : 0)) )
                      << std::endl;
#endif

            index_group_offset_ = index_group_offset;
            is_cells_built = true;
        }

        bool isCellsBuilt() const {
            return is_cells_built;
        }

        /*! Insert particle into the mesh
            @param[in] particle: particle to be inserted
            @param[in] particle_index: index of the particle to be inserted, if it is group index, will be added with index_group_offset_
            @param[in] is_group: flag to indicate if the particle is a group particle
        */
        template <class Tparticle>
        void insertParticleInOrder(const Tparticle& particle, const int particle_index, const bool is_group){
            int index = particle_index;
            if (is_group) index += index_group_offset_;

            if (particle.getRSearch() > r_cell_max){
                cell_indices_[index] = {particle_indices_large_r_search_.size(), -1, -1};
                auto it = std::lower_bound(particle_indices_large_r_search_.begin(), particle_indices_large_r_search_.end(), index);
                ASSERT(it == particle_indices_large_r_search_.end() || *it != index); // index must not exist
                particle_indices_large_r_search_.insert(it, index);
                return;
            }
            std::array<int,3> cell_index;
            getCellIndex(particle.pos, cell_index);
            auto& cell = cells_[cell_index[0]][cell_index[1]][cell_index[2]];
            // insert particle index in order
            auto it = std::lower_bound(cell.indices.begin(), cell.indices.end(), index);
            ASSERT(it == cell.indices.end() || *it != index); // index must not exist
            cell.indices.insert(it, index);
            if (!is_group) cell.group_offset ++;
            ASSERT((!is_group && int(it - cell.indices.begin()) < cell.group_offset) || (is_group && int(it - cell.indices.begin()) >= cell.group_offset));
            if (cell_indices_.size() <= index)
                cell_indices_.resize(index + 1);
            cell_indices_[index] = cell_index;
        }

        /*! Remove particle
            @param[in] particle_index: index of the particle to be removed, if group index, will be added with index_group_offset_
            @param[in] is_group: flag to indicate if the particle is a group particle
        */
        void removeParticle(const int particle_index, const bool is_group){
            int index = particle_index;
            if (is_group) index += index_group_offset_; 
            const auto& cell_index = cell_indices_[index];
            ASSERT(cell_index[0] != -1); // particle must exist
            if (cell_index[1] == -1) {
                // particle was in large r_search list
                auto it = std::lower_bound(particle_indices_large_r_search_.begin(), particle_indices_large_r_search_.end(), index);
                ASSERT(it != particle_indices_large_r_search_.end()); // index must exist
                particle_indices_large_r_search_.erase(it);
                return;
            }
            auto& cell = cells_[cell_index[0]][cell_index[1]][cell_index[2]].indices;
            auto it = std::lower_bound(cell.begin(), cell.end(), index);
            ASSERT(it != cell.end()); // index must exist
            if (!is_group) cells_[cell_index[0]][cell_index[1]][cell_index[2]].group_offset --;
            ASSERT((is_group && (static_cast<int>(it - cell.begin()) < cells_[cell_index[0]][cell_index[1]][cell_index[2]].group_offset)) || (!is_group && (static_cast<int>(it - cell.begin()) >= cells_[cell_index[0]][cell_index[1]][cell_index[2]].group_offset)));
            cell.erase(it);
        }

        /*! Update particle position in the mesh, assume no r_search change
            @param[in] pos: new position array
            @param[in] particle_index: index of the particle to be updated, if group index, will be added with index_group_offset_
            @param[in] is_group: flag to indicate if the particle is a group particle
        */
        void updateParticle(const Float* pos, const int particle_index, const bool is_group){
            int index = particle_index;
            if (is_group) index += index_group_offset_;
            const auto& old_cell_index = cell_indices_[index];
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
                auto it = std::lower_bound(old_cell.begin(), old_cell.end(), index);
                ASSERT(it != old_cell.end()); // index must exist
                old_cell.erase(it);
                // insert into new cell
                auto& cell = cells_[new_cell_index[0]][new_cell_index[1]][new_cell_index[2]].indices;
                it = std::lower_bound(cell.begin(), cell.end(), index);
                ASSERT(it == cell.end() || *it != index); // index must not exist
                cell.insert(it, index);
                cell_indices_[index] = new_cell_index;
            }
        }


        /*! Search for particles within a given radius from a position
            @param[in] particle: target particle
            @param[out] result: vector to store found particle indices
            @param[out] result_group: vector to store found group particle indices
        */
        template <class Tparticle>
        void searchNeighbor(const Tparticle& particle, std::vector<int>& result, std::vector<int>& result_group) const {
            result.clear();
            result_group.clear();

            auto &pos = particle.pos;
            std::array<int,3> cell, cell_min, cell_max;            
            getCellIndex(pos, cell);             

            Float &rsearch = particle.getRSearch();
            // if particle has large r_search, determine search cell range accordingly, otherwise only search neighboring cells
            int dn_cells = (rsearch <= r_cell_max) ? 1 : static_cast<int>(std::ceil(rsearch / r_cell_max));
            for (int i=0; i<3; i++){
                cell_min[i] = std::max(cell[i]-dn_cells, 0);
                cell_max[i] = std::min(cell[i]+dn_cells, n_div_[i]-1);
            }

            for (int i=cell_min[0]; i<=cell_max[0]; i++){
                for (int j=cell_min[1]; j<=cell_max[1]; j++){
                    for (int k=cell_min[2]; k<=cell_max[2]; k++){
                        const auto& cell_particles = cells_[i][j][k].indices;
                        result.insert(result.end(), cell_particles.begin(), cell_particles.begin() + cells_[i][j][k].group_offset);
                        result_group.insert(result_group.end(), cell_particles.begin() + cells_[i][j][k].group_offset, cell_particles.end());
                     }
                }
            }

            // add particles in large r_search list
            result.insert(result.end(), particle_indices_large_r_search_.begin(), particle_indices_large_r_search_.begin() + particle_indices_large_r_search_group_offset_);
            result_group.insert(result_group.end(), particle_indices_large_r_search_.begin() + particle_indices_large_r_search_group_offset_, particle_indices_large_r_search_.end());

            // correct group indices
            for (size_t i=0; i<result_group.size(); i++) result_group[i] -= index_group_offset_;
        }


        //! get n division
        void getNDiv(int n_div[3]) const {
            for (int i=0; i<3; i++){
                n_div[i] = n_div_[i];
            }
        }

        //! get box sizes
        void getBoxSize(Float box_size[3]) const {
            for (int i=0; i<3; i++){
                box_size[i] = box_size_[i];
            }
        }

        //! get box center
        void getBoxCenter(Float box_center[3]) const {
            for (int i=0; i<3; i++){
                box_center[i] = box_center_[i];
            }
        }

        //! for debug, check whether the search neighbor return correct results
        template <class Tparticle, class Tpcm, class TGroup>
        Float checkSearchNeighborForOneParticle(const int particle_index,
                                              const bool is_group,
                                              const ParticleGroup<Tparticle, Tpcm>* particles, 
                                              const COMM::List<int>* particle_indices = nullptr,
                                              const COMM::List<TGroup>* groups = nullptr,
                                              const COMM::List<int>* group_indices = nullptr) const {

            ASSERT(is_cells_built); 
            std::vector<int> result, result_group;
            const Float *pos;
            Float rsearch;
            if (is_group) {
                const auto& particle = (*groups)[particle_index].cm;
                pos = particle.pos;
                rsearch = particle.getRSearch();
                searchNeighbor(particle, particle_index, is_group, result, result_group);
            }
            else {
                const auto& particle = (*particles)[particle_index];
                pos = particle.pos;
                rsearch = particle.getRSearch();
                searchNeighbor(particle, particle_index, is_group, result, result_group);
            }

            int n_found = result.size() + result_group.size();

            // Check if the found particles are correct
            std::vector<int> result_check, result_group_check;

            if (particle_indices != nullptr) {
                for (int i=0; i<particle_indices->getSize(); i++){
                    const int idx = (*particle_indices)[i];
                    if (!is_group && idx == particle_index) continue;
                    const auto& p = (*particles)[idx];
                    Float dist_sq = 0.0;
                    for (int j=0; j<3; j++){
                        Float diff = p.pos[j] - pos[j];
                        dist_sq += diff * diff;
                    }
                    if (std::sqrt(dist_sq) < std::max(rsearch, p.getRSearch())){
                        result_check.push_back(idx);                    
                    }
                }
            }
            else {
                for (int i=0; i<particles->getSize(); i++){
                    if (!is_group && i == particle_index) continue;
                    const auto& p = (*particles)[i];
                    Float dist_sq = 0.0;
                    for (int j=0; j<3; j++){
                        Float diff = p.pos[j] - pos[j];
                        dist_sq += diff * diff;
                    }
                    if (std::sqrt(dist_sq) < std::max(rsearch, p.getRSearch())){
                        result_check.push_back(i);                    
                    }
                }
            }
            // compare result and result_check
            std::sort(result_check.begin(), result_check.end());
            // first generate new result list with distance of two particles less than r_search
            std::vector<int> result_filtered;
            for (size_t i = 0; i < result.size(); i++) {
                const int idx = result[i];
                const auto& p = (*particles)[idx];
                Float dist_sq = 0.0;
                for (int j=0; j<3; j++){
                    Float diff = p.pos[j] - pos[j];
                    dist_sq += diff * diff;
                }
                if (std::sqrt(dist_sq) < std::max(rsearch, p.getRSearch())){
                    result_filtered.push_back(idx);
                }
            }            
            std::sort(result_filtered.begin(), result_filtered.end());
            // check filtered result with result_check, first check size
            if (result_filtered.size() != result_check.size()) {
                std::cerr<< "ParticleMeshForSearchNeighbor::checkSearchNeighborForOneParticle(): size mismatch for particle index "<<particle_index
                         << ", result size = "<<result_filtered.size()
                         << ", result_check size = "<<result_check.size()<<std::endl;
                abort();
            }

            // then check each index
            for (size_t j = 0; j < result_check.size(); j++) {
                if (result_filtered[j] != result_check[j]) {
                    std::cerr<< "ParticleMeshForSearchNeighbor::checkSearchNeighborForOneParticle(): index mismatch for particle index "<<particle_index
                             << ", result["<<j<<"] = "<<result_filtered[j]
                             << ", result_check["<<j<<"] = "<<result_check[j]<<std::endl;
                    abort();
                }
            }

            // now check groups, first get group list with distance less than r_search
            std::vector<int> result_group_filtered;
            for (size_t i = 0; i < result_group.size(); i++) {
                const int idx = result_group[i];
                const auto& g = (*groups)[idx];
                Float dist_sq = 0.0;
                for (int j=0; j<3; j++){
                    Float diff = g.cm.pos[j] - pos[j];
                    dist_sq += diff * diff;
                }
                if (std::sqrt(dist_sq) < std::max(rsearch, g.cm.getRSearch())){
                    result_group_filtered.push_back(idx);
                }
            }

            if (group_indices != nullptr) {
                for (int i=0; i<group_indices->getSize(); i++){
                    const int idx = (*group_indices)[i];
                    if (is_group && idx == particle_index) continue;
                    const auto& g = (*groups)[idx];
                    Float dist_sq = 0.0;
                    for (int j=0; j<3; j++){
                        Float diff = g.cm.pos[j] - pos[j];
                        dist_sq += diff * diff;
                    }
                    if (std::sqrt(dist_sq) < std::max(rsearch, g.cm.getRSearch())){
                        result_group_check.push_back(idx);                    
                    }
                }
            }
            else {
                for (int i=0; i<groups->getSize(); i++){
                    const auto& g = (*groups)[i];
                    if (is_group && i == particle_index) continue;
                    Float dist_sq = 0.0;
                    for (int j=0; j<3; j++){
                        Float diff = g.cm.pos[j] - pos[j];
                        dist_sq += diff * diff;
                    }
                    if (std::sqrt(dist_sq) < std::max(rsearch, g.cm.getRSearch())){
                        result_group_check.push_back(i);                    
                    }
                }
            }

            // compare result_group_filtered and result_group_check
            std::sort(result_group_check.begin(), result_group_check.end());
            // first check size
            if (result_group_filtered.size() != result_group_check.size()) {
                std::cerr<< "ParticleMeshForSearchNeighbor::checkSearchNeighborForOneParticle(): size mismatch for group index "<<particle_index
                         << ", result_group size = "<<result_group_filtered.size()
                         << ", result_group_check size = "<<result_group_check.size()<<std::endl;
                abort();
            }
            // then check each index
            std::sort(result_group_filtered.begin(), result_group_filtered.end());
            for (size_t j = 0; j < result_group_check.size(); j++) {
                if (result_group_filtered[j] != result_group_check[j]) {
                    std::cerr<< "ParticleMeshForSearchNeighbor::checkSearchNeighborForOneParticle(): index mismatch for group index "<<particle_index
                             << ", result_group["<<j<<"] = "<<result_group_filtered[j]
                             << ", result_group_check["<<j<<"] = "<<result_group_check[j]<<std::endl;
                    abort();
                }
            }            

            int n_filtered = result_filtered.size() + result_group_filtered.size();

            return Float(n_filtered)/Float(n_found);
        }

        template <class Tparticle, class Tpcm, class TGroup>
        bool checkSearchNeighborForAllParticles(const ParticleGroup<Tparticle, Tpcm>* particles, 
                                                const COMM::List<int>* particle_indices = nullptr,
                                                const COMM::List<TGroup>* groups = nullptr,
                                                const COMM::List<int>* group_indices = nullptr) const {
            ASSERT(is_cells_built); 
            int n_particles = particle_indices != nullptr ? particle_indices->getSize() : particles->getSize();
            int n_groups = group_indices != nullptr ? group_indices->getSize() : (groups != nullptr ? groups->getSize() : 0);
            Float ratio_match_mean = 0, ratio_match_max = 0, ratio_match_min = 1.0;
            for (int i=0; i<n_particles; i++){
                const int idx = particle_indices != nullptr ? (*particle_indices)[i] : i;
                Float ratio_found = checkSearchNeighborForOneParticle(idx, false, particles, particle_indices, groups, group_indices);
                if (ratio_found > ratio_match_max) ratio_match_max = ratio_found;
                if (ratio_found < ratio_match_min) ratio_match_min = ratio_found;
                ratio_match_mean += ratio_found;
            }
            for (int i=0; i<n_groups; i++){
                const int idx = group_indices != nullptr ? (*group_indices)[i] : i;
                Float ratio_found = checkSearchNeighborForOneParticle(idx, true, particles, particle_indices, groups, group_indices);
                if (ratio_found > ratio_match_max) ratio_match_max = ratio_found;
                if (ratio_found < ratio_match_min) ratio_match_min = ratio_found;
                ratio_match_mean += ratio_found;
            }
            ratio_match_mean /= (n_particles + n_groups);
            std::cout << "ParticleMeshForSearchNeighbor::checkSearchNeighborForAllParticles(): "
                      << " find neighbor match min = " << ratio_match_min
                      << ", mean = " << ratio_match_mean
                      << ", max = " << ratio_match_max << std::endl;

            return true;
        }
        
    };

} // namespace COMM
