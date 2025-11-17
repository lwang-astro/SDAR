#pragma once

#include "Common/Float.h"
#include "particle_group.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

namespace COMM{
    class ParticleMesh{
    private:
        int n_div_[3]; // number of divisions in x,y,z
        Float box_size_[3]; // box sizes in x,y,z
        Float box_center_[3]; // box center in x,y,z
        std::vector<size_t>*** cells_; // cell storage (not used here, but can be used to store particle indices)
        std::vector<std::array<size_t,3>> cell_indices_; // list of all cell indices

    public:
        ParticleMesh(): n_div_{0,0,0},
                        box_size_{0.0, 0.0, 0.0}, 
                        box_center_{0.0, 0.0, 0.0}, 
                        cells_(nullptr), cell_indices_() {}

        template <class Tparticle, class Tpcm>
        ParticleMesh(const ParticleGroup<Tparticle, Tpcm>& particles, const int n_div_tot=1024) {
            // determine box size and center from particles
            Float pos_min[3] = {std::numeric_limits<Float>::max(), std::numeric_limits<Float>::max(), std::numeric_limits<Float>::max()};
            Float pos_max[3] = {std::numeric_limits<Float>::lowest(), std::numeric_limits<Float>::lowest(), std::numeric_limits<Float>::lowest()};
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

            // determine n_div_ based on box_size_ and n_div_tot
            Float box_volume = box_size_[0] * box_size_[1] * box_size[2];
            Float cell_volume_mean = box_volume / static_cast<Float>(n_div_tot);
            Float cell_size_mean = std::cbrt(cell_volume);

            // ensure each dimension have similar box size per cell
            for (int i=0; i<3; i++){
                n_div_[i] = static_cast<int>(std::ceil(box_size_[i] / cell_size_mean));
                if (n_div_[i] < 1) n_div_[i] = 1;
            }

            buildCells(particles);
        }

        template <class Tparticle, class Tpcm>        
        void buildCells(const ParticleGroup<Tparticle, Tpcm>& particles){
            assert(cells_ == nullptr);
            cells_ = new std::vector<size_t>**[n_div_[0]];    
            for (int i=0; i<n_div_[0]; i++){
                cells_[i] = new std::vector<size_t>*[n_div_[1]];
                for (int j=0; j<n_div_[1]; j++){
                    cells_[i][j] = new std::vector<size_t>[n_div_[2]];
                }
            }

            cell_indices_.resize(particles.getSize());
            for (int i=0; i<particles.getSize(); i++){
                const auto& pos = particles[i].pos;
                insertParticle(pos, i);
            }
        }

        void clearCells(){
            assert(cells_ != nullptr);
            for (int i=0; i<n_div_[0]; i++){
                for (int j=0; j<n_div_[1]; j++){
                    delete [] cells_[i][j];
                }
                delete [] cells_[i];
            }
            delete [] cells_;
        }

        ~ParticleMesh(){
            clearCells();
        }        
        
        void getCellIndex(const Float* pos, std::array<size_t,3> & cell_index) const {
            for (int i=0; i<3; i++){
                Float rel_pos = pos[i] - (box_center_[i] - 0.5*box_size_[i]);
                int idx = static_cast<int>(std::floor(rel_pos / box_size_[i] * n_div_[i]));
                if (idx < 0) idx = 0;
                if (idx >= n_div_[i]) idx = n_div_[i] - 1;
                cell_index[i] = idx;
            }
        }

        void insertParticle(const Float* pos, const size_t particle_index){
            std::array<size_t,3> cell_index;
            getCellIndex(pos, cell_index);
            cells_[cell_index[0]][cell_index[1]][cell_index[2]].push_back(particle_index);
            cell_indices_[particle_index] = cell_index;
        }
        
        void updateParticle(const Float* pos, const size_t particle_index){
            const auto& old_cell_index = cell_indices_[particle_index];
            std::array<size_t,3> new_cell_index;
            getCellIndex(pos, new_cell_index);
            if (old_cell_index != new_cell_index){
                // remove from old cell
                auto& old_cell = cells_[old_cell_index[0]][old_cell_index[1]][old_cell_index[2]];
                auto it = std::find(old_cell.begin(), old_cell.end(), particle_index);
                if (it != old_cell.end()){
                    old_cell.erase(it);
                }
                // insert into new cell
                cells_[new_cell_index[0]][new_cell_index[1]][new_cell_index[2]].push_back(particle_index);
                cell_indices_[particle_index] = new_cell_index;
            }
        }

        const std::vector<size_t>& getCellParticles(const std::array<size_t,3>& cell_index) const {
            return cells_[cell_index[0]][cell_index[1]][cell_index[2]];
        }

        void searchPosRange(const Float pos_min[3], const Float pos_max[3], std::vector<size_t>& result) const {
            std::array<size_t,3> cell_min, cell_max;
            getCellIndex(pos_min, cell_min);
            getCellIndex(pos_max, cell_max);
            for (int i=cell_min[0]; i<=cell_max[0]; i++){
                for (int j=cell_min[1]; j<=cell_max[1]; j++){
                    for (int k=cell_min[2]; k<=cell_max[2]; k++){
                        const auto& cell_particles = cells_[i][j][k];
                        result.insert(result.end(), cell_particles.begin(), cell_particles.end());
                    }
                }
            }
        }

        void searchNeighbor(const Float pos[3], const Float radius, std::vector<size_t>& result) const {
            Float pos_min[3], pos_max[3];
            for (int i=0; i<3; i++){
                pos_min[i] = pos[i] - radius;
                pos_max[i] = pos[i] + radius;
            }
            searchPosRange(pos_min, pos_max, result);
        }
    };

} // namespace COMM
