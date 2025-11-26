#pragma once

#include "Common/Float.h"
#include "Common/particle_group.h"
#include "Common/list.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <numeric> // for std::iota
#include <type_traits> 

namespace COMM {

    // Polyfill for std::void_t (C++11/C++14 compatible)
    template<typename... Ts> struct make_void { typedef void type; };
    template<typename... Ts> using void_t = typename make_void<Ts...>::type;

    // Helper to abstract size retrieval for different container types
    template <class T>
    size_t get_container_size(const std::vector<T>& c) { return c.size(); }

    template <class Tptcl, class Tcm>
    size_t get_container_size(const ParticleGroup<Tptcl, Tcm>& c) { return c.getSize(); }

    template <class T>
    size_t get_container_size(const List<T>& c) { return c.getSize(); }

    template <class T>
    int get_container_size(const T* c) { return -1; } // For raw pointers, size must be provided separately

    // ---------------------------------------------------------
    // Unified Accessor using SFINAE
    // ---------------------------------------------------------
    
    // Primary template: Assumes T is a Particle (has .pos and .r_search)
    template <typename T, typename = void>
    struct TargetAccessor {
        static const Float* getPos(const T& p) { return &p.pos[0]; }
        static Float getRSearch(const T& p) { return p.r_search; }
    };

    // Specialization: Assumes T is a Group (has .particles.cm)
    // Detects if T has a member named 'particles'
    template <typename T>
    struct TargetAccessor<T, COMM::void_t<decltype(T::particles)>> {
        static const Float* getPos(const T& p) { return &p.particles.cm.pos[0]; }
        static Float getRSearch(const T& p) { return p.particles.cm.r_search; }
    };

    // ---------------------------------------------------------
    // Core KD-Tree Implementation (Generic)
    // ---------------------------------------------------------
    // CHANGED: Removed Accessor template parameter. 
    // Now uses TargetAccessor<T> internally for all types.
    class KDTreeCore {
    public:
        struct Node {
            int axis;           // Split axis
            int left = -1;
            int right = -1;
            int parent = -1;    // Parent index for bottom-up updates
            
            int index;          // Store index
            
            Float pos[3];       // Cache position
            Float r_search;     // Cache r_search
            
            // Bounding box & Max search radius in subtree (for pruning)
            Float min_box[3];
            Float max_box[3];
            Float max_r_subtree;
            Float min_r_subtree; // Min search radius in subtree
            
            bool is_removed = false;

            // Constructor uses TargetAccessor to get data
            template <typename ParticleType>
            Node(int idx, const ParticleType& p, int ax, int p_idx = -1) 
                : axis(ax), parent(p_idx), index(idx) { 
                // CHANGED: Use TargetAccessor
                const Float* p_pos = TargetAccessor<ParticleType>::getPos(p);
                for(int k=0; k<3; k++) {
                    pos[k] = p_pos[k];
                    min_box[k] = p_pos[k];
                    max_box[k] = p_pos[k];
                }
                r_search = TargetAccessor<ParticleType>::getRSearch(p);
                max_r_subtree = r_search;
                min_r_subtree = r_search;
            }
        };

    private:
        std::vector<Node> nodes_;
        std::vector<int> index_to_node_; 
        int root_ = -1;
        int active_count_ = 0;
        int removed_count_ = 0;

        //! Update bounding box and max_r_subtree for a node
        /*!
          @param[in] node_idx: index of the node to update
        */
        void update_node_stats(int node_idx) {
            if (node_idx == -1) return;
            Node& node = nodes_[node_idx];
            
            // Reset to self
            for(int k=0; k<3; k++) {
                node.min_box[k] = node.pos[k];
                node.max_box[k] = node.pos[k];
            }
            node.max_r_subtree = node.r_search;
            node.min_r_subtree = node.r_search;

            int children[] = {node.left, node.right};
            for(int child : children) {
                if(child != -1) {
                    Node& c = nodes_[child];
                    for(int k=0; k<3; k++) {
                        node.min_box[k] = std::min(node.min_box[k], c.min_box[k]);
                        node.max_box[k] = std::max(node.max_box[k], c.max_box[k]);
                    }
                    node.max_r_subtree = std::max(node.max_r_subtree, c.max_r_subtree);
                    node.min_r_subtree = std::min(node.min_r_subtree, c.min_r_subtree);
                }
            }
        }

        //! Recursive function to build the tree
        template <class TContainer>
        int build_recursive(std::vector<int>& indices, const TContainer& particles, 
                            int start, int end, int depth, int parent_idx) {
            if (start >= end) return -1;
            int axis = depth % 3;
            int mid = (start + end) / 2;
            
            // Sort indices based on TargetAccessor::getPos
            // We need to deduce the value type of the container to use TargetAccessor
            // For std::vector<T>, value_type is T. For ParticleGroup, operator[] returns reference to Tptcl.
            // Using auto& p = particles[a] works.
            std::nth_element(indices.begin() + start, indices.begin() + mid, indices.begin() + end,
                [&particles, axis](int a, int b) { 
                    // CHANGED: Use TargetAccessor with decltype
                    using PType = typename std::decay<decltype(particles[0])>::type;
                    return TargetAccessor<PType>::getPos(particles[a])[axis] < TargetAccessor<PType>::getPos(particles[b])[axis]; 
                });

            int idx = indices[mid];
            // Pass parent_idx to constructor
            nodes_.emplace_back(idx, particles[idx], axis, parent_idx);
            int current_node_idx = nodes_.size() - 1;

            // Update map
            if (idx >= static_cast<int>(index_to_node_.size())) index_to_node_.resize(idx + 1, -1);
            index_to_node_[idx] = current_node_idx;

            int left_child = build_recursive(indices, particles, start, mid, depth + 1, current_node_idx);
            nodes_[current_node_idx].left = left_child;

            int right_child = build_recursive(indices, particles, mid + 1, end, depth + 1, current_node_idx);
            nodes_[current_node_idx].right = right_child;

            update_node_stats(current_node_idx);
            return current_node_idx;
        }

        //! Calculate squared distance from a point to a box
        /*!
          @param[in] point: 3D point coordinates
          @param[in] min_b: minimum coordinates of the box
          @param[in] max_b: maximum coordinates of the box
          @return squared distance
        */
        Float dist_sq_point_to_box(const Float point[3], const Float min_b[3], const Float max_b[3]) {
            Float d2 = 0.0;
            for (int i = 0; i < 3; ++i) {
                if (point[i] < min_b[i]) d2 += (min_b[i] - point[i]) * (min_b[i] - point[i]);
                else if (point[i] > max_b[i]) d2 += (point[i] - max_b[i]) * (point[i] - max_b[i]);
            }
            return d2;
        }

        //! Recursive function to search neighbors
        /*!
          @param[in] node_idx: current node index
          @param[in] target: target particle/group
          @param[out] neighbor_list: list to store found neighbor indices
        */
        // CHANGED: Removed ignore_index parameter
        template <class Ttarget>
        void search_recursive(int node_idx, const Ttarget& target, std::vector<int>& neighbor_list) {
            if (node_idx == -1) return;
            Node& node = nodes_[node_idx];

            // Use TargetAccessor to get target properties
            const Float* t_pos = TargetAccessor<Ttarget>::getPos(target);
            Float t_r = TargetAccessor<Ttarget>::getRSearch(target);

            // 1. Pruning
            Float d2_box = dist_sq_point_to_box(t_pos, node.min_box, node.max_box);
            Float max_dist = std::max(t_r, node.max_r_subtree);
            
            // Use a safe tolerance
            if (d2_box > max_dist * max_dist * 1.0001 + 1e-18) return;

            // 2. Check Particle
            if (!node.is_removed) {
                Float d2 = 0.0;
                for(int k=0; k<3; k++) d2 += (t_pos[k] - node.pos[k])*(t_pos[k] - node.pos[k]);
                
                Float r_crit = std::max(t_r, node.r_search);
                if (d2 < r_crit * r_crit) {
                    neighbor_list.push_back(node.index);
                }
            }

            // 3. Recurse
            search_recursive(node.left, target, neighbor_list);
            search_recursive(node.right, target, neighbor_list);
        }

        //! Recursive function to search neighbors and apply function
        template <class Ttarget, typename Func>
        void search_recursive_apply(int node_idx, const Ttarget& target, Func&& func) {
            if (node_idx == -1) return;
            Node& node = nodes_[node_idx];

            const Float* t_pos = TargetAccessor<Ttarget>::getPos(target);
            Float t_r = TargetAccessor<Ttarget>::getRSearch(target);

            // 1. Pruning
            Float d2_box = dist_sq_point_to_box(t_pos, node.min_box, node.max_box);
            Float max_dist = std::max(t_r, node.max_r_subtree);
            
            if (d2_box > max_dist * max_dist * 1.0001 + 1e-18) return;

            // 2. Check Particle
            if (!node.is_removed) {
                Float d2 = 0.0;
                for(int k=0; k<3; k++) d2 += (t_pos[k] - node.pos[k])*(t_pos[k] - node.pos[k]);
                
                Float r_crit = std::max(t_r, node.r_search);
                if (d2 < r_crit * r_crit) {
                    // FOUND NEIGHBOR: Call the function directly
                    func(node.index);
                }
            }

            // 3. Recurse
            search_recursive_apply(node.left, target, func);
            search_recursive_apply(node.right, target, func);
        }

    public:
        KDTreeCore() {}

        //! Clear the tree
        void clear() {
            nodes_.clear();
            // We don't necessarily clear index_to_node_ to save reallocation, 
            // but we should fill it with -1 if we want to be safe. 
            // For performance, we just rely on resize in build.
            std::fill(index_to_node_.begin(), index_to_node_.end(), -1);
            root_ = -1;
            active_count_ = 0;
            removed_count_ = 0;
        }

        //! Build tree from particles (all or subset)
        /*!
          @param[in] particles: container of particles
          @param[in] subset_indices: optional pointer to list of indices to include. If nullptr, use all/n_particle.
          @param[in] n_particle: number of particles to include (only used if subset_indices is nullptr). If -1, use container size.
        */
        template <class TContainer>
        void build(const TContainer& particles, const List<int>* subset_indices = nullptr, int n_particle = -1) {
            clear();
            
            std::vector<int> indices;
            int max_idx = 0;

            if (subset_indices) {
                // --- Subset Mode ---
                size_t n = subset_indices->getSize();
                if (n == 0) return;
                
                indices.resize(n);
                // We need to find max_idx for resizing index_to_node_
                for (size_t i = 0; i < n; ++i) {
                    indices[i] = (*subset_indices)[i];
                    if (indices[i] > max_idx) max_idx = indices[i];
                }
            } else {
                // --- All / First N Mode ---
                size_t n;
                if (n_particle != -1) n = n_particle;
                else {
                    int n_container = get_container_size(particles);
                    assert(n_container >= 0); // Ensure size is known
                    n = static_cast<size_t>(n_container);
                }
                if (n == 0) return;

                indices.resize(n);
                std::iota(indices.begin(), indices.end(), 0);
                max_idx = n > 0 ? static_cast<int>(n) - 1 : 0;
            }

            // Common logic: Reserve memory
            size_t reserve_cap = static_cast<size_t>(indices.size() * 1.25) + 16;
            if (nodes_.capacity() < reserve_cap) nodes_.reserve(reserve_cap);

            // Common logic: Resize map
            if (static_cast<int>(index_to_node_.size()) <= max_idx) index_to_node_.resize(max_idx + 1, -1);

            active_count_ = indices.size();
            root_ = build_recursive(indices, particles, 0, indices.size(), 0, -1);
        }

        //! Search neighbors for a target
        /*!
          @param[in] target: target particle/group
          @param[out] neighbor_list: output list of neighbor indices
        */
        // CHANGED: Removed ignore_index parameter
        template <class Ttarget>
        void search(const Ttarget& target, std::vector<int>& neighbor_list) {
            if (root_ == -1) return;
            search_recursive(root_, target, neighbor_list);
        }

        //! Search neighbors and apply function
        /*!
          @param[in] target: target particle/group
          @param[in] func: function to apply to each found neighbor index
        */
        template <class Ttarget, typename Func>
        void searchApply(const Ttarget& target, Func&& func) {
            if (root_ == -1) return;
            search_recursive_apply(root_, target, std::forward<Func>(func));
        }

        //! Insert a single particle into the tree
        /*!
          @param[in] index: index of the particle in the external container
          @param[in] p: particle object
        */
        // CHANGED: Made insert a template
        template <typename ParticleType>
        void insert(int index, const ParticleType& p) {
            // Ensure map is large enough
            if (index >= static_cast<int>(index_to_node_.size())) index_to_node_.resize(index + 1, -1);

            // CHANGED: Reserve small initial capacity if empty to avoid frequent reallocs on start
            if (nodes_.capacity() == 0) nodes_.reserve(16);

            if (root_ == -1) {
                nodes_.emplace_back(index, p, 0, -1); // Parent -1
                root_ = 0;
                active_count_++;
                index_to_node_[index] = 0; 
                return;
            }
            
            // CHANGED: Use TargetAccessor
            const Float* p_pos = TargetAccessor<ParticleType>::getPos(p);
            Float p_r = TargetAccessor<ParticleType>::getRSearch(p);

            int curr = root_;
            while(true) {
                // Use index access to avoid reference invalidation issues later, 
                // but for updating stats NOW, reference is fine as long as we don't emplace_back yet.
                Node& node = nodes_[curr];
                
                // Expand bbox and max_r
                for(int k=0; k<3; k++) {
                    node.min_box[k] = std::min(node.min_box[k], p_pos[k]);
                    node.max_box[k] = std::max(node.max_box[k], p_pos[k]);
                }
                node.max_r_subtree = std::max(node.max_r_subtree, p_r);
                node.min_r_subtree = std::min(node.min_r_subtree, p_r);

                int axis = node.axis;
                bool go_left = p_pos[axis] < node.pos[axis];
                
                // Check if child exists
                int child_idx = go_left ? node.left : node.right;
                
                if (child_idx == -1) {
                    // We need to insert here.
                    int next_axis = (axis + 1) % 3;
                    
                    // CAUTION: emplace_back may invalidate 'node' reference!
                    // Do not use 'node' or pointers to its members after this line.
                    nodes_.emplace_back(index, p, next_axis, curr);
                    int new_idx = nodes_.size() - 1;
                    
                    // Re-access the parent node using the stable index 'curr'
                    if (go_left) nodes_[curr].left = new_idx;
                    else         nodes_[curr].right = new_idx;
                    
                    active_count_++;
                    index_to_node_[index] = new_idx;
                    break;
                }
                curr = child_idx;
            }
        }

        //! Remove a particle from the tree (lazy deletion)
        /*!
          @param[in] index: index of the particle to remove
        */
        void remove(int index) {
            // Optimized remove using map: O(1)
            if (index < 0 || index >= static_cast<int>(index_to_node_.size())) return;
            
            int node_idx = index_to_node_[index];
            if (node_idx != -1) {
                Node& node = nodes_[node_idx];
                if (!node.is_removed) {
                    node.is_removed = true;
                    active_count_--;
                    removed_count_++;
                }
                // We keep the mapping pointing to the removed node or set to -1?
                // Setting to -1 allows re-insertion logic to know it's gone from active duty.
                index_to_node_[index] = -1; 
            }
        }

        //! Update a particle's position/radius in the tree
        /*!
          @note This performs a lazy remove followed by an insert. 
                Frequent updates will degrade tree quality and increase memory usage.
                Call build() periodically to rebalance.
          @param[in] index: index of the particle
          @param[in] p: new particle object (with updated pos/r_search)
          @param[in] force_reinsert: if true, force remove+insert (default false)
        */
        // CHANGED: Made update a template
        template <typename ParticleType>
        void update(int index, const ParticleType& p, bool force_reinsert = false) {
            if (index < 0 || index >= static_cast<int>(index_to_node_.size())) return;
            int node_idx = index_to_node_[index];
            if (node_idx == -1) {
                insert(index, p);
                return;
            }

            if (force_reinsert) {
                remove(index);
                insert(index, p);
                return;
            }

            Node& node = nodes_[node_idx];
            if (node.is_removed) {
                remove(index);
                insert(index, p);
                return;
            }

            // CHANGED: Use TargetAccessor
            const Float* new_pos = TargetAccessor<ParticleType>::getPos(p);
            Float new_r = TargetAccessor<ParticleType>::getRSearch(p);

            // Check displacement
            Float d2 = 0.0;
            for(int k=0; k<3; k++) d2 += (new_pos[k] - node.pos[k])*(new_pos[k] - node.pos[k]);
            
            // Threshold: if moved less than 10% of search radius (or a small constant), update in-place
            // This is a heuristic. 
            Float threshold = node.r_search * 0.1; 
            if (threshold < 1e-4) threshold = 1e-4; // Min threshold

            if (d2 < threshold * threshold) {
                // --- In-Place Update ---
                // 1. Update Node Data
                for(int k=0; k<3; k++) node.pos[k] = new_pos[k];
                node.r_search = new_r;

                // 2. Bottom-Up Update of Bounding Boxes
                int curr = node_idx;
                while (curr != -1) {
                    Node& c = nodes_[curr];
                    
                    // Re-calculate stats for 'c' based on itself and children
                    // Reset to self (which is now updated)
                    for(int k=0; k<3; k++) {
                        c.min_box[k] = c.pos[k];
                        c.max_box[k] = c.pos[k];
                    }
                    c.max_r_subtree = c.r_search;
                    c.min_r_subtree = c.r_search;

                    // Merge children
                    int children[] = {c.left, c.right};
                    for(int child : children) {
                        if(child != -1) {
                            Node& child_node = nodes_[child];
                            // Important: Child might be removed, but its box is still valid for tree structure?
                            // Actually, if child is removed, we shouldn't count it? 
                            // Standard KDTree usually keeps removed nodes in structure until rebuild.
                            // So we include them to maintain box integrity.
                            for(int k=0; k<3; k++) {
                                c.min_box[k] = std::min(c.min_box[k], child_node.min_box[k]);
                                c.max_box[k] = std::max(c.max_box[k], child_node.max_box[k]);
                            }
                            c.max_r_subtree = std::max(c.max_r_subtree, child_node.max_r_subtree);
                            c.min_r_subtree = std::min(c.min_r_subtree, child_node.min_r_subtree);
                        }
                    }
                    
                    // Move up
                    curr = c.parent;
                }
            } else {
                // Moved too much, re-insert
                remove(index);
                insert(index, p);
            }
        }

        //! Check if the tree is currently built (has active nodes)
        bool isBuilt() const {
            return root_ != -1;
        }

        //! Check if building the tree is beneficial based on particle distribution
        /*!
          @param[in] particles: container of particles
          @param[in] n_limit: minimum number of particles to justify tree build (default 16)
          @param[in] r_ratio_limit: if min_r_search > max_extent * r_ratio_limit, return false (default 1.0)
          @return true if tree should be built
        */
        template <class TContainer>
        bool checkBuildCondition(const TContainer& particles, int n_limit = 32, Float r_ratio_limit = 0.3) {
            size_t n = get_container_size(particles);
            if (n <= static_cast<size_t>(n_limit)) return false;

            Float min_box[3] = {NUMERIC_FLOAT_MAX, NUMERIC_FLOAT_MAX, NUMERIC_FLOAT_MAX};
            Float max_box[3] = {-NUMERIC_FLOAT_MAX, -NUMERIC_FLOAT_MAX, -NUMERIC_FLOAT_MAX};
            Float min_r = NUMERIC_FLOAT_MAX;

            for (size_t i = 0; i < n; ++i) {
                const auto& p = particles[i];
                // CHANGED: Use TargetAccessor with decltype
                using PType = typename std::decay<decltype(p)>::type;
                const Float* pos = TargetAccessor<PType>::getPos(p);
                Float r = TargetAccessor<PType>::getRSearch(p);

                for(int k=0; k<3; k++) {
                    if (pos[k] < min_box[k]) min_box[k] = pos[k];
                    if (pos[k] > max_box[k]) max_box[k] = pos[k];
                }
                if (r < min_r) min_r = r;
            }

            Float max_extent = 0.0;
            for(int k=0; k<3; k++) max_extent = std::max(max_extent, max_box[k] - min_box[k]);

            if (max_extent <= 0.0) return false; 
            
            // If the smallest search radius is larger than the system size * ratio, 
            // then ALL particles have search radius larger than system size.
            // In this case, tree search is effectively O(N^2) with overhead.
            if (min_r > max_extent * r_ratio_limit) return false;

            return true;
        }

        //! Check if building the tree is beneficial (subset version)
        template <class TContainer>
        bool checkBuildCondition(const TContainer& particles, const List<int>& subset_indices, int n_limit = 32, Float r_ratio_limit = 0.3) {
            size_t n = subset_indices.getSize();
            if (n <= static_cast<size_t>(n_limit)) return false;

            Float min_box[3] = {NUMERIC_FLOAT_MAX, NUMERIC_FLOAT_MAX, NUMERIC_FLOAT_MAX};
            Float max_box[3] = {-NUMERIC_FLOAT_MAX, -NUMERIC_FLOAT_MAX, -NUMERIC_FLOAT_MAX};
            Float min_r = NUMERIC_FLOAT_MAX;

            for (size_t k = 0; k < n; ++k) {
                int idx = subset_indices[k];
                const auto& p = particles[idx];
                // CHANGED: Use TargetAccessor with decltype
                using PType = typename std::decay<decltype(p)>::type;
                const Float* pos = TargetAccessor<PType>::getPos(p);
                Float r = TargetAccessor<PType>::getRSearch(p);

                for(int d=0; d<3; d++) {
                    if (pos[d] < min_box[d]) min_box[d] = pos[d];
                    if (pos[d] > max_box[d]) max_box[d] = pos[d];
                }
                if (r < min_r) min_r = r;
            }

            Float max_extent = 0.0;
            for(int d=0; d<3; d++) max_extent = std::max(max_extent, max_box[d] - min_box[d]);

            if (max_extent <= 0.0) return false; 
            if (min_r > max_extent * r_ratio_limit) return false;

            return true;
        }
    };

    // ---------------------------------------------------------
    // ParticleKDTree (Manages both Particles and Groups)
    // ---------------------------------------------------------
    class ParticleKDTree {
    private:
        // CHANGED: KDTreeCore is now generic, no template args needed
        KDTreeCore tree_ptcl_;
        KDTreeCore tree_group_;

    public:
        //! Clear both particle and group trees
        void clear() {
            tree_ptcl_.clear();
            tree_group_.clear();
        }

        // --- Build Methods ---
        
        //! Add particles to the tree (generic)
        /*!
          @param[in] particles: container of particles (ParticleGroup, std::vector, or raw pointer)
          @param[in] indices: optional pointer to list of indices to include. If nullptr, use all/n_particles.
          @param[in] n_particles: number of particles (used if particles is raw pointer or to limit count).
        */
        template <class TContainer>
        void addParticles(const TContainer& particles, const COMM::List<int>* indices = nullptr, int n_particles = -1) {
            tree_ptcl_.build(particles, indices, n_particles);
        }

        //! Add groups to the tree (generic)
        /*!
          @param[in] groups: container of groups (List, std::vector, or raw pointer)
          @param[in] indices: optional pointer to list of indices to include. If nullptr, use all/n_groups.
          @param[in] n_groups: number of groups (used if groups is raw pointer or to limit count).
        */
        template <class TContainer>
        void addGroups(const TContainer& groups, const COMM::List<int>* indices = nullptr, int n_groups = -1) {
            tree_group_.build(groups, indices, n_groups);
        }

        // --- Search Methods ---        
        //! Search neighbors for a target in both particle and group trees
        /*!
          @param[in] target: target particle or group
          @param[out] p_idx_list: output list of neighbor particle indices
          @param[out] g_idx_list: output list of neighbor group indices
        */
        // CHANGED: Removed target_index and target_is_group parameters
        template <class Ttarget>
        void searchNeighbor(const Ttarget& target, 
                            std::vector<int>& p_idx_list, 
                            std::vector<int>& g_idx_list) {
            p_idx_list.clear();
            g_idx_list.clear();

            tree_ptcl_.search(target, p_idx_list);
            tree_group_.search(target, g_idx_list);
        }

        //! Search neighbors only in particle tree
        /*!
          @param[in] target: target particle or group
          @param[out] p_idx_list: output list of neighbor particle indices
        */
        // CHANGED: Removed target_index and target_is_group parameters
        template <class Ttarget>
        void searchNeighborParticles(const Ttarget& target, std::vector<int>& p_idx_list) {
            p_idx_list.clear();
            tree_ptcl_.search(target, p_idx_list);
        }

        //! Search neighbors only in group tree
        /*!
          @param[in] target: target particle or group
          @param[out] g_idx_list: output list of neighbor group indices
        */
        // CHANGED: Removed target_index and target_is_group parameters
        template <class Ttarget>
        void searchNeighborGroups(const Ttarget& target, std::vector<int>& g_idx_list) {
            g_idx_list.clear();
            tree_group_.search(target, g_idx_list);
        }

        //! Search neighbors in particle tree and apply function
        /*!
            @param[in] target: target particle or group
            @param[in] func: function to apply to each found neighbor index
        */
        template <class Ttarget, typename Func>
        void searchNeighborParticlesApply(const Ttarget& target, Func&& func) {
            tree_ptcl_.searchApply(target, std::forward<Func>(func));
        }

        //! Search neighbors in group tree and apply function
        /*!
            @param[in] target: target particle or group
            @param[in] func: function to apply to each found neighbor index
        */
        template <class Ttarget, typename Func>
        void searchNeighborGroupsApply(const Ttarget& target, Func&& func) {
            tree_group_.searchApply(target, std::forward<Func>(func));
        }

        // --- Dynamic Update Methods ---
        
        //! Insert a single particle into the tree
        /*!
          @param[in] index: index of the particle
          @param[in] p: particle object
        */
        // CHANGED: Added template parameter
        template <class Tptcl>
        void InsertParticle(int index, const Tptcl& p) { tree_ptcl_.insert(index, p); }

        //! Insert a single group into the tree
        /*!
          @param[in] index: index of the group
          @param[in] g: group object
        */
        // CHANGED: Added template parameter
        template <class Tgroup>
        void InsertGroup(int index, const Tgroup& g) { tree_group_.insert(index, g); }

        //! Remove a particle from the tree
        /*!
          @param[in] index: index of the particle to remove
        */
        void removeParticle(int index) { tree_ptcl_.remove(index); }

        //! Remove a group from the tree
        /*!
          @param[in] index: index of the group to remove
        */
        void removeGroup(int index) { tree_group_.remove(index); }

        //! Update a particle in the tree
        /*!
          @param[in] index: index of the particle
          @param[in] p: new particle data
        */
        // CHANGED: Added template parameter
        template <class Tptcl>
        void updateParticle(int index, const Tptcl& p) { tree_ptcl_.update(index, p); }

        //! Update a group in the tree
        /*!
          @param[in] index: index of the group
          @param[in] g: new group data
        */
        // CHANGED: Added template parameter
        template <class Tgroup>
        void updateGroup(int index, const Tgroup& g) { tree_group_.update(index, g); }

        // --- Status Check Methods ---

        //! Check if particle tree is built
        bool isParticleTreeBuilt() const { return tree_ptcl_.isBuilt(); }

        //! Check if group tree is built
        bool isGroupTreeBuilt() const { return tree_group_.isBuilt(); }

        //! Check if particle tree should be built
        template <class Tptcl, class Tcm>
        bool checkParticleBuildCondition(const ParticleGroup<Tptcl, Tcm>& particles, int n_limit = 64, Float r_ratio_limit = 0.3) {
            return tree_ptcl_.checkBuildCondition(particles, n_limit, r_ratio_limit);
        }

        //! Check if particle tree should be built (subset)
        template <class Tptcl, class Tcm>
        bool checkParticleBuildCondition(const ParticleGroup<Tptcl, Tcm>& particles, const COMM::List<int>& indices, int n_limit = 64, Float r_ratio_limit = 0.3) {
            return tree_ptcl_.checkBuildCondition(particles, indices, n_limit, r_ratio_limit);
        }

        //! Check if group tree should be built
        template <class Tgroup>
        bool checkGroupBuildCondition(const List<Tgroup>& groups, int n_limit = 64, Float r_ratio_limit = 0.3) {
            return tree_group_.checkBuildCondition(groups, n_limit, r_ratio_limit);
        }

        //! Check if group tree should be built (subset)
        template <class Tgroup>
        bool checkGroupBuildCondition(const List<Tgroup>& groups, const COMM::List<int>& indices, int n_limit = 64, Float r_ratio_limit = 0.3) {
            return tree_group_.checkBuildCondition(groups, indices, n_limit, r_ratio_limit);
        }

    };

} // namespace COMM