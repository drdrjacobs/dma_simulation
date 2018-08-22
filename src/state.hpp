/// @file
///
/// @brief Header for the State class which stores all data that changes as
/// simulation advances.
///

#pragma once

#include <vector>
#include <random>
#include <memory>
#include <string>

#include "constants.hpp"
#include "KDTreeVectorOfVectorsAdaptor.h"
#include "cells.hpp"

/// @brief State class manages data structures that change as simulation
/// advances.
///
/// Works for two or three dimensions. Set during compile time.
///
/// Runs exact dynamics and approximate dynamics depending on whether particle
/// is close to plated, see @ref Simulation.
///
/// All nearest neighbor computations for exact dynamics are run through the
/// nanoflann library, see:
///
/// https://github.com/jlblancoc/nanoflann
///
/// As a result of this, plated (plated particles) are maintained in a kd tree.
///
/// For approximate dynamics, plated are stored in a cells structure, see @ref 
/// Cells.
///
class State {
public:
    
    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vec>, double, 
					 kDims> KDTree;
    typedef Cells::CellIndices CellIndices;

    // Documented in the cpp
    static const int kNoCollision;

    /// @brief blank constructor, must call ethier load_state or 
    /// set_up_new_state.
    ///
    State() {};
    /// @brief empy destructor
    ///
    ~State() {};
    /// @brief Gets current size of cluster.
    ///
    /// @returns cluster_size: the current size of the cluster
    ///
    int get_cluster_size() const {return plated_cloud_.size();}
    /// @brief Gets particle position.
    ///
    /// @returns particle: the particle position
    ///
    Vec get_particle() {return particle_;}
    /// @brief Sets particle position.
    ///
    /// @param new_particle: new position
    ///
    void set_particle(Vec new_particle) {particle_ = new_particle;}
    /// @brief Gets plated_cloud.
    ///
    /// @param returns plated cloud: vector tracking all plated
    ///
    std::vector<Vec>& get_plated_cloud() {return plated_cloud_;}

    void set_up_new_state(double cell_length, int max_leaf_size, int seed);
    void load_state(double cell_length, int max_leaf_size, 
		    std::string load_path);
    void save_state() const;
    // for propogation of dynamics
    void add_new_particle();
    bool has_neighbors() const;
    void stick_particle();
    double check_collisions_loop(Vec jump_unit_vector,
				 double jump_length, Vec jump,
				 CellIndices &bounce_cell_indices,
				 size_t &bounce_plated_index) const;
    bool resolve_jump(Vec jump, double p);
    bool take_small_step(double dt, double jump_cutoff, double p_);
    double find_nearest_neighbor() const;
    void take_large_step();
    void check_for_regeneration();
    // more output
    int check_overlaps(bool verbose = false) const;
    void write_xyz() const;

private:
    /// Indicates outcome of forward step
    enum Status {free, stuck, rejection};
    Status attempt_bounce(double minimum_collision_distance,
			  CellIndices bounce_cell_indices,
			  size_t bounce_plated_index,
			  Vec initial_particle,
			  Vec jump_unit_vector, double jump_length,
			  Vec &jump);
    Status resolve_stick_or_bounce(double p, double minimum_collision_distance,
				   CellIndices bounce_cell_indices,
				   size_t bounce_plated_index,
				   Vec initial_particle,
				   Vec jump_unit_vector, double jump_length,
				   Vec &jump);

    /// Position of diffusing particle
    Vec particle_;
    /// cloud of plated, std::vector of eigen vectors, used for kd_tree
    std::vector<Vec> plated_cloud_;
    /// nanoflann based k dimensional tree for finding nearest neighbor of 
    /// particle, have to use smart pointer here due to how KDTree is set up
    std::unique_ptr<KDTree> kd_tree_;
    /// Object handles storing plated in cell structure
    Cells cells_;
    //// tracks radius of cluster, distance from origin to furthest plated
    double radius_;
    /// random number generator
    std::mt19937 gen_;
    /// uniform [0, 1) distribution for random number generator
    Uniform uniform_;
};

// Free functions
double get_collision_distance(Vec particle, Vec plated_r, 
			      Vec jump_unit_vector);
double calculate_collisions(double minimum_collision_distance,
			    Vec particle, Vec plated_r,
			    Vec jump_unit_vector, double jump_length, 
			    Vec jump);
