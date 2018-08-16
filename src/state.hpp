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

// use eigen for vector operations
#include <Eigen/Dense>

#include "vec.hpp"
#include "KDTreeVectorOfVectorsAdaptor.h"
#include "cells.hpp"

// forward declaration for friend class
class Simulation;

/// @brief State class includes all data structures that change as simulation
/// advances.
///
/// Works for two or three dimensions. Set during compile time.
///
class State {
public:
    
    friend class Simulation;

    State() {};
    ~State() {};
    
    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vec>, double, 
					 DIMENSIONS> KDTree;
    typedef std::uniform_real_distribution<double> Uniform;

    int check_overlaps() const;
    void write_xyz() const;
    double find_nearest_neighbor() const;
    void save_state() const;
    void load_state(std::string load_path, int max_leaf_size);

private:
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

    // for resolving collisions
    /// initial position of particle before step
    Vec initial_particle_;
    /// indices of cell containing plated that particle just bounced off of
    Cells::CellIndices bounce_cell_indices_;
    /// index of plated that particle just bounced off of
    size_t bounce_plated_index_;
};
