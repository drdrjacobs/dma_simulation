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

#include "KDTreeVectorOfVectorsAdaptor.h"

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
    
    typedef Eigen::Matrix<float, DIMENSIONS, 1> Vec;
    typedef std::vector<State::Vec> PlatedCloud;
    typedef KDTreeVectorOfVectorsAdaptor<PlatedCloud, float, 
					 DIMENSIONS> KDTree;
    typedef std::uniform_real_distribution<float> Uniform;

    // Documented in the cpp
    static const int kDims;

    void write_xyz();
    float find_nearest_neighbor();
    void save_state();
    void load_state(std::string load_path, int max_leaf_size);

private:
    /// Position of diffusing particle
    Vec particle_;
    /// cloud of plated, std::vector of eigen vectors, used for kd_tree
    PlatedCloud plated_cloud_;
    /// nanoflann based k dimensional tree for finding nearest neighbor of 
    /// particle, have to use smart pointer here due to how KDTree is set up
    std::unique_ptr<KDTree> kd_tree_;
    //// tracks radius of cluster, distance from origin to furthest plated
    float radius_;
    /// random number generator
    std::mt19937 gen_;
    /// uniform [0, 1) distribution for random number generator
    Uniform uniform_;
};
