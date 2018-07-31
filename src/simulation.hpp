/// @file
///
/// @brief Header for the Simulation class which coordinates all calculations.
///

#pragma once

#include <string>
#include <vector>
#include <map>
#include <random>

#include "KDTreeVectorOfVectorsAdaptor.h"

/// @brief Simulation class which runs all calculations.
///
/// The simulation uses reduced units in terms of the particle radius, a, for 
/// length and the diffusion constant, D. This means the time unit is a^2 / D.
/// Time is not usually a feature of generalized aggregation, but including it
/// is actually the most convenient way to demonstrate the model's 
/// universality.
///
/// Uses a two zone approximation. When particles are far from plated, large 
/// jumps are used and the dynamics are exact. When particles are close to 
/// plated, a switch is made to approximate dynamics and time becomes 
/// important. In this second zone, the particle displacement in each 
/// dimension is drawn from a gaussian distrubution with a standard deviation
/// of sqrt(2 dt), D = 1.
///
/// Working with floats since the process is inherently Brownian so precision 
/// is not important.
///
/// All distance computations are run through the nanoflann library, see:
///
/// https://github.com/jlblancoc/nanoflann
///
/// As a result of this, plated (plated particles) are maintained in a kd tree.
///
/// Works for two or three dimensions. Set during compile time.
///
class Simulation {
public:
    Simulation();
    ~Simulation();
    
    // Documented in the cpp
    static const int kDims;
    static const std::string kParamsFilename;
    static const float kSpatialEpsilon;
   
private:
    typedef KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<float>>,
					 float, DIMENSIONS> KDTree;

    void run_simulation();
    std::map<std::string, std::string> read_params_file() const;
    std::string initialize_params();
    int set_up_structures(std::string restart_path,
			  std::vector<std::vector<float>> & plated,
			  KDTree & kd_tree);
    std::vector<float> generate_point_on_ball(int kDims, float radius);

    // Essential physical variables.
    
    /// Timestep, controls size of gaussian particle jumps when close to plated
    /// otherwise dynamics are exactly Brownian
    float dt_;
    /// The max distance a particle can jump in one timestep in one dimension,
    /// in terms of the standard deviations of the non-cutoff parent normal
    /// distribution, applies when close to plated, otherwise dynamics are 
    /// exactly Brownian
    float jump_cutoff_;

    /// The sticking probability at each contact
    float p_;
    /// size of the cluster to generate
    int cluster_size_;

    // Essential IO
    
    /// Writes simulation xyz frame after every interval in terms of the number
    /// of particles that have been added
    int write_frame_interval_;

    // Internal simulation parameters 
    
    /// Cells are cubic, side length of each cube is the max length particle 
    /// can jump in one dt_ plus one diameter + epsilon
    /// ensures all collisions can be resolved
    float cell_length_;
    /// optimization parameter for kd tree
    int max_leaf_size_;
    /// seed for rng
    int seed_;
    /// Rng on cpu
    std::mt19937 cpu_gen_;
};

