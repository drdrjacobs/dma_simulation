/// @file
///
/// @brief Header for the Simulation class which coordinates all calculations.
///

#pragma once

#include <string>
#include <vector>
#include <map>
#include <random>

/// @brief Simulation class which runs all calculations.
///
/// The simulation uses reduced units in terms of the particle radius, a, for 
/// length and the diffusion constant, D. This means the time unit is a^2 / D. 
/// Time is not usually a feature of generalized aggregation, but including it
/// when particles move over their step size is actually the most convenient 
/// way to demonstrate the model's universality.
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
class Simulation {
public:
    Simulation();
    ~Simulation();
    
    // Documented in the cpp
    static const std::string kParamsFilename;
   
private:
    void run_simulation();
    std::map<std::string, std::string> read_params_file() const;
    std::string initialize_params();

    // Essential physical variables.
    
    /// dimension of cluster
    int dimension_;
    /// step size that particle travels over
    float step_size_;
    /// The sticking probability at each contact
    float p_;
    /// size of the cluster to generate
    int cluster_size_;

    // Essential IO
    
    /// Writes simulation xyz frame after every interval in terms of the number
    /// of particles that have been added
    int write_frame_interval_;

    // Internal simulation parameters 
    
    /// seed for rng
    int seed_;
    /// Rng on cpu
    std::mt19937 cpu_gen_;
};

