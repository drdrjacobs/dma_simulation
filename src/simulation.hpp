/// @file
///
/// @brief Header for the Simulation class which coordinates all calculations.
///

#pragma once

#include <string>
#include <map>
#include <utility>

#include "constants.hpp"
#include "state.hpp"

/// @brief Simulation class which runs all calculations by interfacing with the
/// system @ref State.
///
/// Works for two or three dimensions. Set during compile time.
///
/// The simulation uses reduced units in terms of the particle radius, a, for 
/// length and the diffusion constant, D. This means the time unit is a^2 / D.
/// Time is not usually a feature of generalized aggregation, but including it
/// is actually the most convenient way to demonstrate the model's 
/// universality.
///
/// Uses a two zone approximation. When free particle is far from plated, 
/// large jumps are used and the dynamics are exact. When particles are close 
/// to plated, a switch is made to approximate dynamics and time becomes 
/// important. In this second zone, the particle displacement in each 
/// dimension is drawn from a truncated normal distrubution with a standard 
/// deviation of sqrt(2 dt), D = 1.
///
/// For more info on these dynamics see the @ref State class.
///
class Simulation {
public:
    // Documented in the cpp
    static const std::string kParamsFilename;

    Simulation();
    /// @brief Blank destructor.
    ///
    ~Simulation() {};

    void run_simulation();
   
private:
    void save_parallel_state() const;
    void load_parallel_state(std::string parallel_path);
    void resolve_sticking();
    std::map<std::string, std::string> read_params_file() const;
    std::pair<std::string, std::string> initialize_params();

    // Essential physical variables.
    
    /// Timestep, controls size of gaussian particle jumps when close to plated
    /// otherwise dynamics are exactly Brownian
    double dt_;
    /// The max distance a particle can jump in one timestep in one dimension,
    /// in terms of the standard deviations of the non-cutoff parent normal
    /// distribution, applies when close to plated, otherwise dynamics are 
    /// exactly Brownian
    double jump_cutoff_;

    /// The sticking probability at each contact
    double p_;
    /// size of the cluster to generate
    int cluster_size_;

    // Essential IO
    
    /// Writes restart file after every interval in terms of the number
    /// of particles that have been added, setting to 0 means no write outs
    int write_restart_interval_;
    /// Writes progress every interval in terms of the number of particles 
    /// that have been added
    int progress_interval_;

    // Internal simulation parameters 
    
    /// Cells are square/cubic, side length of each square/cube is the max 
    /// length particle can jump in one dt_ plus one diameter + epsilon
    /// ensures all collisions can be resolved
    double cell_length_;
    /// optimization parameter for kd tree
    int max_leaf_size_;
    /// if true, rejection moves are used instead of bounces
    bool rejection_only_;
    /// seed for rng
    int seed_;    

    /// random number generator
    std::mt19937 gen_;
    /// uniform [0, 1) distribution for random number generator
    Uniform uniform_;
    /// State object manages some data structures that change as simulation 
    /// advances
    State state_;

    // For parallelism

    /// How many steps to run before threads communicate and update cluster
    int parallel_interval_;
    int N_threads_;
    /// Tracks active particles for each thread
    std::vector<Vec> particles_;
    /// Tracks whether given particle is stuck for each thread
    std::vector<int> stuck_statuses_;
    /// random number generator for each thread
    std::vector<std::mt19937> gens_;
    /// uniform [0, 1) distribution for random number generator for each thread
    std::vector<Uniform> uniforms_;
};
