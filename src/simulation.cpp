/// @file
///
/// @brief Implementation of the Simulation class which coordinates all
/// calculations.
///

#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <experimental/filesystem>

#include "simulation.hpp"
#include "sampling.hpp"

// define globals
const std::string Simulation::kParamsFilename = "params.txt";

/// @brief Constructor setus up simulation for run.
///                                   
Simulation::Simulation() {
    // get simulattion parameters
    std::cout << "DIMENSIONS = " << kDims << std::endl;
    std::string restart_path = initialize_params();
    if (restart_path.empty()) {
	state_.set_up_new_state(cell_length_, max_leaf_size_, seed_, 
				rejection_only_);
    }
    else {
	state_.load_state(cell_length_, max_leaf_size_, restart_path, 
			  rejection_only_);
	// make copy of restart file so that it is not deleted later
	auto option = \
	    std::experimental::filesystem::copy_options::overwrite_existing;
	std::experimental::filesystem::copy_file(restart_path, 
						 "save_" + restart_path, 
						 option);
						 
    }
}

/// @brief Runs simulation starting from current state.
///
void Simulation::run_simulation() { 
    // i is number of particle that is being added
    for (int i = state_.get_cluster_size() + 1; i <= cluster_size_; i++) {
	// tracks whether current particle has stuck
	bool stuck = false;
	state_.add_new_particle();
	while (!stuck) {
	    // check to see if there are any plated in current cell or 
	    // surrounding cells
	    if (state_.has_neighbors()) {
		// take a step forward in time since a collision is 
		// possible
		stuck = state_.take_small_step(dt_, jump_cutoff_, p_);
	    }
	    else {
		// nearest plated is far away, take large step
		state_.take_large_step();
	    }
	    // regenerate particle from first-hit distribution if it gets too 
	    // far away from cluster
	    if (!stuck) {
		state_.check_for_regeneration(dt_);
	    }
	}
	if (write_restart_interval_ > 0 && i % write_restart_interval_ == 0) {
	    state_.save_state();
	    // remove file from two intervals ago
	    int remove = i - 2 * write_restart_interval_;
	    std::string remove_path = ("restart_" + std::to_string(remove) + 
				       ".ser");
	    std::remove(remove_path.c_str());
	}
	if (i % progress_interval_ == 0) {
	    std::cout << "N_plated = " << i << std::endl;
	}
    }
    state_.write_xyz();
    state_.save_state();
    std::cout << "Done." << std::endl;
}

/// @brief Extracts parameters from params file.
///
/// @returns params_map: map that contains params std::string key value pairs
///
std::map<std::string, std::string> Simulation::read_params_file() const {
    std::map<std::string, std::string> params_map;
    
    std::ifstream params_file;
    params_file.open(kParamsFilename);
    if (params_file.fail()) {
	std::cout << "File " << kParamsFilename << " not found!" << std::endl;
	exit(-1);
    }
    std::string line;
    std::string type;
    std::string name;
    std::string dummy;
    std::string value;
    // param format is "type name = value"
    // ignore type, just for reference, parse as strings for now
    // order in file does not matter
    while (std::getline(params_file, line)) {
	if (!line.empty() && line[0] != '#') {
	    std::istringstream iss(line);
	    iss >> type >> name >> dummy >> value;
	    params_map[name] = value;
	}
    }
    std::cout << "# Params" << std::endl;
    for (const auto &k_v_pair : params_map ) {
        std::cout << "# " << k_v_pair.first << " = " << k_v_pair.second <<
	    std::endl;
    }
    std::cout << std::endl;
    return params_map;
}

/// @brief Defines simulation parameters in member variables. 
/// 
/// @returns initial_N_ions: number of ions initially in box (sei + bath)
///
std::string Simulation::initialize_params() {
    // parse everything from file as std::strings
    std::map<std::string, std::string> params_map = read_params_file();

    std::string restart_path;
    if (params_map.find("restart_path") != params_map.end()) {
	restart_path = params_map["restart_path"];
    }
    else {
	restart_path = "";
    }

    write_restart_interval_ = \
	static_cast <int> (std::stod(params_map["write_restart_interval"]));
    progress_interval_ = \
	static_cast <int> (std::stod(params_map["progress_interval"]));
    cluster_size_ = static_cast <int> (std::stod(params_map["cluster_size"]));
    max_leaf_size_ = std::stoi(params_map["max_leaf_size"]);
    rejection_only_ = std::stoi(params_map["rejection_only"]);
    seed_ = std::stoi(params_map["seed"]);

    // set expected length of particle movement in 2d or 3d
    double rms_jump_size = std::stod(params_map["rms_jump_size"]);
    // in one dimensional Brownian motion, the rms length of paticle step is 
    // sqrt(2 * D * dt_), in reduced units the D disappears 
    // so
    // rms_jump_size = sqrt(2 * kDims * dt_)
    dt_ = std::pow(rms_jump_size, 2) / (2 * kDims);

    // set jump_cutoff_ for cell based spatial parallelism
    jump_cutoff_ = std::stod(params_map["jump_cutoff"]);
    std::cout << "variance ratio (truncated to non-truncated normal) = " <<
	Sampling::calculate_variance_ratio(jump_cutoff_) << std::endl;
    // cell length must be larger than maximum jump size + diameter + epsilon 
    // so that all collisions can be resolved
    double max_jump_length = std::sqrt(2 * kDims * dt_) * jump_cutoff_;
    std::cout << "max_jump_length = " << max_jump_length << std::endl;
    cell_length_ = max_jump_length + kDiameter + kSpatialEpsilon;
    std::cout << "cell_length_ = " << cell_length_ << std::endl;

    if (params_map.find("kappa") != params_map.end()) {
	double kappa = std::stod(params_map["kappa"]);
	p_ = kappa * std::sqrt(dt_);
	if (p_ > 1.0) {
	    std::string error = "sticking probability, ";
	    error += "p = kappa * sqrt(dt) = " + std::to_string(p_);
	    error += " > 1.0, use lower kappa or dt (rms_jump_size).";
	    throw std::invalid_argument(error);	       
	}
    }
    else {
	double fraction_max_kappa = \
	    std::stod(params_map["fraction_max_kappa"]);
	// this line is confusing,
	// 
	// p_max = kappa_max * sqrt(dt) = 1.0 
	// kappa_max = p_max / sqrt(dt) = 1.0 / sqrt(dt)
	// p_ = kappa * sqrt(dt) = fraction_max_kappa * max_kappa * sqrt(dt)
	// p_ = fraction_kappa_max * 1.0 / sqrt(dt) * sqrt(dt)
	//
	// so finally
	p_ = fraction_max_kappa;
    }
    std::cout << "p_ = fraction_max_kappa = " << p_ << std::endl;
    
    // kappa = p_ / sqrt(dt) 
    double kappa = p_ / std::sqrt(dt_);
    std::cout << "da = kappa = " << kappa << std::endl << std::endl;
    
    return restart_path;
}
