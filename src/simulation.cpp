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
#include "omp.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "simulation.hpp"
#include "sampling.hpp"

// define globals
const std::string Simulation::kParamsFilename = "params.txt";

/// @brief Constructor setus up simulation for run.
///                                   
Simulation::Simulation() {
    // get simulattion parameters
    std::cout << "DIMENSIONS = " << kDims << std::endl;
    auto paths = initialize_params();
    std::string restart_path = paths.first;
    std::string parallel_path = paths.second;

    /// parallel infrastructure
    #pragma omp parallel 
    {
	int tid = omp_get_thread_num();
	if (tid == 0) {
	    N_threads_ = omp_get_num_threads();
	}
    }    
    std::cout << "N_threads_ = " << N_threads_ << std::endl << std::endl;

    if (restart_path.empty()) {
	std::cout << "Must always provide serial file to start from!" << 
	    std::endl;
	std::exit(0);
    }
    else {
	// load serial state
	state_.load_state(cell_length_, max_leaf_size_, restart_path, 
			  rejection_only_);
	// make copy of restart file so that it is not deleted later
        auto option = \
            std::experimental::filesystem::copy_options::overwrite_existing;
	std::experimental::filesystem::copy_file(restart_path,
                                                 "save_" + restart_path,
                                                 option);
	if (parallel_path.empty()) {
	    // restarting from serial state
	    gen_ = std::mt19937(seed_);
	    uniform_ = Uniform(0.0, 1.0);
 
	    bool stuck = false;
	    for (int i = 0; i < N_threads_; i++) {
		particles_.push_back(state_.create_new_particle(gen_, 
								uniform_));
		stuck_statuses_.push_back(stuck);
		int seed = gen_();
		gens_.push_back(std::mt19937(seed));
		uniforms_.push_back(Uniform(0.0, 1.0));
	    }
	}
	else {
	    load_parallel_state(parallel_path);
	}
    }
}

/// @brief Serialization writes out exact state of parallel simulation.
///
/// Serialization is run through std::vector and boost binary_oarchive.
///
/// Requires the params file and number of threads to be set exactly the same 
/// to work upon loading.
///
/// Note that running two trajectories starting from the same serialized state
/// should produce exactly the same results.
///
void Simulation::save_parallel_state() const {
    std::string restart_string = ("restart_" +
                                  std::to_string(state_.get_cluster_size()) +
                                  ".prl");
    std::ofstream stream(restart_string);
    boost::archive::binary_oarchive archive(stream);
    std::vector<std::vector<double>> serialize_particles;
    for (auto p : particles_) {
	std::vector<double> v(p.data(), p.data() + kDims);
        serialize_particles.push_back(v);
    }
    archive & serialize_particles;
    archive & stuck_statuses_;

    stream << N_threads_ << std::endl;
    stream << gen_ << std::endl;
    stream << uniform_ << std::endl;
    for (auto g : gens_) {
	stream << g << std::endl;
    }
    for (auto u : uniforms_) {
	stream << u << std::endl;
    }
}

/// @brief Loads exact serialized parallel state of a simulation.
///
/// Serialization is run through std::vector and boost binary_iarchives.
///
/// Requires the simulation parameters and number of threads to be set exactly 
/// the same as when the state was saved to work!
///
/// Note that running two trajectories starting from the same serialized state
/// should produce exactly the same results.
///
/// @param parallel_path: file to load from
///
void Simulation::load_parallel_state(std::string parallel_path) {
    std::cout << "Restarting from parallel state." << std::endl;

    std::ifstream stream(parallel_path);
    boost::archive::binary_iarchive archive(stream);
    std::vector<std::vector<double>> serialize_particles;
    archive & serialize_particles;
    for (auto p : serialize_particles) {
        Vec v(p.data());
	particles_.push_back(v);
    }
    archive & stuck_statuses_;

    int saved_N_threads;
    stream >> saved_N_threads;
    if (N_threads_ != saved_N_threads) {
	std::cout << "Must use " << saved_N_threads << 
	    " threads for parallel restart!" << std::endl;
	std::exit(0);
    }
    stream >> gen_;
    stream >> uniform_;
    std::mt19937 tmp_gen;
    Uniform tmp_uniform;
    for (int i = 0; i < N_threads_; i++) {
	stream >> tmp_gen;
	gens_.push_back(tmp_gen);
    }
    for (int i = 0; i < N_threads_; i++) {
	stream >> tmp_uniform;
	uniforms_.push_back(tmp_uniform);
    }
    // make copy of restart file so that it is not deleted later
    auto option = \
	std::experimental::filesystem::copy_options::overwrite_existing;
    std::experimental::filesystem::copy_file(parallel_path,
					     "save_" + parallel_path,
					     option);
}

/// @brief Runs simulation starting from current state.
///
void Simulation::run_simulation() { 
    int progress_counter = state_.get_cluster_size() / progress_interval_;
    int restart_counter = -1;
    int one_ago_save = state_.get_cluster_size();
    int two_ago_save = -1;
    if (write_restart_interval_ > 0) {
	restart_counter = state_.get_cluster_size() / write_restart_interval_;
    }
    while (state_.get_cluster_size() < cluster_size_) {
	// Run in parallel for a while, if cluster is large, sticking 
	// probability is small and time not too long, bias should be minimal
        #pragma omp parallel 
	{	
	    int id = omp_get_thread_num();
	    // counter tracks number of steps (large or small) taken
	    int counter = 0;
	    while (!stuck_statuses_[id] && counter < parallel_interval_) {
		// check to see if there are any plated in current cell or 
		// surrounding cells
		if (state_.has_neighbors(particles_[id])) {
		    // take a step forward in time since a collision is 
		    // possible
		    stuck_statuses_[id] =				\
			state_.take_small_step(dt_, jump_cutoff_, p_, 
					       particles_[id], gens_[id], 
					       uniforms_[id]);
		}
		else {
		    // nearest plated is far away, take large step
		    state_.take_large_step(particles_[id], gens_[id],
					   uniforms_[id]);
		}
		// regenerate particle from first-hit distribution if it gets 
		// too far away from cluster
		if (!stuck_statuses_[id]) {
		    state_.check_for_regeneration(dt_, particles_[id], 
						  gens_[id], uniforms_[id]);
		}
		counter += 1;
	    }
	}
	// add particles to cluster while resolving overlaps
	resolve_sticking();
	int tmp_progress = state_.get_cluster_size() / progress_interval_;
	if (tmp_progress > progress_counter) {
	    progress_counter = tmp_progress;
	    std::cout << "N_plated = " << tmp_progress * progress_interval_ << 
		std::endl;
	}
	if (write_restart_interval_ > 0) {
	    int tmp_restart = (state_.get_cluster_size() / 
			       write_restart_interval_);
	    if (tmp_restart > restart_counter) {
		restart_counter = tmp_restart;
		state_.save_state();
		save_parallel_state();
		// remove files from two intervals ago
		std::string remove_path = ("restart_" + 
					   std::to_string(two_ago_save) + 
					   ".ser");
		std::remove(remove_path.c_str());
		remove_path = ("restart_" + std::to_string(two_ago_save) + 
			       ".prl");
		std::remove(remove_path.c_str());
		two_ago_save = one_ago_save;
		one_ago_save = state_.get_cluster_size();
	    }
	}
    }
    state_.write_xyz();
    std::cout << "Done." << std::endl;
}

/// @brief Adds stuck particles to cluster after removing overlaps.
///
void Simulation::resolve_sticking() { 
    // first remove overlaps between plated particles
    // expected to be very rare
    for (int i = 0; i < N_threads_; i++) {
	if (stuck_statuses_[i]) {
	    for (int j = i + 1; j < N_threads_; j++) {
		if (stuck_statuses_[j]) {
		    double distance = (particles_[i] - particles_[j]).norm();
		    if (distance < kDiameter + kSpatialEpsilon) {
			stuck_statuses_[j] = false;
			//std::cout << "Stuck overlap!" << std::endl;
		    }
		}
	    }
	    // now add to cluster
	    state_.stick_particle(particles_[i]);
	}
    }
    // reset any overlaps with non-plated
    for (int i = 0; i < N_threads_; i++) {
	if (!stuck_statuses_[i]) {
	    for (int j = 0; j < N_threads_; j++) {
		if (stuck_statuses_[j] and i != j) {
		    double distance = (particles_[i] - particles_[j]).norm();
		    if (distance < kDiameter + kSpatialEpsilon) {
			particles_[i] = state_.create_new_particle(gen_, 
								   uniform_);
			//std::cout << "Free overlap!" << std::endl;
		    }
		}
	    }
	}
    }
    // finally switch stuck particles to new particles
    for (int i = 0; i < N_threads_; i++) {
	if (stuck_statuses_[i]) {
	    particles_[i] = state_.create_new_particle(gen_, uniform_);
	    stuck_statuses_[i] = false;
	}
    }
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
std::pair<std::string, std::string> Simulation::initialize_params() {
    // parse everything from file as std::strings
    std::map<std::string, std::string> params_map = read_params_file();

    std::string restart_path;
    if (params_map.find("restart_path") != params_map.end()) {
	restart_path = params_map["restart_path"];
    }
    else {
	restart_path = "";
    }
    std::string parallel_path;
    if (params_map.find("parallel_path") != params_map.end()) {
	parallel_path = params_map["parallel_path"];
    }
    else {
	parallel_path = "";
    }
    std::pair<std::string, std::string> paths = std::make_pair(restart_path, 
							       parallel_path);

    write_restart_interval_ = \
	static_cast <int> (std::stod(params_map["write_restart_interval"]));
    progress_interval_ = \
	static_cast <int> (std::stod(params_map["progress_interval"]));
    cluster_size_ = static_cast <int> (std::stod(params_map["cluster_size"]));
    max_leaf_size_ = std::stoi(params_map["max_leaf_size"]);
    rejection_only_ = std::stoi(params_map["rejection_only"]);
    parallel_interval_ = std::stoi(params_map["parallel_interval"]);
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
    
    return paths;
}
