/// @file
///
/// @brief Implementation of the Simulation class which coordinates all
/// calculations.
///

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <cassert>
#include <sstream>

// nanoflann library for distance computations
#include "nanoflann.hpp"
#include "KDTreeVectorOfVectorsAdaptor.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "simulation.hpp"

// define globals
const int Simulation::kDims = DIMENSIONS;
const std::string Simulation::kParamsFilename = "params.txt";
typedef KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<float>>, 
				     float, Simulation::kDims> KDTree;

/// @brief Constructor starts simulation.
///
Simulation::Simulation() {
    run_simulation();
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

    step_size_ = std::stof(params_map["step_size"]);
    write_frame_interval_ = std::stoi(params_map["write_frame_interval"]);
    cluster_size_ = static_cast <int> (std::stof(params_map["cluster_size"]));
    seed_ = std::stoi(params_map["seed"]);

    float fraction_max_kappa = std::stof(params_map["fraction_max_kappa"]);
    // this line is confusing,
    // 
    // in one dimensional Brownian motion, the length of paticle step is 
    // sqrt(2 dt), where dt is the time in reduced units, (so effectively 
    // D = 1)
    // step_size_ = sqrt(2 * kDims * dt_), 
    // dt = step_size_^2 / 2 * kDims
    //
    // p_max = kappa_max * sqrt(dt) = 1.0 
    // kappa_max = p_max / sqrt(dt) = 1.0 / sqrt(dt)
    // kappa = fraction_max_kappa * kappa_max
    // p_ = kappa * sqrt(dt) = fraction_max_kappa * max_kappa * sqrt(dt)
    // p_ = fraction_kappa_max * 1.0 / sqrt(dt) * sqrt(dt)
    //
    // so finally
    p_ = fraction_max_kappa;
    std::cout << "p_ = fraction_max_kappa = " << p_ << std::endl;
    
    // kappa = p_ / sqrt(dt) 
    float kappa = p_ * std::sqrt(2 * kDims) / step_size_;
    std::cout << "da = kappa^2 = " << std::pow(kappa, 2) << std::endl;
    return restart_path;
}

/// @brief Coordinates running simulation.
///
void Simulation::run_simulation() {
    std::cout << "DIMENSION = " << kDims << std::endl;
    std::string restart_path = initialize_params();

    std::vector<std::vector<float>> plated;
    
    // optimization parameter
    int max_leaf_size = 10;
    KDTree kd_tree(kDims, plated, max_leaf_size);

    if (!restart_path.empty()) {
	std::cout << "Restarting from file " << restart_path << std::endl;
        // load_state(setup_params.restart_path);
    }
    else {
	// add single plated at origin to start
	plated.resize(1);
	plated.at(0).resize(kDims);
	for (int i = 0; i < kDims; i++) {
	    plated.at(0).at(i) = 0;
	}
	(*kd_tree.index).addPoints(0, 0);
    }

    std::vector<float> query_pt(kDims);
    for (int i = 0; i < kDims; i++) {
	query_pt.at(i) = i + 1;
    }

    // do a knn search
    const size_t num_results = 1;
    std::vector<size_t> ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<double> resultSet(num_results);

    resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
    (*kd_tree.index).findNeighbors(resultSet, &query_pt[0], 
				   nanoflann::SearchParams());

    std::cout << "knnSearch(nn="<<num_results<<"): \n";
    for (size_t i = 0; i < num_results; i++)
	std::cout << "ret_index["<<i<<"]=" << ret_indexes[i] << \
	    " out_dist_sqr=" << out_dists_sqr[i] << std::endl;
}

/// @brief Blank destructor.
///
Simulation::~Simulation() {}
