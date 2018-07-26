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
#include "utils.h"

#include "simulation.hpp"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

// include this header to serialize vectors
#include <boost/serialization/vector.hpp>

// define globals
const std::string Simulation::kParamsFilename = "params.txt";

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

    dimension_ = std::stoi(params_map["dimension"]);
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
    // step_size_ = sqrt(2 * dimension_ * dt_), 
    // dt = step_size_^2 / 2 * dimension_
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
    float kappa = p_ * std::sqrt(2 * dimension_) / step_size_;
    std::cout << "da = kappa^2 = " << std::pow(kappa, 2) << std::endl;
    return restart_path;
}

/// @brief Coordinates running simulation.
///
void Simulation::run_simulation() {
    std::string restart_path = initialize_params();

    PointCloud<float> cloud;
    // construct a kd-tree index:
    int max_leaf_size = 10;
    /*
    nanoflann::KDTreeSingleIndexDynamicAdaptor
	<nanoflann::L2_Simple_Adaptor<float, PointCloud<float>>,
	 PointCloud<float>,
	 dimension_> kd_tree(dimension_, cloud, 
			     nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size));
    */
    if (!restart_path.empty()) {
	std::cout << "Restarting from file " << restart_path << std::endl;
        // load_state(setup_params.restart_path);
    }
    else {

    }
    /*
    // Generate points:
    generateRandomPointCloud(cloud, N);

    num_t query_pt[3] = { 0.5, 0.5, 0.5 };
    // add points in chunks at a time
    size_t chunk_size = 100;
    for(size_t i = 0; i < N; i = i + chunk_size) {
	size_t end = min(size_t(i + chunk_size), N - 1);
	// Inserts all points from [i, end]
	index.addPoints(i, end);
    }

    dump_mem_usage();
    {
        // do a knn search
        const size_t num_results = 1;
        size_t ret_index;
        num_t out_dist_sqr;
	nanoflann::KNNResultSet<num_t> resultSet(num_results);
        resultSet.init(&ret_index, &out_dist_sqr );
        index.findNeighbors(resultSet, query_pt, nanoflann::SearchParams(10));

	std::cout << "knnSearch(nn="<<num_results<<"): \n";
	std::cout << "ret_index=" << ret_index << " out_dist_sqr=" << 
	    out_dist_sqr << endl;
    }
    {
        // Unsorted radius search:
        const num_t radius = 1;
	std::vector<std::pair<size_t, num_t> > indices_dists;
        RadiusResultSet<num_t, size_t> resultSet(radius, indices_dists);

        index.findNeighbors(resultSet, query_pt, nanoflann::SearchParams());

        // Get worst (furthest) point, without sorting:
	std::pair<size_t,num_t> worst_pair = resultSet.worst_item();
        cout << "Worst pair: idx=" << worst_pair.first << " dist=" << 
	    worst_pair.second << endl;
    }
    */
}

/// @brief Blank destructor.
///
Simulation::~Simulation() {}
