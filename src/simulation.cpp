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

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "simulation.hpp"

// define globals
const int Simulation::kDims = DIMENSIONS;
const std::string Simulation::kParamsFilename = "params.txt";
/// small epsilon used in distance calculations, float so added f
const float Simulation::kSpatialEpsilon = 0.001f;

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

    write_frame_interval_ = std::stoi(params_map["write_frame_interval"]);
    cluster_size_ = static_cast <int> (std::stof(params_map["cluster_size"]));
    max_leaf_size_ = std::stoi(params_map["max_leaf_size"]);
    seed_ = std::stoi(params_map["seed"]);

    // set expected length of particle movement in 2d or 3d
    float rms_jump_size = std::stof(params_map["rms_jump_size"]);
    // in one dimensional Brownian motion, the rms length of paticle step is 
    // sqrt(2 * D * dt_), in reduced units the D disappears 
    // so
    // rms_jump_size = sqrt(2 * kDims * dt_)
    dt_ = std::pow(rms_jump_size, 2) / (2 * kDims);

    // set jump_cutoff_ for cell based spatial parallelism
    jump_cutoff_ = std::stof(params_map["jump_cutoff"]);
    // cell length must be larger than maximum jump size + 2 * diameter + 
    // epsilon so that all collisions can be resolved
    int diameter = 2;
    float max_jump_length = std::sqrt(2 * kDims * dt_) * jump_cutoff_;
    std::cout << "max_jump_length = " << max_jump_length << std::endl;
    cell_length_ = max_jump_length + diameter + kSpatialEpsilon;
    std::cout << "cell_length_ = " << cell_length_ << std::endl;

    float fraction_max_kappa = std::stof(params_map["fraction_max_kappa"]);
    // this line is confusing,
    // 
    // p_max = kappa_max * sqrt(dt) = 1.0 
    // kappa_max = p_max / sqrt(dt) = 1.0 / sqrt(dt)
    // p_ = kappa * sqrt(dt) = fraction_max_kappa * max_kappa * sqrt(dt)
    // p_ = fraction_kappa_max * 1.0 / sqrt(dt) * sqrt(dt)
    //
    // so finally
    p_ = fraction_max_kappa;
    std::cout << "p_ = fraction_max_kappa = " << p_ << std::endl;
    
    // kappa = p_ / sqrt(dt) 
    float kappa = p_ / std::sqrt(dt_);
    std::cout << "da = kappa^2 = " << std::pow(kappa, 2) << std::endl;
    
    return restart_path;
}

/// @brief Sets up state object of simulation.
///
/// @param restart_path: path of file that may contain saved state
///
void Simulation::set_up_state(std::string restart_path) {
    state_.kd_tree.reset(new State::KDTree(kDims, state_.plated_cloud, 
					   max_leaf_size_));
    if (!restart_path.empty()) {
	std::cout << "Restarting from file " << restart_path << std::endl;
        // load_state(setup_params.restart_path);
    }
    else {
        // add single plated at origin to start
        state_.plated_cloud.resize(1);
        for (int i = 0; i < kDims; i++) {
            state_.plated_cloud.at(0)(i) = 0;
        }
        (*(*state_.kd_tree).index).addPoints(0, 0);
	state_.radius = 0;

	// set up rng
	state_.gen = std::mt19937(seed_);
	state_.distribution = State::Uniform(0.0, 1.0);
    }
}

/// @brief Finds a random point on a ball with given radius and dimension.
///
/// @param kDims: the dimensions of the ball
/// @param radius: the radius of the ball
/// @param[out] gen: random number generator
/// @param[out] distribution: random number uniform distribution over [0, 1)
///
/// @returns result: vector representing random point
///
Simulation::Vec Simulation::generate_point_on_ball(int kDims, float radius,
						   std::mt19937 & gen,
						   State::Uniform & dist) {
    Vec result;
    float theta = 2.0 * M_PI * dist(gen);
    if (kDims == 2) {
	float x = std::cos(theta);
	float y = std::sin(theta);
	result << x, y;
	result *= radius;
    }
    else if (kDims == 3) {
	// need random point on sphere, see: 
	// http://mathworld.wolfram.com/SpherePointPicking.html
	float cos_phi = 2.0 * dist(gen) - 1.0;
	float x = std::sqrt(1 - std::pow(cos_phi, 2)) * std::cos(theta);
	float y = std::sqrt(1 - std::pow(cos_phi, 2)) * std::sin(theta);
	float z = cos_phi;
	result << x, y, z;
	result *= radius;
    }
    return result;
}

/// Transforms jump from uniform [0, 1) into a truncated normal distribution.
///
/// Truncated normal distribution is implemented using the inverse transform
/// method. Truncation occurs after certain number of standard deviations. Mean
/// is 0.
///
/// See: https://en.wikipedia.org/wiki/Truncated_normal_distribution
///
/// @param standard_deviation: standard deviation of parent non-truncated
///     normal distribution
/// @param cutoff: truncate after this many standard deviations of the
///     non-truncated parent normal distribution.
/// @param[out] v: vector of uniform [0, 1) representing particle's jump
/// 
void truncated_normal_transform(float standard_deviation, float cutoff,
				Simulation::Vec & v) {
    boost::math::normal_distribution<float> normal;
    float phi = boost::math::cdf(normal, cutoff);
    float tmp;
    for (int i = 0; i < v.size(); i++) {
	// use property that phi(-x) = 1 - phi(x) 
	tmp = v[i] * (2.0 * phi - 1.0);
	tmp = (1.0 - phi) + tmp;
	v[i] = boost::math::quantile(normal, tmp) * standard_deviation;
    }
}

/// @brief Take a small step forward in time.
///
/// @returns stuck: true if particle sticks to cluster
///
bool Simulation::step_forward() {
    Vec jump;
    for (int i = 0; i < Simulation::kDims; i++) {
	jump(i) = state_.distribution(state_.gen);
    }
    // standard deviation of Brownian jump in one dimension
    float standard_deviation = std::sqrt(2 * dt_);
    truncated_normal_transform(standard_deviation, jump_cutoff_, jump);

    // find neighbors within max interaction length, cell_length_^2
    // neighbors has one pair for each neighbor
    // first element in pair is index of neighbor in plated_cloud
    // second element is squared distance between particle
    std::vector<std::pair<size_t, float>> neighbors;
    typedef nanoflann::RadiusResultSet<float, size_t> RadiusResultSet;
    RadiusResultSet result_set(std::pow(cell_length_, 2), neighbors);
    (*(*state_.kd_tree).index).findNeighbors(result_set, 
					     state_.particle.data(), 
					     nanoflann::SearchParams());
    bool stuck = false;
    // minimum distance along jump vector at which particle contacts a plated
    float minimum_contact_distance = 999;
    // unit vector along jump, account for precision 
    Vec l = ((state_.particle + jump) - state_.particle);
    float jump_length = l.norm();
    l = l / jump_length;
    for (auto n : neighbors) {
	Vec plated_r = state_.plated_cloud[n.first];
	// Going to calculate point at which line defined by 
	// jump_vector makes contact with sphere defined by 
	// plated. See:                           
	//
	// https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
	Vec diff = state_.particle - plated_r;
	// radius is 2, not 1 here since looking at the collision surface
	int squared_radius = 4;
	float discriminant = (std::pow(l.dot(diff), 2) - 
			      diff.squaredNorm() + squared_radius);
	if (discriminant >= 0) {
	    // line pierces sphere twice, so take the closest contact
            // ignore negative d values, means ion has to go the other
	    // way along l
            float d1 = -l.dot(diff) + std::sqrt(discriminant);
            float d2 = -l.dot(diff) - std::sqrt(discriminant);
            if (d1 < 0.0) {
                d1 = 999;
            }
            if (d2 < 0.0) {
                d2 = 999;
            }
            float d = std::min(d1, d2);
            if (d <= jump_length && d < minimum_contact_distance) {
		// collision occured since 0 < d < jump_length
                minimum_contact_distance = d;
	    }
	}
    }
    return stuck;
}

/// @brief Coordinates running simulation.
///
void Simulation::run_simulation() {
    // get simulattion parameters
    std::cout << "DIMENSION = " << kDims << std::endl;
    std::string restart_path = initialize_params();
    set_up_state(restart_path);
    
    // tracks whether current particle has stuck
    bool stuck = false;
    int diameter = 2;
    // main loop
    for (int i = 0; i < cluster_size_; i++) {
	// radius at which new particle should be generated at
	float generation_radius = state_.radius + diameter + kSpatialEpsilon;
	// position vector of diffusing particle
	state_.particle = generate_point_on_ball(kDims, generation_radius, 
						 state_.gen, 
						 state_.distribution);
	std::cout << "particle:\n" << state_.particle << std::endl;
	while (!stuck) {
	    // do a k nearest neighbors search with k = 1
	    const size_t k = 1;
	    size_t index;
	    float squared_distance;
	    int diameter = 2;
	    // cell_length_ = max_jump_distance + diameter + epsilon is
	    // max interaction length within one step
	    float squared_cell_length = std::pow(cell_length_, 2);

	    (*state_.kd_tree).query(state_.particle.data(), k, &index, 
				    &squared_distance);

	    if (squared_distance > squared_cell_length) {
		// particle is farther than max interaction length from nearest
		// plated so take as large a step as possible using exact
		// Brownian dynamics
		float jump_length = (std::sqrt(squared_distance) - diameter 
				     - kSpatialEpsilon);
		Vec jump = generate_point_on_ball(kDims, jump_length, 
						  state_.gen,
						  state_.distribution);
		state_.particle += jump;
		std::cout << "particle:\n" << state_.particle << std::endl;
	    }
	    else {
		// take a step forward in time since a collision is possible
		stuck = step_forward();
	    }
	    // see if particle crossed boundary, if it did regenerate it
	    if (state_.particle.norm() > 3 * generation_radius) {
		state_.particle = generate_point_on_ball(kDims, 
							 generation_radius, 
							 state_.gen, 
							 state_.distribution);
	    }
	    exit(0);
	}
    }
}

/// @brief Blank destructor.
///
Simulation::~Simulation() {}
