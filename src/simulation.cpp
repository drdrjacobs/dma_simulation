/// @file
///
/// @brief Implementation of the Simulation class which coordinates all
/// calculations.
///

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>

#include "simulation.hpp"
#include "sampling.hpp"

// define globals
/// Used in collision detection loop to indicate particle has not hit any 
/// plated
const int Simulation::kNoCollision = 99999;
const std::string Simulation::kParamsFilename = "params.txt";

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
    cluster_size_ = static_cast <int> (std::stod(params_map["cluster_size"]));
    max_leaf_size_ = std::stoi(params_map["max_leaf_size"]);
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
    // cell length must be larger than maximum jump size + diameter + epsilon 
    // so that all collisions can be resolved
    double max_jump_length = std::sqrt(2 * kDims * dt_) * jump_cutoff_;
    std::cout << "max_jump_length = " << max_jump_length << std::endl;
    cell_length_ = max_jump_length + kDiameter + kSpatialEpsilon;
    std::cout << "cell_length_ = " << cell_length_ << std::endl;

    double fraction_max_kappa = std::stod(params_map["fraction_max_kappa"]);
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
    double kappa = p_ / std::sqrt(dt_);
    std::cout << "da = kappa^2 = " << std::pow(kappa, 2) << std::endl << 
	std::endl;
    
    return restart_path;
}

/// @brief Finds all contacts between particle and plated along jump 
/// trajectory.
///
/// @param jump_unit_vector: unit vector pointing along jump, accounts for 
///     precision
/// @param jump_length: length of jump along unit vector, accounts for precsion
/// @param jump: vector that discribes particle jump
///
/// @returns minimum_contact_distance: distance along jump vector at which 
///    particle makes contact with plated, kNoCollision if no contact is made
///
double Simulation::calculate_collisions(Vec jump_unit_vector, 
					double jump_length, Vec jump) {
    // minimum distance along jump vector at which particle contacts a plated
    double minimum_contact_distance = kNoCollision;
    // store details of plated that particle makes contact with
    Cells::CellIndices store_bounce_cell_indices = Cells::CellIndices();
    size_t store_bounce_plated_index = -1;
    for (int offset = 0; offset < Cells::kCellsToLoopOver; offset++) {
	Cells::CellIndices cell_indices = \
	    state_.cells_.offset_get_cell_indices(state_.particle_, offset);
	const std::vector<Vec>& cell = state_.cells_.get_cell(cell_indices);
	for (size_t i = 0; i < cell.size(); i++) {
	    // if plated is not plated that was just bounced against
	    if (cell_indices != state_.bounce_cell_indices_ ||
		i != state_.bounce_plated_index_) {
		Vec plated_r = cell[i];
		
		// check to make sure did not start too close to plated
		if ((state_.particle_ - plated_r).norm() <= 
		    kDiameter - kSpatialEpsilon) {
		    std::cout << "Particle started too close to plated!" << 
			std::endl;
		    std::exit(-1);
		}

		// Going to calculate point at which line defined by 
		// jump_vector makes contact with n-sphere defined by 
		// plated. See:                           
		//
		// https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
		Vec diff = state_.particle_ - plated_r;
		// radius is 2, not 1 here since looking at the collision 
		// surface
		int squared_radius = 4;
		double discriminant = (std::pow(jump_unit_vector.dot(diff), 2)
				       - diff.squaredNorm() + squared_radius);
		if (discriminant >= 0) {
		    // line pierces n-sphere twice, so take the closest contact
		    // ignore negative d values, means ion has to go the other
		    // way along l
		    double d1 = (-jump_unit_vector.dot(diff) + 
				 std::sqrt(discriminant));
		    double d2 = (-jump_unit_vector.dot(diff) - 
				 std::sqrt(discriminant));
		    if (d1 < 0.0) {
			d1 = kNoCollision;
		    }
		    if (d2 < 0.0) {
			d2 = kNoCollision;
		    }
		    double d = std::min(d1, d2);
		    if (d <= jump_length && d < minimum_contact_distance) {
			// potential collision since 0 < d < jump_length
			minimum_contact_distance = d;
			store_bounce_cell_indices = cell_indices;
			store_bounce_plated_index = i;
		    }
		}
		if (minimum_contact_distance == kNoCollision) {
		    // No contact is made according to the caclulation, 
		    // however, this calculation is not perfect due to 
		    // floating point precision and so may miss cases where 
		    // the ion jumps very near the boundary.
		    //
		    // For this reason, need to check the endpoint separately 
		    // to ensure exactly that no ion ever ends up less than or 
		    // equal to diameter distance away from a plated.
		    //
		    Vec new_particle_r = state_.particle_ + jump;
		    double new_center_center_distance = (new_particle_r - 
							 plated_r).norm();
		    // Add epsilon here to create a small skin. I think this 
		    // will resolve all stability issues.
		    //
		    // The skin is needed since there is no strict guarantee 
		    // that the algorithm will find the boundary if the 
		    // center-center distance is > diameter. However, 
		    // presumably if the 
		    // center-center distance > diameter + skin, there won't 
		    // be any problems.
		    if (new_center_center_distance <= 
			kDiameter + kSpatialEpsilon) {
			// contact actually was made
			minimum_contact_distance = jump_length;
			store_bounce_cell_indices = cell_indices;
			store_bounce_plated_index = i;
		    }
		}
	    }
	}
    }
    // update tracking of plated that was bounced off of
    state_.bounce_cell_indices_ = store_bounce_cell_indices;
    state_.bounce_plated_index_ = store_bounce_plated_index;
    return minimum_contact_distance;
}

/// @brief Moves particle along given jump, updates jump if bounce occurs 
/// and updates simulation state.
///
/// If particle hits plated and does not stick, does a specular reflection
/// (bounce) instead.
///
/// @param minimum_contact_distance: distance along jump vector at which 
///    particle makes contact with plated, kNoCollision if no contact is made
/// @param jump_unit_vector: unit vector pointing along jump, accounts for 
///     precision
/// @param jump_length: length of jump along unit vector, accounts for precsion
/// @param jump[out]: vector that discribes particle jump, updated with
///     new direction and length if bounce occurs
///
/// @returns stuck: true if particle sticks to plated, false otherwise
///
Simulation::Status Simulation::resolve_jump(double minimum_contact_distance,
					    Vec jump_unit_vector, 
					    double jump_length,
					    Vec &jump) {
    Status status = free;
    if (minimum_contact_distance == kNoCollision) {
	// no contact is made, take full jump
	state_.particle_ += jump;
    }
    else {
	// collision has occured
	// put particle in contact with plated
	state_.particle_ += minimum_contact_distance * jump_unit_vector;    
	if (state_.uniform_(state_.gen_) < p_) {
	    // particle stuck
	    status = stuck;
	    // add particle to cluster
	    state_.plated_cloud_.push_back(state_.particle_);
	    // add new plated in cluster to tree
	    size_t new_index = state_.plated_cloud_.size() - 1; 
	    (*(*state_.kd_tree_).index).addPoints(new_index, new_index);
	    // add newly plated to cells
	    state_.cells_.add_to_cells(state_.particle_);
	    // update cluster radius if necessary
	    double newly_plated_radius = state_.particle_.norm();
	    if (newly_plated_radius > state_.radius_) {
		state_.radius_ = newly_plated_radius;
	    }
	}
	else {
	    // bounce
	    const std::vector<Vec>& cell = \
		state_.cells_.get_cell(state_.bounce_cell_indices_);
	    Vec plated_r = cell[state_.bounce_plated_index_];
	    // unit normal to n-sphere at collision point
	    Vec unit_normal = state_.particle_ - plated_r;
	    unit_normal = unit_normal / unit_normal.norm();
	    jump = jump_unit_vector * (jump_length - minimum_contact_distance);
	    // minus sign since dot product will be negative
	    jump = jump - 2 * jump.dot(unit_normal) * unit_normal;
	    Vec new_particle = state_.particle_ + jump;
	    double new_center_center_distance =  (new_particle - 
						  plated_r).norm();
	    if (new_center_center_distance <= kDiameter + kSpatialEpsilon) {
		// new jump will not get particle out of contact with plated
		// rejection move instead
		status = rejection;
		state_.particle_ = state_.initial_particle_;    
		std::cout << "rejection!" << std::endl;
	    }
	}
    }
    return status;
}

/// @brief Take a small step forward in time.
///
/// Farms out work to series of helper methods.
///
/// @returns stuck_boolean: true if particle sticks to cluster
///
bool Simulation::step_forward() {
    Vec jump = Sampling::generate_jump(dt_, jump_cutoff_, 
				       state_.gen_, state_.uniform_);
    // set to dummy value so that gets into while loop to start
    double minimum_contact_distance = -1;
    Status status = free;
    state_.initial_particle_ = state_.particle_;
    // defualt indices, won't matter since bounce_plated_index will not match
    state_.bounce_cell_indices_ = Cells::CellIndices();
    // nonsensical index since has not bounced yet
    state_.bounce_plated_index_ = -1;
    // if particle is free, but contact was made (2nd condition) there was a 
    // bounce/reflection
    while (status == free && minimum_contact_distance != kNoCollision) {
	// unit vector along jump and jump length, account for precision 
	Vec jump_unit_vector = ((state_.particle_ + jump) - state_.particle_);
	double jump_length = jump_unit_vector.norm();
	jump_unit_vector = jump_unit_vector / jump_length;
	minimum_contact_distance = calculate_collisions(jump_unit_vector, 
							jump_length, jump);
	status = resolve_jump(minimum_contact_distance, 
			      jump_unit_vector, jump_length,
			      jump);
    }
    bool stuck_boolean = false;
    if (status == stuck) {
	stuck_boolean = true;
    }
    return stuck_boolean; 
}

/// @brief Coordinates running simulation.
///
void Simulation::run_simulation() {
    // get simulattion parameters
    std::cout << "DIMENSIONS = " << kDims << std::endl;
    std::string restart_path = initialize_params();
    if (restart_path.empty()) {
	state_.set_up_new_state(cell_length_, max_leaf_size_, seed_);
    }
    else {
	state_.load_state(cell_length_, max_leaf_size_, restart_path);
    }
    
    // main loop, already start with cluster size of 1
    for (int i = state_.plated_cloud_.size() + 1; i <= cluster_size_; i++) {
	// tracks whether current particle has stuck
	bool stuck = false;
	// radius at which new particle should be generated at
	// need 2 times epsilon, due to epsilon skin around particles
	double generation_radius = (state_.radius_ + kDiameter + 
				   2 * kSpatialEpsilon);
	// position vector of diffusing particle
	state_.particle_ = Sampling::generate_point_on_sphere(generation_radius, 
							      state_.gen_,
							      state_.uniform_);
	while (!stuck) {
	    // check to see if there are any plated in current cell or 
	    // surrounding cells
	    if (state_.cells_.has_neighbors(state_.particle_)) {
		// take a step forward in time since a collision is 
		// possible
		stuck = step_forward();
	    }
	    else {
		double squared_distance = state_.find_nearest_neighbor();
		// need 2 times epsilon, due to epsilon skin around particles
		double jump_length = (std::sqrt(squared_distance) - kDiameter 
				      - 2 * kSpatialEpsilon);
		Vec jump = Sampling::generate_point_on_sphere(jump_length, 
							      state_.gen_,
							      state_.uniform_);
		state_.particle_ += jump;
	    }
	    // see if particle crossed boundary, if it did regenerate it
	    if (!stuck && state_.particle_.norm() > 
		generation_radius + kSpatialEpsilon) {
		state_.particle_ = Sampling::sample_first_hit(state_.particle_,
							      generation_radius,
							      state_.gen_, 
							      state_.uniform_);
	    }
	}
	if (i % write_frame_interval_ == 0) {
	    std::cout << "N_plated = " << i << std::endl;
	    state_.write_xyz();
	}
    }
    state_.write_xyz();
    state_.save_state();
    std::cout << "Done." << std::endl;
}
