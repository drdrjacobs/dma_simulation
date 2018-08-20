/// @file
///
/// @brief Implementation of the State class which stores data that change as 
/// simulation advances.
///

#include <iostream>
#include <cmath>
#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "state.hpp"
#include "sampling.hpp"

// define globals
/// Used in collision detection loop to indicate particle has not hit any
/// plated
const int State::kNoCollision = 99999;

/// @brief Sets up new state with single plated at origin.
///
/// @param cell_length: the length of each cell
/// @param max_leaf_size: optimization parameter for kd tree
/// @param seed: seed for state random number generator
///
void State::set_up_new_state(double cell_length, int max_leaf_size, int seed) {
    // set up cells
    cells_.set_up_cells(cell_length);

    // single plated at origin to start
    plated_cloud_.resize(1);
    for (int i = 0; i < kDims; i++) {
	plated_cloud_.at(0)(i) = 0;
    }
    // assign smart pointer
    kd_tree_.reset(new KDTree(kDims, plated_cloud_, max_leaf_size));
    (*(*kd_tree_).index).addPoints(0, 0);
    radius_ = 0;
    
    cells_.add_to_cells(plated_cloud_.at(0));
    
    // set up rng
    gen_ = std::mt19937(seed);
    uniform_ = Uniform(0.0, 1.0);
}

/// @brief Loads exact serialized state of a simulation. 
///
/// Serialization is run through std::vector and boost binary_iarchives.
///
/// Requires the simulation parameters to be set exactly the same as when the
/// state was saved to work!
///
/// Note that running two trajectories starting from the same serialized state
/// should produce exactly the same results.
///
/// @param cell_length: the length of each cell
/// @param max_leaf_size: optimization parameter for kd tree
/// @param load_path: path to file to load
///
void State::load_state(double cell_length, int max_leaf_size, 
		       std::string load_path) {
    std::ifstream stream(load_path);
    boost::archive::binary_iarchive archive(stream);

    std::vector<std::vector<double>> serialize_plated;
    archive & serialize_plated;

    // set up cells
    cells_.set_up_cells(cell_length);
    plated_cloud_.clear();
    radius_ = 0;
    for (auto p : serialize_plated) {
	Vec v(p.data());
	if (v.norm() > radius_) {
	    radius_ = v.norm();
	}
	plated_cloud_.push_back(v);
	// reconstruct cells
	cells_.add_to_cells(v);
    }

    // load state of random number generator
    stream >> gen_;
    stream >> uniform_;
    stream.close();

    // set up kd_tree
    kd_tree_.reset(new State::KDTree(kDims, plated_cloud_, max_leaf_size));
    (*(*kd_tree_).index).addPoints(0, plated_cloud_.size() - 1);
}

/// @brief Serialization writes out exact state of simulation.
///
/// Serialization is run through std::vector and boost binary_oarchive.
///
/// Requires the params file to be set exactly the same to work upon loading.
///
/// Note that running two trajectories starting from the same serialized state
/// should produce exactly the same results.
///
void State::save_state() const {
    std::string restart_string = ("restart_" + 
				  std::to_string(plated_cloud_.size()) + 
				  ".ser");
    std::ofstream stream(restart_string);
    boost::archive::binary_oarchive archive(stream);

    std::vector<std::vector<double>> serialize_plated;
    for (auto p : plated_cloud_) {
	std::vector<double> v(p.data(), p.data() + kDims);
	serialize_plated.push_back(v);
    }
    archive & serialize_plated;

    // save state of random number generator
    stream << gen_ << std::endl;
    stream << uniform_ << std::endl;
    stream.close();
}

/// @brief Adds new particle to simulation that starts just outside cluster 
/// radius.
///
void State::add_new_particle() {
    // radius at which new particle should be generated at need 2 times 
    // epsilon, due to epsilon skin around particles           
    double generation_radius = (radius_ + kDiameter + 2 * kSpatialEpsilon);
    // position vector of diffusing particle
    particle_ = Sampling::generate_point_on_sphere(generation_radius,
						   gen_, uniform_);
}

/// @brief Checks to see if there are neigboring plated around particle.
///
/// @returns result: true if there are plated in the particle's cell or the 
///     cells around the particle's cell, false otherwise
///
bool State::has_neighbors() const {
    bool result = cells_.has_neighbors(particle_);
    return result;
}

/// @brief Gets distance particle would have to travel along jump_unit_vector
///     to hit plated.
///
/// @param particle: position of particle
/// @param plated_r: the position of the plated
/// @param jump_unit_vector: unit vector pointing along jump, accounts for 
///     precision
/// @param jump_length: length of jump
/// @param jump: vector that discribes particle jump
///
/// @returns collision_distance: distance particle would have to travel along 
///     jump_unit_vector to hit plated, kNoCollision if this is not possible,
///     distance is negative, or distance > jump_length
///
double get_collision_distance(Vec particle, Vec plated_r, 
			      Vec jump_unit_vector, double jump_length,
			      Vec jump) {
    // check to make sure did not start too close to plated
    if ((particle - plated_r).norm() <= kDiameter - kSpatialEpsilon) {
	std::cout << "Particle started too close to plated!" << std::endl;
	std::exit(-1);
    }
    // Going to calculate point at which line defined by jump_vector makes 
    // contact with n-sphere defined by plated. See:                           
    //
    // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    Vec diff = particle - plated_r;
    // radius is 2, not 1 here since looking at the collision surface
    int squared_radius = 4;
    double discriminant = (std::pow(jump_unit_vector.dot(diff), 2) 
			   - diff.squaredNorm() + squared_radius);
    double collision_distance = State::kNoCollision;
    if (discriminant >= 0) {
	// line pierces n-sphere twice, so take the closest contact ignore 
	// negative d values, means particle has to go the other way along
	// jump_unit_vector
	double d1 = (-jump_unit_vector.dot(diff) + 
		     std::sqrt(discriminant));
	double d2 = (-jump_unit_vector.dot(diff) - 
		     std::sqrt(discriminant));
	if (d1 < 0.0) {
	    d1 = State::kNoCollision;
	}
	if (d2 < 0.0) {
	    d2 = State::kNoCollision;
	}
	double d = std::min(d1, d2);
	if (d <= jump_length) {
	    collision_distance = d;
	}
    }
    if (collision_distance == State::kNoCollision) {
	// No contact is made according to the caclulation, however, this 
	// calculation is not perfect due to floating point precision and 
	// may miss cases where the particle jumps very near the boundary.
	//
	// For this reason, need to check the endpoint separately to ensure 
	// exactly that particle never ends up too close to plated.
	Vec new_particle_r = particle + jump;
	double new_center_center_distance = (new_particle_r - plated_r).norm();
	// Add epsilon here to create a small skin. The skin is needed since 
	// there is no strict guarantee that algorithm will find the boundary 
	// if the center-center distance is > diameter. However, presumably if 
	// the center-center distance > diameter + skin, there won't be any 
	// problems.
	if (new_center_center_distance <= kDiameter + kSpatialEpsilon) {
	    // contact actually was made
	    collision_distance = jump_length;
	}
    }
    return collision_distance;
}

/// @brief Finds all contacts between particle and plated along jump.
///
/// @param jump_unit_vector: unit vector pointing along jump, accounts for 
///     precision
/// @param jump_length: length of jump along unit vector, accounts for precsion
/// @param jump: vector that discribes particle jump
/// @param[out] bounce_cell_indices: indices of cell containing plated that 
///     particle just bounced off of
/// @param[out] bounce_plated_index: index of plated that particle just 
///     bounced off of
///
/// @returns minimum_collision_distance: distance along jump vector at which 
///    particle makes contact with plated, kNoCollision if no contact is made
///
double State::check_collisions_loop(Vec jump_unit_vector, 
				    double jump_length, Vec jump,
				    CellIndices &bounce_cell_indices,
				    size_t &bounce_plated_index) const {
    // minimum distance along jump vector at which particle contacts a plated
    double minimum_collision_distance = kNoCollision;
    // store details of plated that particle makes contact with
    CellIndices store_bounce_cell_indices = CellIndices();
    size_t store_bounce_plated_index = -1;
    for (int offset = 0; offset < Cells::kCellsToLoopOver; offset++) {
	CellIndices cell_indices = cells_.offset_get_cell_indices(particle_, 
								  offset);
	const std::vector<Vec>& cell = cells_.get_cell(cell_indices);
	for (size_t i = 0; i < cell.size(); i++) {
	    // if plated is not plated that was just bounced against
	    if (cell_indices != bounce_cell_indices ||
		i != bounce_plated_index) {
		Vec plated_r = cell[i];
		double collision_distance = \
		    get_collision_distance(particle_, plated_r,
					   jump_unit_vector, jump_length, 
					   jump);
		if (collision_distance < minimum_collision_distance) { 
		    minimum_collision_distance = collision_distance;
		    store_bounce_cell_indices = cell_indices;
		    store_bounce_plated_index = i;
		}
	    }
	}
    }
    // update tracking of plated that was bounced off of
    bounce_cell_indices = store_bounce_cell_indices;
    bounce_plated_index = store_bounce_plated_index;
    return minimum_collision_distance;
}

/// @brief Tranforms particle into plated at its current position.
///
void State::stick_particle() {
    // add particle to cluster
    plated_cloud_.push_back(particle_);
    // add new plated in cluster to tree
    size_t new_index = plated_cloud_.size() - 1; 
    (*(*kd_tree_).index).addPoints(new_index, new_index);
    // add newly plated to cells
    cells_.add_to_cells(particle_);
    // update cluster radius if necessary
    double newly_plated_radius = particle_.norm();
    if (newly_plated_radius > radius_) {
	radius_ = newly_plated_radius;
    }
}

/// @brief Attempts to bounce particle off plated and updates jump if succeeds.
///
/// Rejection move is used if bounce fails so particle position does not 
/// change. Failure occurs if endpoint of bounce does not seperate plated and 
/// particle.
///
/// @param minimum_collision_distance: distance along jump vector at which 
///    particle makes contact with plated, kNoCollision if no contact is made
/// @param bounce_cell_indices: indices of cell containing plated that 
///     particle just bounced off of
/// @param bounce_plated_index: index of plated that particle just 
///     bounced off of
/// @param initial_particle: initial position of particle before step
/// @param jump_unit_vector: unit vector pointing along jump, accounts for 
///     precision
/// @param jump_length: length of jump along unit vector, accounts for precsion
/// @param jump[out]: vector that discribes particle jump, updated with
///     new direction and length if bounce occurs
///
/// @returns status: free if bounce succeeds, rejection otherwise
///
State::Status State::attempt_bounce(double minimum_collision_distance,
				    CellIndices bounce_cell_indices,
				    size_t bounce_plated_index,
				    Vec initial_particle,
				    Vec jump_unit_vector, double jump_length,
				    Vec &jump) {
    Status status = free;
    const std::vector<Vec>& cell = cells_.get_cell(bounce_cell_indices);
    Vec plated_r = cell[bounce_plated_index];
    // unit normal to n-sphere at collision point
    Vec unit_normal = particle_ - plated_r;
    unit_normal = unit_normal / unit_normal.norm();
    jump = (jump_unit_vector * 
	    (jump_length - minimum_collision_distance));
    // minus sign since dot product will be negative
    jump = jump - 2 * jump.dot(unit_normal) * unit_normal;
    Vec new_particle = particle_ + jump;
    double new_center_center_distance = (new_particle - plated_r).norm();
    if (new_center_center_distance <= kDiameter + kSpatialEpsilon) {
	// new jump will not get particle out of contact with plated
	// rejection move instead
	status = rejection;
	particle_ = initial_particle;    
	std::cout << "rejection!" << std::endl;
    }
    return status;
}

/// @brief Resolves whether particle sticks or bounces (or rejects) and 
/// updates jump (or particle) appropriately.
///
/// Has side effects!
///
/// If particle hits plated and does not stick, does a specular reflection
/// (bounce) instead. In rare cases rejection move will be used so particle
/// position will not change.
///
/// @param p: the sticking probability
/// @param minimum_collision_distance: distance along jump vector at which 
///    particle makes contact with plated, kNoCollision if no contact is made
/// @param bounce_cell_indices: indices of cell containing plated that 
///     particle just bounced off of
/// @param bounce_plated_index: index of plated that particle just 
///     bounced off of
/// @param initial_particle: initial position of particle before step
/// @param jump_unit_vector: unit vector pointing along jump, accounts for 
///     precision
/// @param jump_length: length of jump along unit vector, accounts for precsion
/// @param jump[out]: vector that discribes particle jump, updated with
///     new direction and length if bounce occurs
///
/// @returns stuck: true if particle sticks to plated, false otherwise
///
State::Status State::resolve_stick_or_bounce(double p, 
					     double minimum_collision_distance,
					     CellIndices bounce_cell_indices,
					     size_t bounce_plated_index,
					     Vec initial_particle,
					     Vec jump_unit_vector, double 
					     jump_length,
					     Vec &jump) {
    Status status = free;
    if (minimum_collision_distance == kNoCollision) {
	// no contact is made, take full jump
	particle_ += jump;
    }
    else {
	// collision has occured
	// put particle in contact with plated
	particle_ += minimum_collision_distance * jump_unit_vector;    
	if (uniform_(gen_) < p) {
	    // particle stuck
	    status = stuck;
	    stick_particle();
	}
	else {
	    // bounce, assigns new jump if successful or turns into rejection
	    // move
	    status = attempt_bounce(minimum_collision_distance, 
				    bounce_cell_indices, bounce_plated_index,
				    initial_particle,
				    jump_unit_vector, jump_length, jump);
	}
    }
    return status;
}

/// @brief Move particle according to jump using approximate dynamics to 
/// resolve collisions.
///
/// Has side effects!
///
/// @param jump: the jump the particle is taking
/// @param p: the sticking probability
///
/// @returns stuck_boolean: true if particle sticks to cluster
///
bool State::resolve_jump(Vec jump, double p) {
    // set to dummy value so that gets into while loop to start
    double minimum_collision_distance = -1;
    Status status = free;
    // initial position of particle before step
    Vec initial_particle = particle_;
    // indices of cell containing plated that particle just bounced off of
    // default indices, won't matter since bounce_plated_index will not match
    CellIndices bounce_cell_indices = CellIndices();
    /// index of plated that particle just bounced off of
    // nonsensical index since has not bounced yet
    size_t bounce_plated_index = -1;

    // if particle is free, but contact was made (2nd condition) there was a 
    // bounce/reflection
    while (status == free && minimum_collision_distance != kNoCollision) {
	// unit vector along jump and jump length, account for precision 
	Vec jump_unit_vector = ((particle_ + jump) - particle_);
	double jump_length = jump_unit_vector.norm();
	jump_unit_vector = jump_unit_vector / jump_length;
	// const method
	minimum_collision_distance = \
	    check_collisions_loop(jump_unit_vector, 
				  jump_length, jump,
				  bounce_cell_indices,
				  bounce_plated_index);
	// updates state so all side effects are in here
	status = resolve_stick_or_bounce(p, minimum_collision_distance, 
					 bounce_cell_indices, 
					 bounce_plated_index,
					 initial_particle,
					 jump_unit_vector, jump_length,
					 jump);
    }
    bool stuck_boolean = false;
    if (status == stuck) {
	stuck_boolean = true;
    }
    return stuck_boolean; 
}

/// @brief Take a small step forward in time using approximate dynamics to
/// resolve collisions.
///
/// @param dt: the time step, determines the jump distribution standard
///     deviation
/// @param jump_cutoff: max distance a particle can jump in one timestep in
///     one dimension, in terms of the standard deviations of the non-cutoff
///     parent normal distribution
/// @param p: the sticking probability
///
/// @returns stuck: true if particle sticks to cluster
///
bool State::take_small_step(double dt, double jump_cutoff, double p) {
    Vec jump = Sampling::generate_jump(dt, jump_cutoff, gen_, uniform_);
    bool stuck = resolve_jump(jump, p);
    return stuck;
}


/// @brief Finds the nearest neighbor plated of particle.
///
/// Uses nanoflann implementation.
///
/// @returns squared_distance: the squared distance between particle and 
///     nearest neighbor plated
///
double State::find_nearest_neighbor() const {
    const size_t k = 1;
    size_t index;
    double squared_distance;
    (*kd_tree_).query(particle_.data(), k, &index, &squared_distance);
    return squared_distance;
}

/// @brief takes large (exact dynamics) step forward.
///
/// Finds nearest neighbor and takes large step on n-sphere such that particle
/// will just barely avoid a collision.
///
void State::take_large_step() {
    double squared_distance = find_nearest_neighbor();
    // need 2 times epsilon, due to epsilon skin around particles   
    double jump_length = (std::sqrt(squared_distance) - kDiameter 
			  - 2 * kSpatialEpsilon);
    Vec jump = Sampling::generate_point_on_sphere(jump_length,
						  gen_, uniform_);
    particle_ += jump;
}

/// @brief Checks to see if particle is too far from cluster and regenerates it
/// from first-hit distribution if it is.
///
void State::check_for_regeneration() {
    double generation_radius = (radius_ + kDiameter + 2 * kSpatialEpsilon);
    if (particle_.norm() > generation_radius + kSpatialEpsilon) {
	particle_ = Sampling::sample_first_hit(particle_, generation_radius,
					       gen_, uniform_);
    }
}

/// @brief Finds overlaps between plated.
///
/// @returns count: the number of pairwise overlaps found
///
int State::check_overlaps() const {
    int count = 0;
    // pairwise distance of less than this threshold is an overlap
    double threshold = kDiameter - kSpatialEpsilon;
    // for each central cell
    for (std::pair<CellIndices, std::vector<Vec>> element : 
	     cells_.get_cell_map()) {
	CellIndices central_cell_indices = element.first;
	const std::vector<Vec>& central_cell = element.second;
	// loop through all cells around central cell
	for (int offset = 0; offset < Cells::kCellsToLoopOver; offset++) {
	    CellIndices other_cell_indices = \
		cells_.offset_get_cell_indices(central_cell_indices, offset);
	    const std::vector<Vec>& other_cell = \
		cells_.get_cell(other_cell_indices);
	    // for each plated in central cell
	    for (size_t i = 0; i < central_cell.size(); i++) {
		Vec r_i = central_cell.at(i);
		int start = 0;
		// avoid counting self as overlap
		if (central_cell_indices == other_cell_indices) {
		    start = i + 1;
		}
		// for each plated in other cell
		for (size_t j = start; j < other_cell.size(); j++) {
		    Vec r_j = other_cell.at(j);
		    double distance = (r_i - r_j).norm();
		    if (distance <= threshold) {
			count += 1;
			std::cout << distance << std::endl;
		    }
		}
	    }
	}
    }
    return count;
}

/// @brief Writes current plated configuartion to a .xyz file.
///
void State::write_xyz() const {
    int count = check_overlaps();
    if (count > 0) {
	std::cout << "Found overlaps between plated!" << std::endl;
	std::exit(-1);
    }

    int N_plated = plated_cloud_.size();
    std::ofstream xyz_file;
    std::string path = ("frame_" + std::to_string(N_plated) + ".xyz");
    xyz_file.open(path);
    // header is total number of plated
    xyz_file << N_plated << "\n\n";
    // ions                                                                 
    for (int i = 0; i < N_plated; i++) {
	// plated are type N
	xyz_file << "N";
	for (int j = 0; j < kDims; j++) {
	    xyz_file  << " " << plated_cloud_[i][j];
	}
	if (kDims == 2) {
	    // always write out in 3d for ovito compatibility
	    xyz_file  << " 0";
	}
	// fourth column is particle number, like an effective plating time
	// included for compatibility with dendrite simulation .xyz format
	xyz_file << " " << i << "\n";
    }
    xyz_file.close();
}
