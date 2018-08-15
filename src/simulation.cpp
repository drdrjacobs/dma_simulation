/// @file
///
/// @brief Implementation of the Simulation class which coordinates all
/// calculations.
///

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <sstream>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>

#include "simulation.hpp"

// define globals
const int Simulation::kDims = DIMENSIONS;
const int Simulation::kDiameter = 2;
const std::string Simulation::kParamsFilename = "params.txt";
/// small epsilon used in distance calculations
const double Simulation::kSpatialEpsilon = 1e-7;

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

/// @brief Sets up state object of simulation.
///
/// @param restart_path: path of file that may contain saved state
///
void Simulation::set_up_state(std::string restart_path) {
    state_.cells_.set_up_cells(cell_length_);
    if (!restart_path.empty()) {
	std::cout << "Restarting from file " << restart_path << std::endl;
        state_.load_state(restart_path, max_leaf_size_);
    }
    else {
	// start new simulation,
        // add single plated at origin to start
        state_.plated_cloud_.resize(1);
        for (int i = 0; i < kDims; i++) {
            state_.plated_cloud_.at(0)(i) = 0;
        }
	state_.kd_tree_.reset(new State::KDTree(kDims, state_.plated_cloud_, 
						max_leaf_size_));
        (*(*state_.kd_tree_).index).addPoints(0, 0);
	state_.radius_ = 0;

	state_.cells_.add_to_cells(state_.plated_cloud_.at(0));

	// set up rng
	state_.gen_ = std::mt19937(seed_);
	state_.uniform_ = State::Uniform(0.0, 1.0);
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
Vec Simulation::generate_point_on_ball(int kDims, double radius, 
				       std::mt19937 &gen, 
				       State::Uniform &uniform) {
    Vec result;
    double theta = 2.0 * M_PI * uniform(gen);
    if (kDims == 2) {
	double x = std::cos(theta);
	double y = std::sin(theta);
	result << x, y;
	result *= radius;
    }
    else if (kDims == 3) {
	// need random point on sphere, see: 
	// http://mathworld.wolfram.com/SpherePointPicking.html
	double cos_phi = 2.0 * uniform(gen) - 1.0;
	double x = std::sqrt(1 - std::pow(cos_phi, 2)) * std::cos(theta);
	double y = std::sqrt(1 - std::pow(cos_phi, 2)) * std::sin(theta);
	double z = cos_phi;
	result << x, y, z;
	result *= radius;
    }
    return result;
}

/// @brief Draws Brownian jump from truncated normal distribution with 
///     appropriate width.
///
/// Truncated normal distribution is implemented using the inverse transform
/// method. Truncation occurs after certain number of standard deviations.
///
/// See: https://en.wikipedia.org/wiki/Truncated_normal_distribution
///
/// @returns jump: jump drawn from truncated normal distribution 
///
Vec Simulation::generate_jump() {
    Vec jump;
    for (int i = 0; i < Simulation::kDims; i++) {
	jump(i) = state_.uniform_(state_.gen_);
    }
    // standard deviation of parent Brownian non-truncated normal distribution
    // in one dimension
    double standard_deviation = std::sqrt(2 * dt_);
    // truncate after jump_cutoff_ standard deviations of the
    // non-truncated parent normal distribution.
    //
    // See: https://en.wikipedia.org/wiki/Truncated_normal_distribution
    boost::math::normal_distribution<double> normal;
    double phi = boost::math::cdf(normal, jump_cutoff_);
    double tmp;
    for (int i = 0; i < kDims; i++) {
	// use property that phi(-x) = 1 - phi(x) 
	tmp = jump[i] * (2.0 * phi - 1.0);
	tmp = (1.0 - phi) + tmp;
	jump[i] = boost::math::quantile(normal, tmp) * standard_deviation;
    }
    return jump;
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
///    particle makes contact with plated, inf if no contact is made
///
double Simulation::calculate_collisions(Vec jump_unit_vector, 
					double jump_length, Vec jump) {
    // minimum distance along jump vector at which particle contacts a plated
    double minimum_contact_distance = std::numeric_limits<double>::infinity();
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
		// jump_vector makes contact with sphere defined by 
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
		    // line pierces sphere twice, so take the closest contact
		    // ignore negative d values, means ion has to go the other
		    // way along l
		    double d1 = (-jump_unit_vector.dot(diff) + 
				 std::sqrt(discriminant));
		    double d2 = (-jump_unit_vector.dot(diff) - 
				 std::sqrt(discriminant));
		    if (d1 < 0.0) {
			d1 = std::numeric_limits<double>::infinity();
		    }
		    if (d2 < 0.0) {
			d2 = std::numeric_limits<double>::infinity();
		    }
		    double d = std::min(d1, d2);
		    if (d <= jump_length && d < minimum_contact_distance) {
			// potential collision since 0 < d < jump_length
			minimum_contact_distance = d;
			store_bounce_cell_indices = cell_indices;
			store_bounce_plated_index = i;
		    }
		}
		if (std::isinf(minimum_contact_distance)) {
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
///    particle makes contact with plated, inf if no contact is made
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
    if (std::isinf(minimum_contact_distance)) {
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
	    // unit normal to sphere at collision point
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
    Vec jump = generate_jump();
    // set to dummy value so that gets into while loop to start
    double minimum_contact_distance = -1;
    Status status = free;
    state_.initial_particle_ = state_.particle_;
    // defualt indices, won't matter since bounce_plated_index will not match
    state_.bounce_cell_indices_ = Cells::CellIndices();
    // nonsensical index since has not bounced yet
    state_.bounce_plated_index_ = -1;
    // minimum_contact_distance is inf if no contact is made
    // so !isinf means contact was made
    while (status == free && !std::isinf(minimum_contact_distance)) {
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

/// @brief Analytically samples position where particle on the outside of a
/// ball hits the ball.
///
/// Reference for 2d is:
///
/// Sander, E., Sander, L. M. & Ziff, R. M. Fractals and Fractal Correlations.
/// Comput. Phys. 8, 420 (1994).
///
/// See first-hit_distribution_in_3d.pdf that derives these formulas for 3d.
///
/// @param kDims: the dimensions of the ball
/// @param particle: the position vector of the particle
/// @param radius: the radius of the ball that particle will hit
/// @param[out] gen: random number generator
/// @param[out] distribution: random number uniform distribution over [0, 1)
///
/// @returns result: vector on ball which is a sample 
///
Vec Simulation::sample_first_hit(int kDims, Vec particle, double radius,
				 std::mt19937 &gen, 
				 State::Uniform &uniform) {
    Vec result;
    if (kDims == 2) {
	const int X = 0;
	const int Y = 1;
	// see references in header
	double a = particle.norm();
	double rho = uniform(gen);
	double nu = (a - radius) / (a + radius) * std::tan(M_PI * rho);
	result(X) = (1 - std::pow(nu, 2)) * particle(X) - 2 * nu * particle(Y);
	result(X) = radius / a * result(X) / (1 + std::pow(nu, 2));
	result(Y) = (1 - std::pow(nu, 2)) * particle(Y) + 2 * nu * particle(X);
	result(Y) = radius / a * result(Y) / (1 + std::pow(nu, 2));
    }
    else if (kDims == 3) {
	// see references in header
	double alpha = particle.norm() / radius;
	if (uniform(gen) > 1 / alpha) {
	    // particle escapes to infinity, draw new particle from uniform
	    result = generate_point_on_ball(kDims, radius, 
					    state_.gen_, 
					    state_.uniform_);
	}
	else {
	    // sample first-hit distribution in 3d
	    result = sample_first_hit_3d(particle, radius, state_.gen_, 
					 state_.uniform_);
	}
    }
    return result;
}

/// @brief Analytically samples position where the particle hits sphere given 
/// that it does hit the sphere
///
/// See first-hit_distribution_in_3d.pdf that derives these formulas for 3d.
///
/// This assumes that the particle will hit the sphere. There is also a chance 
/// the particle will escape off to infinity. This is not accounted for here.
///
/// @param particle: the position vector of the particle
/// @param radius: the radius of the sphere that particle will hit
/// @param[out] gen: random number generator
/// @param[out] distribution: random number uniform distribution over [0, 1)
///
/// @returns result: vector on sphere which is sample from first-hit 
///     distribution, given that particle hits the sphere
///
Vec Simulation::sample_first_hit_3d(Vec particle, double radius,
				    std::mt19937 &gen, 
				    State::Uniform &uniform) {
    Vec result;
    const int X = 0;
    const int Y = 1;
    const int Z = 2;
    // see references in header
    double alpha = particle.norm() / radius;
    // Q is cdf
    double Q = uniform(gen);
    // eta = cos(theta)
    double eta = (std::pow(alpha, 2) + 1) / (2 * alpha);
    eta -= (std::pow((std::pow(alpha, 2) - 1), 2) / 
	    (2 * alpha * std::pow(alpha - 1 + 2 * Q, 2)));
    double theta = std::acos(eta);
    double phi = uniform(gen) * 2 * M_PI;
    // vector based on z axis
    result(X) = radius * std::sin(theta) * std::cos(phi);
    result(Y) = radius * std::sin(theta) * std::sin(phi);
    result(Z) = radius * std::cos(theta);
    result = first_hit_3d_rotation(result, particle);
    return result;
}

/// @brief Rotates hit vector to be in correct orientation with respect to 
/// the particle.
///
/// @param hit_vector: where the particle hits the sphere in the reference 
///     frame where the particle lies on the z axis
/// @param particle: the position of the particle
///
/// @returns result: the rotated hit vector
///
Vec Simulation::first_hit_3d_rotation(Vec hit_vector, Vec particle) {
    Vec result = hit_vector;
    // need to rotate result so that z axis points along
    // particle vector
    Vec particle_unit_vector = particle / particle.norm();
    Vec z_unit_vector;
    z_unit_vector << 0, 0, 1;
    double angle = std::acos(z_unit_vector.dot(particle_unit_vector));
    double epsilon = 1e-8;
    if (angle > epsilon) {
	// compiler macro prevents Eigen static assertions from
	// triggering for undefined cross product in 2d even though
	// code would never actually run
        # if DIMENSIONS == 3
	    Vec axis =  z_unit_vector.cross(particle_unit_vector);
	    // cross product is not normalized
	    axis = axis / axis.norm();
	    Eigen::AngleAxis<double> rotation(angle, axis);
	    result = rotation * result;
        # endif
    }
    // otherwise, particle is on z axis and we're done
    return result;
}

/// @brief Coordinates running simulation.
///
void Simulation::run_simulation() {
    // get simulattion parameters
    std::cout << "DIMENSION = " << kDims << std::endl;
    std::string restart_path = initialize_params();
    set_up_state(restart_path);
    
    // main loop, already start with cluster size of 1
    for (int i = state_.plated_cloud_.size() + 1; i <= cluster_size_; i++) {
	// tracks whether current particle has stuck
	bool stuck = false;
	// radius at which new particle should be generated at
	// need 2 times epsilon, due to epsilon skin around particles
	double generation_radius = (state_.radius_ + kDiameter + 
				   2 * kSpatialEpsilon);
	// position vector of diffusing particle
	state_.particle_ = generate_point_on_ball(kDims, generation_radius, 
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
		Vec jump = generate_point_on_ball(kDims, jump_length, 
						  state_.gen_,
						  state_.uniform_);
		state_.particle_ += jump;
	    }
	    // see if particle crossed boundary, if it did regenerate it
	    if (!stuck && state_.particle_.norm() > 
		generation_radius + kSpatialEpsilon) {
		state_.particle_ = sample_first_hit(kDims, state_.particle_,
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

/// @brief Blank destructor.
///
Simulation::~Simulation() {}
