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

// define globals
const int State::kDims = DIMENSIONS;

/// @brief Writes current plated configuartion to a .xyz file.
///
void State::write_xyz() {
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

/// @brief Finds the nearest neighbor plated of particle.
///
/// Uses nanoflann implementation.
///
/// @returns squared_distance: the squared distance between particle and 
///     nearest neighbor plated
///
float State::find_nearest_neighbor() {
    const size_t k = 1;
    size_t index;
    float squared_distance;
    (*kd_tree_).query(particle_.data(), k, &index, &squared_distance);
    return squared_distance;
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
void State::save_state() {
    std::string restart_string = ("restart_" + 
				  std::to_string(plated_cloud_.size()) + 
				  ".ser");
    std::ofstream stream(restart_string);
    boost::archive::binary_oarchive archive(stream);

    std::vector<std::vector<float>> serialize_plated;
    for (auto p : plated_cloud_) {
	std::vector<float> v(p.data(), p.data() + p.size());
	serialize_plated.push_back(v);
    }
    archive & serialize_plated;

    // save state of random number generator
    stream << gen_ << std::endl;
    stream << uniform_ << std::endl;
    stream.close();
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
/// @param load_path: path to file to load
/// @param max_leaf_size: optimization parameter for kd tree
///
void State::load_state(std::string load_path, int max_leaf_size) {
    std::ifstream stream(load_path);
    boost::archive::binary_iarchive archive(stream);

    std::vector<std::vector<float>> serialize_plated;
    archive & serialize_plated;

    plated_cloud_.clear();
    radius_ = 0;
    for (auto p : serialize_plated) {
	Vec v(p.data());
	if (v.norm() > radius_) {
	    radius_ = v.norm();
	}
	plated_cloud_.push_back(v);
    }

    // load state of random number generator
    stream >> gen_;
    stream >> uniform_;
    stream.close();

    // set up kd_tree
    kd_tree_.reset(new State::KDTree(kDims, plated_cloud_, max_leaf_size));
    (*(*kd_tree_).index).addPoints(0, plated_cloud_.size() - 1);
}
