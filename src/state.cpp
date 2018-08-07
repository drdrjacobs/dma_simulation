/// @file
///
/// @brief Implementation of the State class which stores data that change as 
/// simulation advances.
///

#include <iostream>
#include <cmath>
#include <fstream>

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
