/// @file
///
/// @brief Header that defines global constants and typedefs used throughout 
/// the simulation.
///

#pragma once

// use eigen for vector operations
#include <Eigen/Dense>

#include <string>
#include <random>

// define globals
const int kDims = DIMENSIONS;
const int kDiameter = 2;
/// small epsilon used in distance calculations
const double kSpatialEpsilon = 1e-8;

// typedefs
/// Global typdef of vector with correct dimensions
typedef Eigen::Matrix<double, kDims, 1> Vec;
typedef std::uniform_real_distribution<double> Uniform;
