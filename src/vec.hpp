/// @file
///
/// @brief Header that defines Eigen Vec type with appropriate dimension 
/// defined at compile time.
///

#pragma once

// use eigen for vector operations
#include <Eigen/Dense>

#include "constants.hpp"

/// Global typdef of vector with correct dimensions 
typedef Eigen::Matrix<double, kDims, 1> Vec;
