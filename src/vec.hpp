/// @file
///
/// @brief Header that defines Eigen Vec type with appropriate dimension 
/// defined at compile time.
///

// use eigen for vector operations
#include <Eigen/Dense>

/// Global typdef of vector with correct dimensions 
typedef Eigen::Matrix<double, DIMENSIONS, 1> Vec;
