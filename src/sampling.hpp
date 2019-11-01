/// @file
///
/// @brief Header for the Sampling namespace that includes functionality to 
/// draw from different probability distributions.
///

#pragma once

#include <random>

#include "constants.hpp"

/// @brief Sampling namespace includes functionality to draw from various 
/// probability distributions.
namespace Sampling {
    Vec generate_point_on_sphere(double radius,
				 std::mt19937 &gen, Uniform &uniform);
    Vec generate_point_on_plane(double height, double L,
				std::mt19937 &gen, Uniform &uniform);    
    Vec generate_jump(double dt, double jump_cutoff,
		      std::mt19937 &gen, Uniform &uniform);
    double calculate_variance_ratio(double jump_cutoff);
    Vec sample_first_hit(Vec particle, double height,
			 std::mt19937 &gen, Uniform &uniform);
}
