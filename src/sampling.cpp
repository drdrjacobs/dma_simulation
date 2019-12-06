/// @file
///
/// @brief Implementation of the Sampling namespace for sampling different 
/// probability distributions.
///

#include <iostream>
#include <cmath>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>

#include "sampling.hpp"

/// @brief Finds a random point on an n-sphere with given radius and dimension.
///
/// @param radius: the radius of the n-sphere
/// @param[out] gen: random number generator
/// @param[out] uniform: random number uniform distribution over [0, 1)
///
/// @returns result: vector representing random point
///
Vec Sampling::generate_point_on_sphere(double radius,
                                       std::mt19937 &gen, Uniform &uniform) {
    Vec result;
    double theta = 2.0 * M_PI * uniform(gen);
    if (kDims == 2) {
        // random point on circle
        double x = std::cos(theta);
        double y = std::sin(theta);
        result << x, y;
        result *= radius;
    }
    else if (kDims == 3) {
        // need random point on 2-sphere, see:
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

/// @brief Finds a random point on plane (or line) within pbcs
///
/// @param height: the height of the plane (y coordinate)
/// @param L: dimensions of the pbcs in the x (and z if 3d) coordinate(s)
///     pbcs run from [-L/2, L/2)
/// @param[out] gen: random number generator
/// @param[out] uniform: random number uniform distribution over [0, 1)
///
/// @returns result: vector representing random point
///
Vec Sampling::generate_point_on_plane(double height, double L,
				      std::mt19937 &gen, Uniform &uniform) {
    Vec result;
    const int X = 0;
    const int Y = 1;
    result[X] = L * uniform(gen) - L / 2.0;
    result[Y] = height;
    if (kDims == 3) {
	const int Z = 2;
        result[Z] = L * uniform(gen) - L / 2.0;
    }
    return result;
}

/// @brief Draws Brownian jump from truncated normal distribution with 
///     given width.
///
/// Truncated normal distribution is implemented using the inverse transform
/// method. Truncation occurs after certain number of standard deviations.
///
/// See: https://en.wikipedia.org/wiki/Truncated_normal_distribution
///
/// @param dt: the time step, determines the jump distribution standard 
///     deviation
/// @param jump_cutoff: max distance a particle can jump in one timestep in 
///     one dimension, in terms of the standard deviations of the non-cutoff 
///     parent normal distribution
/// @param[out] gen: random number generator
/// @param[out] uniform: random number uniform distribution over [0, 1)
///
/// @returns jump: jump drawn from truncated normal distribution 
///
Vec Sampling::generate_jump(double dt, double jump_cutoff, 
			     std::mt19937 &gen, Uniform &uniform) {
    Vec jump;
    for (int i = 0; i < kDims; i++) {
	jump(i) = uniform(gen);
    }
    // standard deviation of parent Brownian non-truncated normal distribution
    // in one dimension
    double standard_deviation = std::sqrt(2 * dt);
    // truncate after jump_cutoff standard deviations of the non-truncated 
    // parent normal distribution.
    //
    // See: https://en.wikipedia.org/wiki/Truncated_normal_distribution
    boost::math::normal_distribution<double> normal;
    double phi = boost::math::cdf(normal, jump_cutoff);
    double tmp;
    for (int i = 0; i < kDims; i++) {
	// use property that phi(-x) = 1 - phi(x) 
	tmp = jump[i] * (2.0 * phi - 1.0);
	tmp = (1.0 - phi) + tmp;
	jump[i] = boost::math::quantile(normal, tmp) * standard_deviation;
    }
    return jump;
}

/// @brief Calculates ratio of truncated to non-truncated normal distribution
/// variances.
///
/// @param jump_cutoff: max distance a particle can jump in one timestep in 
///     one dimension, in terms of the standard deviations of the non-cutoff 
///     parent normal distribution
///
/// @returns ratio: ratio of truncated to non-truncated normal distribution 
///     variances
///
double Sampling::calculate_variance_ratio(double jump_cutoff) {
    // truncate after jump_cutoff standard deviations of the non-truncated 
    // parent normal distribution.
    //
    // See: https://en.wikipedia.org/wiki/Truncated_normal_distribution
    double beta = jump_cutoff;
    double alpha = -jump_cutoff;
    boost::math::normal_distribution<double> normal;
    // note capital letter cdf
    double Phi_beta = boost::math::cdf(normal, jump_cutoff);
    // use property that Phi(-x) = 1 - Phi(x) 
    double Z = 2.0 * Phi_beta - 1.0;
    // lower case letter pdf
    double phi_beta = boost::math::pdf(normal, beta);
    double phi_alpha = boost::math::pdf(normal, alpha);
    double ratio = 1 + (alpha * phi_alpha - beta * phi_beta) / Z;
    ratio += -std::pow((phi_alpha - phi_beta) / Z, 2);
    return ratio;
}

/// @brief Analytically samples position where particle above absorbing line or
/// plane hits.
///
/// See first-hit_distributions.pdf that derives these formulas.
///
/// @param particle: the position vector of the particle
/// @param height: the height (y coordinate) of the line or plane that 
///     particle will hit
/// @param[out] gen: random number generator
/// @param[out] uniform: random number uniform distribution over [0, 1)
///
/// @returns result: vector to point on plane or line where particle hits
///
Vec Sampling::sample_first_hit(Vec particle, double height,
			       std::mt19937 &gen, Uniform &uniform) {
    Vec result;
    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    double height_above = particle[Y] - height;
    result[Y] = height;
    // inverse transform method
    double cdf = uniform(gen);    
    if (kDims == 2) {
	double nondimensional_delta_x = std::tan(M_PI * (cdf - 0.5));
	double delta_x = nondimensional_delta_x * height_above;
	result[X] = particle[X] + delta_x;
    }
    else if (kDims == 3) {
	double nondimensional_radius = cdf / std::sqrt(1 - std::pow(cdf, 2));
	double radius = nondimensional_radius * height_above;
	// polar angle sampled uniformly
	double theta = 2 * M_PI * uniform(gen);
	result[X] = particle[X] + std::cos(theta) * radius;
	result[Z] = particle[Z] + std::sin(theta) * radius;
    }
    return result;
}
