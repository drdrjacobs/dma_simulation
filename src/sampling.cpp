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

/// @brief Analytically samples position where particle on the outside of a
/// n-sphere hits the n-sphere.
///
/// Reference for 2d is:
///
/// Sander, E., Sander, L. M. & Ziff, R. M. Fractals and Fractal Correlations.
/// Comput. Phys. 8, 420 (1994).
///
/// See first-hit_distribution_in_3d.pdf that derives these formulas for 3d.
///
/// @param particle: the position vector of the particle
/// @param radius: the radius of the n-sphere that particle will hit
/// @param[out] gen: random number generator
/// @param[out] uniform: random number uniform distribution over [0, 1)
///
/// @returns result: vector on n-sphere which is a sample 
///
Vec Sampling::sample_first_hit(Vec particle, double radius,
			       std::mt19937 &gen, Uniform &uniform) {
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
	result = sample_first_hit_3d(particle, radius, gen, uniform);
    }
    return result;
}

/// @brief Analytically samples position where the particle hits 2-sphere 
/// given that it does hit the 2-sphere
///
/// See first-hit_distribution_in_3d.pdf that derives these formulas for 3d.
///
/// This assumes that the particle will hit the 2-sphere. There is also a 
/// chance the particle will escape off to infinity. This is not accounted 
/// for here.
///
/// @param particle: the position vector of the particle
/// @param radius: the radius of the 2-sphere that particle will hit
/// @param[out] gen: random number generator
/// @param[out] uniform: random number uniform distribution over [0, 1)
///
/// @returns result: vector on 2-sphere which is sample from first-hit 
///     distribution, given that particle hits the sphere
///
Vec Sampling::sample_first_hit_3d(Vec particle, double radius,
				  std::mt19937 &gen, Uniform &uniform) {
    Vec result;
    double alpha = particle.norm() / radius;
    if (uniform(gen) > 1 / alpha) {
	// particle escapes to infinity, draw new particle from uniform
	result = generate_point_on_sphere(radius, gen, uniform);
    }
    else {
	// sample first-hit distribution in 3d given that particle does hit
	// 2-sphere
	const int X = 0;
	const int Y = 1;
	const int Z = 2;
	// see references in header
	// Q is cdf
	double Q = uniform(gen);
	// eta = cos(theta)
	double eta = (std::pow(alpha, 2) + 1) / (2 * alpha);
	eta -= (std::pow((std::pow(alpha, 2) - 1), 2) / 
		(2 * alpha * std::pow(alpha - 1 + 2 * Q, 2)));
	// sometimes when alpha is small, eta ends up outside of [-1, 1]
	// fix this
	if (eta > 1) {
	    eta = 1;
	}
	else if (eta < -1) {
	    eta = -1;
	}
	double theta = std::acos(eta);
	double phi = uniform(gen) * 2 * M_PI;
	// vector based on z axis
	result(X) = radius * std::sin(theta) * std::cos(phi);
	result(Y) = radius * std::sin(theta) * std::sin(phi);
	result(Z) = radius * std::cos(theta);
	result = first_hit_3d_rotation(result, particle);
    }
    return result;
}

/// @brief Rotates hit vector to be in correct orientation with respect to 
/// the particle.
///
/// @param hit_vector: where the particle hits the 2-sphere in the reference 
///     frame where the particle lies on the z axis
/// @param particle: the position of the particle
///
/// @returns result: the rotated hit vector
///
Vec Sampling::first_hit_3d_rotation(Vec hit_vector, Vec particle) {
    Vec result = hit_vector;
    // need to rotate result so that z axis points along
    // particle vector
    Vec particle_unit_vector = particle / particle.norm();
    Vec z_unit_vector;
    z_unit_vector << 0, 0, 1;
    // in range [0, pi]
    double angle = std::acos(z_unit_vector.dot(particle_unit_vector));
    double epsilon = 1e-8;
    if (angle > epsilon and M_PI - angle > epsilon) {
        // if constexpr new c++17 feature provides compile time if 
        // this prevents Eigen static assertions from triggering for 
        // undefined operations in 2d (cross product, rotation about 
        // arbitrary axis) even though code would never actually run 
        if constexpr(kDims == 3) {
	    Vec axis = z_unit_vector.cross(particle_unit_vector);
	    // cross product is not normalized
	    axis = axis / axis.norm();
	    Eigen::AngleAxis<double> rotation(angle, axis);
	    result = rotation * result;
	}
    }
    else if (M_PI - angle < epsilon) {
        // particle is on -z axis, rotation is degenerate, choose rotation that
	// reflects across xy plane, this preserves uniform distribution in phi
        const int Z = 2;
        result(Z) = -result(Z);
    }
    // otherwise, particle is on z axis and we're done
    return result;
}
