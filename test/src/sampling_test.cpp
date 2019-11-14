/// @file
///
/// @brief Tests cells methods.
///

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/filesystem.hpp>

#include "../lib/catch.hpp"
#include "../../src/constants.hpp"
#include "../../src/sampling.hpp"

TEST_CASE("test generate_point_on_sphere", "[sampling]") {
    int seed = 0;
    std::mt19937 gen(seed);
    Uniform uniform(0.0, 1.0);
    int samples;
    std::vector<double> radii = {1, 4.5, 20};
    for (auto radius : radii) {
	std::string tag = ", radius = " + std::to_string(radius);
	SECTION( "test radius is what was specified" + tag) {
	    double epsilon = 1e-8;
	    Approx approx_zero = Approx(0).margin(epsilon);
	    samples = 1000;
	    bool error = false;

	    for (int i = 0; i < samples; i++) {
		Vec sample = Sampling::generate_point_on_sphere(radius, 
								gen, uniform);
		if (!(sample.norm() - radius == approx_zero)) {
		    error = true;
		}
	    }
	    REQUIRE(!error);
	}
    }
    // spit out samples to be analyzed with python
    double radius = 1.0;
    samples = 100000;
    std::ofstream file;
    std::string dir = "test_distributions/"; 
    if (!boost::filesystem::is_directory(dir)) {
	std::cout << "Can only run sampling tests from path test/" << 
	    std::endl;
	exit(-1);
    }
    std::string path = (dir + "/generate_point_on_sphere_" + 
			std::to_string(kDims) + "d.txt");
    file.open(path);
    for (int i = 0; i < samples; i++) {
	Vec sample = Sampling::generate_point_on_sphere(radius, 
							gen, uniform);
	for (int j = 0; j < kDims; j++) {
	    file << sample(j) << " ";
	}
	file << std::endl;
    }
    file.close();
}

TEST_CASE("test generate_point_on_plane", "[sampling]") {
    int seed = 0;
    std::mt19937 gen(seed);
    Uniform uniform(0.0, 1.0);
    int samples;
    std::vector<double> heights = {10, 4.5, 20};
    std::vector<double> Ls = {3.0, 8.5, 16.5};
    const int Y = 1;
    for (auto height : heights) {
	for (auto L : Ls) {
	    std::string tag = (", height = " + std::to_string(height) + 
			       ", L = "  + std::to_string(L));
	    SECTION( "test height is what was specified" + tag) {
		double epsilon = 1e-8;
		Approx approx_zero = Approx(0).margin(epsilon);
		samples = 1000;
		bool error = false;

		for (int i = 0; i < samples; i++) {
		    Vec sample = Sampling::generate_point_on_plane(height, L,
								   gen, 
								   uniform);
		    if (!(sample(Y) - height == approx_zero)) {
			error = true;
		    }
		    for (int j = 0; j < kDims; j++) {
			if (j != Y && (sample(j) < -L / 2.0 || 
				       sample(j) >= L / 2.0)) {
			    error = true;
			}
		    }
		}
		REQUIRE(!error);
	    }
	}
    }

    // spit out samples to be analyzed with python
    std::ofstream file;
    std::string dir = "test_distributions/"; 
    if (!boost::filesystem::is_directory(dir)) {
	std::cout << "Can only run sampling tests from path test/" << 
	    std::endl;
	exit(-1);
    }
    file << std::scientific << std::setprecision(12);
    heights = {10, 4.5};
    Ls = {3.0, 8.5};
    for (auto height : heights) {
	for (auto L : Ls) {
	    samples = 100000;
	    std::string path = (dir + "/generate_point_on_plane_" + 
				std::to_string(kDims) + "d_h_" + 
				std::to_string(height) + "_L_" + 
				std::to_string(L) + ".txt");
	    file.open(path);
	    for (int i = 0; i < samples; i++) {
		Vec sample = Sampling::generate_point_on_plane(height, L, 
							       gen, uniform);
		for (int j = 0; j < kDims; j++) {
		    file << sample(j) << " ";
		}
		file << std::endl;
	    }
	    file.close();
	}
    }
}

TEST_CASE("test generate_jump", "[sampling]") {
    int seed = 0;
    std::mt19937 gen(seed);
    Uniform uniform(0.0, 1.0);
    double dt = 0.1;
    int samples;
    std::vector<double> jump_cutoffs = {0.5, 1.0, 2};
    for (auto jump_cutoff : jump_cutoffs) {
	std::string tag = ", jump_cutoff = " + std::to_string(jump_cutoff);
	bool error = false;
	SECTION( "test distribution is properly bounded" + tag) {
	    samples = 1000;
	    for (int i = 0; i < samples; i++) {
		Vec sample = Sampling::generate_jump(dt, jump_cutoff,
						     gen, uniform);
		for (int j = 0; j < kDims; j++) {
		    if (sample[j] / std::sqrt(2 * dt) > jump_cutoff ||
			sample[j] / std::sqrt(2 * dt) < -jump_cutoff) {
			error = true;
		    }
		}
	    }
	    REQUIRE(!error);
	}

	// spit out samples to be analyzed with python
	samples = 100000;
	std::ofstream file;
	std::string dir = "test_distributions/"; 
	if (!boost::filesystem::is_directory(dir)) {
	    std::cout << "Can only run sampling tests from path test/" << 
		std::endl;
	    exit(-1);
	}
	std::string path = (dir + "/generate_jump_" + 
			    std::to_string(jump_cutoff) + "_" +
			    std::to_string(kDims) + "d.txt");
	file.open(path);
	file << "variance = " << 
	    Sampling::calculate_variance_ratio(jump_cutoff) << std::endl;
	for (int i = 0; i < samples; i++) {
	    Vec sample = Sampling::generate_jump(dt, jump_cutoff, 
						 gen, uniform);
	    for (int j = 0; j < kDims; j++) {
		file << sample(j) / std::sqrt(2 * dt) << " ";
	    }
	    file << std::endl;
	}
	file.close();
    }
}

TEST_CASE("test sample_first_hit", "[sampling]") {
    int seed = 0;
    std::mt19937 gen(seed);
    Uniform uniform(0.0, 1.0);
    int samples;
    std::vector<Vec> particles;
    Vec v;
    if (kDims == 2) {
	v << 2, 2;
	particles.push_back(v);
	v << -10.5, 11.5;
	particles.push_back(v);
    }
    else {
	v << 2, 2, -6.7;
	particles.push_back(v);
	v << -10.5, 11.5, 22.0;
	particles.push_back(v);
    }
    std::vector<double> height_fractions = {0.2, 0.5, 0.9};
    const int X = 0;
    const int Y = 1;
    for (auto particle : particles) {
	for (auto height_fraction : height_fractions) {
	    std::string tag = (", particle = " + 
			       std::to_string(particle(X)) + ", " + 
			       std::to_string(particle(Y)) + 
			       ", height_fraction = " + 
			       std::to_string(height_fraction));
	    SECTION( "test height is what was specified" + tag) {
		double epsilon = 1e-8;
		Approx approx_zero = Approx(0).margin(epsilon);
		samples = 1000;
		bool error = false;
		double height = height_fraction * particle(Y);
		for (int i = 0; i < samples; i++) {
		    Vec sample = Sampling::sample_first_hit(particle, 
							    height, 
							    gen, uniform);
		    if (!(sample(Y) - height == approx_zero)) {
			error = true;
		    }
		}
		REQUIRE(!error);
	    }
	}
    }
    // spit out samples to be analyzed with python
    std::ofstream file;
    std::string dir = "test_distributions/"; 
    if (!boost::filesystem::is_directory(dir)) {
	std::cout << "Can only run sampling tests from path test/" << 
	    std::endl;
	exit(-1);
    }
    file << std::scientific << std::setprecision(12);
    int count = 0;
    for (auto particle : particles) {
	for (auto height_fraction : height_fractions) {
	    double height = height_fraction * particle(Y);
	    samples = 100000;
	    std::string path = (dir + "/sample_first_hit_" + 
				std::to_string(kDims) + "d_" + 
				std::to_string(count) + ".txt");
	    file.open(path);
	    file << "# height = " << height << std::endl;
	    file << "# partitle = ";
	    for (int i = 0; i < kDims; i++) file << particle(i) << ", ";
	    file << std::endl;
	    for (int i = 0; i < samples; i++) {
		Vec sample = Sampling::sample_first_hit(particle, height, 
							gen, uniform);
		for (int j = 0; j < kDims; j++) {
		    file << sample(j) << " ";
		}
		file << std::endl;
	    }
	    file.close();
	    count++;
	}
    }
}
