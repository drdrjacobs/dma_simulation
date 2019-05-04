/// @file
///
/// @brief Tests cells methods.
///

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

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
	}
	REQUIRE(!error);

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
    std::vector<double> radii = {0.5, 10.0, 17.3};
    Vec particle;
    for (auto radius : radii) {
	std::string tag = ", radius = " + std::to_string(radius);
	SECTION( "test radius is what was specified" + tag) {
	    for (int i = 0; i < kDims; i++) {
		particle(i) = std::pow(-1, i) * radius * i;
	    }
	    double epsilon = 1e-8;
	    Approx approx_zero = Approx(0).margin(epsilon);
	    samples = 1000;
	    bool error = false;
	    for (int i = 0; i < samples; i++) {
		Vec sample = Sampling::sample_first_hit(particle, radius, 
							gen, uniform);
		if (!(sample.norm() - radius == approx_zero)) {
		    error = true;
		}
	    }
	    REQUIRE(!error);
	}
    }
    // spit out samples to be analyzed with python
    samples = 100000;
    double radius = 1.0;
    std::ofstream file;
    for (int option = 0; option < 4; option++) {
	std::string path = ("test_distributions/sample_first_hit_" + 
			    std::to_string(option) + "_" + 
			    std::to_string(kDims) + "d.txt");
	file.open(path);
	double alpha = 2;
	if (option > 1) {
	    alpha = 4;
	}
	if (option % 2 == 0) {
	    for (int i = 0; i < kDims; i++) {
		particle(i) = 1;
	    }
	}
	else {
	    for (int i = 0; i < kDims; i++) {
		particle(i) = std::pow(-1, i + 1);
	    }
	}
	particle = alpha * particle / particle.norm();
	file << "particle = ";
	for (int i = 0; i < kDims; i++) {
	    file << particle(i) << " ";
	}
	file << std::endl;
	for (int i = 0; i < samples; i++) {
	    Vec sample = Sampling::sample_first_hit(particle, radius, 
						    gen, uniform);
	    for (int j = 0; j < kDims; j++) {
		file << sample(j) << " ";
	    }
	    file << std::endl;
	}
	file.close();
    }
}

TEST_CASE("test first_hit_3d_rotation", "[sampling]") {
    if (kDims == 3) {
	Vec hit_vector;
	Vec particle;
	Vec answer;
	double epsilon = 1e-8;
	SECTION( "test no rotation if pointing along same axis") {
	    hit_vector << 3, 0, 3;
	    particle << 0, 0, 5;
	    Vec result = Sampling::first_hit_3d_rotation(hit_vector, particle);
	    REQUIRE(result.cwiseEqual(hit_vector).count() == kDims);
	}
	SECTION( "test no rotation if almost pointing along same axis") {
	    hit_vector << 3, 0, 3;
	    particle << epsilon, epsilon, 5;
	    Vec result = Sampling::first_hit_3d_rotation(hit_vector, particle);
	    REQUIRE(result.cwiseEqual(hit_vector).count() == kDims);
	}
	SECTION( "test rotation") {
	    hit_vector << 1, 2, 3;
	    particle << 5, 0, 0;
	    answer << 3, 2, -1;
	    Vec result = Sampling::first_hit_3d_rotation(hit_vector, particle);
	    REQUIRE(result.isApprox(answer, epsilon));
	}
	SECTION( "test reflection across xy plane if pointing along -z axis") {
	    hit_vector << 3, 0, 3;
	    answer << 3, 0, -3;
	    particle << 0, 0, -5;
	    Vec result = Sampling::first_hit_3d_rotation(hit_vector, particle);
	    REQUIRE(result.cwiseEqual(answer).count() == kDims);
	}
	SECTION( "test reflection across xy plane if almost pointing along "
		 "-z axis") {
	    hit_vector << 3, 0, 3;
	    answer << 3, 0, -3;
	    particle << epsilon, epsilon, -5;
	    Vec result = Sampling::first_hit_3d_rotation(hit_vector, particle);
	    REQUIRE(result.cwiseEqual(answer).count() == kDims);
	}
    }
}
