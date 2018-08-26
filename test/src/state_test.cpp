/// @file
///
/// @brief Tests cells methods.
///

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "../lib/catch.hpp"
#include "../../src/constants.hpp"
#include "../../src/state.hpp"

void test_resolve_jump_sticking_2d() {
    State state;
    double cell_length = 10;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(cell_length, max_leaf_size, seed);
    
    Vec particle;
    Vec jump;
    double p = 1.0;
    double epsilon = 1e-8;
    Vec answer;
    const int X = 0;
    const int Y = 1;
    // test cases where everything is within a cell and where particle moves
    // across cells
    std::vector<Vec> offsets;
    Vec v;
    v << 2, 2;
    offsets.push_back(v);
    v << 7, 7;
    offsets.push_back(v);
    v << 97, 14;
    offsets.push_back(v);
    v << 94, -13;
    offsets.push_back(v);
    v << -9, -9;
    offsets.push_back(v);
    v << -9, -33;
    offsets.push_back(v);

    for (auto offset : offsets) {
	std::string tag = ", offset = " + std::to_string(offset(X)) + ", " + 
	    std::to_string(offset(Y));
	SECTION( "test sticking" + tag) {
	    particle = offset;
	    state.set_particle(particle);
	    state.stick_particle();
	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 2;
	    }
	    state.set_particle(particle);
	    state.stick_particle();

	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 4;
	    }
	    state.set_particle(particle);
	    jump << -4, -4;
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(stuck);
	    for (int i = 0; i < kDims; i++) {
		answer(i) = offset(i) + 2 + std::sqrt(2);
	    }
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));

	    particle = answer;
	    particle(X) -= 2;
	    particle(Y) += 2;
	    state.set_particle(particle);
	    jump << 4, -4;
	    stuck = state.resolve_jump(jump, p);
	    REQUIRE(stuck);
	    answer(X) = offset(X) + 2;
	    answer(Y) = offset(Y) + 2 + 2 * std::sqrt(2);
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));
	}
    }
}

void test_resolve_jump_sticking_3d() {
    State state;
    double cell_length = 10;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(cell_length, max_leaf_size, seed);
    
    Vec particle;
    Vec jump;
    double p = 1.0;
    double epsilon = 1e-8;
    Vec answer;
    const int X = 0;
    const int Y = 1;
    const int Z = 2;
    // test cases where everything is within a cell and where particle moves
    // across cells
    std::vector<Vec> offsets;
    Vec v;
    v << 2, 2, 2;
    offsets.push_back(v);
    v << 7, 7, 7;
    offsets.push_back(v);
    v << 97, 14, 56;
    offsets.push_back(v);
    v << 94, -14, -33;
    offsets.push_back(v);
    v << -9, -9, -9;
    offsets.push_back(v);
    v << -9, -36, -77;
    offsets.push_back(v);
    
    for (auto offset : offsets) {
	std::string tag = ", offset = " + std::to_string(offset(X)) + ", " + 
	    std::to_string(offset(Y)) + ", " + std::to_string(offset(Z));
	SECTION( "test sticking within same cell" + tag) {
	    particle << offset;
	    state.set_particle(particle);
	    state.stick_particle();
	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 2;
	    }
	    state.set_particle(particle);
	    state.stick_particle();

	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 4;
	    }
	    state.set_particle(particle);
	    jump << -4, -4, -4;
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(stuck);
	    for (int i = 0; i < kDims; i++) {
		answer(i) = offset(i) + 2 + 2 * std::sqrt(3) / 3;
	    }
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));

	    particle = answer;
	    particle(X) -= 2;
	    particle(Y) += 2;
	    particle(Z) += 2;
	    state.set_particle(particle);
	    jump << 4, -4, -4;
	    stuck = state.resolve_jump(jump, p);
	    REQUIRE(stuck);
	    answer(X) = offset(X) + 2;
	    answer(Y) = offset(Y) + 2 + 4 * std::sqrt(3) / 3;
	    answer(Z) = offset(Z) + 2 + 4 * std::sqrt(3) / 3;
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));
	}
    }
}

TEST_CASE("test resolve_jump sticking") {
    if (kDims == 2) {
	test_resolve_jump_sticking_2d();
    }
    else {
	test_resolve_jump_sticking_3d();
    }
}

void test_resolve_jump_bouncing_2d() {
    State state;
    double cell_length = 10;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(cell_length, max_leaf_size, seed);
    
    Vec particle;
    Vec jump;
    double p = 0.0;
    double epsilon = 1e-8;
    Vec answer;
    const int X = 0;
    const int Y = 1;
    // test cases where everything is within a cell and where particle moves
    // across cells
    std::vector<Vec> offsets;
    Vec v;
    v << 2, 2;
    offsets.push_back(v);
    v << 7, 7;
    offsets.push_back(v);
    v << 97, 14;
    offsets.push_back(v);
    v << 94, -13;
    offsets.push_back(v);
    v << -9, -9;
    offsets.push_back(v);
    v << -9, -33;
    offsets.push_back(v);

    for (auto offset : offsets) {
	std::string tag = ", offset = " + std::to_string(offset(X)) + ", " + 
	    std::to_string(offset(Y));
	SECTION( "test bounce once straight" + tag) {
	    particle = offset;
	    state.set_particle(particle);
	    state.stick_particle();
	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 2;
	    }
	    state.set_particle(particle);
	    state.stick_particle();

	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 4;
	    }
	    state.set_particle(particle);
	    jump << -1, -1;
	    jump = 2 * (2 * std::sqrt(2) - 2) * (jump / jump.norm());
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(!stuck);
	    for (int i = 0; i < kDims; i++) {
		answer(i) = offset(i) + 2;
	    }
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));
	}
	SECTION( "test bounce once angled" + tag) {
	    particle = offset;
	    state.set_particle(particle);
	    state.stick_particle();
	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 2;
	    }
	    state.set_particle(particle);
	    state.stick_particle();

	    particle(X) = offset(X) + 3;
	    particle(Y) = offset(Y) + 5;
	    state.set_particle(particle);
	    jump << -1, -1;
	    jump = 2 * std::sqrt(2) * (jump / jump.norm());
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(!stuck);
	    answer(X) = offset(X) + 1;
	    answer(Y) = offset(Y) + 5;
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
	}
	SECTION( "test bounce twice angled" + tag) {
	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 2;
	    }
	    state.set_particle(particle);
	    state.stick_particle();
	    particle(X) = offset(X) - 1;
	    particle(Y) = offset(Y) + 5;
	    state.set_particle(particle);
	    state.stick_particle();

	    particle(X) = offset(X) + 3;
	    particle(Y) = offset(Y) + 5;
	    state.set_particle(particle);
	    jump << -1, -1;
	    jump = 3 * std::sqrt(2) * (jump / jump.norm());
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(!stuck);
	    answer(X) = offset(X) + 2;
	    answer(Y) = offset(Y) + 6;
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
	}
	SECTION( "test rejection" + tag) {
	    particle = offset;
	    state.set_particle(particle);
	    state.stick_particle();

	    particle(X) = offset(X) + 1;
	    particle(Y) = offset(Y) + 3;
	    state.set_particle(particle);
	    jump << -1, -1;
	    jump = (std::sqrt(2) * (1 + kSpatialEpsilon / 2) * 
		    (jump / jump.norm()));
	    REQUIRE(jump.norm() > 1.0);
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(!stuck);
	    answer = particle;
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
	}
    }
}

void test_resolve_jump_bouncing_3d() {
    State state;
    double cell_length = 10;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(cell_length, max_leaf_size, seed);
    
    Vec particle;
    Vec jump;
    double p = 0.0;
    double epsilon = 1e-8;
    Vec answer;
    const int X = 0;
    const int Y = 1;
    const int Z = 2;
    // test cases where everything is within a cell and where particle moves
    // across cells
    std::vector<Vec> offsets;
    Vec v;
    v << 2, 2, 2;
    offsets.push_back(v);
    v << 7, 7, 7;
    offsets.push_back(v);
    v << 97, 14, 56;
    offsets.push_back(v);
    v << 94, -14, -33;
    offsets.push_back(v);
    v << -9, -9, -9;
    offsets.push_back(v);
    v << -9, -36, -77;
    offsets.push_back(v);

    for (auto offset : offsets) {
	std::string tag = ", offset = " + std::to_string(offset(X)) + ", " + 
	    std::to_string(offset(Y)) + ", " + std::to_string(offset(Z));
	SECTION( "test bounce once straight" + tag) {
	    particle = offset;
	    state.set_particle(particle);
	    state.stick_particle();
	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 2;
	    }
	    state.set_particle(particle);
	    state.stick_particle();

	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 4;
	    }
	    state.set_particle(particle);
	    jump << -1, -1, -1;
	    jump = 2 * (2 * std::sqrt(3) - 2) * (jump / jump.norm());
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(!stuck);
	    for (int i = 0; i < kDims; i++) {
		answer(i) = offset(i) + 2;
	    }
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));
	}
	SECTION( "test bounce once angled" + tag) {
	    particle = offset;
	    state.set_particle(particle);
	    state.stick_particle();
	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 2;
	    }
	    state.set_particle(particle);
	    state.stick_particle();

	    particle(X) = offset(X) + 3;
	    particle(Y) = offset(Y) + 5;
	    particle(Z) = offset(Z) + 3;
	    state.set_particle(particle);
	    jump << -1, -1, -1;
	    jump = 2 * std::sqrt(3) * (jump / jump.norm());
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(!stuck);
	    answer(X) = offset(X) + 1;
	    answer(Y) = offset(Y) + 5;
	    answer(Z) = offset(Z) + 1;
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
	}
	SECTION( "test bounce twice angled" + tag) {
	    for (int i = 0; i < kDims; i++) {
		particle(i) = offset(i) + 2;
	    }
	    state.set_particle(particle);
	    state.stick_particle();
	    particle(X) = offset(X) - 1;
	    particle(Y) = offset(Y) + 5;
	    particle(Z) = offset(Z) + 1;
	    state.set_particle(particle);
	    state.stick_particle();

	    particle(X) = offset(X) + 3;
	    particle(Y) = offset(Y) + 5;
	    particle(Z) = offset(Z) + 3;
	    state.set_particle(particle);
	    jump << -1, -1, -1;
	    jump = 3 * std::sqrt(3) * (jump / jump.norm());
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(!stuck);
	    answer(X) = offset(X) + 2;
	    answer(Y) = offset(Y) + 6;
	    answer(Z) = offset(Z);
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
	}
	SECTION( "test rejection" + tag) {
	    particle = offset;
	    state.set_particle(particle);
	    state.stick_particle();

	    particle(X) = offset(X) + 1;
	    particle(Y) = offset(Y) + 3;
	    particle(Z) = offset(Z) + 1;
	    state.set_particle(particle);
	    jump << -1, -1, -1;
	    jump = std::sqrt(3) * (jump / jump.norm());
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(!stuck);
	    answer = particle;
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
	}
    }
}

TEST_CASE("test resolve_jump bouncing") {
    if (kDims == 2) {
	test_resolve_jump_bouncing_2d();
    }
    else {
	test_resolve_jump_bouncing_3d();
    }
}

TEST_CASE("test find_nearest_neighbor") {
    State state;
    double cell_length = 1;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(cell_length, max_leaf_size, seed);
    
    Vec particle;
    double epsilon = 1e-8;
    double squared_distance;
    double answer;
    Approx approx_zero = Approx(0).margin(epsilon);
    SECTION( "test with positive numbers") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = ((i + 1) * (i + 1)) / 2;
	}
	state.set_particle(particle);
	state.stick_particle();
	
	for (int i = 0; i < kDims; i++) {
	    particle(i) += (i + 1);
	}
	state.set_particle(particle);
	squared_distance = state.find_nearest_neighbor();
	// gives right distance for 2 and 3 dimensions
	answer = 9 * (kDims - 2) + 5;
	REQUIRE(squared_distance - answer == approx_zero);
    }
    SECTION( "test with negative numbers") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = -((i + 2) * (i + 2)) / 2;
	}
	state.set_particle(particle);
	state.stick_particle();
	
	for (int i = 0; i < kDims; i++) {
	    particle(i) -= (i + 1);
	}
	state.set_particle(particle);
	squared_distance = state.find_nearest_neighbor();
	// gives right distance for 2 and 3 dimensions
	answer = 9 * (kDims - 2) + 5;
	REQUIRE(squared_distance - answer == approx_zero);
    }
}

TEST_CASE("test take_large_step") {
    State state;
    double cell_length = 1;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(cell_length, max_leaf_size, seed);
    
    Vec particle;
    double epsilon = 1e-8;
    double answer;
    Approx approx_zero = Approx(0).margin(epsilon);
    SECTION("test with positive numbers") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = ((i + 1) * (i + 1)) / 2;
	}
	state.set_particle(particle);
	state.stick_particle();
	
	for (int i = 0; i < kDims; i++) {
	    particle(i) += (i + 1);
	}
	state.set_particle(particle);
	
	answer = (std::sqrt(state.find_nearest_neighbor()) - kDiameter - 
		  2 * kSpatialEpsilon);
	state.take_large_step();
	REQUIRE((particle - state.get_particle()).norm() - answer == 
		approx_zero);
    }
    SECTION("test with negative numbers") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = -((i + 2) * (i + 3)) / 2;
	}
	state.set_particle(particle);
	state.stick_particle();
	
	for (int i = 0; i < kDims; i++) {
	    particle(i) -= (i + 1);
	}
	state.set_particle(particle);
	
	answer = (std::sqrt(state.find_nearest_neighbor()) - kDiameter - 
		  2 * kSpatialEpsilon);
	state.take_large_step();
	REQUIRE((particle - state.get_particle()).norm() - answer == 
		approx_zero);
    }
}

TEST_CASE("test check_for_regeneration") {
    State state;
    double cell_length = 1;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(cell_length, max_leaf_size, seed);
    
    Vec particle;
    double epsilon = 1e-8;
    double answer;
    Approx approx_zero = Approx(0).margin(epsilon);
    for (int i = 0; i < kDims; i++) {
	particle(i) = ((i + 1) * (i + 1)) / 2;
    }
    state.set_particle(particle);
    state.stick_particle();
    answer = particle.norm() + kDiameter + 2 * kSpatialEpsilon;
    
    for (int i = 0; i < kDims; i++) {
	particle(i) += (i + 10);
    }
    state.set_particle(particle);
    
    state.check_for_regeneration();
    REQUIRE((particle - state.get_particle()).norm() > epsilon);
    REQUIRE(state.get_particle().norm() - answer == approx_zero);
}

TEST_CASE("test check_overlaps") {
    State state;
    double cell_length = 1;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(cell_length, max_leaf_size, seed);
    
    Vec particle;
    double epsilon = 5e-9;
    SECTION("test with no overlaps") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = -2 * std::sqrt(kDims) / kDims;
	}
	state.set_particle(particle);
	state.stick_particle();

	// small overlap should be ignored
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 2 * std::sqrt(kDims) / kDims - epsilon;
	}
	state.set_particle(particle);
	state.stick_particle();
	
	REQUIRE(state.check_overlaps() == 0);
    }
    SECTION("test with large overlap") {
	for (int i = 0; i < kDims; i++) {
	    // float division
	    particle(i) = kDiameter / 4.0;
	}
	state.set_particle(particle);
	state.stick_particle();
	for (int i = 0; i < kDims; i++) {
	    // float division
	    particle(i) = -kDiameter / 4.0;
	}
	state.set_particle(particle);
	state.stick_particle();
	REQUIRE(state.check_overlaps() == 3);
    }
    SECTION("test with small overlap") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 2 * std::sqrt(kDims) / kDims - kSpatialEpsilon;
	}
	state.set_particle(particle);
	state.stick_particle();
	REQUIRE(state.check_overlaps() == 1);
    }
}


