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
    double L = 100;
    double cell_length = 10;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(L, cell_length, max_leaf_size, seed);
    
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
    v << -13, 17;
    offsets.push_back(v);
    // will check pbcs
    v << 47, 37;
    offsets.push_back(v);
    for (auto offset : offsets) {
	std::string tag = ", offset = " + std::to_string(offset(X)) + ", " + 
	    std::to_string(offset(Y));
	SECTION("test sticking" + tag) {
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
	    answer = enforce_pbcs(answer, L);
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));
	    
	    // test sticking to what just stuck
	    particle = answer;
	    particle(X) -= 2;
	    particle(Y) += 2;
	    state.set_particle(particle);
	    jump << 4, -4;
	    stuck = state.resolve_jump(jump, p);
	    REQUIRE(stuck);
	    answer(X) = offset(X) + 2;
	    answer(Y) = offset(Y) + 2 + 2 * std::sqrt(2);
	    answer = enforce_pbcs(answer, L);
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));
	}
    }

    offsets.clear();
    v << 2, 2;
    offsets.push_back(v);
    // tests pbcs
    v << 47, 4;
    offsets.push_back(v);
    for (auto offset : offsets) {
	std::string tag = ", offset = " + std::to_string(offset(X)) + ", " + 
	    std::to_string(offset(Y));
	SECTION("test sticking bottom" + tag) {
	    particle = offset;
	    state.set_particle(particle);
	    jump << 4, -4;
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(stuck);
	    answer = offset;
	    answer(X) += offset(Y);
	    answer(Y) = 0.0;	    
	    answer = enforce_pbcs(answer, L);
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));
	}
    }
}

void test_resolve_jump_sticking_3d() {
    State state;
    double L = 100;
    double cell_length = 10;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(L, cell_length, max_leaf_size, seed);
    
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
    v << -13, 17, -13;
    offsets.push_back(v);
    // will check pbcs
    v << 47, 37, 47;
    offsets.push_back(v);
    for (auto offset : offsets) {
	std::string tag = ", offset = " + std::to_string(offset(X)) + ", " + 
	    std::to_string(offset(Y)) + ", " + std::to_string(offset(Z));
	SECTION("test sticking" + tag) {
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
	    jump << -4, -4, -4;
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(stuck);
	    for (int i = 0; i < kDims; i++) {
		answer(i) = offset(i) + 2 + 2 * std::sqrt(3) / 3;
	    }
	    answer = enforce_pbcs(answer, L);
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));
	    
	    // test sticking to what just stuck
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
	    answer = enforce_pbcs(answer, L);
	    REQUIRE(state.get_plated_cloud().back().isApprox(answer, epsilon));
	}
    }

    offsets.clear();
    v << 2, 2, 2;
    offsets.push_back(v);
    // tests pbcs
    v << 47, 4, 47;
    offsets.push_back(v);
    for (auto offset : offsets) {
	std::string tag = ", offset = " + std::to_string(offset(X)) + ", " + 
	    std::to_string(offset(Y));
	SECTION("test sticking bottom" + tag) {
	    particle = offset;
	    state.set_particle(particle);
	    jump << 5, -5, 5;
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(stuck);
	    answer = offset;
	    answer(X) += offset(Y);
	    answer(Y) = 0.0;	    
	    answer(Z) += offset(Y);
	    answer = enforce_pbcs(answer, L);
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
    double L = 100;
    double cell_length = 10;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(L, cell_length, max_leaf_size, seed);
    
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
    v << -13, 17;
    offsets.push_back(v);
    // will check pbcs
    v << 47, 37;
    offsets.push_back(v);
    v << -51.5, 37;
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
	    answer = particle;
	    answer = enforce_pbcs(answer, L);
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
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
	    answer = enforce_pbcs(answer, L);
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
	    answer = enforce_pbcs(answer, L);
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
	    answer = enforce_pbcs(particle, L);
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
	}
	SECTION( "test other rejection" + tag) {
	    // This occurs when bounce happens too close border of 
	    // another plated
	    particle = offset;
	    state.set_particle(particle);
	    state.stick_particle();

	    particle = offset;
	    particle(X) = offset(X) - 2 - kSpatialEpsilon / 2.0;
	    particle(Y) = offset(Y) + 2;
	    state.set_particle(particle);
	    state.stick_particle();

	    particle(X) = offset(X) + 1;
	    particle(Y) = offset(Y) + 3;
	    state.set_particle(particle);
	    jump << -3, -3;
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(!stuck);
	    answer = enforce_pbcs(particle, L);
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
	}
    }
}

void test_resolve_jump_bouncing_3d() {
    State state;
    double L = 100;
    double cell_length = 10;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(L, cell_length, max_leaf_size, seed);
    
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
    v << -13, 17, -13;
    offsets.push_back(v);
    // will check pbcs
    v << 47, 37, 47;
    offsets.push_back(v);
    v << -51.5, 37, -51.5;
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
	    answer = particle;
            answer = enforce_pbcs(answer, L);
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
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
            answer = enforce_pbcs(answer, L);
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
	    answer = enforce_pbcs(answer, L);
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
	    answer = enforce_pbcs(answer, L);
	    REQUIRE(state.get_particle().isApprox(answer, epsilon));
	}
	SECTION( "test other rejection" + tag) {
	    // This occurs when bounce happens too close border of 
	    // another plated
	    particle = offset;
	    state.set_particle(particle);
	    state.stick_particle();

	    particle = offset;
	    particle(X) = offset(X) - 2 - kSpatialEpsilon / 2.0;
	    particle(Y) = offset(Y) + 2;
	    state.set_particle(particle);
	    state.stick_particle();

	    particle(X) = offset(X) + 1;
	    particle(Y) = offset(Y) + 3;
	    particle(Z) = offset(Z) + 1;
	    state.set_particle(particle);
	    jump << -3, -3, -3;
	    bool stuck = state.resolve_jump(jump, p);
	    REQUIRE(!stuck);
	    answer = particle;
	    answer = enforce_pbcs(answer, L);
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
    double L = 10;
    double cell_length = 1;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(L, cell_length, max_leaf_size, seed);
    
    Vec particle;
    double epsilon = 1e-8;
    double squared_distance;
    double answer;
    Approx approx_zero = Approx(0).margin(epsilon);
    const int X = 0;
    const int Y = 1;
    SECTION( "test base case") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 0;
	}
	state.set_particle(particle);
	state.stick_particle();
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 1;
	}
	state.set_particle(particle);
	state.stick_particle();
	
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 4;
	}
	state.set_particle(particle);
	squared_distance = state.find_nearest_neighbor();
	// gives right distance for 2 and 3 dimensions
	answer = std::sqrt(9 * kDims);
	REQUIRE(squared_distance - answer == approx_zero);
    }
    SECTION( "test base case with negatives") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = -2;
	}
	particle(Y) = 0;
	state.set_particle(particle);
	state.stick_particle();
	for (int i = 0; i < kDims; i++) {
	    particle(i) = -1;
	}
	particle(Y) = 1;
	state.set_particle(particle);
	state.stick_particle();
	
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 2;
	}
	state.set_particle(particle);
	squared_distance = state.find_nearest_neighbor();
	if (kDims == 2) answer = std::sqrt(10);
	else answer = std::sqrt(19);
	REQUIRE(squared_distance - answer == approx_zero);
    }
    SECTION( "test with pbcs") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 3;
	}
	state.set_particle(particle);
	state.stick_particle();       
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 4.5;
	}
	state.set_particle(particle);
	state.stick_particle();
	
	for (int i = 0; i < kDims; i++) {
	    particle(i) = -4.5;
	}
	particle(Y) = 3.5;
	state.set_particle(particle);
	squared_distance = state.find_nearest_neighbor();
	if (kDims == 2) answer = std::sqrt(2);
	else answer = std::sqrt(3);
	REQUIRE(squared_distance - answer == approx_zero);
    }
    SECTION( "test with alt pbcs") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 3;
	}
	state.set_particle(particle);
	state.stick_particle();       
	for (int i = 0; i < kDims; i++) {
	    particle(i) = -4.5;
	}
	particle(Y) = 4.5;
	state.set_particle(particle);
	state.stick_particle();
	
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 4.5;
	}
	particle(Y) = 3.5;
	if (kDims == 3) particle(X) = -4.5;
	// test flip across Z only in 3d
	state.set_particle(particle);
	squared_distance = state.find_nearest_neighbor();
	answer = std::sqrt(2);
	REQUIRE(squared_distance - answer == approx_zero);
    }
}

TEST_CASE("test take_large_step") {
    State state;
    double L = 100;
    double cell_length = 10;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(L, cell_length, max_leaf_size, seed);
    
    Vec particle;
    double epsilon = 1e-8;
    double answer;
    Approx approx_zero = Approx(0).margin(epsilon);
    SECTION("test base case") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = ((i + 1) * (i + 1)) / 2;
	}
	state.set_particle(particle);
	state.stick_particle();
	
	for (int i = 0; i < kDims; i++) {
	    particle(i) += (i + 1);
	}
	state.set_particle(particle);
	
	answer = (state.find_nearest_neighbor() - kDiameter - 
		  2 * kSpatialEpsilon);
	state.take_large_step();
	REQUIRE((particle - state.get_particle()).norm() - answer == 
		approx_zero);
    }
    SECTION("test detection of bottom of box") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 3;
	}
	state.set_particle(particle);
	state.stick_particle();
	
	for (int i = 0; i < kDims; i++) {
	    particle(i) = ((i + 1) * (i + 1)) / 2 - 1;
	}
	state.set_particle(particle);
	
	answer = 1 - 2 * kSpatialEpsilon;
	state.take_large_step();
	REQUIRE((particle - state.get_particle()).norm() - answer == 
		approx_zero);
    }
}

TEST_CASE("test check_for_regeneration") {
    State state;
    double L = 100;
    double cell_length = 10;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(L, cell_length, max_leaf_size, seed);
    
    Vec particle;
    double epsilon = 1e-8;
    double answer;
    const int Y = 1;
    Approx approx_zero = Approx(0).margin(epsilon);
    for (int i = 0; i < kDims; i++) {
	particle(i) = ((i + 1) * (i + 1)) / 2;
    }
    state.set_particle(particle);
    state.stick_particle();
    SECTION("test regeneration") {
	answer = particle(Y) + kDiameter + 2 * kSpatialEpsilon;
    
	for (int i = 0; i < kDims; i++) {
	    particle(i) += (i + 10);
	}
	state.set_particle(particle);
	
	double dt = 1;
	state.check_for_regeneration(dt);
	REQUIRE((particle - state.get_particle()).norm() > epsilon);
	REQUIRE(state.get_particle()(Y) - answer == approx_zero);
    }
    SECTION("test no regeneration") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) += i + 1.1;
	}
	state.set_particle(particle);
	
	double dt = 1;
	state.check_for_regeneration(dt);
	REQUIRE(state.get_particle().isApprox(particle, epsilon));
    }
}

TEST_CASE("test check_overlaps") {
    State state;
    double L = 10;
    double cell_length = 1;
    int max_leaf_size = 10;
    int seed = 0;
    state.set_up_new_state(L, cell_length, max_leaf_size, seed);

    Vec particle;
    double epsilon = 5e-9;
    const int Y = 1;
    SECTION("test with no overlaps") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 0;
	}
	particle(Y) += 3;
	state.set_particle(particle);
	state.stick_particle();
	for (int i = 0; i < kDims; i++) {
	    particle(i) = -2 * std::sqrt(kDims) / kDims;
	}
	particle(Y) += 3;
	state.set_particle(particle);
	state.stick_particle();

	// small overlap should be ignored
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 2 * std::sqrt(kDims) / kDims - epsilon;
	}
	particle(Y) += 3;
	state.set_particle(particle);
	state.stick_particle();
	
	REQUIRE(state.check_overlaps() == 0);
    }
    SECTION("test with large overlap") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 0;
	}
	particle(Y) += 3;
	state.set_particle(particle);
	state.stick_particle();
	for (int i = 0; i < kDims; i++) {
	    // float division
	    particle(i) = kDiameter / 4.0;
	}
	particle(Y) += 3;
	state.set_particle(particle);
	state.stick_particle();
	for (int i = 0; i < kDims; i++) {
	    // float division
	    particle(i) = -kDiameter / 4.0;
	}
	particle(Y) += 3;
	state.set_particle(particle);
	state.stick_particle();
	REQUIRE(state.check_overlaps() == 3);
    }
    SECTION("test with small overlap") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 0;
	}
	particle(Y) += 3;
	state.set_particle(particle);
	state.stick_particle();
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 2 * std::sqrt(kDims) / kDims - kSpatialEpsilon;
	}
	particle(Y) += 3;
	state.set_particle(particle);
	state.stick_particle();
	REQUIRE(state.check_overlaps() == 1);
    }
    SECTION("test pbcs") {
	for (int i = 0; i < kDims; i++) {
	    particle(i) = 4.75;
	}
	state.set_particle(particle);
	state.stick_particle();
	for (int i = 0; i < kDims; i++) {
	    particle(i) = -4.75;
	}
	particle(Y) = 4.75;
	state.set_particle(particle);
	state.stick_particle();
	REQUIRE(state.check_overlaps() == 1);
    }
}

