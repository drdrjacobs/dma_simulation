/// @file
///
/// @brief Tests cells methods.
///

#include <vector>
#include <string>
#include <set>
#include <iostream>

#include "../lib/catch.hpp"
#include "../../src/constants.hpp"
#include "../../src/sampling.hpp"

TEST_CASE("test first_hit_3d_rotation") {
    if (kDims == 3) {
	Vec hit_vector;
	Vec particle;
	Vec answer;
	double epsilon = 1e-8;
	SECTION( "test no rotation if point along same axis") {
	    hit_vector << 3, 0, 3;
	    particle << 0, 0, 5;
	    Vec result = Sampling::first_hit_3d_rotation(hit_vector, particle);
	    REQUIRE(result.cwiseEqual(hit_vector).count() == kDims);
	}
	SECTION( "test no rotation if almost point along same axis") {
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
    }
}
