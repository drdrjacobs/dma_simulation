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
#include "../../src/cells.hpp"
#include "../../src/state.hpp"

TEST_CASE("test offset_get_cell_indices") {
    Cells cells;  
    double cell_length = 1;
    double L = 20;
    cells.set_up_cells(cell_length, L);

    std::vector<Cells::CellIndices> central_cells;
    Cells::CellIndices c;
    // inculdes tests of pbcs
    for (int i = 0; i < kDims; i++) {
	c[i] = 2;
    }
    central_cells.push_back(c);
    for (int i = 0; i < kDims; i++) {
	c[i] = std::pow(-1, i + 1) * (L / 2);
    }
    central_cells.push_back(c);
    for (int i = 0; i < kDims; i++) {
	c[i] = (L / 2) - 1;
    }
    central_cells.push_back(c);
    for (int i = 0; i < kDims; i++) {
	c[i] = (L / 2) - 1;
    }
    const int X = 0;
    c[X] = 0;
    central_cells.push_back(c);
    for (auto central_cell : central_cells) {
	std::string tag = ", central cell = ";
	for (int i = 0; i < kDims; i++) {
	    tag += std::to_string(central_cell[i]) + ", ";
	}
	Cells::CellIndices offset_cell;
	std::set<Cells::CellIndices> store_cells; 
	for (int offset = 0; offset < Cells::kCellsToLoopOver; offset++) {
	    offset_cell = cells.offset_get_cell_indices(central_cell, offset);
	    store_cells.insert(offset_cell);
	}
	SECTION("test total number of cells that get looped through" + tag) {
	    REQUIRE(store_cells.size() == std::pow(3, kDims));
	}
	SECTION( "test individual cells that get looped through" + tag) {
	    for (auto cell : store_cells) {
		for (int i = 0; i < kDims; i++) {
		    int smaller = central_cell[i] - 1;
		    if (i != 1 && smaller < -L / 2) smaller = (L / 2) - 1;
		    int larger = central_cell[i] + 1;
		    if (i != 1 && larger > (L / 2) - 1) larger = -L / 2;
		    REQUIRE((cell[i] == central_cell[i] ||  
			     cell[i] == smaller ||
			     cell[i] == larger));
		}
	    }
	}
    }
}

TEST_CASE("test has_neighbors") {
    Cells cells;  
    double cell_length = 1;
    double L = 20;
    cells.set_up_cells(cell_length, L);

    std::vector<Vec> particles;
    Vec v;
    // inculdes tests of pbcs
    for (int i = 0; i < kDims; i++) {
	v[i] = 1.5;
    }
    particles.push_back(v);
    for (int i = 0; i < kDims; i++) {
	v[i] = std::pow(-1, i + 1) * (L / 2) + 0.5;
    }
    particles.push_back(v);
    for (int i = 0; i < kDims; i++) {
	v[i] = (L / 2) - 0.5;
    }
    particles.push_back(v);
    for (int i = 0; i < kDims; i++) {
	v[i] = (L / 2) - 0.5;
    }
    const int X = 0;
    v[X] = 0;
    particles.push_back(v);

    for (auto particle : particles) {
	std::string tag = ", particle = ";
	for (int i = 0; i < kDims; i++) {
	    tag += std::to_string(particle[i]) + ", ";
	}
	Vec plated;
	for (int i = 0; i < kDims; i++) {
	    plated(i) = particle(i);
	}
	const int X = 0;
	const int Y = 1;
	const int Z = 2;
	SECTION( "test no neighbors with empty cells" + tag) {
	    REQUIRE(!cells.has_neighbors(particle));
	}
	SECTION( "test no neighbors if nearest neighbor two cells out X" + 
		 tag) {
	    plated(X) = particle(X) + 2 * cell_length;
	    plated = enforce_pbcs(plated, L);
	    cells.add_to_cells(plated);
	    REQUIRE(!cells.has_neighbors(particle));
	}
	SECTION( "test has neighbors if nearest neighbor one cell out X" + 
		 tag) {
	    plated(X) = particle(X) + cell_length ;
	    plated = enforce_pbcs(plated, L);
	    cells.add_to_cells(plated);
	    REQUIRE(cells.has_neighbors(particle));
	}
	SECTION( "test has neighbors if nearest neighbor one cell out -X" + 
		 tag) {
	    plated(X) = particle(X) - cell_length ;
	    plated = enforce_pbcs(plated, L);
	    cells.add_to_cells(plated);
	    REQUIRE(cells.has_neighbors(particle));
	}
	SECTION( "test has neighbors if nearest neighbor one cell out XY" + 
		 tag) {
	    plated(X) = particle(X) + cell_length ;
	    plated(Y) = particle(Y) + cell_length ;
	    plated = enforce_pbcs(plated, L);
	    cells.add_to_cells(plated);
	    REQUIRE(cells.has_neighbors(particle));
	}
	if (kDims == 3) {
	    SECTION("test has neighbors if nearest neighbor one cell out XZ" + 
		    tag) {
		plated(X) = particle(X) + cell_length ;
		plated(Z) = particle(Z) + cell_length ;
		plated = enforce_pbcs(plated, L);
		cells.add_to_cells(plated);
		REQUIRE(cells.has_neighbors(particle));
	    }
	}
	SECTION( "test has neighbors if nearest neighbor in same cell" + tag) {
	    for (int i = 0; i < kDims; i++) {
		plated(i) = particle(i) + cell_length / 4.0;
	    }
	    cells.add_to_cells(plated);
	    REQUIRE(cells.has_neighbors(particle));
	}
    }
}

TEST_CASE("test get_cell_indices") {
    Cells cells;    
    Vec v;
    Cells::CellIndices answer;
    double epsilon = 1e-8;
    double cell_length;

    double L = 100;
    std::vector<double> options = {0.111, 2, 5.424, 8.999};
    for (auto option : options) {
	cell_length = option;
	cells.set_up_cells(cell_length, L);
	std::string tag = ", cell_length = " + std::to_string(cell_length);
	SECTION( "base test" + tag) {
	    for (int i = 0; i < kDims; i++) {
		v(i) = cell_length / 2.0;
		answer[i] = 0;
	    }
	    REQUIRE(cells.get_cell_indices(v) == answer);
	}
	SECTION( "test positive edge goes up" + tag) {
	    for (int i = 0; i < kDims; i++) {
		v(i) = cell_length;
		answer[i] = 1;
	    }
	    REQUIRE(cells.get_cell_indices(v) == answer);
	}
	SECTION( "test negative edge goes up" + tag) {
	    for (int i = 0; i < kDims; i++) {
		v(i) = -cell_length;
		answer[i] = -1;
	    }
	    REQUIRE(cells.get_cell_indices(v) == answer);
	}
	SECTION( "test below positive edge goes down" + tag) {
	    for (int i = 0; i < kDims; i++) {
		v(i) = cell_length - epsilon;
		answer[i] = 0;
	    }
	    REQUIRE(cells.get_cell_indices(v) == answer);
	}
	SECTION( "test below negative edge goes down" + tag) {
	    for (int i = 0; i < kDims; i++) {
		v(i) = -cell_length - epsilon;
		answer[i] = -2;
	    }
	    REQUIRE(cells.get_cell_indices(v) == answer);
	}
    }
}
