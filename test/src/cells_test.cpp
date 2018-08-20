/// @file
///
/// @brief Tests cells methods.
///

#include <vector>
#include <string>
#include <set>

#include "../lib/catch.hpp"
#include "../../src/constants.hpp"
#include "../../src/cells.hpp"

TEST_CASE("test offset_get_cell_indices") {
    Cells cells;  
    double cell_length = 1;
    cells.set_up_cells(cell_length);
    Cells::CellIndices central_cell;
    for (int i = 0; i < kDims; i++) {
	central_cell[i] = i * i;
    }
    Cells::CellIndices offset_cell;

    std::set<Cells::CellIndices> store_cells; 
    for (int offset = 0; offset < Cells::kCellsToLoopOver; offset++) {
	offset_cell = cells.offset_get_cell_indices(central_cell, offset);
	store_cells.insert(offset_cell);
    }

    SECTION( "test total number of cells that get looped through") {
	REQUIRE(store_cells.size() == std::pow(3, kDims));
    }
    SECTION( "test individual cells that get looped through") {
	for (auto cell : store_cells) {
	    for (int i = 0; i < kDims; i++) {
		REQUIRE((cell[i] == i * i || cell[i] == i * i + 1 || 
			 cell[i] == i * i - 1));
	    }
	}
    }
}

TEST_CASE("test has_neighbors") {
    Cells cells;  
    double cell_length = 1;
    cells.set_up_cells(cell_length);
    Vec particle;
    std::vector<double> offsets = {0, 10000, -513, -5000};
    for (double offset : offsets) {
	std::string tag = ", offset = " + std::to_string(offset);
	for (int i = 0; i < kDims; i++) {
	    particle(i) = cell_length * (0.5 + offset);
	}
	Vec plated;
	for (int i = 0; i < kDims; i++) {
	    plated(i) = particle(i);
	}
	const int X = 0;
	const int Y = 1;

	SECTION( "test no neighbors with empty cells" + tag) {
	    REQUIRE(cells.has_neighbors(particle) == false);
	}
	SECTION( "test no neighbors if nearest neighbor two cells out" + tag) {
	    plated(X) = particle(X) + 2 * cell_length;
	    cells.add_to_cells(plated);
	    REQUIRE(cells.has_neighbors(particle) == false);
	}
	SECTION( "test has neighbors if nearest neighbor one cell out" + tag) {
	    plated(Y) = particle(Y) + cell_length ;
	    cells.add_to_cells(plated);
	    REQUIRE(cells.has_neighbors(particle) == true);
	}
	SECTION( "test has neighbors if nearest neighbor in same cell" + tag) {
	    for (int i = 0; i < kDims; i++) {
		plated(i) = particle(i) + cell_length / 4.0;
	    }
	    cells.add_to_cells(plated);
	    REQUIRE(cells.has_neighbors(particle) == true);
	}
    }
}

TEST_CASE("test get_cell_indices") {
    Cells cells;    
    Vec v;
    Cells::CellIndices answer;
    double epsilon = 1e-8;
    double cell_length;

    std::vector<double> options = {0.111, 2, 5.424, 8.999};
    for (auto option : options) {
	cell_length = option;
	cells.set_up_cells(cell_length);
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
