/// @file
///
/// @brief Header for the Cells class which stores plated in cell structure.
///

#include <iostream>
#include <cmath>

#include "cells.hpp"

// define globals
const int Cells::kDims = DIMENSIONS;

/// @brief Constructor, must also run set_up_cells to properly initialize.
///
Cells::Cells() {
    cell_length_ = 0;
}

/// @brief Returns reference to cell_map_, for debugging.
///
Cells::CellMap& Cells::get_cell_map() {
    return cell_map_;
}

/// @brief Sets cell_length
///
/// @param cell_length: the length of each cell
///
void Cells::set_up_cells(float cell_length) {
    cell_length_ = cell_length;
}

/// @brief Gets the cell indicies for a given position vector.
///
/// @param v: the position to get the cell indices for
///
/// @returns result: the cell indices associated with v
///
Cells::CellIndices Cells::get_cell_indices(Vec v) const {
    CellIndices result;
    for (int i = 0; i < kDims; i++) {
	result[i] = std::floor(v[i] / cell_length_);
    }
    return result;
}

/// @brief Adds plated to appropriate cells.
///
/// @param plated: position vector of plated to add
///
void Cells::add_to_cells(Vec plated) {
    const int X = 0;
    const int Y = 1;
    const int Z = 2;
    CellIndices indices = get_cell_indices(plated);
    // Need to add plated to its actual cell and all surrounding 3^kDims cells
    for (int i = -1; i < 2; i++) {
	for (int j = -1; j < 2; j++) {
	    CellIndices add_indices = indices;
	    add_indices[X] += i;
	    add_indices[Y] += j;
	    if (kDims == 3) {
		for (int k = -1; k < 2; k++) {
		    add_indices[Z] = indices[Z];
		    add_indices[Z] += k;
		    cell_map_[add_indices].push_back(plated);
		}
	    }
	    else {
		cell_map_[add_indices].push_back(plated);
	    }
	}
    }
}

/// @brief Gets all the plated in the cell of the given vector.
///
/// @param v: the vector to get the surrounding plated for
///
std::vector<Vec>& Cells::get_neighbors(Vec v) {
    CellIndices indices = get_cell_indices(v);
    return cell_map_[indices];
}
