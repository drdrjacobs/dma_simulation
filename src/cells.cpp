/// @file
///
/// @brief Header for the Cells class which stores plated in cell structure.
///

#include <iostream>
#include <cmath>

#include <boost/range/adaptors.hpp>

#include "cells.hpp"

// define globals
const int Cells::kDims = DIMENSIONS;
// number of cells that must be looped over when looking for collisions
const int Cells::kCellsToLoopOver = std::pow(3, DIMENSIONS);

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
void Cells::set_up_cells(double cell_length) {
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
    CellIndices indices = get_cell_indices(plated);
    cell_map_[indices].push_back(plated);
}

/// @brief Gets cell for collision detection loop.
///
/// The cell the given particle is in is considered the central cell. Looping 
/// through offset = 0, offset < kCellsToLoopOver, offset++ will end up
/// looping through all of the cells surrounding the central cell and the 
/// central cell itself.
///
/// @param particle: position vector of particle, defines central cell
/// @param offset: determines which cell, in relation to the central cell is 
///     examined
///
/// @returns cell: std::vector of plated in cell with offset compared to 
///     central cell
///
const std::vector<Vec>& Cells::get_cell(Vec particle, int offset) {
    const int X = 0;
    const int Y = 1;
    const int Z = 2;
    CellIndices indices = get_cell_indices(particle);
    // -1 shifts adjustment in X/Y/Z to range from -1 to 1 instead of 0     
    // to 2
    indices[X] = indices[X] + ((offset % 3) - 1);
    indices[Y] = indices[Y] + ((offset % 9) / 3 - 1);
    if (kDims == 3) {
	indices[Z] = indices[Z] + ((offset / 9) - 1);
    }
    std::vector<Vec>& cell = cell_map_[indices];
    return cell;
}

/// @brief Checks to see if there are neigboring plated around position vector.
///
/// @param v: checks for neighbors at this position
///
/// @returns result: true if there are plated in v's cell or the cells around 
///     v's cell, false otherwise
///
bool Cells::has_neighbors(Vec v) {
    bool result = false;
    for (int offset = 0; offset < kCellsToLoopOver; offset++) {
	if (!get_cell(v, offset).empty()) {
	    result = true;
	}
    }
    return result;
}
