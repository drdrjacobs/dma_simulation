/// @file
///
/// @brief Header for the Cells class which stores plated in cell structure.
///

#pragma once

#include <vector>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "constants.hpp"

/// @brief Cells class stores plated in cell structure.
///
/// Works for two or three dimensions. Set during compile time.
///
class Cells {
public:
    
    typedef std::array<int, kDims> CellIndices;
    typedef std::unordered_map<CellIndices, std::vector<Vec>, 
			       boost::hash<CellIndices>> CellMap;

    // Documented in the cpp
    static const int kCellsToLoopOver;
    
    Cells();
    ~Cells() {};
    const CellMap& get_cell_map() const;
    void set_up_cells(double cell_length);
    void add_to_cells(Vec plated);
    CellIndices offset_get_cell_indices(Vec particle, int offset) const;
    CellIndices offset_get_cell_indices(CellIndices central_cell, 
					int offset) const;
    const std::vector<Vec>& get_cell(CellIndices cell_indices) const;
    bool has_neighbors(Vec v) const;

private:
    /// Plated are stored in an unordered_map based cells_map_ structure.
    CellMap cell_map_;
    /// Cells are square/cubic, side length of each square/cube is the max
    /// length particle can jump in one dt_ plus one diameter + epsilon
    /// ensures all collisions can be resolved
    double cell_length_;
    /// Empty cell is returned if indices not in cell map
    std::vector<Vec> empty_cell_;

    CellIndices get_cell_indices(Vec v) const;
};
