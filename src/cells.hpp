/// @file
///
/// @brief Header for the Cells class which stores plated in cell structure.
///

#pragma once

#include <vector>
#include <unordered_map>

// use eigen for vector operations
#include <Eigen/Dense>

#include <boost/functional/hash.hpp>

#include "vec.hpp"

/// @brief Cells class stores plated in cell structure.
///
/// Works for two or three dimensions. Set during compile time.
///
class Cells {
public:
    
    typedef std::array<int, DIMENSIONS> CellIndices;
    typedef std::unordered_map<CellIndices, std::vector<Vec>, 
			       boost::hash<CellIndices>> CellMap;

    // Documented in the cpp
    static const int kDims;
    
    Cells();
    ~Cells() {};
    CellMap& get_cell_map();
    void set_up_cells(double cell_length);
    void add_to_cells(Vec plated);
    std::vector<Vec>& get_neighbors(Vec v);

private:
    /// Plated are stored in an unordered_map based cells_map_ structure.
    /// Cells are not actually cells in the regular sense, but super cells.
    /// Each cell stores all the plated within a central spatial square/cube 
    /// (standard cell) and all the plated in the spatial squares/cubes around 
    /// that central cube. This creates redudancy, but that is useful for 
    /// getting a parallel speedup. So you can replace every intance of cell 
    /// with "super_cell" in this code.
    CellMap cell_map_;
    /// Cells are square/cubic, side length of each square/cube is the max
    /// length particle can jump in one dt_ plus one diameter + epsilon
    /// ensures all collisions can be resolved
    double cell_length_;

    CellIndices get_cell_indices(Vec v) const;
};
