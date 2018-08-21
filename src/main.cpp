/// @file
///
/// @brief Implementation of the main method for the simualtion.
///

#include <iostream>
#include <iomanip>

#include "simulation.hpp"

/// @brief Main for running simulation.
///
/// @param argc: number of command line arguments
/// @param argv: command line arguments
///
int main(int argc, char **argv) {
    std::cout << std::scientific << std::setprecision(12) << std::endl;
    // Simulation object runs everything
    Simulation simulation;
    simulation.run_simulation();
    return 0;
}
