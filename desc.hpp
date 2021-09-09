#ifndef DESC_HPP
#define DESC_HPP

#include "common.hpp"
#include <string>
#include <vector>

class DescConfig
{
    public:
    
    int                           L,       M,       N;       // Resolution of lambda, and R & Z polynomials
    std::vector<double>           L_lmn,   R_lmn,   Z_lmn;   // Store mode coefficients as flattened 1D vectors
    std::vector<std::vector<int>> L_modes, R_modes, Z_modes; // 2D vector relating 1D to 3D mode indices
    std::vector<double>           iota, pressure;            // Profiles of rotational transform and pressure

    int N_L, N_R, N_Z; // Total number of modes for lambda and surface polynomials
    int N_iota, N_pressure; // Number of modes in the rot. transform and pressure polynomials

    DescConfig(std::string input_file, const int iteration);

    //void construct_equilibrium(const real_t r[3], real_t U[9]);
};

#endif
