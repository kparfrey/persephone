#ifndef DESC_HPP
#define DESC_HPP

#include "common.hpp"
#include <string>
#include <vector>

class DescConfig
{
    private:

    void surface_polynomial_expansion(real_t& f, const real_t r[3], const std::vector<double>& f_lmn, 
                                      const std::vector<std::vector<int>>& f_modes);

    public:
    
    int                           L,       M,       N;       // Resolution of lambda, and R & Z polynomials
    std::vector<double>           L_lmn,   R_lmn,   Z_lmn;   // Store mode coefficients as flattened 1D vectors
    std::vector<std::vector<int>> L_modes, R_modes, Z_modes; // 2D vector relating 1D to 3D mode indices
    std::vector<double>           iota, pressure;            // Profiles of rotational transform and pressure

    double Nfp; // Number of field periods --- stored as a double

    // Should remove these, not needed
    int N_L, N_R, N_Z; // Total number of modes for lambda and surface polynomials
    int N_iota, N_pressure; // Number of modes in the rot. transform and pressure polynomials

    DescConfig(std::string input_file, const int iteration);
    
    void unit_disc_to_physical_space(real_t r[3]);

    //void construct_equilibrium(const real_t r[3], real_t U[9]);
};

#endif
