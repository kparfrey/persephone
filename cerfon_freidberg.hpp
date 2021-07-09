#ifndef CERFON_FREIDBERG_HPP
#define CERFON_FREIDBERG_HPP

#include "common.hpp"
#include "domain_map.hpp"

class CerfonFreidbergConfig
{
    public:

    enum Machine {tftr, jet, iter, mast, nstx, spheromak, frc};

    const Machine machine = tftr;

    const real_t A  = 0.0;
    const real_t B0 = 1.0; // Set directly for now...
                           // Should set q ~ B0 / B_poloidal?
    
    real_t epsilon, kappa, delta, alpha;

    real_t c[8] = {}; // Coeffs in the potential expansion, initialised to zero


    CerfonFreidbergConfig();

    void construct_equilibrium(const real_t r[3], real_t U[9], const real_t gamma);
};

#endif
