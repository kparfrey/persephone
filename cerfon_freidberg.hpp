#ifndef CERFON_FREIDBERG_HPP
#define CERFON_FREIDBERG_HPP

#include "common.hpp"
#include "domain_map.hpp"

class CerfonFreidbergConfig
{
    public:

    enum Machine {tftr, jet, iter, mast, nstx, spheromak, frc};

    const Machine machine = mast;

    const real_t A  = 0.0;
    const real_t B0 = 1e-3; //1e-2; // Entirely determines B_phi if A = 0
                            
    real_t epsilon, kappa, delta, alpha;

    real_t c[8] = {}; // Coeffs in the potential expansion, initialised to zero


    CerfonFreidbergConfig();

    void construct_equilibrium(const real_t r[3], real_t U[9], const real_t gamma);
};

#endif
