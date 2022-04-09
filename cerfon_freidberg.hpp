#ifndef CERFON_FREIDBERG_HPP
#define CERFON_FREIDBERG_HPP

#include "common.hpp"
#include "torus_config.hpp"

class CerfonFreidbergConfig : public TorusConfig
{
    public:

    enum Machine {tftr, jet, iter, mast, nstx, spheromak, frc};

    const Machine machine = jet;

    const real_t A  = 0.0;
    const real_t B0 = 1e-3; // Entirely determines B_phi if A = 0
                            
    real_t epsilon, kappa, delta, alpha;

    real_t c[8] = {}; // Coeffs in the potential expansion, initialised to zero


    CerfonFreidbergConfig();

    void unit_disc_to_physical_space(real_t r[3]) const override;
    void construct_equilibrium(const real_t r_uds[3], const real_t r_phys[3], real_t U[11]) const override;
};

#endif
