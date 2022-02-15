#ifndef TORUS_CONFIG_HPP
#define TORUS_CONFIG_HPP

#include "common.hpp"

class TorusConfig
{
    public:
    real_t gamma;    // Adiabatic index - will be copied in from MHD object
    real_t sqrt_mu0; // Copied in from MHD object

    virtual void unit_disc_to_physical_space(real_t r[3]) const = 0;

    /* construct_equilibrium takes both UDS and physical coords, because different config types will use
     * different coords: CF uses physical, DESC needs UDS. */
    virtual void construct_equilibrium(const real_t r_uds[3], const real_t r_phys[3], real_t U[9]) const = 0;
};

#endif
