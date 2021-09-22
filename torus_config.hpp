#ifndef TORUS_CONFIG_HPP
#define TORUS_CONFIG_HPP

#include "common.hpp"

class TorusConfig
{
    public:

    virtual void unit_disc_to_physical_space(real_t r[3]) const = 0;
    virtual void construct_equilibrium(const real_t r[3], real_t U[9]) const = 0;
};

#endif
