#ifndef PHYSICS_FUNCTORS_HPP
#define PHYSICS_FUNCTORS_HPP

#include "common.hpp"


/* Abstract base classes for functors for finding fluxes, primitive variables,
 * wave speeds etc. at a single point, for different equation systems. */

class ConservedToPrimitive
{
    public:
    inline virtual void operator()(const real_t* const U, real_t* const P) const = 0;
};


class WaveSpeedsFromPrimitive
{
    public:
    inline virtual void operator()(const real_t* const P, real_t* const c,
                                   const int dir) const = 0;
};


class FluxesFromPrimitive
{
    public:
    inline virtual void operator()(const real_t* const P, real_t (*F)[3]) const = 0;
};
#endif
