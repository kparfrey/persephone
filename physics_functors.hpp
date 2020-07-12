#ifndef PHYSICS_FUNCTORS_HPP
#define PHYSICS_FUNCTORS_HPP

#include "common.hpp"


/* Abstract base classes for functors for finding fluxes, primitive variables,
 * wave speeds etc. at a single point, for different equation systems. */

class ConservedToPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ U, 
                                         real_t* const __restrict__ P) const = 0;
};


/* Returns the signed fastest wave speeds in this direction, moving 
 * in the positive direction (c[0]) and negative direction (c[1]) */
class WaveSpeedsFromPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ P, 
                                         real_t* const __restrict__ c,
                                   const int dir) const = 0;
};


/* To speed up the simple straight-element case with no directional mixing, 
 * could overload operator() with a version taking an extra int dir argument 
 * for use in XYZ_straight numerical fluxes */
class FluxesFromPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ P, 
                                         real_t (*__restrict__ F)[3]) const = 0;
};
#endif
