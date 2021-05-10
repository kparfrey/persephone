#ifndef PHYSICS_FUNCTORS_HPP
#define PHYSICS_FUNCTORS_HPP

#include "common.hpp"
#include <string>
using std::string;


class SystemData
{
    public:
    int Nfield;
    string* variables;

    bool diffusive;
    real_t viscosity;   // dynamic viscosity ("mu")
    real_t resistivity; // really magnetic diffusivity...

    /* Only used by MHD */
    real_t c_h; // Wave/transport speed for hyperbolic part of div cleaning
    real_t psi_damping_const; // Mignone & Tzeferacos's alpha, 0 < const < 1
    real_t psi_damping_rate;  // d psi/dt ... = - rate * psi
                              // p_d_rate = CFL * p_d_const / dt
    real_t psi_damping_exp;   // exp(-damping_rate * dt)


    SystemData()
    {
        /* Default values unless overriden in derived classes */
        diffusive = false;
    }
};


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


class DiffusiveFluxes
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ U, 
                                   const real_t (* const __restrict__ dU)[3],
                                         real_t (*__restrict__ F)[3],
                                   const real_t* const __restrict__ args) const = 0;
};
#endif
