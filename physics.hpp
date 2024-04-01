#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include "common.hpp"
#include "spatial_metric.hpp"
#include <string>
using std::string;


class Physics
{
    public:
    /* Defined directly in derived-class constructors */
    EqnSystem system;
    int Nfield;
    string* variables;

    /* The metric now lives here */
    SpatialMetric* metric;

    inline static bool diffusive;
    inline static real_t viscosity;    // kinematic viscosity ("nu")
    inline static real_t resistivity;  // really magnetic diffusivity...
    inline static real_t conductivity;  
    inline static real_t diffusive_timestep_const; // effective "CFL for diffusion"
    inline static bool apply_floors; // Whether to use a flooring procedure

    /* Only used by MHD */
    inline static real_t ch; 
    inline static real_t psi_damping_const; // Dedner's cr; Mignone & Tzeferacos's alpha, 0 < const < 1
    inline static real_t psi_damping_rate;  // Derigs' alpha = ch^2/cp^2 = ch/cr; d psi/dt ... = - rate * psi
                                            // p_d_rate = CFL * p_d_const / dt
    //inline static real_t psi_damping_exp;   // exp(-damping_rate * dt)


    Physics()
    {
        /* Default values unless overriden in derived classes */
        diffusive = false;
        viscosity = 0.0;
        resistivity = 0.0;
        conductivity = 0.0;
        diffusive_timestep_const = 1.0/3.0; // Appears safe even for WaveRect
        apply_floors = false;
    }


    /* Methods 
     * These all take an int mem, which is passed to SpatialMetric
     * and indicates the active memory location. */
    virtual void ConservedToPrimitive(const real_t* const __restrict__ U, 
                                            real_t* const __restrict__ P,
                                      const int mem) const = 0;

    virtual void WaveSpeeds(const real_t* const __restrict__ P, 
                                  real_t* const __restrict__ c,
                            const int dir,
                            const int mem) const = 0;

    virtual void Fluxes(const real_t* const __restrict__ P, 
                              real_t (*__restrict__ F)[3],
                        const int mem) const = 0;

    virtual void DiffusiveFluxes(const real_t* const __restrict__ U, 
                                 const real_t (* const __restrict__ dU)[3],
                                       real_t (*__restrict__ F)[3],
                                 const int mem) const = 0;
   
    virtual void OrthonormaliseVectors(real_t* const P, const int mem) const = 0;

    virtual void Floors(real_t* const __restrict__ U, const int mem) const = 0;

    /* Only used by MHD, and potential future MHD extensions */ 
    //virtual void Fluxes_divB_subsystem(const real_t* const __restrict__ P, 
    //                                         real_t (*__restrict__ F)[3],
    //                                   const int mem) const {}
};


/**************************************************************************/


#if 0
class SystemData
{
    public:
    int Nfield;
    string* variables;

    bool diffusive;
    real_t viscosity;   // kinematic viscosity ("nu")
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
        viscosity = 0.0;
        resistivity = 0.0;
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
    real_t ch_sq; // MHD only: square of the div-cleaning wavespeed

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
#endif
