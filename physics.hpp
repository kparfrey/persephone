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
    int Nfield; // Total number of fields
    int Ncons;  // Number of conserved variables
    int Nflux;  // Number of conserved variables to calculate fluxes for
    string* variables;

    /* The metric now lives here */
    SpatialMetric* metric;

    inline static bool diffusive;
    inline static real_t viscosity;   // kinematic viscosity ("nu")
    inline static real_t resistivity; // really magnetic diffusivity...
    inline static real_t diffusive_timestep_const; // effective "CFL for diffusion"
    inline static bool apply_floors; // Whether to use a flooring procedure

    Physics()
    {
        /* Default values unless overriden in derived classes */
        diffusive = false;
        viscosity = 0.0;
        resistivity = 0.0;
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
};

#endif
