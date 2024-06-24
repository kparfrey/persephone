#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include <string>
#include "common.hpp"
#include "spatial_metric.hpp"
using std::string;


template <class T>
class Physics
{
    public:
    
    /* Replace with MetricType once the full metric has been implemented */
    DiagonalSpatialMetric metric;

    /* Defined directly in derived-class constructors */
    inline static EqnSystem system;
    inline static int Nfield;
    inline static string* variables;

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
    void ConservedToPrimitive(const real_t* const __restrict__ U, 
                                    real_t* const __restrict__ P,
                              const int mem) const
    {
        static_cast<T const*>(this)->ConservedToPrimitive_(U, P, mem);
        return;
    }

    void PrimitiveToConserved(const real_t* const __restrict__ P, 
                                    real_t* const __restrict__ U,
                              const int mem) const
    {
        static_cast<T const*>(this)->PrimitiveToConserved_(P, U, mem);
        return;
    }

    void WaveSpeeds(const real_t* const __restrict__ P, 
                          real_t* const __restrict__ c,
                    const int dir,
                    const int mem) const
    {
        static_cast<T const*>(this)->WaveSpeeds_(P, c, dir, mem);
        return;
    }

    void Fluxes(const real_t* const __restrict__ P, 
                      real_t (*__restrict__ F)[3],
                const int mem) const
    {
        static_cast<T const*>(this)->Fluxes_(P, F, mem);
        return;
    }

    void DiffusiveFluxes(const real_t* const __restrict__ U, 
                         const real_t (* const __restrict__ dU)[3],
                               real_t (*__restrict__ F)[3],
                         const int mem) const
    {
        static_cast<T const*>(this)->DiffusiveFluxes_(U, dU, F, mem);
        return;
    }
   
    void OrthonormaliseVectors(real_t* const P, const int mem) const
    {
        static_cast<T const*>(this)->OrthonormaliseVectors_(P, mem);
        return;
    }

    void Floors(real_t* const __restrict__ U, const int mem) const 
    {
        static_cast<T const*>(this)->Floors_(U, mem);
        return;
    }

    /* Only used by MHD, and potential future MHD extensions */ 
    //void Fluxes_divB_subsystem(const real_t* const __restrict__ P, 
    //                                 real_t (*__restrict__ F)[3],
    //                           const int mem) const {}
};

#endif
