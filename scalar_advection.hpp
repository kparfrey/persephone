#ifndef SCALAR_ADVECTION_HPP
#define SCALAR_ADVECTION_HPP

#include "physics_functors.hpp"

#ifdef __CUDACC__
#define ACCEL_DECORATOR __host__ __device__
#else
#define ACCEL_DECORATOR
#endif


class UtoP_scalar_advection : public ConservedToPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const U, real_t* const P)
    {
        P[0] = U[0];
        return;
    }
};


class WaveSpeeds_scalar_advection: public WaveSpeedsFromPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const P, real_t* const c)
    {
        return;
    }
};


class  Fluxes_scalar_advection: public FluxesFromPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const P, real_t* const F)
    {
        return;
    }
};
#endif
