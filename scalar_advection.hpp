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
    inline virtual void operator()(const real_t* const U, real_t* const P) const
    {
        P[0] = U[0];
        return;
    }
};


class WaveSpeeds_scalar_advection: public WaveSpeedsFromPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const P, real_t* const c) const
    {
        return;
    }
};


class  Fluxes_scalar_advection: public FluxesFromPrimitive
{
    public:
    real_t wave_speed[3] = {0.7, 0.3, 0.0};

    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const P, real_t (*F)[3]) const
    {
        for (int i: dirs)
            F[0][i] = wave_speed[i] * P[0];

        return;
    }
};
#endif
