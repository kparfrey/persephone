#ifndef SCALAR_ADVECTION_HPP
#define SCALAR_ADVECTION_HPP

#include "physics_functors.hpp"

//#include <string>


static real_t wave_speed_scalar_advection[3] = {0.7, -0.3, 0.2};



class SystemData_scalar_advection : public SystemData
{
    public:
    SystemData_scalar_advection()
    {
        Nfield = 1;
        variables = new string [1];
        variables[0] = "phi";
    }
};


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


/* Returns the signed fastest wave speeds in this direction, moving 
 * in the positive direction (c[0]) and negative direction (c[1]) */
class WaveSpeeds_scalar_advection: public WaveSpeedsFromPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const P, real_t* const c,
                                   const int dir) const
    {
        c[0] = MAX(0.0, wave_speed_scalar_advection[dir]);
        c[1] = MIN(0.0, wave_speed_scalar_advection[dir]);

        return;
    }
};


class  Fluxes_scalar_advection: public FluxesFromPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const P, real_t (*F)[3]) const
    {
        for (int i: dirs)
            F[0][i] = wave_speed_scalar_advection[i] * P[0];

        return;
    }
};
#endif
