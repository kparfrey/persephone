#ifndef SCALAR_ADVECTION_HPP
#define SCALAR_ADVECTION_HPP

#include "physics.hpp"


class ScalarAdvection : public Physics
{
    private:
    const real_t wave_speed[3] = {0.7, -0.3, 0.2};

    public:
    ScalarAdvection()
    {
        system = scalar_advection;
        Nfield = 1;
        variables = new string [1];
        variables[0] = "phi";
    }


    void ConservedToPrimitive(const real_t* const __restrict__ U, 
                                    real_t* const __restrict__ P) const override
    {
        P[0] = U[0];
        return;
    }


    void WaveSpeeds(const real_t* const P, real_t* const c,
                    const int dir) const override
    {
        c[0] = MAX(0.0, wave_speed[dir]);
        c[1] = MIN(0.0, wave_speed[dir]);

        return;
    }


    void Fluxes(const real_t* const P, real_t (*F)[3]) const override
    {
        for (int i: dirs)
            F[0][i] = wave_speed[i] * P[0];

        return;
    }


    void DiffusiveFluxes(const real_t* const __restrict__ U, 
                         const real_t (* const __restrict__ dU)[3],
                               real_t (*__restrict__ F)[3],
                         const real_t* const __restrict__ args) const override
    {
        return;
    }
};


/****************************************************************/

#if 0 
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
#endif
