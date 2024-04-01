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
                                    real_t* const __restrict__ P,
                              const int) const override
    {
        P[0] = U[0];
        return;
    }


    void WaveSpeeds(const real_t* const P, real_t* const c,
                    const int dir, const int) const override
    {
        c[0] = MAX(0.0, wave_speed[dir]);
        c[1] = MIN(0.0, wave_speed[dir]);

        return;
    }


    void Fluxes(const real_t* const P, real_t (*F)[3], const int) const override
    {
        for (int i: dirs)
            F[0][i] = wave_speed[i] * P[0];

        return;
    }


    void DiffusiveFluxes(const real_t* const __restrict__ U, 
                         const real_t (* const __restrict__ dU)[3],
                               real_t (*__restrict__ F)[3],
                         const int) const override
    {
        return;
    }
    

    void OrthonormaliseVectors(real_t* const P, const int mem) const override
    {
        return;
    }


    void Floors(real_t* const __restrict__ U, const int mem) const override
    {
        return;
    }
};
#endif
