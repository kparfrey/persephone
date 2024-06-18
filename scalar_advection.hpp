#ifndef SCALAR_ADVECTION_HPP
#define SCALAR_ADVECTION_HPP

#include "physics.hpp"


class ScalarAdvection : public Physics<ScalarAdvection>
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


    void ConservedToPrimitive_(const real_t* const __restrict__ U, 
                                     real_t* const __restrict__ P,
                               const int) const 
    {
        P[0] = U[0];
        return;
    }


    void PrimitiveToConserved_(const real_t* const __restrict__ P, 
                                     real_t* const __restrict__ U,
                               const int) const 
    {
        U[0] = P[0];
        return;
    }


    void WaveSpeeds_(const real_t* const P, real_t* const c,
                     const int dir, const int) const 
    {
        c[0] = MAX(0.0, wave_speed[dir]);
        c[1] = MIN(0.0, wave_speed[dir]);

        return;
    }


    void Fluxes_(const real_t* const P, real_t (*F)[3], const int) const 
    {
        for (int i: dirs)
            F[0][i] = wave_speed[i] * P[0];

        return;
    }


    void DiffusiveFluxes_(const real_t* const __restrict__ U, 
                          const real_t (* const __restrict__ dU)[3],
                                real_t (*__restrict__ F)[3],
                          const int) const 
    {
        return;
    }
    

    void OrthonormaliseVectors_(real_t* const P, const int mem) const 
    {
        return;
    }


    void Floors_(real_t* const __restrict__ U, const int mem) const 
    {
        return;
    }
};
#endif
