#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include <cmath>
#include "physics_functors.hpp"

/* U --- conserved
 * 0 rho
 * 1 rho v_0
 * 2 rho v_1
 * 3 rho v_2
 * 4 E  ---  total energy density
 *
 * P --- primitive
 * 0 rho
 * 1 v_0
 * 2 v_1
 * 3 v_2
 * 4 p  ---  pressure
 */

constexpr real_t gamma_navstokes = 5.0/3.0;
constexpr real_t gm1_navstokes = gamma_navstokes - 1.0;


class SystemData_navstokes : public SystemData
{
    public:
    SystemData_navstokes()
    {
        Nfield = 5;
        variables = new string [5];
        variables[0] = "rho";
        variables[1] = "v0";
        variables[2] = "v1";
        variables[3] = "v2";
        variables[4] = "p";

        diffusive = true;
        viscosity = 1e-2;
    }
};


class UtoP_navstokes : public ConservedToPrimitive
{
    public:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure  };

    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ U, 
                                         real_t* const __restrict__ P) const
    {
        P[density] = U[density]; 

        P[v0] = U[mom0] / U[density];
        P[v1] = U[mom1] / U[density];
        P[v2] = U[mom2] / U[density];

        const real_t KE_density = 0.5 * P[density] 
                                   * (P[v0]*P[v0] + P[v1]*P[v1] + P[v2]*P[v2]);

        P[pressure] = gm1_navstokes * (U[tot_energy] - KE_density);

        return;
    }
};


class WaveSpeeds_navstokes : public WaveSpeedsFromPrimitive
{
    public:
    enum primitive {density, v0, v1, v2, pressure};

    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ P, 
                                         real_t* const __restrict__ c,
                                   const int dir) const 
    {
        const real_t sound_speed = std::sqrt(gamma_navstokes * P[pressure] / P[density]);

        c[0] = MAX(0.0, P[v0+dir] + sound_speed);
        c[1] = MIN(0.0, P[v0+dir] - sound_speed);

        return;
    }
};


class Fluxes_navstokes : public FluxesFromPrimitive
{
    public:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure  };

    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ P, 
                                         real_t (*__restrict__ F)[3]) const
    {
        const real_t KE_density = 0.5 * P[density] 
                                   * (P[v0]*P[v0] + P[v1]*P[v1] + P[v2]*P[v2]);
        const real_t E = KE_density + P[pressure] / gm1_navstokes;

        real_t v;

        for (int d: dirs)
        {
            v = P[v0+d]; // velocity in this direction

            F[density][d]    = v * P[density];
            F[tot_energy][d] = v * (E + P[pressure]);

            F[mom0][d] = P[density] * v * P[v0];
            F[mom1][d] = P[density] * v * P[v1];
            F[mom2][d] = P[density] * v * P[v2];

            F[mom0+d][d] += P[pressure];
        }

        return;
    }
};


/* This is the flux that should go into the flux divergence on the
 * "left hand side": i.e. dU/dt + divF = 0, F = F_advective + F_diffusive 
 * Many references define it on the RHS, so the negative of this. */
class DiffusiveFluxes_navstokes : public DiffusiveFluxes
{
    public:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure  };

    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ U, 
                                   const real_t (* const __restrict__ dU)[3],
                                         real_t (*__restrict__ F)[3],
                                   const real_t* const __restrict__ args) const 
    {
        const real_t mu = args[0]; // dynamic viscosity
        const real_t lambda = - (2.0/3.0) * mu; // bulk viscosity
        
        real_t v[3];
        real_t tau[3][3];
        real_t dv[3][3]; // dv[component][deriv dir]
        real_t divv;

        for (int d: dirs)
            v[d] = U[mom0+d] / U[Density];

        for (int comp: dirs)
            for (int deriv: dirs)
                dv[comp][deriv] = (dU[mom0+comp][deriv] - v[comp]*dU[Density][deriv])
                                                                             / U[Density];
        divv = dv[0][0] + dv[1][1] + dv[2][2]; // Assuming Cartesian...

        for (int d: dirs)
            tau[d][d] = 2*mu*dv[d][d] + lambda*divv;

        tau[0][1] = tau[1][0] = mu * (dv[0][1] + dv[1][0]);
        tau[0][2] = tau[2][0] = mu * (dv[0][2] + dv[2][0]);
        tau[1][2] = tau[2][1] = mu * (dv[1][2] + dv[2][1]);

        for (int d: dirs) // flux direction
        {
            F[Density][d] = 0.0;

            for (int i: dirs)
                F[mom0+i][d] = - tau[d][i];

            F[tot_energy][d] = - (v[0]*tau[0][d] + v[1]*tau[1][d] + v[2]*tau[2][d]);
        }

        return;
    }
};
#endif
