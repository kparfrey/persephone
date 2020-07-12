#ifndef EULER_HPP
#define EULER_HPP

#include <cmath>
#include "physics_functors.hpp"

/* U --- conserved
 * 0 rho
 * 1 E  ---  total energy density
 * 2 rho v_0
 * 3 rho v_1
 * 4 rho v_2
 *
 * P --- primitive
 * 0 rho
 * 1 p  ---  pressure
 * 2 v_0
 * 3 v_1
 * 4 v_2
 */

constexpr real_t gamma_euler = 5.0/3.0;
constexpr real_t gm1_euler = gamma_euler - 1.0;

constexpr static int density    = 0;
constexpr static int tot_energy = 1;
constexpr static int pressure   = 1;
constexpr static int v0   = 2;
constexpr static int v1   = 3;
constexpr static int v2   = 4;
constexpr static int mom0 = 2;
constexpr static int mom1 = 3;
constexpr static int mom2 = 4;


class UtoP_euler : public ConservedToPrimitive
{
    public:
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

        P[pressure] = gm1_euler * (U[tot_energy] - KE_density);

        return;
    }
};


class WaveSpeeds_euler : public WaveSpeedsFromPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ P, 
                                         real_t* const __restrict__ c,
                                   const int dir) const 
    {
        const real_t sound_speed = std::sqrt(gamma_euler * P[pressure] / P[density]);

        c[0] = MAX(0.0, P[dir+v0] + sound_speed);
        c[1] = MIN(0.0, P[dir+v0] - sound_speed);

        return;
    }
};


class Fluxes_euler : public FluxesFromPrimitive
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ P, 
                                         real_t (*__restrict__ F)[3]) const
    {
        const real_t KE_density = 0.5 * P[density] 
                                   * (P[v0]*P[v0] + P[v1]*P[v1] + P[v2]*P[v2]);
        const real_t E = KE_density + P[pressure] / gm1_euler;

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
#endif
