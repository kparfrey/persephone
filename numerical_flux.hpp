#ifndef NUMERICAL_FLUX_HPP
#define NUMERICAL_FLUX_HPP

#include "common.hpp"
#include "physics_includes.hpp"

/* In general, create the physical-coordinate numerical flux in all
 * three directions. The dir argument is only used when using straight
 * elements, in which cases only the flux in that direction is nonzero. */
class NumericalFlux
{
    public:
    ConservedToPrimitive*    U_to_P;
    WaveSpeedsFromPrimitive* c_from_P;
    FluxesFromPrimitive*     F_from_P;

    int Nfield;

    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ UL, 
                                   const real_t* const __restrict__ UR,
                                         real_t (*Fnum)[3],
                                   const int dir) const = 0;
};


/* Doesn't include the projection along the normal required
 * when using curved elements. Only calculates the flux in the
 * dir direction. Not fully optimized, since still calculating
 * all three components of FL and FR. */
class HLL_straight : public NumericalFlux
{
    public:
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ UL, 
                                   const real_t* const __restrict__ UR,
                                         real_t (*Fnum)[3],
                                   const int dir) const
    {
        /* Temporary variables --- use a fixed size */
        constexpr int Nv = 10; // since Nfield isn't constexpr
        real_t PL[Nv],    PR[Nv];
        real_t FL[Nv][3], FR[Nv][3];
        real_t cL[2], cR[2];

        (*U_to_P)(UL, PL);
        (*U_to_P)(UR, PR);

        (*c_from_P)(PL, cL, dir);
        (*c_from_P)(PR, cR, dir);

        (*F_from_P)(PL, FL);
        (*F_from_P)(PR, FR);

        //int cmax =   MAX(0.0, MAX(cL[0], cR[0])); // +ve moving waves
        //int cmin = - MIN(0.0, MIN(cL[1], cR[1])); // -ve moving waves

        for (int i = 0; i < Nfield; ++i)
        {
            /*
            Fnum[i][dir] = (cmin*FR[i][dir] + cmax*FL[i][dir]
                             - cmax*cmin*(UR[i] - UL[i])) / (cmax + cmin);
            Fnum[i][dir_plus_one[dir]] = 0.0;
            Fnum[i][dir_plus_two[dir]] = 0.0;
            */
            
            for (int j: dirs)
                Fnum[i][j] = FL[i][j];
                //Fnum[i][j] = 0.5 * (FL[i][j] + FR[i][j]);
        }


        return;
    }
};
#endif
