#ifndef NUMERICAL_FLUX_HPP
#define NUMERICAL_FLUX_HPP

#include "common.hpp"
#include "physics_includes.hpp"

#include "write_screen.hpp"

/* In general, create the physical-coordinate numerical flux in all
 * three directions. The dir argument is only used when using straight
 * elements, in which cases only the flux in that direction is nonzero. */
class NumericalFlux
{
    public:
    Physics* physics;
    //ConservedToPrimitive*    U_to_P;
    //WaveSpeedsFromPrimitive* c_from_P;
    //FluxesFromPrimitive*     F_from_P;

    int Nfield;

    /* Unique flux for both sides of the interface */
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ UL, 
                                   const real_t* const __restrict__ UR,
                                   const real_t* const __restrict__ n,
                                         real_t (*Fnum)[3],
                                   const int dir,
                                   const int mem) const {}

    /* Only creates a unique normal flux --- tangential flux is discontinuous */
    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ UL, 
                                   const real_t* const __restrict__ UR,
                                   const real_t* const __restrict__ n,
                                         real_t (*Fnum_L)[3],
                                         real_t (*Fnum_R)[3],
                                   const int dir,
                                   const int mem) const {}
};


/* This currently averages the tangential fluxes, to give a unique 
 * flux in all three physical components for both sides. May want to
 * reorganise so that each side uses its own tangential flux, since
 * only the normal flux needs to be the same at the interface. */
class HLL : public NumericalFlux
{
    public:

    ACCEL_DECORATOR
    inline void operator()(const real_t* const __restrict__ UL, 
                           const real_t* const __restrict__ UR,
                           const real_t* const __restrict__ n,
                                 real_t (*Fnum)[3],
                           const int, // dir not used
                           const int mem) const override
    {
        /* Temporary variables --- use a fixed size */
        constexpr int Nv = 10; // since Nfield isn't constexpr
        real_t PL[Nv],    PR[Nv];
        real_t FL[Nv][3], FR[Nv][3];
        real_t cL[3][2], cR[3][2];
        real_t cpos[3], cneg[3];
        real_t cp, cn;
        
        real_t Fstar;
        real_t F0; // Here (FL + FR)/2 for now

        //(*U_to_P)(UL, PL);
        //(*U_to_P)(UR, PR);
        physics->ConservedToPrimitive(UL, PL, mem);
        physics->ConservedToPrimitive(UR, PR, mem);

        for (int d: dirs)
        {
            //(*c_from_P)(PL, cL[i], i);
            //(*c_from_P)(PR, cR[i], i);
            physics->WaveSpeeds(PL, cL[d], d, mem);
            physics->WaveSpeeds(PR, cR[d], d, mem);
        }

        /* Signed max wavespeeds in each direction.
         * cpos is max speed parallel to the coord axis (L to R); ie Toro's S_R.
         * cneg is max speed anti-parallel to the coord axis (R to L),
         * keeping the negative sign; ie Toro's S_L.
         * The outer MAX(0,..), MIN(0,..) is to allow the use of the general
         * expression in supersonic cases. Otherwise would need an explicit
         * test for eg if (cneg >= 0.0) Fnum = FL etc. */
        for (int d: dirs)
        {
            cpos[d] = MAX(0.0, MAX(cL[d][0], cR[d][0])); // fastest +ve moving wave
            cneg[d] = MIN(0.0, MIN(cL[d][1], cR[d][1])); // fastest -ve moving wave
        }


        physics->Fluxes(PL, FL, mem);
        physics->Fluxes(PR, FR, mem);

        for (int field = 0; field < Nfield; ++field)
            for (int i: dirs)
            {
                cp = cpos[i];
                cn = cneg[i];

                Fstar = (cp*FL[field][i] - cn*FR[field][i]
                                 + cp*cn*(UR[field] - UL[field])) / (cp - cn);

                F0 = 0.5 * (FL[field][i] + FR[field][i]);


                Fnum[field][i] = F0 * (1.0 - n[i]) + Fstar * n[i];
            }

        return;
    }


    /* This only makes the normal flux consistent across the interface --- the tangential
     * fluxes are left unchanged. */
    ACCEL_DECORATOR
    inline void operator()(const real_t* const __restrict__ UL, 
                           const real_t* const __restrict__ UR,
                           const real_t* const __restrict__ n,
                                 real_t (*Fnum_L)[3],
                                 real_t (*Fnum_R)[3],
                           const int, // dir not used
                           const int mem) const override
    {
        /* Temporary variables --- use a fixed size */
        constexpr int Nv = 10; // since Nfield isn't constexpr
        real_t PL[Nv],    PR[Nv];
        real_t FL[Nv][3], FR[Nv][3];
        real_t cL[3][2], cR[3][2];
        real_t cpos[3], cneg[3];
        real_t cp, cn;
        
        real_t Fstar;

        physics->ConservedToPrimitive(UL, PL, mem);
        physics->ConservedToPrimitive(UR, PR, mem);

        for (int d: dirs)
        {
            physics->WaveSpeeds(PL, cL[d], d, mem);
            physics->WaveSpeeds(PR, cR[d], d, mem);
        }

        physics->Fluxes(PL, FL, mem);
        physics->Fluxes(PR, FR, mem);

        /* Signed max wavespeeds in each direction.
         * cpos is max speed parallel to the coord axis (L to R); ie Toro's S_R.
         * cneg is max speed anti-parallel to the coord axis (R to L),
         * keeping the negative sign; ie Toro's S_L.
         * The outer MAX(0,..), MIN(0,..) is to allow the use of the general
         * expression in supersonic cases. Otherwise would need an explicit
         * test for eg if (cneg >= 0.0) Fnum = FL etc. */
        for (int d: dirs)
        {
            cpos[d] = MAX(0.0, MAX(cL[d][0], cR[d][0])); // fastest +ve moving wave
            cneg[d] = MIN(0.0, MIN(cL[d][1], cR[d][1])); // fastest -ve moving wave
        }


        for (int field = 0; field < Nfield; ++field)
            for (int i: dirs)
            {
                cp = cpos[i];
                cn = cneg[i];

                Fstar = (cp*FL[field][i] - cn*FR[field][i]
                                 + cp*cn*(UR[field] - UL[field])) / (cp - cn);

                Fnum_L[field][i] = FL[field][i] * (1.0 - n[i]) + Fstar * n[i];
                Fnum_R[field][i] = FR[field][i] * (1.0 - n[i]) + Fstar * n[i];
            }

        return;
    }

};


/***
class HLL_divB_subsystem : public NumericalFlux
{
    public:

    ACCEL_DECORATOR
    inline void operator()(const real_t* const __restrict__ UL, 
                           const real_t* const __restrict__ UR,
                           const real_t* const __restrict__ n,
                                 real_t (*Fnum)[3],
                           const int, // dir not used
                           const int mem) const override
    {
        // Temporary variables --- use a fixed size 
        constexpr int Nv = 10; // since Nfield isn't constexpr
        real_t FL[Nv][3], FR[Nv][3];
        real_t PL[Nv],    PR[Nv];
        real_t cL[3][2], cR[3][2];
        real_t cpos[3], cneg[3];
        real_t cp, cn;
        
        real_t Fstar;
        real_t F0; // Here (FL + FR)/2 for now

        // P = U, so use U wherever P would usually go

        physics->ConservedToPrimitive(UL, PL, mem);
        physics->ConservedToPrimitive(UR, PR, mem);

        physics->Fluxes_divB_subsystem(PL, FL, mem);
        physics->Fluxes_divB_subsystem(PR, FR, mem);

        for (int d: dirs)
        {
            physics->WaveSpeeds(PL, cL[d], d, mem);
            physics->WaveSpeeds(PR, cR[d], d, mem);
        }
        
        for (int d: dirs)
        {
            cpos[d] = MAX(0.0, MAX(cL[d][0], cR[d][0])); // fastest +ve moving wave
            cneg[d] = MIN(0.0, MIN(cL[d][1], cR[d][1])); // fastest -ve moving wave
        }


        for (int field = 5; field <= 8; ++field)
            for (int i: dirs)
            {
                cp = cpos[i];
                cn = cneg[i];

                Fstar = (cp*FL[field][i] - cn*FR[field][i]
                                 + cp*cn*(UR[field] - UL[field])) / (cp - cn);

                F0 = 0.5 * (FL[field][i] + FR[field][i]);

                Fnum[field][i] = F0 * (1.0 - n[i]) + Fstar * n[i];
            }

        return;
    }
};
***/




/* Doesn't include the projection along the normal required
 * when using curved elements. Only calculates the flux in the
 * dir direction. Not fully optimized, since still calculating
 * all three components of FL and FR. */
class HLL_straight : public NumericalFlux
{
    public:
    ACCEL_DECORATOR
    inline void operator()(const real_t* const __restrict__ UL, 
                           const real_t* const __restrict__ UR,
                           const real_t* const __restrict__ n, // not used
                                 real_t (*Fnum)[3],
                           const int dir,
                           const int mem) const override
    {
        /* Temporary variables --- use a fixed size */
        constexpr int Nv = 10; // since Nfield isn't constexpr
        real_t PL[Nv],    PR[Nv];
        real_t FL[Nv][3], FR[Nv][3];
        real_t cL[2], cR[2];

#if 0
        (*U_to_P)(UL, PL);
        (*U_to_P)(UR, PR);
        (*c_from_P)(PL, cL, dir);
        (*c_from_P)(PR, cR, dir);
        (*F_from_P)(PL, FL);
        (*F_from_P)(PR, FR);
#endif
        physics->ConservedToPrimitive(UL, PL, mem);
        physics->ConservedToPrimitive(UR, PR, mem);
        physics->WaveSpeeds(PL, cL, dir, mem);
        physics->WaveSpeeds(PR, cR, dir, mem);
        physics->Fluxes(PL, FL, mem);
        physics->Fluxes(PR, FR, mem);

        /* Signed max wavespeeds in each direction.
         * cpos is max speed parallel to the coord axis (L to R); ie Toro's S_R.
         * cneg is max speed anti-parallel to the coord axis (R to L),
         * keeping the negative sign; ie Toro's S_L.
         * The outer MAX(0,..), MIN(0,..) is to allow the use of the general
         * expression in supersonic cases. Otherwise would need an explicit
         * test for eg if (cneg >= 0.0) Fnum = FL etc. */
        real_t cpos = MAX(0.0, MAX(cL[0], cR[0])); // fastest +ve moving wave
        real_t cneg = MIN(0.0, MIN(cL[1], cR[1])); // fastest -ve moving wave

        for (int i = 0; i < Nfield; ++i)
        {
            Fnum[i][dir] = (cpos*FL[i][dir] - cneg*FR[i][dir]
                          + cpos*cneg*(UR[i] - UL[i])) / (cpos - cneg + TINY);
            
            /* Fully straight elements (ie all angles pi/2) --- ref-space
             * fluxes are just a rescaling of the physical fluxes in the same
             * direction, and the other directions are unnecessary. */
            Fnum[i][dir_plus_one[dir]] = 0.0;
            Fnum[i][dir_plus_two[dir]] = 0.0;
        }

        return;
    }
};
#endif
