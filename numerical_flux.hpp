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


#if 0
        const real_t Bdotn_L = n[0]*UL[5] + n[1]*UL[6] + n[2]*UL[7];        
        const real_t Bdotn_R = n[0]*UR[5] + n[1]*UR[6] + n[2]*UR[7];        
        real_t nu[3]; // contravariant components of normal vector
        physics->metric->raise(n, nu, mem);

        real_t Bn_L, Bn_R, Bt_L, Bt_R;
        real_t Bn_m;
        real_t B_L[3], B_R[3]; 
        
        for (int i: dirs)
        {
            Bn_L = Bdotn_L * nu[i];
            Bn_R = Bdotn_R * nu[i];
            Bt_L = UL[5+i] - Bdotn_L*nu[i];
            Bt_R = UR[5+i] - Bdotn_R*nu[i];

            Bn_m = 0.5 * (Bn_L + Bn_R);

            B_L[i] = Bn_m + Bt_L;
            B_R[i] = Bn_m + Bt_R;
        }
        
        const real_t psi_m = 0.5 * (UR[8] + UL[8]);
        //const real_t psi_m = 0.5 * (UR[8] + UL[8]) - 0.5 * Physics::ch * (Bdotn_R - Bdotn_L);

        /**/
        for (int i: dirs)
        {
            PL[5+i] = B_L[i];
            PR[5+i] = B_R[i];
        }
        /**/

        PL[8] = PR[8] = psi_m;
#endif

        physics->Fluxes(PL, FL, mem);
        physics->Fluxes(PR, FR, mem);

        /*
        const real_t Mdotn_L = n[0]*UL[1] + n[1]*UL[2] + n[2]*UL[3];        
        const real_t Mdotn_R = n[0]*UR[1] + n[1]*UR[2] + n[2]*UR[3];        

        if ((Mdotn_L > 0) && (Mdotn_R > 0))
        {
            for (int field = 0; field < Nfield; ++field)
                for (int i: dirs)
                    Fnum[field][i] = FL[field][i];
        }
        else if ((Mdotn_L < 0) && (Mdotn_R < 0))
        {
            for (int field = 0; field < Nfield; ++field)
                for (int i: dirs)
                    Fnum[field][i] = FR[field][i];
        }
        else
        { // Use HLL is the momenta are of opposite signs?
        */
        for (int field = 0; field < Nfield; ++field)
            for (int i: dirs)
            {
                cp = cpos[i];
                cn = cneg[i];

                Fstar = (cp*FL[field][i] - cn*FR[field][i]
                                 + cp*cn*(UR[field] - UL[field])) / (cp - cn);

                F0 = 0.5 * (FL[field][i] + FR[field][i]);

                /*
                if ((Mdotn_L > 0) && (Mdotn_R > 0))
                    Fstar = FL[field][i];
                else if ((Mdotn_L < 0) && (Mdotn_R < 0))
                    Fstar = FR[field][i];
                    */

                Fnum[field][i] = F0 * (1.0 - n[i]) + Fstar * n[i];
            }
        //}

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
        /* Temporary variables --- use a fixed size */
        constexpr int Nv = 10; // since Nfield isn't constexpr
        real_t FL[Nv][3], FR[Nv][3];
        const real_t cp =   Physics::ch;
        const real_t cn = - Physics::ch;
        
        real_t Fstar;
        real_t F0; // Here (FL + FR)/2 for now

        /* P = U, so use U wherever P would usually go */

        physics->Fluxes_divB_subsystem(UL, FL, mem);
        physics->Fluxes_divB_subsystem(UR, FR, mem);

        for (int field = 5; field <= 8; ++field)
            for (int i: dirs)
            {
                Fstar = (cp*FL[field][i] - cn*FR[field][i]
                                 + cp*cn*(UR[field] - UL[field])) / (cp - cn + TINY);

                F0 = 0.5 * (FL[field][i] + FR[field][i]);

                Fnum[field][i] = F0 * (1.0 - n[i]) + Fstar * n[i];
            }

        return;
    }
};




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
