#ifndef DOMAIN_MAP_TORUS_HPP
#define DOMAIN_MAP_TORUS_HPP

#include <cmath>
#include "common.hpp"
#include "domain_map.hpp"
#include "torus_mode_pack.hpp"
#include "write_screen.hpp"



/* All versions of the torus map will use this transformation from the origin-centred
 * unit disc to physical space, applying the Fourier modes stored in boundary_modes. */
static void unit_disc_to_physical_space(real_t r[3], TorusModePack &modes)
{
    /* Unit disc coords were stored as Cartesian/cylindrical for ease of testing.
     * Convert back to polar coords... Probably should just do all of the first
     * division of the disc directly in polar coords... */
    real_t r_uds = std::sqrt(r[0]*r[0] + r[1]*r[1]);
    real_t t_uds = 0.0;

    using std::atan;

    if (r_uds > TINY)
    {
        if (r[0] >= 0.0 )
        {
            if (r[1] >= 0.0)
                t_uds = atan(r[0]/r[1]);
            else
                t_uds = pi_2 + atan(-r[1]/r[0]);
        }
        else
        {
            if (r[1] >= 0)
                t_uds = 3.*pi_2 + atan(-r[1]/r[0]);
            else
                t_uds = pi + atan(r[0]/r[1]);
        }
    }


    /* Apply transformation with Fourier modes */
    /* NB: for this choice (R -> sin, Z -> cos) the unit-disc and physical spaces
     * have the same ordering arrangements. In both, theta = 0 is in the positive
     * Z direction, and theta increases clockwise. This may be different to e.g. VMEC */
    real_t R0 = 0.0;
    real_t Z0 = 0.0;
    for (int k = 0; k < modes.Nk; ++k)
    {
        R0 += modes.Rmk[0][k] * std::cos(k*r[2]);
        Z0 += modes.Zmk[0][k] * std::cos(k*r[2]);
    }

    real_t R = 0.0;
    real_t Z = 0.0;
    for (int m = 1; m < modes.Nm; ++m)
        for (int k = 0; k < modes.Nk; ++k)
        {
            R += modes.Rmk[m][k] * std::sin(m*t_uds) * std::cos(k*r[2]);
            Z += modes.Zmk[m][k] * std::cos(m*t_uds) * std::cos(k*r[2]);
        }

    r[0] = r_uds * R + R0;
    r[1] = r_uds * Z + Z0;

    return;
}



/* r[3] = (R, Z, phi) 
 * This map only goes into unit-disc space, but it holds the modes for the subsequent
 * pointwise map. */
class BasicSquareTorusMap : public DomainMap
{
    public:

    const real_t b = 1.0 / std::sqrt(2.0); // For centred square and pi/4 outward divisions
    const real_t a = 0.5 * b; // radius of the circle circumscribing the inner square / sqrt(2)

    /* The 8 labelled corners dividing up the disc face into 5 groups */
    const real_t global_corners[8][2] = {{-a,  a},  // 0
                                         { a,  a},  // 1
                                         { a, -a},  // 2
                                         {-a, -a},  // 3
                                         {-b,  b},  // 4
                                         { b,  b},  // 5
                                         { b, -b},  // 6
                                         {-b, -b}}; // 7

    int group; 
    real_t lc[4][2];     // local corners: the 4 poloidal corners of this group's quad in the R-Z plane
    real_t theta_offset; // The starting "colatitude" angle for the outer curved edges (group > 0)

    TorusModePack boundary_modes; // Used in the second, pointwise_transformation() to physical space

    /* Set the domain_depth using the constructor */
    BasicSquareTorusMap(const int group, TorusModePack boundary_modes)
    : boundary_modes(boundary_modes)
    {
        this->group = group;

        /* The local corners for the disc division are only in the poloida/disc plane */
        for (int i = 0; i < 2; ++i)
        {
            switch (group)
            {
                case 0:
                    lc[0][i] = global_corners[0][i];
                    lc[1][i] = global_corners[1][i];
                    lc[2][i] = global_corners[2][i];
                    lc[3][i] = global_corners[3][i];
                    break;
                case 1:
                    lc[0][i] = global_corners[0][i];
                    lc[1][i] = global_corners[4][i];
                    lc[2][i] = global_corners[5][i];
                    lc[3][i] = global_corners[1][i];
                    theta_offset = - pi_4;
                    break;
                case 2:
                    lc[0][i] = global_corners[1][i];
                    lc[1][i] = global_corners[5][i];
                    lc[2][i] = global_corners[6][i];
                    lc[3][i] = global_corners[2][i];
                    theta_offset = pi_4;
                    break;
                case 3:
                    lc[0][i] = global_corners[2][i];
                    lc[1][i] = global_corners[6][i];
                    lc[2][i] = global_corners[7][i];
                    lc[3][i] = global_corners[3][i];
                    theta_offset = 3 * pi_4;
                    break;
                case 4:
                    lc[0][i] = global_corners[3][i];
                    lc[1][i] = global_corners[7][i];
                    lc[2][i] = global_corners[4][i];
                    lc[3][i] = global_corners[0][i];
                    theta_offset = 5 * pi_4;
                    break;
            }
        }
    }


    void operator()(const int n, const real_t x, real_t r[3]) override
    {
        real_t start = 0.0;
        real_t end   = 1.0;

        /* Do poloidal coordinates first, r[0] and r[1] */
        /* The straight edges in unit-disc space */
        for (int i = 0; i < 2; ++i) // for each coord direction
        {
            switch (n) // switch on edge label
            {
                case 0:
                case 4:
                    start = lc[0][i];
                    end   = lc[1][i];
                    break;
                case 1:
                case 5:
                    start = lc[1][i];
                    end   = lc[2][i];
                    break;
                case 2:
                case 6:
                    start = lc[3][i];
                    end   = lc[2][i];
                    break;
                case 3:
                case 7:
                    start = lc[0][i];
                    end   = lc[3][i];
                    break;
                case 8:
                    start = end = lc[0][i];
                    break;
                case 9:
                    start = end = lc[1][i];
                    break;
                case 10:
                    start = end = lc[2][i];
                    break;
                case 11:
                    start = end = lc[3][i];
                    break;
                default:
                    write::error("Inappropriate edge number passed to BasicSquareTorusMap", destroy);
            }

            r[i] = start * (1.0 - x) + end * x;
        }

        /* Overwrite for the curved outer edges (always edge 1) */
        if ((group > 0) && ((n == 1) || (n == 5)))
        {
            /* edge-1 is always aligned with increasing theta */
            real_t theta = x * pi_2 + theta_offset;
            r[0] = std::sin(theta); // effectively R or X in unit-disc space
            r[1] = std::cos(theta); // effectively Z or Y in unit-disc space
        }


        /* Toroidal coordinate directly from x[2] */
        if (n < 4)
            r[2] = 0.0;      // "front" face
        else if (n < 8)
            r[2] = 2.0 * pi; // "back" face
        else // send groupwise x[2] linearly into phi
            r[2] = x * 2.0 * pi;
        

        return;
    }


    void pointwise_transformation(real_t r[3]) override
    {
        unit_disc_to_physical_space(r, boundary_modes);
    };
};



#endif
