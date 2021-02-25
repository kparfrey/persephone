#ifndef DOMAIN_MAP_TORUS_HPP
#define DOMAIN_MAP_TORUS_HPP

#include "common.hpp"
#include "domain_map.hpp"

#include <cmath>
#include "write_screen.hpp"



/* r[3] = (R, Z, phi) 
 * This map only goes into unit-disc space. Really a cylinder for now. */
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
    real_t lc[4][2]; // local corners: the 4 corners of this group's quad
    real_t theta_offset; // The starting "colatitude" angle for the outer curved edges (group > 0)

    /* Set the domain_depth using the constructor */
    //BasicSquareTorusMap(): DomainMap(0.5) {}

    /* Replace by a constructor? */
    virtual void fill_local_data(const int group)
    {
        this->group = group;

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

        return;
    }


    // Do I actually need n > 3 here?? Since am only ever doing 2D TFI,
    // should only need four edges
    virtual void operator()(const int n, const real_t x, real_t r[3])
    {
        real_t start;
        real_t end;

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
        if ((group > 0) && (n == 1))
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
};



#endif
