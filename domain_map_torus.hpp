#ifndef DOMAIN_MAP_TORUS_HPP
#define DOMAIN_MAP_TORUS_HPP

#include "common.hpp"
#include "domain_map.hpp"

#include <cmath>



/* r[3] = (R, Z, phi) */
class BasicSquareTorusMap : public DomainMap
{
    public:

    const real_t b = 1.0 / std::sqrt(2.0); // For centred square and pi/4 outward divisions
    const real_t a = 0.5 * b; // radius of the circle circumscribing the inner square / sqrt(2)
    const real_t delta_phi = 0.2;

    /* The 8 labelled corners dividing up the disc face into 5 groups */
    const real_t global_corners[8][2] = {{-a,  a},  // 0
                                         { a,  a},  // 1
                                         { a, -a},  // 2
                                         {-a, -a},  // 3
                                         {-b,  b},  // 4
                                         { b,  b},  // 5
                                         { b, -b},  // 6
                                         {-b, -b}}; // 7

    real_t lc[4][2]; // local corners: the 4 corners of this group's quad


    void fill_local_data(const int group)
    {
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
                    break;
                case 2:
                    lc[0][i] = global_corners[1][i];
                    lc[1][i] = global_corners[5][i];
                    lc[2][i] = global_corners[6][i];
                    lc[3][i] = global_corners[2][i];
                    break;
                case 3:
                    lc[0][i] = global_corners[2][i];
                    lc[1][i] = global_corners[6][i];
                    lc[2][i] = global_corners[7][i];
                    lc[3][i] = global_corners[3][i];
                    break;
                case 4:
                    lc[0][i] = global_corners[3][i];
                    lc[1][i] = global_corners[7][i];
                    lc[2][i] = global_corners[4][i];
                    lc[3][i] = global_corners[0][i];
                    break;
            }
        }

        return;
    }


    // Do I actually need n > 3 here??
    virtual void operator()(const int n, const real_t x, real_t r[3])
    {
        real_t start;
        real_t end;

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
            }

            r[i] = start * (1.0 - x) + end * x;
        }

        if (n < 4)
            r[2] = 0.0;
        else if (n < 8)
            r[2] = delta_phi;
        else
            r[2] = x * delta_phi;

        return;
    }
};



#endif
