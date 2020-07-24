#ifndef DOMAIN_MAP_HPP
#define DOMAIN_MAP_HPP

#include "common.hpp"

#include <cmath>


/* Abstract base class for domain mapping functor. 
 * 0 <= n < 12 : index of edge to return map for
 * 0 <= x <= 1 : groupwise reference coordinate along edge
 * r[3] : physical coordinates in 3-space */
class DomainMap
{
    public:
    virtual void operator()(const int n, const real_t x, real_t r[3]) = 0;
};



/*********************************
 * Specific domain mappings
 *********************************/


/* Based on M1 in Kopriva's book, p231 (243) */
class QuarterAnnulusMap : public DomainMap
{
    public:
    const real_t d2 = 0.5; // r2 = [-d2, d2]

    virtual void operator()(const int n, const real_t x, real_t r[3])
    {
        switch (n)
        {
            case 0:
            case 4:
                r[0] = 1 + 2 * x;
                r[1] = 0.0;
                break;
            case 1:
            case 5:
                r[0] = 3 * std::cos(pi_2 * x);
                r[1] = 3 * std::sin(pi_2 * x);
                break;
            case 2:
            case 6:
                r[0] = 0.0;
                r[1] = 1 + 2 * x;
                break;
            case 3:
            case 7:
                r[0] = std::cos(pi_2 * x);
                r[1] = std::sin(pi_2 * x);
                break;
            case 8:
                r[0] = 1.0;
                r[1] = 0.0;
                break;
            case 9:
                r[0] = 3.0;
                r[1] = 0.0;
                break;
            case 10:
                r[0] = 0.0;
                r[1] = 3.0;
                break;
            case 11:
                r[0] = 0.0;
                r[1] = 1.0;
                break;
        }

        if (n < 4)
            r[2] = -d2;
        else if (n < 8)
            r[2] = d2;
        else
            r[2] = d2 * (2 * x - 1);

        return;
    }
};

#endif
