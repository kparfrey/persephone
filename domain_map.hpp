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

    virtual void pointwise_transformation(real_t r[3]){};
};



/************************************************
 * Specific domain mappings for Cartesian domains
 ************************************************/


/* Based on M1 in Kopriva's book, p231 (243) */
class QuarterAnnulusMap : public DomainMap
{
    public:
    const real_t d2 = 0.5; // r2 = [-d2, d2]

    void operator()(const int n, const real_t x, real_t r[3]) override
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


class BasicRect2D : public DomainMap
{
    public:
    const real_t d2 = 0.5; // r2 = [-d2, d2]
    real_t c[4][2]; // x,y coords of the rectangle's 4 corners

    const real_t xmin = 0.0;
    const real_t xmax = 1.0/std::cos(pi/6.0);
    const real_t ymin = 0.0;
    const real_t ymax = 1.0/std::sin(pi/6.0);


    BasicRect2D()
    {
        c[0][0] =   xmin;
        c[0][1] =   ymin;

        c[1][0] =   xmax;
        c[1][1] =   ymin; 

        c[2][0] =   xmax; 
        c[2][1] =   ymax; 

        c[3][0] =   xmin; 
        c[3][1] =   ymax;
    }


    void operator()(const int n, const real_t x, real_t r[3]) override
    {
        for (int i = 0; i < 2; ++i)
        {
            switch (n)
            {
                case 0:
                case 4:
                    r[i] = c[0][i] + x * (c[1][i] - c[0][i]);
                    break;
                case 1:
                case 5:
                    r[i] = c[1][i] + x * (c[2][i] - c[1][i]);
                    break;
                case 2:
                case 6:
                    r[i] = c[3][i] + x * (c[2][i] - c[3][i]);
                    break;
                case 3:
                case 7:
                    r[i] = c[0][i] + x * (c[3][i] - c[0][i]);
                    break;
                case 8:
                    r[i] = c[0][i];
                    break;
                case 9:
                    r[i] = c[1][i];
                    break;
                case 10:
                    r[i] = c[2][i];
                    break;
                case 11:
                    r[i] = c[3][i];
                    break;
            }
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


class ObliqueRect2D : public DomainMap
{
    public:
    const real_t d2 = 0.5; // r2 = [-d2, d2]
    real_t c[4][2]; // x,y coords of the rectangle's 4 corners


    ObliqueRect2D()
    {
        c[0][0] =   0.0;
        c[0][1] =   0.0;

        c[1][0] =   1.0;
        c[1][1] =   0.0; //4.0;

        c[2][0] =   1.0; //9.0;
        c[2][1] =   1.0; //6.0;

        c[3][0] =   0.0; //1.0;
        c[3][1] =   1.0;
    }


    void operator()(const int n, const real_t x, real_t r[3]) override
    {
        for (int i = 0; i < 2; ++i)
        {
            switch (n)
            {
                case 0:
                case 4:
                    r[i] = c[0][i] + x * (c[1][i] - c[0][i]);
                    break;
                case 1:
                case 5:
                    r[i] = c[1][i] + x * (c[2][i] - c[1][i]);
                    break;
                case 2:
                case 6:
                    r[i] = c[3][i] + x * (c[2][i] - c[3][i]);
                    break;
                case 3:
                case 7:
                    r[i] = c[0][i] + x * (c[3][i] - c[0][i]);
                    break;
                case 8:
                    r[i] = c[0][i];
                    break;
                case 9:
                    r[i] = c[1][i];
                    break;
                case 10:
                    r[i] = c[2][i];
                    break;
                case 11:
                    r[i] = c[3][i];
                    break;
            }
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


class WaveRect2D : public DomainMap
{
    public:
    const real_t d2 = 0.5; // r2 = [-d2, d2]
    real_t c[4][2]; // x,y coords of the rectangle's 4 corners
    real_t A[4][2] = {}; // Amplitude of the sine wave in x & y for the 4 edges
    real_t nwave[4] = {}; // No. of wavelengths for each of the four edges 


    WaveRect2D()
    {
        /* Corners */
        c[0][0] = - 5.0;
        c[0][1] = - 5.0;

        c[1][0] =   5.0;
        c[1][1] = - 4.0;

        c[2][0] =   9.0;
        c[2][1] =   6.0;

        c[3][0] = - 1.0;
        c[3][1] =   5.0;

        /* Wave amplitudes */
        A[0][0] = 0.5;
        A[0][1] = 1.0;
        A[2][0] = 0.5;
        A[2][1] = 1.0;

        A[1][0] = -0.25;
        A[1][1] = 0.25;
        A[3][0] = -0.25;
        A[3][1] = 0.25;

        /* Numbers of wavelengths */
        nwave[0] = 1.0;
        nwave[2] = 1.0;

        nwave[1] = 2.0;
        nwave[3] = 2.0;
    }


    void operator()(const int n, const real_t x, real_t r[3]) override
    {
        for (int i = 0; i < 2; ++i)
        {
            switch (n)
            {
                case 0:
                case 4:
                    r[i] = c[0][i] + x * (c[1][i] - c[0][i])
                             + A[0][i]*std::sin(x*nwave[0]*2*pi);
                    break;
                case 1:
                case 5:
                    r[i] = c[1][i] + x * (c[2][i] - c[1][i])
                             + A[1][i]*std::sin(x*nwave[1]*2*pi);
                    break;
                case 2:
                case 6:
                    r[i] = c[3][i] + x * (c[2][i] - c[3][i])
                             + A[2][i]*std::sin(x*nwave[2]*2*pi);
                    break;
                case 3:
                case 7:
                    r[i] = c[0][i] + x * (c[3][i] - c[0][i])
                             + A[3][i]*std::sin(x*nwave[3]*2*pi);
                    break;
                case 8:
                    r[i] = c[0][i];
                    break;
                case 9:
                    r[i] = c[1][i];
                    break;
                case 10:
                    r[i] = c[2][i];
                    break;
                case 11:
                    r[i] = c[3][i];
                    break;
            }
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
