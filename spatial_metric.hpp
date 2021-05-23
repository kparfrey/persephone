#ifndef SPATIAL_METRIC_HPP
#define SPATIAL_METRIC_HPP

#include "common.hpp"

/* ABC, to be specified to diagonal or full metrics */
class SpatialMetric
{
    public:

    MetricCoords metric_coords;

    
    SpatialMetric(MetricCoords metric_coords) 
        : metric_coords(metric_coords) {}

    virtual void allocate_on_host(const int N) = 0;

    //virtual void move_to_device();

    virtual void fill(const real_t* const __restrict__ r,
                      const int i) = 0;

    virtual void lower(const real_t* const __restrict__ Vu,
                             real_t* const __restrict__ Vl,
                       const int i) = 0;

    /* Convert from an orthonormal basis to contravariant
     * components in a coordinate basis. */
    //virtual void ON_to_upper(const real_t* const __restrict__ Von,
    //                               real_t* const __restrict__ Vu,
    //                         const int i);
};


class DiagonalSpatialMetric : public SpatialMetric
{
    public:

    real_t* g[3]; // g[i] == g_{ii}

    DiagonalSpatialMetric(MetricCoords metric_coords) 
        : SpatialMetric(metric_coords) {}
    
    void allocate_on_host(const int N)
    {
        for (int d: dirs)
            g[d] = new real_t [N]();

        return;
    }

        
    void fill(const real_t* const __restrict__ r,
              const int i)
    {
        switch(metric_coords)
        {
            case Cartesian:
                g[0][i] = 1.0;
                g[1][i] = 1.0;
                g[2][i] = 1.0;
                break;
            case cylindrical:
                g[0][i] = 1.0;
                g[1][i] = 1.0;
                g[2][i] = r[0] * r[0]; // Phi direction
                break;
        }

        return;
    }


    void lower(const real_t* const __restrict__ Vu,
                     real_t* const __restrict__ Vl,
               const int i)
    {
        Vl[0] = g[0][i] * Vu[0];
        Vl[1] = g[1][i] * Vu[1];
        Vl[2] = g[2][i] * Vu[2];

        return;
    }
};

#endif
