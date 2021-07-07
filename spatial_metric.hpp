#ifndef SPATIAL_METRIC_HPP
#define SPATIAL_METRIC_HPP

#include "common.hpp"

/* ABC, to be specified to diagonal or full metrics */
class SpatialMetric
{
    public:

    PhysicalCoords physical_coords;

    real_t* rdetg;
    real_t* rdetg_deriv[3]; // (1/rdetg) * d_j (rdetg) for each j
                            // Only used in finding stress tensor for viscosity

    /* The active memory location when accessing the metric arrays.
     * Need to explicitly set this before using the metric at each location.
     * This method avoids having to copy into temporary pointwise arrays,
     * and in particular having to know whether the metric is diagonal
     * or full. */
    //int mem;
    
    SpatialMetric(PhysicalCoords physical_coords) 
        : physical_coords(physical_coords) {}

    virtual void allocate_on_host(const int N) = 0;

    //virtual void move_to_device();

    virtual void fill(const real_t* const __restrict__ r,
                      const int mem) = 0;

    virtual void lower(const real_t* const __restrict__ Vu,
                             real_t* const __restrict__ Vl,
                       const int mem) = 0;

    virtual void raise(const real_t* const __restrict__ Vl,
                             real_t* const __restrict__ Vu,
                       const int mem) = 0;

    virtual real_t dot(const real_t* const __restrict__ Au,
                       const real_t* const __restrict__ Bu,
                       const int mem) = 0;

    virtual real_t square(const real_t* const __restrict__ Vu,
                          const int mem) = 0;

    virtual real_t square_cov(const real_t* const __restrict__ Vl,
                              const int mem) = 0;

    /* Convert from an orthonormal basis to contravariant
     * components in a coordinate basis. */
    //virtual void ON_to_upper(const real_t* const __restrict__ Von,
    //                               real_t* const __restrict__ Vu,
    //                         const int i);
};


class DiagonalSpatialMetric : public SpatialMetric
{
    public:

    real_t* g[3];    // g[i] == g_{ii}
    real_t* ginv[3]; // ginv[i] == g^{ii}

    DiagonalSpatialMetric(PhysicalCoords physical_coords) 
        : SpatialMetric(physical_coords) {}
    
    void allocate_on_host(const int N) override
    {
        rdetg = new real_t [N]();

        for (int d: dirs)
        {
            g[d]           = new real_t [N]();
            ginv[d]        = new real_t [N]();
            rdetg_deriv[d] = new real_t [N]();
        }

        return;
    }

        
    void fill(const real_t* const __restrict__ r,
              const int mem) override
    {
        switch(physical_coords)
        {
            /* Should I eventually make a separate Cartesian object
             * that doesn't actually store or use the metric at all? */
            case cartesian:
                g[0][mem] = 1.0;
                g[1][mem] = 1.0;
                g[2][mem] = 1.0;
                rdetg[mem] = 1.0;
                rdetg_deriv[0][mem] = 0.0;
                rdetg_deriv[1][mem] = 0.0;
                rdetg_deriv[2][mem] = 0.0;
                break;
            case cylindrical:
                g[0][mem] = 1.0;
                g[1][mem] = 1.0;
                g[2][mem] = r[0] * r[0]; // Phi direction
                rdetg[mem] = r[0];
                rdetg_deriv[0][mem] = 1.0 / r[0];
                rdetg_deriv[1][mem] = 0.0;
                rdetg_deriv[2][mem] = 0.0;
                break;
        }

        /* Since this is a diagonal metric... */
        for (int d: dirs)
            ginv[d][mem] = 1.0 / g[d][mem];

        return;
    }


    void lower(const real_t* const __restrict__ Vu,
                     real_t* const __restrict__ Vl,
               const int mem) override
    {
        Vl[0] = g[0][mem] * Vu[0];
        Vl[1] = g[1][mem] * Vu[1];
        Vl[2] = g[2][mem] * Vu[2];

        return;
    }


    void raise(const real_t* const __restrict__ Vl,
                     real_t* const __restrict__ Vu,
               const int mem) override
    {
        Vu[0] = ginv[0][mem] * Vl[0];
        Vu[1] = ginv[1][mem] * Vl[1];
        Vu[2] = ginv[2][mem] * Vl[2];

        return;
    }


    real_t dot(const real_t* const __restrict__ Au,
               const real_t* const __restrict__ Bu,
               const int mem) override
    {
        return g[0][mem]*Au[0]*Bu[0] + g[1][mem]*Au[1]*Bu[1] + g[2][mem]*Au[2]*Bu[2];
    }


    /* Takes contravariant/upper/vector components */
    real_t square(const real_t* const __restrict__ Vu,
                  const int mem) override
    {
        return g[0][mem]*Vu[0]*Vu[0] + g[1][mem]*Vu[1]*Vu[1] + g[2][mem]*Vu[2]*Vu[2];
    }


    /* Takes covariant/lower/dual components */
    real_t square_cov(const real_t* const __restrict__ Vl,
                      const int mem) override
    {
        return ginv[0][mem]*Vl[0]*Vl[0] + ginv[1][mem]*Vl[1]*Vl[1] + ginv[2][mem]*Vl[2]*Vl[2];
    }

};

#endif
