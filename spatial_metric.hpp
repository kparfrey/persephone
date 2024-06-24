#ifndef SPATIAL_METRIC_HPP
#define SPATIAL_METRIC_HPP

#include <cmath>


template <class T>
class SpatialMetric
{
    public:

    PhysicalCoords* physical_coords; // Create and set in the specialized params classes

    /* Actual metric arrays live in the derived classes */
    real_t* rdetg;
    real_t* rdetg_deriv[3]; // (1/rdetg) * d_j (rdetg) for each j
                            // Only used in finding stress tensor for viscosity

    /* The active memory location when accessing the metric arrays.
     * Need to explicitly set this before using the metric at each location.
     * This method avoids having to copy into temporary pointwise arrays,
     * and in particular having to know whether the metric is diagonal
     * or full. */
    //int mem;
    
    //SpatialMetric(PhysicalCoords physical_coords) 
    //    : physical_coords(physical_coords) {}

    void allocate_on_host(const int N)
    {
        rdetg = new real_t [N]();
        for (int d: dirs)
            rdetg_deriv[d] = new real_t [N]();

        static_cast<T *>(this)->allocate_on_host_(N);
        return;
    }

    void fill(const real_t* const __restrict__ r, const int mem) const
    {
        static_cast<T const*>(this)->fill_(r, mem);
        return;
    }

    void lower(const real_t* const __restrict__ Vu,
                     real_t* const __restrict__ Vl,
               const int mem) const
    {
        static_cast<T const*>(this)->lower_(Vu, Vl, mem);
        return;
    }

    void raise(const real_t* const __restrict__ Vl,
                     real_t* const __restrict__ Vu,
               const int mem) const
    {
        static_cast<T const*>(this)->raise_(Vl, Vu, mem);
        return;
    }

    real_t dot(const real_t* const __restrict__ Au,
               const real_t* const __restrict__ Bu,
               const int mem) const
    {
        return static_cast<T const*>(this)->dot_(Au, Bu, mem);
    }

    real_t square(const real_t* const __restrict__ Vu,
                  const int mem) const
    {
        return static_cast<T const*>(this)->square_(Vu, mem);
    }

    real_t square_cov(const real_t* const __restrict__ Vl,
                      const int mem) const
    {
        return static_cast<T const*>(this)->square_cov_(Vl, mem);
    }

    void orthonormals(real_t* const Vu, const int mem) const
    {
        static_cast<T const*>(this)->orthonormals_(Vu, mem);
        return;
    }

    /* Convert from an orthonormal basis to contravariant
     * components in a coordinate basis. */
    //virtual void ON_to_upper(const real_t* const __restrict__ Von,
    //                               real_t* const __restrict__ Vu,
    //                         const int i);
};


class DiagonalSpatialMetric : public SpatialMetric<DiagonalSpatialMetric>
{
    public:
    
    /* Actual metric arrays live in the specialized classes, since the
     * non-diagonal form would need e.g. g[3][3] */
    real_t* g[3];    // g[i] == g_{ii}
    real_t* ginv[3]; // ginv[i] == g^{ii}


    //DiagonalSpatialMetric(PhysicalCoords physical_coords) 
    //    : SpatialMetric(physical_coords) {}
    
    void allocate_on_host_(const int N)
    {
        for (int d: dirs)
        {
            g[d]           = new real_t [N]();
            ginv[d]        = new real_t [N]();
        }

        return;
    }

        
    void fill_(const real_t* const __restrict__ r, const int mem) const
    {
        switch(*physical_coords)
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


    void lower_(const real_t* const __restrict__ Vu,
                      real_t* const __restrict__ Vl,
                const int mem) const
    {
        Vl[0] = g[0][mem] * Vu[0];
        Vl[1] = g[1][mem] * Vu[1];
        Vl[2] = g[2][mem] * Vu[2];

        return;
    }


    void raise_(const real_t* const __restrict__ Vl,
                      real_t* const __restrict__ Vu,
                const int mem) const
    {
        Vu[0] = ginv[0][mem] * Vl[0];
        Vu[1] = ginv[1][mem] * Vl[1];
        Vu[2] = ginv[2][mem] * Vl[2];

        return;
    }


    real_t dot_(const real_t* const __restrict__ Au,
                const real_t* const __restrict__ Bu,
                const int mem) const
    {
        return g[0][mem]*Au[0]*Bu[0] + g[1][mem]*Au[1]*Bu[1] + g[2][mem]*Au[2]*Bu[2];
    }


    /* Takes contravariant/upper/vector components */
    real_t square_(const real_t* const __restrict__ Vu,
                   const int mem) const
    {
        return g[0][mem]*Vu[0]*Vu[0] + g[1][mem]*Vu[1]*Vu[1] + g[2][mem]*Vu[2]*Vu[2];
    }


    /* Takes covariant/lower/dual components */
    real_t square_cov_(const real_t* const __restrict__ Vl,
                       const int mem) const
    {
        return ginv[0][mem]*Vl[0]*Vl[0] + ginv[1][mem]*Vl[1]*Vl[1] + ginv[2][mem]*Vl[2]*Vl[2];
    }


    /* For writing orthonormal components to disk 
     * Save o.n. components back into input array */
    void orthonormals_(real_t* const Vu, const int mem) const
    {
        for (int d: dirs)
            Vu[d] *= std::sqrt(g[d][mem]);

        return;
    }

};

#endif
