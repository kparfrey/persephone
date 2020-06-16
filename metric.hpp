#ifndef METRIC_HPP
#define METRIC_HPP

#include "common.hpp"
#include "kernels.hpp"

/* REFERENCE coords are the 3D Cartesian coordinates local to the
 * reference element --- i.e. x[0-2] in [0,1]
 * Synonymous with COMPUTATIONAL/LOGICAL/ELEMENTAL coordinates */

/* For now, all metric arrays are on solution points only 
 * --- assume that all metric operations can be performed before
 *     the transform to flux points 
 * If it turns out we need the metric at flux points we'll need arrays
 * like real_t* g_flux[3][3][3], where the FIRST array index picks out
 * the transform direction, like with ElementBlock.rf[3][3] */

class Metric
{
    private:
    void allocate(int N);


    public:
    /* Transforms between physical and reference coords */
    real_t*    J[3][3]; // d(PHYSICAL)/d(REFERENCE)
    real_t* Jinv[3][3]; // d(REFERENCE)/d(PHYSICAL)

    /* Reference coords */
    real_t*    g[3][3]; // Spatial 3-metric (dual/covariant)
    real_t* ginv[3][3]; // Inverse 3-metric (vector/contravariant)
    real_t* rdetg;      // Square root of 3-metric determinant

    /* Physical coords */
    real_t* gphys[3][3]; // Spatial 3-metric, in chosen physical coords


    /* Methods */
    void setup(int Nelem[3], int Ns_block, real_t corners[8][3]);
    void transform_twoTensor(real_t* T_in[3][3], real_t* T_out[3][3], int N,
                             CoordTransDir ctd=phys2ref, Components c=covariant);
    
};


#endif
