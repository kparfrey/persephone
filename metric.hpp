#ifndef METRIC_HPP
#define METRIC_HPP

#include "common.hpp"
#include "tensor_field.hpp"


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
    void allocate(const int N);


    public:
    /* Transforms between physical and reference coords */
    TensorField    J; // d(PHYSICAL)/d(REFERENCE)
    TensorField Jinv; // d(REFERENCE)/d(PHYSICAL)

    /* Reference coords */
    TensorField g;     // Spatial 3-metric (dual/covariant)
    TensorField ginv;  // Inverse 3-metric (vector/contravariant)
    ScalarField rdetg; // Square root of 3-metric determinant

    /* Physical coords */
    TensorField gphys; // Spatial 3-metric, in chosen physical coords


    /* Methods */
    void setup(const int Nelem[3], const int Ns_block, const real_t corners[8][3]);
    void transform_twoTensor(const TensorField& T_in, TensorField& T_out, 
                             const int N,
                             const CoordTransDir ctd=phys2ref, 
                             const Components c=covariant);
    
};


#endif
