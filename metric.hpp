#ifndef METRIC_HPP
#define METRIC_HPP

#include "common.hpp"
#include "tensor_field.hpp"

class Edge;
class ElementBlock;

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
    void allocate_on_host(const int Ns, const int Nf[3]);


    public:
    /* On solution points... */
    ScalarField Jrdetg; // |J| sqrt{det(g)}


    /* ...and on flux points */
    VectorField S[3];  // Physical->Reference vector transform
                       // S[i] = Jrdetg d ref^i/d phys^j 
    //TensorField g[3];  // Physical components
    // Will need others eventually



    /* Transforms between physical and reference coords */
    //TensorField    J; // d(PHYSICAL)/d(REFERENCE)
    //TensorField Jinv; // d(REFERENCE)/d(PHYSICAL)

    /* Reference coords */ 
    //TensorField g;     // Spatial 3-metric (dual/covariant)
    //TensorField ginv;  // Inverse 3-metric (vector/contravariant)
    //ScalarField rdetg; // Square root of 3-metric determinant
    
    /* ...and on flux points */
    //TensorField     J_f[3];
    //TensorField  Jinv_f[3];
    //TensorField     g_f[3];
    //TensorField  ginv_f[3];
    //ScalarField rdetg_f[3]; 


    /* Physical coords */
    //TensorField gphys; // Spatial 3-metric, in chosen physical coords
    //TensorField gphys_f[3];


    /* Methods */
    void setup_simple(const int Nelem[3], const int Ns_block, 
                      const int Nf_dir_block[3], const real_t corners[8][3]);
    void setup_full(ElementBlock& eb);
    void move_to_device();
    
};


#endif
