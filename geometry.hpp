#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "common.hpp"
#include "tensor_field.hpp"

class Edge;
class ElementBlock;

/* REFERENCE coords are the 3D Cartesian coordinates local to the
 * reference element --- i.e. x[0-2] in [0,1]
 * Synonymous with COMPUTATIONAL/LOGICAL/ELEMENTAL coordinates */

/* This class stores the computational geometry, arising from the 
 * deformation of the reference elements. The physical geometry,
 * ie the spatial metric, is located in the Physics class. */
class Geometry
{
    private:
    void allocate_on_host(const int Ns, const int Nf[3]);


    public:
    /* On solution points... */
    ScalarField Jrdetg;  // |J| sqrt{det(g)}
    VectorField dxdr[3]; // dxdr[i](j,...) = d x^i / d r^j 
                         // Only used to find grad(U) for diffusive flux


    /* ...and on flux points */
    VectorField S[3];  // Physical->Reference vector transform
                       // Has rdetg because used to calculate divF
                       // S[i] = |J| rdetg d ref^i/d phys^j 
    VectorField timestep_transform[3]; // Transforms physical v^i to 
                                       // ref v^j and includes dx^j 
                                       // tt[i] = (1/Dx^i) d x^i/ d r^j
 

    /* At element faces */
    VectorField normal[3]; // One VectorField for each transform dir


    /* Methods */
    void setup_full(ElementBlock& eb); // Rename without _full
    void move_to_device();
    
};
#endif
