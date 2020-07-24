#ifndef EDGE_HPP
#define EDGE_HPP

#include "common.hpp"
#include "tensor_field.hpp"

/* Just a placeholder class for now for generally curved edge object */
/* Think you should only need real Edge objects around processes/element blocks
 * --- within a block you can find all the physical locations and derivatives
 *     by linear blending inwards from the outer edges. */
class Edge
{
    public:
    int dir; // The reference-space direction this edge is aligned with
    int N;                  // Number of interpolation nodes
    real_t* __restrict__ x; // length parameter along the curve, x in [0,1]
    double* __restrict__ w; // Weights for Lagrange interpolation
    VectorField r;          // 3D physical coords at the nodes

    void setup(const int Ns, const real_t* const __restrict__ xs);
    void eval(const real_t s, real_t* const __restrict__ rs);
    void diff(const real_t s, real_t* const __restrict__ drdx);
};

#endif
