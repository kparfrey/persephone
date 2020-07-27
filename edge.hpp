#ifndef EDGE_HPP
#define EDGE_HPP

#include "common.hpp"
#include "tensor_field.hpp"

class Edge
{
    public:
    int id;   // This edge's index, 0 <= id < 12
    int dir;  // The reference-space direction this edge is aligned with
    real_t offset[3]; // Offset of this edge in each dir (unit cube coords)

    int N;                  // Number of interpolation nodes
    real_t* __restrict__ x; // length parameter along the curve, x in [0,1]
    double* __restrict__ w; // Weights for Lagrange interpolation
    VectorField r;          // 3D physical coords at the nodes
    real_t endpoints[2][3]; // Physical 3-space coords of the two endpoints
                            // 0: x = 0, 1: x = 1

    void setup(const int id, const int Nf[3], const VectorField xf);
    void free();
    void eval(const real_t s, real_t* const __restrict__ rs) const;
    void diff(const real_t s, real_t* const __restrict__ drdx) const;
};

#endif
