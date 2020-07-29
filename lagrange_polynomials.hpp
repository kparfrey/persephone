#ifndef LAGRANGE_POLYNOMIALS_HPP
#define LAGRANGE_POLYNOMIALS_HPP

#include "common.hpp"
#include "tensor_field.hpp"

namespace lagrange
{
    double* barycentric_weights(const real_t* const x, const int N);

    real_t interpolate(const real_t* const __restrict__ x,
                       const real_t* const __restrict__ w,
                       const real_t* const __restrict__ f,
                       const real_t s, const int N);

    void interpolate_vector(const real_t* const __restrict__ x,
                            const real_t* const __restrict__ w,
                            const VectorField V,
                            const real_t s, const int N,
                                  real_t* const __restrict__ Vs);

    real_t differentiate(const real_t* const __restrict__ x,
                         const real_t* const __restrict__ w,
                         const real_t* const __restrict__ f,
                         const real_t s, const int N);

    real_t diff_at_node(const real_t* const __restrict__ x,
                        const real_t* const __restrict__ w,
                        const real_t* const __restrict__ f,
                        const int i, const int N);

    real_t* barycentric_interpolation_matrix(const real_t* const x0, const int N0,
                                             const real_t* const x1, const int N1);

    real_t* differentiation_matrix(const real_t* const x, const int N);

    real_t* matrix_matrix_product(const real_t* const B, const real_t* const C,
                                  const int Nrows_B, const int Ncols_C, const int Nsum);
}

#endif
