#ifndef LAGRANGE_POLYNOMIALS_HPP
#define LAGRANGE_POLYNOMIALS_HPP

#include "common.hpp"

namespace lagrange
{
    real_t* barycentric_interpolation_matrix(const real_t* const x0, const int N0,
                                             const real_t* const x1, const int N1);

    real_t* differentiation_matrix(const real_t* const x, const int N);

    real_t* matrix_matrix_product(const real_t* const B, const real_t* const C,
                                  const int Nrows_B, const int Ncols_C, const int Nsum);
}

#endif
