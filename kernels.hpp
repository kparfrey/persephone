#include "common.hpp"
#include "tensor_field.hpp"

namespace kernels
{
    /* Memory allocation on device --- don't use for host-side only memory */
    real_t* alloc(int n);     // Zero-initialize

    real_t* alloc_raw(int n); // Don't initialize



    void free(real_t* a);

    void add_2_vectors(real_t* __restrict__ v1, real_t* __restrict__ v2, 
                       real_t               a1, real_t               a2, 
                       real_t* __restrict__ result, const int N);

    void soln_to_flux(const VectorField soln2flux, 
                      const real_t* const __restrict__ Q, 
                            VectorField Qf, 
                      const LengthBucket lb);

    /*
    void matrix_vector_product(const real_t* const __restrict__ matrix, 
                               const real_t* const __restrict__ vector, 
                                     real_t* const __restrict__ product, 
                               const LengthBucket lb, const Operation op,
                               const int dir);
    */

}
