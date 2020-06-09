#include "common.hpp"

namespace kernels
{
    real_t* alloc(int n);

    real_t* alloc_raw(int n);

    void free(real_t* a);

    void add_2_vectors(real_t* __restrict__ v1, real_t* __restrict__ v2, 
                       real_t               a1, real_t               a2, 
                       real_t* __restrict__ result, int N);

}
