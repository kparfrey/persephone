#include "common.hpp"

namespace kernels
{
    real_t* alloc(int n);

    real_t* alloc_raw(int n);

    void free(real_t* a);
}
