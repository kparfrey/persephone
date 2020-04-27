// Strictly speaking this file should import the absolute minimum
// --- it probably shouldn't need anything except the definition of real_t ?
#include "common.hpp"
#include "kernels.hpp"

namespace kernels
{

    real_t* alloc(int n)
    {
        real_t *ptr = new real_t [n];
        return ptr;
    }

}
