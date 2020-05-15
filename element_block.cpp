#include "element_block.hpp"
#include "kernels.hpp"

void ElementBlock::allocate()
{
    for (int i: dirs)
    {
        xs[i] = kernels::alloc(Ns[i]);
        xf[i] = kernels::alloc(Nf[i]);
    }

    rs = kernels::alloc(Nelem_block * Ns_elem);
    rf = kernels::alloc(Nelem_block * Nf_elem);

    fields = kernels::alloc(Nfield * Nelem_block * Ns_elem);

    return;
}


void ElementBlock::free()
{
    return;
}


void ElementBlock::set_computational_coords()
{
    /* Points are defined on [-1, 1], so that xf[0] = -1 */

    /* Solution / Gauss points */
    for (int d: dirs)
        for (int i = 0; i < Ns[d]; ++i)
            xs[d][i] = - cos(pi * (i + 0.5) / Ns[d]);

    /* Flux / Lobatto points */
    for (int d: dirs)
        for (int i = 0; i < Nf[d]; ++i)
            xf[d][i] = - cos(pi * i  / (Nf[d] - 1.0));

    return;
}
