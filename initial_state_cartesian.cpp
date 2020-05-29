#include "initial_state_cartesian.hpp"
#include "element_block.hpp"

static real_t scalar_function(real_t r[3])
{
    real_t sigma[3] = {0.5, 0.5, 0.5};
    real_t r0[3]    = {0.0, 0.0, 0.0};

    real_t rf;
    real_t f;
    real_t arg = 0.0;
    for (int i: dirs)
        rf   = (r[i] - r0[i]) / sigma[i];
        arg += 0.5 * rf * rf ;
    
    f = exp(-arg);

    return f;
}


void set_scalar(ElementBlock& eb)
{
    int mem_offset;
    int loc;
    real_t r[3];

    for (int ie = 0; ie < eb.Nelem[0]; ++ie)
    for (int je = 0; je < eb.Nelem[1]; ++je)
    for (int ke = 0; ke < eb.Nelem[2]; ++ke)
    {
        mem_offset = eb.id_elem(ie, je, ke) * eb.Ns_elem;

        for (int i = 0; i < eb.Ns[0]; ++i)
        for (int j = 0; j < eb.Ns[1]; ++j)
        for (int k = 0; k < eb.Ns[2]; ++k)
        {
            loc = eb.ids(i,j,k) + mem_offset;
            for (int d: dirs)
                r[d] = eb.rs[d][loc];

            eb.fields[loc] = scalar_function(r);
        }
    }

    return;
}
