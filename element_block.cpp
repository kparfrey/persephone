#include "element_block.hpp"
#include "kernels.hpp"


void ElementBlock::setup()
{
    /* No. of flux points in each transform line, in each direction */
    Nf[0] = Ns[0] + 1;
    Nf[1] = Ns[1] + 1;
    Nf[2] = Ns[2] + 1;

    for (int i: dirs)
        Ns_tot[i] = Ns[i] * Nelem[i];

    /* Cyclic ordering for indexing flux-point arrays */
    Nsf[0] = Ns[1];
    Nsf[1] = Ns[2];
    Nsf[2] = Ns[0];

    /* No. of solution points per element */
    Ns_elem = Ns[0] * Ns[1] * Ns[2];

    /* No. of solution points in the block */
    Ns_block = Ns_elem * Nelem_block;
    
    /* No. of flux points in each direction, in an element */
    Nf_dir[0] = Ns[1] * Ns[2] * Nf[0];
    Nf_dir[1] = Ns[2] * Ns[0] * Nf[1];
    Nf_dir[2] = Ns[0] * Ns[1] * Nf[2];

    Nf_elem = Nf_dir[0] + Nf_dir[1] + Nf_dir[2]; // 3 Nf Ns^2 if all directions equal

    allocate();
    set_computational_coords();
    set_physical_coords();

    /* Set up the metric */
    metric.setup(Nelem, Ns_block, corners);

    return;
}


void ElementBlock::allocate()
{
    for (int i: dirs)
    {
        xs[i] = kernels::alloc(Ns[i]);
        xf[i] = kernels::alloc(Nf[i]);
        rs[i] = kernels::alloc(Ns_block);
    }

    for (int itrans: dirs)
        for (int d: dirs)
            rf[itrans][d] = kernels::alloc(Nelem_block * Nf_dir[itrans]);

    fields = kernels::alloc(Nfield * Ns_block);

    return;
}


void ElementBlock::free()
{
    return;
}


/* This needs to be adapted to work on multiple devices
 * --- presumably should just set on the host and copy to device? */
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


void ElementBlock::set_physical_coords()
{
    /* For now assume that everything is trivially Cartesian.
     * This will all need to be replaced for more general mappings */
    real_t dr_elemblock[3];
    real_t dr_elem[3];
    real_t elem_origin[3]; // i.e. corner-0 coordinates for a given element
    int mem_offset;
    int mem_loc;
    int elem_idx_1D;

    dr_elemblock[0] = corners[1][0] - corners[0][0];
    dr_elemblock[1] = corners[3][1] - corners[0][1];
    dr_elemblock[2] = corners[4][2] - corners[0][2];

    for (int i: dirs)
        dr_elem[i] = dr_elemblock[i] / Nelem[i];

    for (int ie = 0; ie < Nelem[0]; ++ie)
    for (int je = 0; je < Nelem[1]; ++je)
    for (int ke = 0; ke < Nelem[2]; ++ke)
    {
        elem_origin[0] = corners[0][0] + ie*dr_elem[0]; 
        elem_origin[1] = corners[0][1] + je*dr_elem[1]; 
        elem_origin[2] = corners[0][2] + ke*dr_elem[2]; 
        elem_idx_1D = id_elem(ie, je, ke);

        /* Solution points */
        mem_offset = elem_idx_1D * Ns_elem;
        for (int i = 0; i < Ns[0]; ++i)
        for (int j = 0; j < Ns[1]; ++j)
        for (int k = 0; k < Ns[2]; ++k)
        {
            mem_loc = ids(i,j,k) + mem_offset;
            rs[0][mem_loc] = elem_origin[0] + 0.5 * (xs[0][i] + 1.0) * dr_elem[0];
            rs[1][mem_loc] = elem_origin[1] + 0.5 * (xs[1][j] + 1.0) * dr_elem[1];
            rs[2][mem_loc] = elem_origin[2] + 0.5 * (xs[2][k] + 1.0) * dr_elem[2];
        }

        /* Flux points in each direction */

        mem_offset = elem_idx_1D * Nf_dir[0];
        for (int k = 0; k < Ns[2]; ++k)
        for (int j = 0; j < Ns[1]; ++j)
        for (int i = 0; i < Nf[0]; ++i)
        {
            mem_loc = idf(i,j,k,0) + mem_offset;
            rf[0][0][mem_loc] = elem_origin[0] + 0.5 * (xf[0][i] + 1.0) * dr_elem[0];
            rf[0][1][mem_loc] = elem_origin[1] + 0.5 * (xs[1][j] + 1.0) * dr_elem[1];
            rf[0][2][mem_loc] = elem_origin[2] + 0.5 * (xs[2][k] + 1.0) * dr_elem[2];
        }

        mem_offset = elem_idx_1D * Nf_dir[1];
        for (int i = 0; i < Ns[0]; ++i)
        for (int k = 0; k < Ns[2]; ++k)
        for (int j = 0; j < Nf[1]; ++j)
        {
            mem_loc = idf(j,k,i,1) + mem_offset;
            rf[1][0][mem_loc] = elem_origin[0] + 0.5 * (xs[0][i] + 1.0) * dr_elem[0];
            rf[1][1][mem_loc] = elem_origin[1] + 0.5 * (xf[1][j] + 1.0) * dr_elem[1];
            rf[1][2][mem_loc] = elem_origin[2] + 0.5 * (xs[2][k] + 1.0) * dr_elem[2];
        }

        mem_offset = elem_idx_1D * Nf_dir[2];
        for (int j = 0; j < Ns[1]; ++j)
        for (int i = 0; i < Ns[0]; ++i)
        for (int k = 0; k < Nf[2]; ++k)
        {
            mem_loc = idf(k,i,j,2) + mem_offset;
            rf[2][0][mem_loc] = elem_origin[0] + 0.5 * (xs[0][i] + 1.0) * dr_elem[0];
            rf[2][1][mem_loc] = elem_origin[1] + 0.5 * (xs[1][j] + 1.0) * dr_elem[1];
            rf[2][2][mem_loc] = elem_origin[2] + 0.5 * (xf[2][k] + 1.0) * dr_elem[2];
        }
    }


    return;
}
