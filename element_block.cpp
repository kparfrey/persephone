#include <cmath>
#include "element_block.hpp"
#include "kernels.hpp"
#include "write_screen.hpp"


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

    allocate_on_host();
    set_computational_coords();
    set_physical_coords();
    fill_spectral_difference_matrices();

    metric.setup(Nelem, Ns_block, corners);

    return;
}


void ElementBlock::allocate_on_host()
{
    for (int i: dirs)
    {
        xs(i) = new real_t [Ns[i]];
        xf(i) = new real_t [Nf[i]];
        rs(i) = new real_t [Ns_block];
    }

    for (int itrans: dirs)
        for (int d: dirs)
            rf[itrans](d) = new real_t [Nelem_block * Nf_dir[itrans]];

    fields = new real_t [Nfield * Ns_block];

    for (int i: dirs)
    {
        int matrix_size = Ns[i] * Nf[i]; 
        soln2flux[i]      = new real_t [matrix_size]; // Nf x Ns matrices
        fluxDeriv2soln[i] = new real_t [matrix_size]; // Ns x Nf matrices
    }

    return;
}


void ElementBlock::move_to_device()
{
    write::message("Copying ElementBlock's data to device");

    metric.move_to_device();

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
    /* Points are defined on [0, 1], so that xf[0] = 0 */

    /* Solution / Gauss points */
    for (int d: dirs)
        for (int i = 0; i < Ns[d]; ++i)
            xs(d,i) = 0.5 * (1.0 - std::cos(pi * (i + 0.5) / Ns[d]));

    /* Flux / Lobatto points */
    for (int d: dirs)
        for (int i = 0; i < Nf[d]; ++i)
            xf(d,i) = 0.5 * (1.0 - std::cos(pi * i  / (Nf[d] - 1.0)));

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
            rs(0,mem_loc) = elem_origin[0] + xs(0,i) * dr_elem[0];
            rs(1,mem_loc) = elem_origin[1] + xs(1,j) * dr_elem[1];
            rs(2,mem_loc) = elem_origin[2] + xs(2,k) * dr_elem[2];
        }

        /* Flux points in each direction */

        mem_offset = elem_idx_1D * Nf_dir[0];
        for (int k = 0; k < Ns[2]; ++k)
        for (int j = 0; j < Ns[1]; ++j)
        for (int i = 0; i < Nf[0]; ++i)
        {
            mem_loc = idf(i,j,k,0) + mem_offset;
            rf[0](0,mem_loc) = elem_origin[0] + xf(0,i) * dr_elem[0];
            rf[0](1,mem_loc) = elem_origin[1] + xs(1,j) * dr_elem[1];
            rf[0](2,mem_loc) = elem_origin[2] + xs(2,k) * dr_elem[2];
        }

        mem_offset = elem_idx_1D * Nf_dir[1];
        for (int i = 0; i < Ns[0]; ++i)
        for (int k = 0; k < Ns[2]; ++k)
        for (int j = 0; j < Nf[1]; ++j)
        {
            mem_loc = idf(j,k,i,1) + mem_offset;
            rf[1](0,mem_loc) = elem_origin[0] + xs(0,i) * dr_elem[0];
            rf[1](1,mem_loc) = elem_origin[1] + xf(1,j) * dr_elem[1];
            rf[1](2,mem_loc) = elem_origin[2] + xs(2,k) * dr_elem[2];
        }

        mem_offset = elem_idx_1D * Nf_dir[2];
        for (int j = 0; j < Ns[1]; ++j)
        for (int i = 0; i < Ns[0]; ++i)
        for (int k = 0; k < Nf[2]; ++k)
        {
            mem_loc = idf(k,i,j,2) + mem_offset;
            rf[2](0,mem_loc) = elem_origin[0] + xs(0,i) * dr_elem[0];
            rf[2](1,mem_loc) = elem_origin[1] + xs(1,j) * dr_elem[1];
            rf[2](2,mem_loc) = elem_origin[2] + xf(2,k) * dr_elem[2];
        }
    }


    return;
}


double* ElementBlock::barycentric_weights(const real_t* const x, const int N)
{
    double* weights = new double [N];
    
    double p;
    for (int j = 0; j < N; ++j)
    {
        p = 1.0;
        for (int i = 0; i < N; ++i)
        {
            if (i != j)
                p *= x[j] - x[i];
        }

        weights[j] = 1.0 / p;
    }


    return weights;
}


/* Single interpolation matrix: from x0 points --> x1 points */
real_t* ElementBlock::barycentric_interpolation_matrix(const real_t* const x0, const int N0,
                                                       const real_t* const x1, const int N1)
{
    /* Store like N1 x N0 matrix: the 0/original-point index moves fastest */
    real_t* matrix = new real_t [N1*N0];

    double* w = barycentric_weights(x0, N0);

    double s;
    for (int k = 0; k < N1; ++k)
    {
        s = 0.0;
        for (int i = 0; i < N0; ++i)
            s += w[i] / (x1[k] - x0[i]);

        for (int j = 0; j < N0; ++j)
            matrix[k*N0 + j] = w[j] / (s * (x1[k] - x0[j]));
    }

    delete[] w;

    return matrix;
}


void ElementBlock::fill_spectral_difference_matrices()
{
    /*** solution -> flux interpolation ***/
    write::message("Filling solution -> flux point matrices");
    for (int i: dirs)
        soln2flux[i] = barycentric_interpolation_matrix(xs(i), Ns[i], xf(i), Nf[i]);


    /*** Derivative from flux-point values, interpolation back to solution points */
    write::message("Filling flux-derivative -> solution matrices");
    for (int d: dirs)
    {
        int nf = Nf[d];
        int ns = Ns[d];

        double* w = barycentric_weights(xf(d), nf);

        /* Matrix for derivative using flux points, at flux points */
        double* flux_deriv = new double [nf*nf]();
        
        /* i != j entries */
        for (int j = 0; j < nf; ++j)
            for (int i = 0; i < nf; ++i)
                if (i != j)
                    flux_deriv[i*nf + j] = (w[j]/w[i]) / (xf(d,i) - xf(d,j));

        /* i == j entries */
        double s;
        for (int j = 0; j < nf; ++j)
        {
            s = 0.0;
            for (int i = 0; i < nf; ++i)
                if (i != j)
                    s += flux_deriv[j*nf + i];
        
            flux_deriv[j*nf + j] = - s;
        }

        /* Matrix for interpolating from flux to solution points */
        real_t* flux2soln = barycentric_interpolation_matrix(xf(d), nf, xs(d), ns);

        /* Compose the operations to give a single matrix */
        //fluxDeriv2soln[d] = matrix_matrix_product(....)

        delete[] flux_deriv;
        delete[] flux2soln;
    }

    return;
}
