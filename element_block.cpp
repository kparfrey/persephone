#include <cmath>
#include "element_block.hpp"
#include "kernels.hpp"
#include "write_screen.hpp"
#include "lagrange_polynomials.hpp"
#include "transfinite_map.hpp"
#include "domain_map.hpp"
#include "edge.hpp"


void ElementBlock::setup()
{
    /* No. of flux points in each transform line, in each direction */
    Nf[0] = Ns[0] + 1;
    Nf[1] = Ns[1] + 1;
    Nf[2] = Ns[2] + 1;

    for (int i: dirs)
        Ns_tot[i] = Ns[i] * Nelem[i];

    /* No. of solution points per element */
    Ns_elem = Ns[0] * Ns[1] * Ns[2];

    /* No. of solution points in the block */
    Ns_block = Ns_elem * Nelem_block;
    
    /* No. of flux points in each direction, in an element */
    Nf_dir[0] = Ns[1] * Ns[2] * Nf[0];
    Nf_dir[1] = Ns[2] * Ns[0] * Nf[1];
    Nf_dir[2] = Ns[0] * Ns[1] * Nf[2];

    Nf_elem = Nf_dir[0] + Nf_dir[1] + Nf_dir[2]; // 3 Nf Ns^2 if all directions equal

    /* No. of i-direction flux points in the block */
    for (int i: dirs)
        Nf_dir_block[i] = Nf_dir[i] * Nelem_block;


    /* Setup LengthBucket... */
    lengths.Nfield  = Nfield;
    lengths.Ns_elem = Ns_elem;
    lengths.Ns_block = Ns_block;
    for (int i: dirs)
    {
        lengths.Nelem[i]  = Nelem[i];
        lengths.Ns[i]     = Ns[i];
        lengths.Nf[i]     = Nf[i];
        lengths.Nf_dir[i] = Nf_dir[i];
        lengths.Nf_dir_block[i] = Nf_dir_block[i];
    }


    allocate_on_host();

    set_computational_coords();

    if (geometry == simple_geometry)
    {
        write::message("Setting up coords & metric in simple geometry");
        set_physical_coords_simple();
        metric.setup_simple(Nelem, Ns_block, Nf_dir_block, corners);
    }
    else if (geometry == full_geometry)
    {
        /* Set those parts of the edges which don't change
         * element to element */
        for (int i = 0; i < Nelem_block; ++i)
            for (int j = 0; j < 4; ++j)
                edges[i][j].setup(j, Nf, xf); // Use Lobatto points

        write::message("Setting up coords & metric in full geometry");
        set_physical_coords_full();
        metric.setup_full(*this);
    }
    else
        write::error("Chosen GeometryClass not recognized.");

    fill_spectral_difference_matrices();

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
        soln2flux(i)      = new real_t [matrix_size]; // Nf x Ns matrices
        fluxDeriv2soln(i) = new real_t [matrix_size]; // Ns x Nf matrices
    }

    if (geometry == full_geometry)
    {
        edges = new Edge* [Nelem_block];
        for (int i = 0; i < Nelem_block; ++i)
            edges[i] = new Edge [4]; // 4 edges in 2D
    }

    return;
}


void ElementBlock::free_setup_memory()
{
    if (geometry == full_geometry)
    {
        for (int i = 0; i < Nelem_block; ++i)
        {
            for (int j = 0; j < 4; ++j)
                edges[i][j].free();

            delete[] edges[i];
        }

        delete[] edges;
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


/* Assume that everything is trivially Cartesian. */
void ElementBlock::set_physical_coords_simple()
{
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


    /* Solution points */
    for (int ke = 0; ke < Nelem[2]; ++ke)
    for (int je = 0; je < Nelem[1]; ++je)
    for (int ie = 0; ie < Nelem[0]; ++ie)
    {
        elem_origin[0] = corners[0][0] + ie*dr_elem[0]; 
        elem_origin[1] = corners[0][1] + je*dr_elem[1]; 
        elem_origin[2] = corners[0][2] + ke*dr_elem[2]; 

        elem_idx_1D = id_elem(ie, je, ke);
        mem_offset = elem_idx_1D * Ns_elem;

        for (int k = 0; k < Ns[2]; ++k)
        for (int j = 0; j < Ns[1]; ++j)
        for (int i = 0; i < Ns[0]; ++i)
        {
            mem_loc = ids(i,j,k) + mem_offset;
            rs(0,mem_loc) = elem_origin[0] + xs(0,i) * dr_elem[0];
            rs(1,mem_loc) = elem_origin[1] + xs(1,j) * dr_elem[1];
            rs(2,mem_loc) = elem_origin[2] + xs(2,k) * dr_elem[2];
        }
    }

    
    /* Flux points in each direction 
     * Use transform-direction-relative indices -- d is the transform
     * direction, d1 and d2 are the indices cyclically one and two
     * slots above. 
     * These arrays aren't actually used (yet?) --- Lagrange polynomials
     * etc. are calculated using the reference-space x arrays. Possibly
     * should remove / comment-out these if won't have explicit spatial
     * dependence in the flux calculation. */
    for (int d: dirs)
    {
        int d1 = dir_plus_one[d];
        int d2 = dir_plus_two[d];
        
        for (int ne2 = 0; ne2 < Nelem[d2]; ++ne2)
        for (int ne1 = 0; ne1 < Nelem[d1]; ++ne1)
        for (int ne0 = 0; ne0 < Nelem[d];  ++ne0)
        {
            elem_idx_1D = (ne2*Nelem[d1] + ne1)*Nelem[d] + ne0;
            mem_offset  = elem_idx_1D * Nf_dir[d];

            for (int k = 0; k < Ns[d2]; ++k)
            for (int j = 0; j < Ns[d1]; ++j)
            for (int i = 0; i < Nf[d];  ++i)
            {
                mem_loc = idf(i,j,k,d) + mem_offset;
                rf[d](d, mem_loc) = elem_origin[d]  + xf(d, i) * dr_elem[d];
                rf[d](d1,mem_loc) = elem_origin[d1] + xs(d1,j) * dr_elem[d1];
                rf[d](d2,mem_loc) = elem_origin[d2] + xs(d2,k) * dr_elem[d2];
            }
        }
    }
    

    return;
}




/* Not currently using physical coordinates of flux points, rf.
 * Add code to fill this array if end up needing it. */
void ElementBlock::set_physical_coords_full()
{
    int mem_offset;
    int mem_loc;
    int elem_idx_1D;

    /* Groupwise constructs */
    real_t group_width[3]; // Total width in units where each element is a unit cube
    for (int d: dirs)
        group_width[d] = (real_t) Nproc_group[d]*Nelem[d];

    real_t group_corners[8][3];       // Group corner coords in physical space
    (*map)(0, 0.0, group_corners[0]); // Corner 0 is at groupwise-x=0.0 for edge 0 etc.
    (*map)(1, 0.0, group_corners[1]); 
    (*map)(2, 1.0, group_corners[2]); 
    (*map)(3, 1.0, group_corners[3]); 
    (*map)(4, 0.0, group_corners[4]); 
    (*map)(5, 0.0, group_corners[5]); 
    (*map)(6, 1.0, group_corners[6]); 
    (*map)(7, 1.0, group_corners[7]); 

    /* Elementwise constructs */
    real_t elem_corners[4][3]; // For 2D TF map: r0, r1 coords of 4 corners

    /* Point constructs */
    real_t xg[2];     // Groupwise reference-space coord
    real_t xe[3];     // Elementwise reference-space coord
    //int point_idx[2]; // Elementwise i,j indices of this point
    real_t rp[3];     // Physical coordinates of this point



    /* Solution points */
    for (int ke = 0; ke < Nelem[2]; ++ke)
    for (int je = 0; je < Nelem[1]; ++je)
    for (int ie = 0; ie < Nelem[0]; ++ie)
    {
        elem_idx_1D = id_elem(ie, je, ke);
        mem_offset = elem_idx_1D * Ns_elem;

        Edge* elem_edges = edges[elem_idx_1D];

        /* Use "analytic" transfinite interpolation to interpolate from
         * the group mapping to this element's edge. */
        for (int iedge = 0; iedge < 4; ++iedge)
        {
            Edge& edge = elem_edges[iedge];
            real_t* offset = edge.offset;

            real_t along[3];
            for (int d: dirs)
                along[d] = (edge.dir == d) ? 1.0 : 0.0;

            for (int i = 0; i < edge.N; ++i)
            {
                xg[0] = (group_idx[0]*Nelem[0] + ie + along[0]*edge.x[i] + offset[0]) 
                                                                       / group_width[0];
                xg[1] = (group_idx[1]*Nelem[1] + je + along[1]*edge.x[i] + offset[1]) 
                                                                       / group_width[1];
                analytic_transfinite_map_2D(xg, map, group_corners, rp);
                for (int j: dirs) edge.r(j, i) = rp[j]; 
            }

            /* Storing at Gauss points, so evaluate interpolant to find endpoints */
            //edge.eval(0.0, edge.endpoints[0]); // r(s=0) --> endpoints[0]
            //edge.eval(1.0, edge.endpoints[1]); // r(s=1) --> endpoints[1]

            /* Storing at Lobatto points, directly read off endpoints */
            for (int i = 0; i < 3; ++i)
            {
                edge.endpoints[0][i] = edge.r(i, 0);
                edge.endpoints[1][i] = edge.r(i, edge.N-1);
            }
        }

        for (int d: dirs)
        {
            elem_corners[0][d] = elem_edges[0].endpoints[0][d];
            elem_corners[1][d] = elem_edges[1].endpoints[0][d];
            elem_corners[2][d] = elem_edges[2].endpoints[1][d];
            elem_corners[3][d] = elem_edges[3].endpoints[1][d];
        }
        
        for (int k = 0; k < Ns[2]; ++k)
        for (int j = 0; j < Ns[1]; ++j)
        for (int i = 0; i < Ns[0]; ++i)
        {
            mem_loc = ids(i,j,k) + mem_offset;

            /* Do 2D transfinite map using polynomial interpolation */
            xe[0] = xs(0,i);
            xe[1] = xs(1,j);
            //point_idx[0] = i;
            //point_idx[1] = j;

            polynomial_transfinite_map_2D(xe, elem_edges, elem_corners, rp);

            rs(0,mem_loc) = rp[0];
            rs(1,mem_loc) = rp[1];

            /* For now just set r[2] = 0 --- 2D domain*/
            rs(2,mem_loc) = 0.0;
        }
    }

    return;
}


void ElementBlock::fill_spectral_difference_matrices()
{
    /*** solution -> flux interpolation ***/
    write::message("Filling solution -> flux point matrices");
    for (int i: dirs)
        soln2flux(i) = lagrange::barycentric_interpolation_matrix(xs(i), Ns[i], xf(i), Nf[i]);


    /*** Derivative from flux-point values, interpolation back to solution points */
    write::message("Filling flux-derivative -> solution matrices");
    for (int d: dirs)
    {
        int nf = Nf[d];
        int ns = Ns[d];

        /* Matrix for derivative using flux points, at flux points */
        real_t* D = lagrange::differentiation_matrix(xf(d), nf);

        /* Matrix for interpolating from flux to solution points */
        real_t* flux2soln = lagrange::barycentric_interpolation_matrix(xf(d), nf, xs(d), ns);

        /* Compose the operations to give a single matrix */
        fluxDeriv2soln(d) = lagrange::matrix_matrix_product(flux2soln, D, ns, nf, nf);

        delete[] D;
        delete[] flux2soln;
    }

    return;
}
