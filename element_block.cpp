#include <cmath>
#include "element_block.hpp"
#include "kernels.hpp"
#include "write_screen.hpp"
#include "lagrange_polynomials.hpp"

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
        set_physical_coords_simple();
    else // full_geometry
        set_physical_coords_full();

    fill_spectral_difference_matrices();

    metric.setup(Nelem, Ns_block, Nf_dir_block, corners);

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

/* x is in groupwise reference space, x = [0,1] */
static void analytic_transfinite_map_2D(const real_t x[2], DomainMap* const map,
                                              real_t r[3])
{
    real_t Gx[4][3]; 
    real_t G00[3]; // Gamma_0(0)
    real_t G01[3]; // Gamma_0(1)
    real_t G20[3];
    real_t G21[3];
    
    /* Gamma(xi) or Gamma(eta) */
    (*map)(0, x[0], Gx[0]);
    (*map)(1, x[1], Gx[1]);
    (*map)(2, x[0], Gx[2]);
    (*map)(3, x[1], Gx[3]);

    /* Eventually replace with group corners, passed to this function */
    (*map)(0, 0.0, G00);
    (*map)(0, 1.0, G01);
    (*map)(2, 0.0, G20);
    (*map)(2, 1.0, G21);

    for (int i = 0; i < 2; ++i) // Iterate over coord-vector components
    {
        r[i] =   (1-x[0])*Gx[3][i] + x[0]*Gx[1][i] 
               + (1-x[1])*Gx[0][i] + x[1]*Gx[2][i]
               - (1-x[0])*((1-x[1])*G00[i] + x[1]*G20[i])
               -    x[0] *((1-x[1])*G01[i] + x[1]*G21[i]);
    }

    r[2] = 0.0;

    return;
}


/* x is in elementwise reference space, x = [0,1] */
static void polynomial_transfinite_map_2D(const real_t x[2], const int point_idx[2],
                                          Edge edges[4],
                                          const real_t corners[4][3], real_t r[3])
{
    real_t Gx[4][3]; 
    
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 2; ++j)
            Gx[i][j] = edges[i].r(j, point_idx[edges[i].dir]);
        //edges[i].eval(x[edges[i].dir], Gx[i]);

#if 0
    write::variable<int>("Edge 0 dir:  ", edges[0].dir);
    write::variable<int>("Edge 1 dir:  ", edges[1].dir);
    write::variable<int>("Edge 2 dir:  ", edges[2].dir);
    write::variable<int>("Edge 3 dir:  ", edges[3].dir);
    write::variable<real_t>("Edge-0 length coord:  ", x[edges[0].dir]);
    write::variable<real_t>("Edge-1 length coord:  ", x[edges[1].dir]);
    write::variable<real_t>("Gamma-0 x coord:  ", Gx[0][0]);
    write::variable<real_t>("Gamma-1 x coord:  ", Gx[1][0]);
#endif


    for (int i = 0; i < 2; ++i) // Iterate over coord-vector components
    {
        r[i] =   (1-x[0])*Gx[3][i] + x[0]*Gx[1][i] 
               + (1-x[1])*Gx[0][i] + x[1]*Gx[2][i]
               - (1-x[0])*((1-x[1])*corners[0][i] + x[1]*corners[3][i])
               -    x[0] *((1-x[1])*corners[1][i] + x[1]*corners[2][i]);
    }

    r[2] = 0.0;

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
    (*map)(0, 0.0, group_corners[0]); // Corner 0 is at x=0.0 for edge 0 etc.
    (*map)(1, 0.0, group_corners[1]); 
    (*map)(2, 1.0, group_corners[2]); 
    (*map)(3, 1.0, group_corners[3]); 
    (*map)(4, 0.0, group_corners[4]); 
    (*map)(5, 0.0, group_corners[5]); 
    (*map)(6, 1.0, group_corners[6]); 
    (*map)(7, 1.0, group_corners[7]); 

    /* Elementwise constructs */
    real_t elem_corners[4][3]; // For 2D TF map: r0, r1 coords of 4 corners
    Edge edges[4];             // Need 4 element edges for a 2D map

    /* Set up those parts of edges that don't change between elements */
    for (int i = 0; i < 4; ++i)
        edges[i].setup(i, Ns, xs);

    /* Point constructs */
    real_t xg[2];     // Groupwise reference-space coord
    real_t xe[2];     // Elementwise reference-space coord
    int point_idx[2]; // Elementwise i,j indices of this point
    real_t rp[3];     // Physical coordinates of this point



    /* Solution points */
    for (int ke = 0; ke < Nelem[2]; ++ke)
    for (int je = 0; je < Nelem[1]; ++je)
    for (int ie = 0; ie < Nelem[0]; ++ie)
    {
        elem_idx_1D = id_elem(ie, je, ke);
        mem_offset = elem_idx_1D * Ns_elem;

        /* Use "analytic" transfinite interpolation to interpolate from
         * the group mapping to this element's edge. */
        for (int iedge = 0; iedge < 4; ++iedge)
        {
            Edge& edge = edges[iedge];
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
                analytic_transfinite_map_2D(xg, map, rp);
                for (int j: dirs) edge.r(j, i) = rp[j]; 
            }
        }

        edges[0].eval(0.0, elem_corners[0]);
        edges[1].eval(0.0, elem_corners[1]);
        edges[2].eval(1.0, elem_corners[2]);
        edges[3].eval(1.0, elem_corners[3]);
        
        for (int k = 0; k < Ns[2]; ++k)
        for (int j = 0; j < Ns[1]; ++j)
        for (int i = 0; i < Ns[0]; ++i)
        {
            mem_loc = ids(i,j,k) + mem_offset;

            /* Do 2D transfinite map using polynomial interpolation */
            xe[0] = xs(0,i);
            xe[1] = xs(1,j);
            point_idx[0] = i;
            point_idx[1] = j;
            polynomial_transfinite_map_2D(xe, point_idx, edges, elem_corners, rp);

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
