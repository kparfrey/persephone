#include <cmath>
#include "common.hpp"
#include "element_block.hpp"
#include "kernels.hpp"
#include "write_screen.hpp"
#include "lagrange_polynomials.hpp"
#include "transfinite_map.hpp"
#include "domain_map.hpp"
#include "domain_map_torus.hpp"
#include "edge.hpp"
#include "geometry_labels.hpp"
#include "legendre_roots.hpp"


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

    /* Set those parts of the edges which don't change
     * element to element */
    for (int i = 0; i < Nelem_block; ++i)
        for (int j = 0; j < 12; ++j)
            edges[i][j].setup(j, Nf, xf); // Use Lobatto points

    write::message("Setting up coords and metric");
    set_physical_coords_full();

    write::message("Setting up the grid geometry");
    geometry.setup_full(*this);

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
        
        rs_pre_transform(i) = new real_t [Ns_block];
    }

    //for (int itrans: dirs)
    //    for (int d: dirs)
    //        rf[itrans](d) = new real_t [Nelem_block * Nf_dir[itrans]];

    fields = new real_t [Nfield * Ns_block];

    physics_soln.metric->allocate_on_host(Ns_block);

    if (physics_soln.system == mhd)
    {
        divB = new real_t [Ns_block];
        //divB_init = new real_t [Ns_block];
    }

    for (int i: dirs)
    {
        int matrix_size = Ns[i] * Nf[i]; 
        soln2flux(i)      = new real_t [matrix_size]; // Nf x Ns matrices
        fluxDeriv2soln(i) = new real_t [matrix_size]; // Ns x Nf matrices
    
        physics[i].metric->allocate_on_host(Nf_dir_block[i]);
    }

    edges = new Edge* [Nelem_block];
    for (int i = 0; i < Nelem_block; ++i)
        edges[i] = new Edge [12]; // Each 3D element has 12 edges 

    return;
}


void ElementBlock::free_setup_memory()
{
    for (int i = 0; i < Nelem_block; ++i)
    {
        for (int j = 0; j < 12; ++j)
            edges[i][j].free();

        delete[] edges[i];
    }

    delete[] edges;


    for (int i: dirs)
        delete[] rs_pre_transform(i);

    return;
}


void ElementBlock::move_to_device()
{
    write::message("Copying ElementBlock's data to device");

    geometry.move_to_device();

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

    /* Solution points --- Gauss-Chebyshev */
    for (int d: dirs)
        for (int i = 0; i < Ns[d]; ++i)
            xs(d,i) = 0.5 * (1.0 - std::cos(pi * (i + 0.5) / Ns[d]));

    /* Flux points 
     * Lobatto-Chebyshev points --- weakly unstable */
    //for (int d: dirs)
    //    for (int i = 0; i < Nf[d]; ++i)
    //        xf(d,i) = 0.5 * (1.0 - std::cos(pi * i  / (Nf[d] - 1.0)));

    /* Flux points --- Gauss-Legendre */
    for (int d: dirs)
    {
        int n = Ns[d] - 1; // Order of Legendre polynomial & no. of interior points needed
        xf(d,       0) = 0.0;
        xf(d, Nf[d]-1) = 1.0;
        legendre::find_roots(xf(d), n);
    }

    return;
}


/* Assume that everything is trivially Cartesian. */
/* Save for a reminder of how to set the flux point locations */
#if 0
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
#endif // 0



/* Not currently using physical coordinates of flux points, rf.
 * Add code to fill this array if end up needing it. 
 * This also fills the SpatialMetric on both solution and flux points */
void ElementBlock::set_physical_coords_full()
{
    int mem_offset;
    int mem_loc;
    int elem_idx_1D;

    /* Groupwise constructs */
    real_t group_width[3]; // Total width in units where each element is a unit cube
    for (int d: dirs)
        group_width[d] = (real_t) Nproc_group[d]*Nelem[d];

    /* Group corner coords in physical space
     * Use edge_to_corner[8] to find corner n's coord (0 or 1) along edge n */
    real_t group_corners[8][3];       
    for (int i: icorners)
        (*map)(i, double(edge_to_corner[i]), group_corners[i]); 
    
    /* Elementwise constructs */
    Edge* elem_edges; // Pointer to the first of this element's edges
    real_t elem_corners[8][3]; // r[3] coords of an element's 8 corners

    /* Point constructs */
    real_t xg[3];     // Groupwise reference-space coord
    real_t xe[3];     // Elementwise reference-space coord
    real_t rp[3];     // Physical coordinates of this point


    /* Solution points */
    for (int ke = 0; ke < Nelem[2]; ++ke)
    for (int je = 0; je < Nelem[1]; ++je)
    for (int ie = 0; ie < Nelem[0]; ++ie)
    {
        elem_idx_1D = id_elem(ie, je, ke);
        mem_offset = elem_idx_1D * Ns_elem;

        /* Pointer to this element's edges */
        elem_edges = edges[elem_idx_1D];

        /* Use "analytic" transfinite interpolation to interpolate from
         * the group mapping to this element's edge. */
        for (int iedge: iedges)
        {
            /* Reference to this particular edge */
            Edge& edge = elem_edges[iedge];
            real_t* offset = edge.offset;

            real_t along[3];
            for (int d: dirs)
                along[d] = (edge.dir == d) ? 1.0 : 0.0;

            /* Find the coordinates along each edge */
            for (int i = 0; i < edge.N; ++i)
            {
                xg[0] = (group_idx[0]*Nelem[0] + ie + along[0]*edge.x[i] + offset[0]) 
                                                                       / group_width[0];
                xg[1] = (group_idx[1]*Nelem[1] + je + along[1]*edge.x[i] + offset[1]) 
                                                                       / group_width[1];
                xg[2] = (group_idx[2]*Nelem[2] + ke + along[2]*edge.x[i] + offset[2]) 
                                                                       / group_width[2];
                
                /* First do the group-boundary to element-boundary interpolation */
                /* "Analytic" because the group's boundary is given by functions */
                /* For the torus this maps reference space to unit-disc space    */
                analytic_transfinite_map_3D(xg, map, group_corners, rp);

                /* Then apply an additional coordinate map to each point, on a
                 * point-by-point basis. For the torus, this maps unit-disc space
                 * to physical space. */
                map->pointwise_transformation(rp);
                
                for (int j: dirs) edge.r(j, i) = rp[j]; 
            }

            /* Storing at Lobatto points, directly read off endpoints */
            for (int j: dirs)
            {
                edge.endpoints[0][j] = edge.r(j, 0);
                edge.endpoints[1][j] = edge.r(j, edge.N-1);
            }
        }

        for (int i: icorners)
        for (int d: dirs)
            elem_corners[i][d] = elem_edges[i].endpoints[edge_to_corner[i]][d];

        for (int k = 0; k < Ns[2]; ++k)
        for (int j = 0; j < Ns[1]; ++j)
        for (int i = 0; i < Ns[0]; ++i)
        {
            mem_loc = ids(i,j,k) + mem_offset;
            
            xe[0] = xs(0,i);
            xe[1] = xs(1,j);
            xe[2] = xs(2,k);

            /* Do 3D transfinite map using polynomial interpolation */
            polynomial_transfinite_map_3D(xe, elem_edges, elem_corners, rp);

            for (int d: dirs)
                rs(d,mem_loc) = rp[d];

            /* Fill in the spatial metric arrays */
            physics_soln.metric->fill(rp, mem_loc);
        }


        /* Additionally, do an analytic map to the "pre-transform" "physical" coords,
         * since these are sometimes needed for setup. For the torus, this is equivalent
         * to storing the unit-disc-space coords of each solution point */
        for (int k = 0; k < Ns[2]; ++k)
        for (int j = 0; j < Ns[1]; ++j)
        for (int i = 0; i < Ns[0]; ++i)
        {
            mem_loc = ids(i,j,k) + mem_offset;
            
            /* Elementwise computational-space coords */
            xe[0] = xs(0,i);
            xe[1] = xs(1,j);
            xe[2] = xs(2,k);

            /* Groupwise computational-space coords */
            xg[0] = (group_idx[0]*Nelem[0] + ie + xe[0]) / group_width[0];
            xg[1] = (group_idx[1]*Nelem[1] + je + xe[1]) / group_width[1];
            xg[2] = (group_idx[2]*Nelem[2] + ke + xe[2]) / group_width[2];

            analytic_transfinite_map_3D(xg, map, group_corners, rp);

            cartesian_to_polar(rp); // Store (rho, theta, phi) into rs_pre_transform

            for (int d: dirs)
                rs_pre_transform(d,mem_loc) = rp[d];
        }
    }


    /* Flux points */
    int id_elem_f, id_elem_s;
    int mem_offset_f; 

    for (int d: dirs)
    {
        int ne0, ne1, ne2; 
        int* ie = nullptr;
        int* je = nullptr;
        int* ke = nullptr;

        int d1 = dir_plus_one[d];
        int d2 = dir_plus_two[d];
        
        
        /* Identical to kernels::relative_to_fixed_indices(). Reproduced here so don't
         * need to include kernels. Split off to a function if I end up needing it 
         * somewhere else. --- Duplicated from Geometry::setup(), should put 
         * somewhere else to avoid so much duplication */
        switch(d)
        {
            case 0:
                ie  = &ne0;
                je  = &ne1;
                ke  = &ne2;
                break;
            case 1:
                ie  = &ne2;
                je  = &ne0;
                ke  = &ne1;
                break;
            case 2:
                ie  = &ne1;
                je  = &ne2;
                ke  = &ne0;
                break;
        }

        for (ne2 = 0; ne2 < Nelem[d2]; ++ne2)
        for (ne1 = 0; ne1 < Nelem[d1]; ++ne1)
        for (ne0 = 0; ne0 < Nelem[d];  ++ne0)
        {
            id_elem_f    = (ne2*Nelem[d1] + ne1)*Nelem[d] + ne0;
            mem_offset_f = id_elem_f * Nf_dir[d];

            id_elem_s    = ((*ke)*Nelem[1] + *je)*Nelem[0] + *ie;

            elem_edges = edges[id_elem_s];

            for (int i: icorners)
                for (int m: dirs)
                    elem_corners[i][m] = elem_edges[i].endpoints[edge_to_corner[i]][m];

            for (int n2 = 0; n2 < Ns[d2]; ++n2)
            for (int n1 = 0; n1 < Ns[d1]; ++n1)
            for (int n0 = 0; n0 < Nf[d];  ++n0)
            {
                mem_loc = idf(n0,n1,n2,d) + mem_offset_f;

                xe[d]  = xf(d, n0);
                xe[d1] = xs(d1,n1);
                xe[d2] = xs(d2,n2);

                /* Fill spatial metric --- need to interpolate to find r at flux
                 * point first, since this hasn't been stored.                   */
                polynomial_transfinite_map_3D(xe, elem_edges, elem_corners, rp);
                physics[d].metric->fill(rp, mem_loc);
            } // end loop over flux points
        } // end loop over elements
    } // end loop over direction of flux-point set


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

    write::message("Filling Chebyshev filtering matrices");
    /* For filtering data lying on the Gauss/solution points */
    for (int d: dirs)
        chebyshev_filter(d) = lagrange::chebyshev_filtering_matrix(Ns[d]);

    write::message("Finished filling spectral-difference matrices");

    return;
}
