#include <cmath>
#include "geometry.hpp"
#include "kernels.hpp"
#include "matrix.hpp"
#include "write_screen.hpp"
#include "element_block.hpp"
#include "edge.hpp"
#include "transfinite_map.hpp"
#include "geometry_labels.hpp"
#include "physics.hpp"


void Geometry::allocate_on_host(const int Ns, const int Nf[3])
{
    /* Ns and Nf[3] are the total quantities in the block */

    /* Solution points */
    Jrdetg = new real_t [Ns]();
    Qinteg = new real_t [Ns]();
    
    for(int dref: dirs)
        for(int dphys: dirs)
            dxdr[dref](dphys) = new real_t [Ns]();


    /* Flux points */
    for (int d: dirs)
    {
        for (int j: dirs)
        {
            S[d](j)  = new real_t [Nf[d]]();
            timestep_transform[d](j) = new real_t [Nf[d]]();
        }
    }
    
    return;
}


void Geometry::move_to_device()
{
    write::message("Copying Geometry's data to device");

    return;
}


void Geometry::setup_full(ElementBlock& eb)
{
    using write::str;

    allocate_on_host(eb.Ns_block, eb.Nf_dir_block);

    /* Elementwise constructs */
    Edge* elem_edges;
    real_t elem_corners[8][3]; // For 3D TF map
    int elem_idx_1D;
    int mem_offset;

    /* Pointwise constructs */
    real_t xe[3];
    real_t drdx[3][3]; 
    //real_t drdx[3]; // For 2D method, which is ref-dir by ref-dir
    
    Matrix J; // dPHYS/dREF -- J.arr[phys][ref]
    real_t detJ;
    //real_t rdetg = 1.0; //Assume physical coords are flat-spacetime Cartesian
    int mem_loc;

    /* For finding the timestep_transform */
    real_t dx, dx_l, dx_r;
    real_t tt;
    eb.timestep_transform_max = 0.0;

    /* Face normals only used in full geometry, allocate here */
    for (int dir: dirs)
    {
        int d1 = dir_plus_one[dir];
        int d2 = dir_plus_two[dir];
        int Nfaces = eb.Nelem[d1] * eb.Nelem[d2] * (eb.Nelem[dir]+1);
        int points_per_face = eb.Ns[d1]*eb.Ns[d2];

        for (int dphys: dirs)
            normal[dir](dphys) = new real_t [Nfaces * points_per_face]();
    }


    /* Solution points */
    //real_t Q_prefac = 8.0 * pi*pi*pi /(eb.Ns[0]*eb.Ns[1]*eb.Ns[2]); // for Qinteg

    for (int ke = 0; ke < eb.Nelem[2]; ++ke)
    for (int je = 0; je < eb.Nelem[1]; ++je)
    for (int ie = 0; ie < eb.Nelem[0]; ++ie)
    {
        elem_idx_1D = eb.id_elem(ie, je, ke);
        mem_offset = elem_idx_1D * eb.Ns_elem;

        elem_edges = eb.edges[elem_idx_1D];

        for (int i: icorners)
            for (int m: dirs)
                elem_corners[i][m] = elem_edges[i].endpoints[edge_to_corner[i]][m];

        for (int k = 0; k < eb.Ns[2]; ++k)
        for (int j = 0; j < eb.Ns[1]; ++j)
        for (int i = 0; i < eb.Ns[0]; ++i)
        {
            real_t Jarr[3][3] = {}; 
            mem_loc = eb.ids(i,j,k) + mem_offset;

            xe[0] = eb.xs(0,i);
            xe[1] = eb.xs(1,j);
            xe[2] = eb.xs(2,k);

            /***  2D Method --- leave here for now
            for (int dref: dirs)
            {
                drdx_transfinite_map_2D(dref, xe, elem_edges, elem_corners, drdx);
                for (int dphys: dirs)
                    Jarr[dphys][dref] = drdx[dphys];
            }
             ***/

            drdx_transfinite_map_3D(xe, elem_edges, elem_corners, drdx);
            
            for (int dphys: dirs)
            for (int dref: dirs)
                Jarr[dphys][dref] = drdx[dref][dphys];
            
            J.fill(Jarr);
            J.find_determinant();
            detJ = std::abs(J.det);
            
            Jrdetg(mem_loc) = detJ * eb.physics_soln->metric->rdetg[mem_loc];

            /* sqrt(1-x^2) for x in [-1,1] --> 2 sqrt(x-x^2) for x in [0,1] */
            //Qinteg(mem_loc) = Q_prefac * Jrdetg(mem_loc) * 
            //                    std::sqrt((xe[0]-xe[0]*xe[0]) * (xe[1]-xe[1]*xe[1]) * (xe[2]-xe[2]*xe[2])); 

            /* Using Legendre weights */
            real_t wl[5] = {0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689};
            Qinteg(mem_loc) = Jrdetg(mem_loc) * wl[i] * wl[j]; // Ignore k direction

            /* Used to find grad(U) for diffusive terms */
            J.find_inverse();
            
            /* J.inv[ref][phys] */
            for (int dref: dirs)
                for (int dphys: dirs)
                    dxdr[dref](dphys,mem_loc) = J.inv[dref][dphys];
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
         * somewhere else. */
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

        for (ne2 = 0; ne2 < eb.Nelem[d2]; ++ne2)
        for (ne1 = 0; ne1 < eb.Nelem[d1]; ++ne1)
        for (ne0 = 0; ne0 < eb.Nelem[d];  ++ne0)
        {
            id_elem_f    = (ne2*eb.Nelem[d1] + ne1)*eb.Nelem[d] + ne0;
            mem_offset_f = id_elem_f * eb.Nf_dir[d];

            id_elem_s    = ((*ke)*eb.Nelem[1] + *je)*eb.Nelem[0] + *ie;

            elem_edges = eb.edges[id_elem_s];

            for (int i: icorners)
                for (int m: dirs)
                    elem_corners[i][m] = elem_edges[i].endpoints[edge_to_corner[i]][m];

            /* S matrix for this element */
            for (int n2 = 0; n2 < eb.Ns[d2]; ++n2)
            for (int n1 = 0; n1 < eb.Ns[d1]; ++n1)
            for (int n0 = 0; n0 < eb.Nf[d];  ++n0)
            {
                real_t Jarr[3][3] = {}; 
                mem_loc = eb.idf(n0,n1,n2,d) + mem_offset_f;

                xe[d]  = eb.xf(d, n0);
                xe[d1] = eb.xs(d1,n1);
                xe[d2] = eb.xs(d2,n2);

                /*** 2D Method --- leave here for now
                for (int dref: dirs)
                {
                    drdx_transfinite_map_2D(dref, xe, elem_edges, elem_corners, drdx);
                    for (int dphys: dirs)
                        Jarr[dphys][dref] = drdx[dphys];
                }
                 ***/
                
                drdx_transfinite_map_3D(xe, elem_edges, elem_corners, drdx);
                
                for (int dphys: dirs)
                for (int dref: dirs)
                    Jarr[dphys][dref] = drdx[dref][dphys];

                J.fill(Jarr);
                J.find_determinant();
                J.find_inverse();
                detJ = std::abs(J.det);

                /* Assume J.inv[ref][phys]... */
                for (int dphys: dirs)
                    S[d](dphys,mem_loc) = detJ * eb.physics[d]->metric->rdetg[mem_loc] 
                                               * J.inv[d][dphys];


                /* Calculate helper array for finding timestep limit */
                /* Double check that xf[0] = 0, xf[-1] = 1 */
                if (n0 == 0)
                    dx = eb.xf(d, 1);
                else if (n0 == eb.Nf[d]-1)
                    dx = 1.0 - eb.xf(d,eb.Nf[d]-2);
                else
                {
                    dx_l = eb.xf(d,n0)   - eb.xf(d,n0-1);
                    dx_r = eb.xf(d,n0+1) - eb.xf(d,n0);
                    //dx = std::min(dx_l, dx_r);  // Cautious
                    dx = 0.5 * (dx_l + dx_r);  // Slightly less cautious...
                }

                /* tt[ref](phys,mem) = (1/delta x^ref) * dx^ref/dr^phys */
                for (int dphys: dirs)
                {
                    //timestep_transform[d](dphys,mem_loc) = (1.0/dx) * J.inv[d][dphys];
                    tt = (1.0/dx) * J.inv[d][dphys];
                    timestep_transform[d](dphys,mem_loc) = tt;

                    if (tt > eb.timestep_transform_max)
                        eb.timestep_transform_max = tt;
                }
            } // end loop over every flux point



            /* This needs to be tidied and rationalised... */
            /* Face normals in this transform direction for this element  */
            int points_per_face = eb.Ns[d1]*eb.Ns[d2];
            int mem_loc_normals;
            real_t s[3];
            real_t smag;
            int id_elem_normals = (ne2*eb.Nelem[d1] + ne1)*(eb.Nelem[d]+1) + ne0;
            SpatialMetric* metric = eb.physics[d]->metric;

            for (int n2 = 0; n2 < eb.Ns[d2]; ++n2)
            for (int n1 = 0; n1 < eb.Ns[d1]; ++n1)
            {
                /* Have n0 = 0 face for all but the last element */
                mem_loc         = mem_offset_f + eb.idf(0,n1,n2,d);
                mem_loc_normals = id_elem_normals * points_per_face + n2*eb.Ns[d1] + n1;
                
                for (int dphys: dirs)
                    s[dphys] = S[d](dphys,mem_loc);
               
                /* Was assuming physical coord system is Cartesian... */
                smag = std::sqrt(metric->square_cov(s, mem_loc)); // s is phys-space dual vector
                //smag = std::sqrt(metric->square(s, mem_loc));
                //smag = std::sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);

                for (int dphys: dirs)
                    normal[d](dphys,mem_loc_normals) = s[dphys]/smag;

                /* Last element in this line, also store the rightmost face */
                if (ne0 == eb.Nelem[d]-1)
                {
                    mem_loc         = mem_offset_f + eb.idf(eb.Nf[d]-1,n1,n2,d);
                    mem_loc_normals = (id_elem_normals+1) * points_per_face + n2*eb.Ns[d1] + n1;
                    
                    for (int dphys: dirs)
                        s[dphys] = S[d](dphys,mem_loc);
                    
                    smag = std::sqrt(metric->square_cov(s, mem_loc));
                    //smag = std::sqrt(metric->square(s, mem_loc));
                    //smag = std::sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);

                    for (int dphys: dirs)
                        normal[d](dphys,mem_loc_normals) = s[dphys]/smag;
                }
            }
        }
    }


    return;
}



/* Do a general coordinate transformation for a two-tensor between reference
 * and physical coordinates */
/* This function may no longer be needed */
/***
void Geometry::transform_twoTensor(const TensorField& T_in, TensorField& T_out, const int N,
                                 const CoordTransDir ctd, const Components c)
{
    TensorField* V;

    switch(c)
    {
        case covariant:
            if (ctd == phys2ref)
                V = &J;
            else
                V = &Jinv;
            break;
        case contravariant:
            if (ctd == phys2ref)
                V = &Jinv;
            else
                V = &J;
            break;
    }
    
    int lsum;
    for (int n = 0; n < N; ++n) 
        for (int a: dirs)
        for (int b: dirs)
        {
            lsum = 0.0;
            for (int i: dirs)
            for (int j: dirs)
                    lsum += (*V)(a,i,n) * (*V)(b,j,n) * T_in(i,j,n);

            T_out(a,b,n) = lsum;
        }

    return;
}
***/
