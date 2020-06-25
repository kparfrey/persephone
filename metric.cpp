#include <cmath>
#include "metric.hpp"
#include "kernels.hpp"
#include "matrix.hpp"
#include "write_screen.hpp"


void Metric::allocate_on_host(const int Ns, const int Nf[3])
{
    /* Ns and Nf[3] are the total quantities in the block */

    /* Solution points */
    /***
    for (int i: dirs)
        for (int j: dirs)
        {
            //   J(i,j)  = new real_t [Ns]();
            //Jinv(i,j)  = new real_t [Ns]();
            //   g(i,j)  = new real_t [Ns]();
            //ginv(i,j)  = new real_t [Ns]();
            //gphys(i,j) = new real_t [Ns]();
        }
     ***/

    Jrdetg = new real_t [Ns]();


    /* Flux points */
    for (int d: dirs)
    {
        /***
        for (int i: dirs)
            for (int j: dirs)
            {
                //    J_f[d](i,j) = new real_t [Nf[d]]();
                // Jinv_f[d](i,j) = new real_t [Nf[d]]();
                //    g_f[d](i,j) = new real_t [Nf[d]]();
                // ginv_f[d](i,j) = new real_t [Nf[d]]();
                //gphys_f[d](i,j) = new real_t [Nf[d]]();
            }
         ***/

        for (int j: dirs)
            S[d](j) = new real_t [Nf[d]]();

        //rdetg_f[d] = new real_t [Nf[d]]();
    }

    return;
}


void Metric::move_to_device()
{
    write::message("Copying Metric's data to device");

    return;
}


/* Simplified version that assumes basic Cartesian shape */
void Metric::setup(const int Nelem[3], const int Ns_block,
                   const int Nf_dir_block[3], const real_t corners[8][3])
{
    allocate_on_host(Ns_block, Nf_dir_block);


    real_t dr_elemblock[3];
    real_t dr_elem[3];

    dr_elemblock[0] = corners[1][0] - corners[0][0];
    dr_elemblock[1] = corners[3][1] - corners[0][1];
    dr_elemblock[2] = corners[4][2] - corners[0][2];

    for (int i: dirs)
        dr_elem[i] = dr_elemblock[i] / Nelem[i];

    //real_t identity[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
    Matrix J; // dPHYS/dREF -- J.arr[phys][ref]
    real_t detJ;
    real_t rdetg = 1.0; //Assume for now -- flat spacetime Cartesian

    /* Solution points */
    for (int n = 0; n < Ns_block; ++n) 
    {
        real_t Jarr[3][3] = {}; 
        for (int i: dirs)
            Jarr[i][i] = dr_elem[i]; // All elems identical

        J.fill(Jarr);
        J.find_determinant();
        detJ = std::abs(J.det);

        Jrdetg(n) = detJ * rdetg;
    }
    

    /* Flux points */
    for (int d: dirs)
        for (int n = 0; n < Nf_dir_block[d]; ++n)
        {
            real_t Jarr[3][3] = {}; 
            for (int i: dirs)
                Jarr[i][i] = dr_elem[i]; // All elems identical

            J.fill(Jarr);
            J.find_determinant();
            J.find_inverse();
            detJ = std::abs(J.det);

            /* Assume J.inv[ref][phys]... */
            for (int j: dirs)
                S[d](j,n) = detJ * rdetg * J.inv[d][j];
        }



    /* Find reference-space metric by performing a general coordinate
     * transformation on the physical metric */
    //transform_twoTensor(gphys, g, Ns_block, phys2ref, covariant);


    /* Find the inverse metric and sqrt(g) in reference space automatically */
    /*
    Matrix gmat;       // Matrix object to be reused
    real_t garr[3][3]; // Array to be reused
    for (int n = 0; n < Ns_block; ++n) // For every point in the block...
    {
        for (int i: dirs)
        for (int j: dirs)
            garr[i][j] = g(i,j,n); // Put this point's metric into the array

        gmat.fill(garr);
        gmat.find_inverse();
        
        for (int i: dirs)
        for (int j: dirs)
            ginv(i,j,n) = gmat.inv[i][j];

        rdetg(n) = std::sqrt(gmat.det);
    }
    */

    return;
}


/* Do a general coordinate transformation for a two-tensor between reference
 * and physical coordinates */
/* This function may no longer be needed */
/***
void Metric::transform_twoTensor(const TensorField& T_in, TensorField& T_out, const int N,
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
