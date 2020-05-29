#include "metric.hpp"
#include "kernels.hpp"
#include "matrix.hpp"

void Metric::allocate(int N)
{
    /* For now, assume that N is the no. of solution points per block,
     * and that no metric quantities on flux points are required. */

    for (int i: dirs)
        for (int j: dirs)
        {
               J[i][j]  = kernels::alloc(N);
            Jinv[i][j]  = kernels::alloc(N);
               g[i][j]  = kernels::alloc(N);
            ginv[i][j]  = kernels::alloc(N);
            gphys[i][j] = kernels::alloc(N);
        }

    rdetg = kernels::alloc(N);

    return;
}


/* Simplified version that assumes basic Cartesian shape */
void Metric::setup(int Nelem[3], int Ns_block, real_t corners[8][3])
{
    allocate(Ns_block);


    real_t dr_elemblock[3];
    real_t dr_elem[3];

    dr_elemblock[0] = corners[1][0] - corners[0][0];
    dr_elemblock[1] = corners[3][1] - corners[0][1];
    dr_elemblock[2] = corners[4][2] - corners[0][2];

    for (int i: dirs)
        dr_elem[i] = dr_elemblock[i] / Nelem[i];

    /* Can leave all off-diagonal components as zero */
    for (int d: dirs)
        for (int i = 0; i < Ns_block; ++i) // All elems identical
        {
            J[d][d][i]    = dr_elem[d]/2.0; // Since length of ref element is 2
            Jinv[d][d][i] = 1.0 / J[d][d][i];
            gphys[d][d][i] = 1.0;
        }


    /* Find reference-space metric by performing a general coordinate
     * transformation on the physical metric */
    transform_twoTensor(gphys, g, Ns_block, phys2ref, covariant);


    /* Find the inverse metric and sqrt(g) in reference space automatically */
    Matrix gmat;       // Matrix object to be reused
    real_t garr[3][3]; // Array to be reused
    for (int n = 0; n < Ns_block; ++n) // For every point in the block...
    {
        for (int i: dirs)
        for (int j: dirs)
            garr[i][j] = g[i][j][n]; // Put this point's metric into the array

        gmat.fill(garr);
        gmat.find_inverse();
        
        for (int i: dirs)
        for (int j: dirs)
            ginv[i][j][n] = gmat.inv[i][j];

        rdetg[n] = sqrt( gmat.det );
    }

    return;
}


/* Do a general coordinate transformation for a two-tensor between reference
 * and physical coordinates */
void Metric::transform_twoTensor(real_t* T_in[3][3], real_t* T_out[3][3], int N,
                                 CoordTransDir ctd, Components c)
{
    /* Pointer to a 3x3 array of pointers to real_t... */
    real_t* (*V)[3][3];

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
    
    for (int a: dirs)
    for (int b: dirs)
        for (int i: dirs)
        for (int j: dirs)
            for (int n = 0; n < N; ++n) 
                T_out[a][b][n] = (*V)[a][i][n] * (*V)[b][j][n] * T_in[i][j][n];

    return;
}
