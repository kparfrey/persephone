#include "transfinite_map.hpp"

#include "domain_map.hpp"
#include "edge.hpp"


static void transfinite_map_2D(const real_t Gx[4][3], const real_t corners[4][3], 
                               const real_t x[2], real_t r[3])
{
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


/* x is in groupwise reference space, x = [0,1] */
void analytic_transfinite_map_2D(const real_t x[2], DomainMap* const map,
                                 const real_t corners[4][3], real_t r[3])
{
    real_t Gx[4][3]; /* Gamma(xi) or Gamma(eta) */
    
    (*map)(0, x[0], Gx[0]);
    (*map)(1, x[1], Gx[1]);
    (*map)(2, x[0], Gx[2]);
    (*map)(3, x[1], Gx[3]);

    transfinite_map_2D(Gx, corners, x, r);

    return;
}


/* x is in elementwise reference space, x = [0,1] */
void polynomial_transfinite_map_2D(const real_t x[2], const int point_idx[2],
                                   const Edge edges[4], const real_t corners[4][3], 
                                   real_t r[3])
{
    real_t Gx[4][3]; 
    
    /* Don't need to interpolate, since only need Gamma(x) at the 
     * existing locations along the edge (if have stored the Gauss points) */
    for (int i = 0; i < 4; ++i) // for each edge...
        for (int j = 0; j < 2; ++j) // for each coord component...
            Gx[i][j] = edges[i].r(j, point_idx[edges[i].dir]);

    transfinite_map_2D(Gx, corners, x, r);

    return;
}


/* dir is the reference-space direction along which to take the derivative */
void drdx_transfinite_map_2D(const int dir, const real_t x[2], const int point_idx[2],
                             const Edge edges[4], const real_t corners[4][3], 
                                   real_t dr[3])
{
    real_t Gx[4][3]; // Keep all four rows for clarity, even though will only use 
    real_t dG[4][3]; // two rows in each of these for a given direction
    
    switch (dir)
    {
        case 0: // xi
            edges[0].diff(x[0], dG[0]);
            edges[2].diff(x[0], dG[2]);

            for (int i = 0; i < 2; ++i) // For each component r^i...
            {
                Gx[1][i] = edges[1].r(i, point_idx[1]);
                Gx[3][i] = edges[3].r(i, point_idx[1]);
                
                dr[i] = (1-x[1]) * (dG[0][i] + corners[0][i] - corners[1][i])
                         + x[1]  * (dG[2][i] + corners[3][i] - corners[2][i])
                                 +  Gx[1][i] - Gx[3][i];
            }

            break;

        case 1: // eta
            edges[1].diff(x[1], dG[1]);
            edges[3].diff(x[1], dG[3]);

            for (int i = 0; i < 2; ++i)
            {
                Gx[0][i] = edges[0].r(i, point_idx[0]);
                Gx[2][i] = edges[2].r(i, point_idx[0]);
                
                dr[i] = (1-x[0]) * (dG[3][i] + corners[0][i] - corners[3][i])
                         + x[0]  * (dG[1][i] + corners[1][i] - corners[2][i])
                                 +  Gx[2][i] - Gx[0][i];
            }

            break;
    }

    dr[2] = 0.0;

    return;
}
