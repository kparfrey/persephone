#include "transfinite_map.hpp"

#include "domain_map.hpp"
#include "domain_map_torus.hpp"
#include "edge.hpp"

/* 3D label of each of the 6 faces' 4 corners. Ordering is via cyclic permutation from
 * the face-normal direction; ie face 0 has normal-dir = 2, so its face coords are 
 * {0, 1} with the 0 direction moving first. */
constexpr int corner_map[6][4] = {{0, 1, 2, 3},  // face 0
                                  {4, 5, 6, 7},  //      1
                                  {0, 3, 7, 4},  //      2
                                  {1, 2, 6, 5},  //      3
                                  {0, 4, 5, 1},  //      4
                                  {3, 7, 6, 2}}; //      5


/* 3D label of each of the 6 faces' 4 edges. Ordering is as for corner_map ---
 * these two much have identical ordering systems so that the edges and corners
 * create the correct 2D edge-corner diagram. */
constexpr int edge_map[6][4] = {{ 0,  1,  2,  3},  // face 0
                                { 4,  5,  6,  7},  //      1
                                { 3, 11,  7,  8},  //      2
                                { 1, 10,  5,  9},  //      3
                                { 8,  4,  9,  0},  //      4
                                {11,  6, 10,  2}}; //      5



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


static void transfinite_map_3D(const real_t Gx[12][3], const real_t corners[8][3], 
                               const real_t Fx[6][3],  const real_t x[3], 
                                     real_t r[3])
{
    /* u, v, w etc. language aligns with Eriksson (1985) */
    const real_t u = x[0]; 
    const real_t v = x[1];
    const real_t w = x[2];

    real_t pu, pv, pw;    // Eriksson's pi_u etc.
    real_t puv, puw, pvw; // Eriksson's pi_u pi_v etc.
    real_t puvw;          // Eriksson's pi_u pi_v pi_w

    for (int i: dirs) // Iterate over coord-vector components
    {
        /* Face terms */
        pu = (1 - u) * Fx[2][i] + u * Fx[3][i];
        pv = (1 - v) * Fx[4][i] + v * Fx[5][i];
        pw = (1 - w) * Fx[0][i] + w * Fx[1][i];

        /* Edge terms */
        puv =   (1 - u) * (1 - v) * Gx[8][i]  + u  * (1 - v) * Gx[9][i]
              + (1 - u) *      v  * Gx[11][i] + u  *      v  * Gx[10][i];

        puw =   (1 - u) * (1 - w) * Gx[3][i]  + u  * (1 - w) * Gx[1][i]
              + (1 - u) *      w  * Gx[7][i]  + u  *      w  * Gx[5][i];
        
        pvw =   (1 - v) * (1 - w) * Gx[0][i]  + v  * (1 - w) * Gx[2][i]
              + (1 - v) *      w  * Gx[4][i]  + v  *      w  * Gx[6][i];

        /* Corner terms */
        puvw =   (1 - u) * (1 - v) * (1 - w) * corners[0][i]
               +      u  * (1 - v) * (1 - w) * corners[1][i]
               + (1 - u) *      v  * (1 - w) * corners[3][i]
               +      u  *      v  * (1 - w) * corners[2][i]
               + (1 - u) * (1 - v) *      w  * corners[4][i]
               +      u  * (1 - v) *      w  * corners[5][i]
               + (1 - u) *      v  *      w  * corners[7][i]
               +      u  *      v  *      w  * corners[6][i];

        /* Compose the final Boolean sum */
        r[i] = pu + pv + pw - puv - puw - pvw + puvw;
    }

    return;
}


/* x is in groupwise reference space, x = [0,1] 
 * Used to map from the Group boundary edges, given in map, to an element's edges */
void analytic_transfinite_map_2D(const real_t x[2], DomainMap* const map,
                                 const real_t corners[4][3], real_t r[3])
{
    /* These are the values at the Group edges corresponding to x[0] and x[1] */
    real_t Gx[4][3]; /* Gamma(xi) or Gamma(eta) */
    
    (*map)(0, x[0], Gx[0]);
    (*map)(1, x[1], Gx[1]);
    (*map)(2, x[0], Gx[2]);
    (*map)(3, x[1], Gx[3]);

    transfinite_map_2D(Gx, corners, x, r);

    return;
}


/* For group-boundary to element-boundary mapping. Should rename. */
void analytic_transfinite_map_3D(const real_t x[3], DomainMap* const map,
                                 const real_t corners[8][3], real_t r[3])
{
    real_t Gx[12][3]; /* Edge values: e.g. f(u,0,1) or f(1,v,1) */
    real_t Fx[6][3];  /* Face values: e.g. f(u,v,0) or f(u,1,w) */
    
    for (int i: iedges)
        (*map)(i, x[edge_dir[i]], Gx[i]);

    /* Do the 2D TFI subproblem on each of the group's faces */
    /* Possibly split into function for calling by polynomial TFI too? */
    real_t fx[2];    // 2D groupwise coord on this face
    real_t fc[4][3]; // This face's 4 corners
    real_t fe[4][3]; // Values evaluated at this face's 4 edges
    for (int f: ifaces)
    {
        fx[0] = dir_plus_one[face_normal[f]];
        fx[1] = dir_plus_two[face_normal[f]];

        /* Iterate over this face's 4 corners and edges */
        for (int j = 0; j < 4; ++j)
        for (int d: dirs)
        {
            fc[j][d] = corners[corner_map[f][j]][d]; 
            fe[j][d] = Gx[edge_map[f][j]][d]; 
        }

        transfinite_map_2D(fe, fc, fx, Fx[f]);
    }

    transfinite_map_3D(Gx, corners, Fx, x, r);

    return;
}




/* x is in elementwise reference space, x = [0,1] */
void polynomial_transfinite_map_2D(const real_t x[3],
                                   const Edge edges[4], const real_t corners[4][3], 
                                   real_t r[3])
{
    real_t Gx[4][3]; 
    
    /* Storing at Lobatto points -- need to interpolate to solution points */
    for (int i = 0; i < 4; ++i) // for each edge...
        edges[i].eval(x[edges[i].dir], Gx[i]);

    transfinite_map_2D(Gx, corners, x, r);

    return;
}


/* dir is the reference-space direction along which to take the derivative */
void drdx_transfinite_map_2D(const int dir, const real_t x[3], 
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
            edges[1].eval(x[1], Gx[1]);
            edges[3].eval(x[1], Gx[3]);

            for (int i = 0; i < 2; ++i) // For each component r^i...
                dr[i] = (1-x[1]) * (dG[0][i] + corners[0][i] - corners[1][i])
                         + x[1]  * (dG[2][i] + corners[3][i] - corners[2][i])
                                 +  Gx[1][i] - Gx[3][i];

            dr[2] = 0.0;

            break;

        case 1: // eta
            edges[1].diff(x[1], dG[1]);
            edges[3].diff(x[1], dG[3]);
            edges[0].eval(x[0], Gx[0]);
            edges[2].eval(x[0], Gx[2]);

            for (int i = 0; i < 2; ++i)
                dr[i] = (1-x[0]) * (dG[3][i] + corners[0][i] - corners[3][i])
                         + x[0]  * (dG[1][i] + corners[1][i] - corners[2][i])
                                 +  Gx[2][i] - Gx[0][i];

            dr[2] = 0.0;

            break;

        case 2:
            dr[0] = 0.0;
            dr[1] = 0.0;
            dr[2] = 1.0;
            break;
    }

    return;
}
