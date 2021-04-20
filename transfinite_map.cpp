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



/* A 2D mapping in that it uses 4 edges. Apply to all 3 coordinates, so are effectively
 * doing an interpolation on a warped plane in 3D space. Used to find the required values
 * on the hex element's 6 faces when doing 3D (volumetric) interpolation. */
static void transfinite_map_2D(const real_t G[4][3], const real_t corners[4][3], 
                               const real_t x[2],          real_t r[3])
{
    for (int i: dirs) // Iterate over coord-vector components
    {
        r[i] =   (1-x[0])*G[3][i] + x[0]*G[1][i] 
               + (1-x[1])*G[0][i] + x[1]*G[2][i]
               - (1-x[0])*((1-x[1])*corners[0][i] + x[1]*corners[3][i])
               -    x[0] *((1-x[1])*corners[1][i] + x[1]*corners[2][i]);
    }

    return;
}


static void transfinite_map_3D(const real_t G[12][3], const real_t corners[8][3], 
                               const real_t F[6][3],  const real_t x[3], 
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
        pu = (1 - u) * F[2][i] + u * F[3][i];
        pv = (1 - v) * F[4][i] + v * F[5][i];
        pw = (1 - w) * F[0][i] + w * F[1][i];

        /* Edge terms */
        puv =   (1 - u) * (1 - v) * G[8][i]  + u  * (1 - v) * G[9][i]
              + (1 - u) *      v  * G[11][i] + u  *      v  * G[10][i];

        puw =   (1 - u) * (1 - w) * G[3][i]  + u  * (1 - w) * G[1][i]
              + (1 - u) *      w  * G[7][i]  + u  *      w  * G[5][i];
        
        pvw =   (1 - v) * (1 - w) * G[0][i]  + v  * (1 - w) * G[2][i]
              + (1 - v) *      w  * G[4][i]  + v  *      w  * G[6][i];

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


/* 3D TFI needs values at arbitrary locations on the hex's 6 faces, eg f(u,v,0),
 * but we define the hex via its 12 edges. Find the required 6 face values with
 * 2D TFI using each face's 4 corners and 4 edges */
static void face_values_via_2D_TFI(const real_t G[12][3],
                                   const real_t corners[8][3],
                                   const real_t x[3],
                                         real_t F[6][3])
{
    real_t xf[2];    // 2D groupwise coord on this face
    real_t cf[4][3]; // This face's 4 corners
    real_t Gf[4][3]; // Values evaluated at this face's 4 edges
    for (int f: ifaces)
    {
        xf[0] = x[dir_plus_one[face_normal[f]]];
        xf[1] = x[dir_plus_two[face_normal[f]]];

        /* Iterate over this face's 4 corners and edges */
        for (int j = 0; j < 4; ++j)
        for (int d: dirs)
        {
            cf[j][d] = corners[corner_map[f][j]][d]; 
            Gf[j][d] = G[edge_map[f][j]][d]; 
        }

        transfinite_map_2D(Gf, cf, xf, F[f]);
    }

    return;
}


/* x is in groupwise reference space, x = [0,1] 
 * Used to map from the Group boundary edges, given in map, to an element's edges */
void analytic_transfinite_map_2D(const real_t x[2], DomainMap* const map,
                                 const real_t corners[4][3], real_t r[3])
{
    /* These are the values at the Group edges corresponding to x[0] and x[1] */
    real_t G[4][3]; /* Gamma(xi) or Gamma(eta) */
    
    (*map)(0, x[0], G[0]);
    (*map)(1, x[1], G[1]);
    (*map)(2, x[0], G[2]);
    (*map)(3, x[1], G[3]);

    transfinite_map_2D(G, corners, x, r);

    return;
}


/* For group-boundary to element-boundary mapping. Should rename. */
void analytic_transfinite_map_3D(const real_t x[3], DomainMap* const map,
                                 const real_t corners[8][3], real_t r[3])
{
    real_t G[12][3]; /* Edge values: e.g. f(u,0,1) or f(1,v,1) */
    real_t F[6][3];  /* Face values: e.g. f(u,v,0) or f(u,1,w) */
    
    for (int i: iedges)
        (*map)(i, x[edge_dir[i]], G[i]);

    face_values_via_2D_TFI(G, corners, x, F);

    transfinite_map_3D(G, corners, F, x, r);

    return;
}




/* x is in elementwise reference space, x = [0,1] */
void polynomial_transfinite_map_2D(const real_t x[2],
                                   const Edge edges[4], 
                                   const real_t corners[4][3], 
                                         real_t r[3])
{
    real_t G[4][3]; 
    
    /* Storing at Lobatto points -- need to interpolate to solution points */
    for (int i = 0; i < 4; ++i) // for each edge...
        edges[i].eval(x[edges[i].dir], G[i]);

    transfinite_map_2D(G, corners, x, r);

    return;
}


void polynomial_transfinite_map_3D(const real_t x[3],
                                   const Edge edges[12], 
                                   const real_t corners[8][3], 
                                         real_t r[3])
{
    real_t G[12][3]; /* Edge values: e.g. f(u,0,1) or f(1,v,1) */
    real_t F[6][3];  /* Face values: e.g. f(u,v,0) or f(u,1,w) */
    
    /* Storing at Lobatto points -- need to interpolate to solution points */
    for (int i: iedges)
        edges[i].eval(x[edges[i].dir], G[i]);

    face_values_via_2D_TFI(G, corners, x, F);

    transfinite_map_3D(G, corners, F, x, r);

    return;
}


/* dir is the reference-space direction along which to take the derivative */
void drdx_transfinite_map_2D(const int dir, const real_t x[3], 
                             const Edge edges[4], const real_t corners[4][3], 
                                   real_t drdx[3])
{
    real_t G[4][3];  // Keep all four rows for clarity, even though will only use 
    real_t dG[4][3]; // two rows in each of these for a given direction
    
    switch (dir)
    {
        case 0: // xi
            edges[0].diff(x[0], dG[0]);
            edges[2].diff(x[0], dG[2]);
            edges[1].eval(x[1], G[1]);
            edges[3].eval(x[1], G[3]);

            for (int i = 0; i < 2; ++i) // For each component r^i...
                drdx[i] = (1-x[1]) * (dG[0][i] + corners[0][i] - corners[1][i])
                           + x[1]  * (dG[2][i] + corners[3][i] - corners[2][i])
                                   +  G[1][i] - G[3][i];

            drdx[2] = 0.0;

            break;

        case 1: // eta
            edges[1].diff(x[1], dG[1]);
            edges[3].diff(x[1], dG[3]);
            edges[0].eval(x[0], G[0]);
            edges[2].eval(x[0], G[2]);

            for (int i = 0; i < 2; ++i)
                drdx[i] = (1-x[0]) * (dG[3][i] + corners[0][i] - corners[3][i])
                           + x[0]  * (dG[1][i] + corners[1][i] - corners[2][i])
                                   +  G[2][i] - G[0][i];

            drdx[2] = 0.0;

            break;

        case 2:
            drdx[0] = 0.0;
            drdx[1] = 0.0;
            drdx[2] = 1.0;
            break;
    }

    return;
}


/* dF[face-id][derivative-dir][phys-vector-component] */
static void face_derivs_via_2D_TFI(const real_t G[12][3],
                                   const real_t dG[12][3],
                                   const real_t corners[8][3],
                                   const real_t x[3],
                                         real_t dF[6][3][3])
{
    int d0, d1, d2;   // face-normal-dir, ... + 1, and ... + 2
    real_t xf[2];     // 2D groupwise coord on this face
    real_t cf[4][3];  // This face's 4 corners
    real_t Gf[4][3];  // Values evaluated at this face's 4 edges
    real_t dGf[4][3]; // Derivative along this face's 4 edges

    for (int f: ifaces)
    {
        d0    = face_normal[f];
        d1    = dir_plus_one[d0];
        d2    = dir_plus_two[d0];
        xf[0] = x[d1];
        xf[1] = x[d2];

        /* Iterate over this face's 4 corners and edges */
        /* Express everything in 2D language on the face, so can use
         * the known 2D TFI method directly. */
        for (int j = 0; j < 4; ++j) // corner and edge indices in 2D
        for (int i: dirs) // physical vector components
        {
            cf[j][i] = corners[corner_map[f][j]][i]; 
            Gf[j][i] = G[edge_map[f][j]][i]; 
            dGf[j][i] = dG[edge_map[f][j]][i]; 
        }

        // Should reimplement 2D drdx here since the evaluation
        // and differentiation has been taken care of already

        for (int i: dirs) // physical vector components
        {
            /* differentiate along d1 = dir-0 on the 2D face */
            dF[f][d1][i] = (1-xf[1]) * (dGf[0][i] + cf[0][i] - cf[1][i])
                            + xf[1]  * (dGf[2][i] + cf[3][i] - cf[2][i])
                            + Gf[1][i] - Gf[3][i];
            
            /* differentiate along d2 = dir-1 on the 2D face */
            dF[f][d2][i] = (1-xf[0]) * (dGf[3][i] + cf[0][i] - cf[3][i])
                            + xf[0]  * (dGf[1][i] + cf[1][i] - cf[2][i])
                            + Gf[2][i] - Gf[0][i];

            dF[f][d0][i] = 0.0; // derivative normal to the face is zero
        }

    }

    return;
}


/* Only used for polynomial-based interpolation, so include all the Edge
 * and Edge-derivative evaluation in this function. 
 * Store dr as dr[reference][physical] --- note the reverse of Jarr in Metric */
void drdx_transfinite_map_3D(const real_t x[3], const Edge edges[12],
                             const real_t corners[8][3], 
                                   real_t drdx[3][3])
{
    /* u, v, w etc. language aligns with Eriksson (1985) */
    const real_t u = x[0]; 
    const real_t v = x[1];
    const real_t w = x[2];

    real_t pu, pv, pw;    // Eriksson's pi_u etc.
    real_t puv, puw, pvw; // Eriksson's pi_u pi_v etc.
    real_t puvw;          // Eriksson's pi_u pi_v pi_w

    real_t  G[12][3]; /* Edge values: e.g. f(u,0,1) or f(1,v,1) */
    real_t  F[6][3];  /* Face values: e.g. f(u,v,0) or f(u,1,w) */
    real_t dG[12][3]; /* Derivative along each edge             */
    real_t dF[6][3][3]; /* Derivatives on each face             */

    /* Evaluate values on edges and faces. Duplicates work done
     * when calculating the coords in the first place...        */
    for (int i: iedges)
        edges[i].eval(x[edges[i].dir], G[i]);

    face_values_via_2D_TFI(G, corners, x, F);

    /* Calculate derivatives along edges and on faces */
    for (int i: iedges)
        edges[i].diff(x[edges[i].dir], dG[i]);

    face_derivs_via_2D_TFI(G, dG, corners, x, dF);

    /* Deliberately make this section very simple to avoid confusion.
     * Find the derivative along each reference direction by direct
     * differentiation of the transformation in transfinite_map_3D().
     * Don't rearrange or group terms, to make the correspondence to 
     * the original transformation obvious. */
    for (int i: dirs) // Iterate over coord-vector components
    {
        /*** Differentiate wrt reference-direction-0, aka u ***/
        /* Face terms */
        pu = - F[2][i] + F[3][i];
        pv = (1 - v) * dF[4][0][i] + v * dF[5][0][i];
        pw = (1 - w) * dF[0][0][i] + w * dF[1][0][i];

        /* Edge terms */
        puv = - (1 - v) * G[8][i]  + (1 - v) * G[9][i]
              -      v  * G[11][i] +      v  * G[10][i];

        puw = - (1 - w) * G[3][i]  + (1 - w) * G[1][i]
              -      w  * G[7][i]  +      w  * G[5][i];
        
        pvw =   (1 - v) * (1 - w) * dG[0][i]  + v  * (1 - w) * dG[2][i]
              + (1 - v) *      w  * dG[4][i]  + v  *      w  * dG[6][i];

        /* Corner terms */
        puvw = - (1 - v) * (1 - w) * corners[0][i]
               + (1 - v) * (1 - w) * corners[1][i]
               -      v  * (1 - w) * corners[3][i]
               +      v  * (1 - w) * corners[2][i]
               - (1 - v) *      w  * corners[4][i]
               + (1 - v) *      w  * corners[5][i]
               -      v  *      w  * corners[7][i]
               +      v  *      w  * corners[6][i];

        drdx[0][i] = pu + pv + pw - puv - puw - pvw + puvw;
        

        /*** Differentiate wrt reference-direction-1, aka v ***/
        /* Face terms */
        pu = (1 - u) * dF[2][1][i] + u * dF[3][1][i];
        pv = - F[4][i] + F[5][i];
        pw = (1 - w) * dF[0][1][i] + w * dF[1][1][i];

        /* Edge terms */
        puv = - (1 - u) * G[8][i]  - u * G[9][i]
              + (1 - u) * G[11][i] + u * G[10][i];

        puw =   (1 - u) * (1 - w) * dG[3][i]  + u  * (1 - w) * dG[1][i]
              + (1 - u) *      w  * dG[7][i]  + u  *      w  * dG[5][i];
        
        pvw = - (1 - w) * G[0][i]  + (1 - w) * G[2][i]
              -      w  * G[4][i]  +      w  * G[6][i];

        /* Corner terms */
        puvw = - (1 - u) * (1 - w) * corners[0][i]
               -      u  * (1 - w) * corners[1][i]
               + (1 - u) * (1 - w) * corners[3][i]
               +      u  * (1 - w) * corners[2][i]
               - (1 - u) *      w  * corners[4][i]
               -      u  *      w  * corners[5][i]
               + (1 - u) *      w  * corners[7][i]
               +      u  *      w  * corners[6][i];

        drdx[1][i] = pu + pv + pw - puv - puw - pvw + puvw;
        

        /*** Differentiate wrt reference-direction-2, aka w ***/
        /* Face terms */
        pu = (1 - u) * dF[2][2][i] + u * dF[3][2][i];
        pv = (1 - v) * dF[4][2][i] + v * dF[5][2][i];
        pw = - F[0][i] + F[1][i];

        /* Edge terms */
        puv =   (1 - u) * (1 - v) * dG[8][i]  + u  * (1 - v) * dG[9][i]
              + (1 - u) *      v  * dG[11][i] + u  *      v  * dG[10][i];

        puw = - (1 - u) * G[3][i] - u * G[1][i]
              + (1 - u) * G[7][i] + u * G[5][i];
        
        pvw = - (1 - v) * G[0][i] - v * G[2][i]
              + (1 - v) * G[4][i] + v * G[6][i];

        /* Corner terms */
        puvw = - (1 - u) * (1 - v) * corners[0][i]
               -      u  * (1 - v) * corners[1][i]
               - (1 - u) *      v  * corners[3][i]
               -      u  *      v  * corners[2][i]
               + (1 - u) * (1 - v) * corners[4][i]
               +      u  * (1 - v) * corners[5][i]
               + (1 - u) *      v  * corners[7][i]
               +      u  *      v  * corners[6][i];

        drdx[2][i] = pu + pv + pw - puv - puw - pvw + puvw;
    }

    return;
}
