#ifndef ELEMENT_BLOCK_HPP
#define ELEMENT_BLOCK_HPP

#include "common.hpp"
#include "edge.hpp"

class ElementBlock
{
    private:
    /* These only called by the setup() method */
    void allocate();
    void set_computational_coords();
    void set_physical_coords();


    public:

    /* Data */
    int Nelem[3]; // No. of elements in each direction
    int Ns[3];    // No. of soln points in each element, in each direction
    int Nf[3];    // No. of flux points "   "      "     "   "       "
    int Nsf[3];   // For indexing the flux-point arrays; see ElementBlock::idf()
    int Ns_tot[3]; // Total number of soln points in each direction
    int Ns_elem;   // Total no. of soln points in each element
    int Nf_elem;   // Total no. of flux points in each element
    int Nf_dir[3]; // No. of flux points in each direction, in each element
    int Nelem_block;
    int Nfield;

    /* Geometrical information */
    /* For now these are degenerate with the corresponding Process data */
    real_t corners[8][3]; // 3 physical-space coordinates for each of 8 corners
    Edge edges[12]; 

    /* Computational space locations: one element's worth only */
    real_t *xs[3]; // Soln point locations in computational space, in each direction 
    real_t *xf[3]; // Flux point locations... 

    /* Physical/real space: need the whole 3D grid for each r^i, i in [0,1,2] */
    real_t *rs[3]; // Soln points, stored as a 1D array for each vector component
    real_t *rf[3][3]; // Flux points: 3 3D arrays, one for each transform direction
                      // i.e. rf[transform-dir][coord][3D-array]

    real_t *fields; // Field data on solution points

    /* Methods */
    inline int ids(int i, int j, int k); // soln-point index map
    inline int idf(int idir, int j1, int j2, int dir); // flux-point index map
    inline int id_elem(int i, int j, int k);
    inline int ids_full(int i, int j, int k, int ie, int je, int ke); 
    inline int cyclic_add(int i, int add);
    void setup();
    void free();
};


/* C-order, i slowest index, k fastest */
/* The i,j,k index solution points within an element */
inline int ElementBlock::ids(int i, int j, int k)
{
    //return Ns[0]*Ns[1]*Ns[2]*var + (i*Ns[1] + j)*Ns[2] + k;
    return (i*Ns[1] + j)*Ns[2] + k;
}


/* 3D --> 1D indices for flux points within an element */
/* For dir = 0: idir -> i/0, j1 -> j/1, j2 -> k/2 etc.
 * ie j2 is "two above"; for dir = 2, j2 -> j/1 */
/* This assumes CPU ordering (transform direction is fastest) ---
 * will need to use conditional compilation for reversed GPU version */
inline int ElementBlock::idf(int idir, int j1, int j2, int dir)
{
    return (j2*Nsf[dir] + j1)*Nf[dir] + idir;
}


/* Turn 3D indices of elements in the block into a 1D index */
inline int ElementBlock::id_elem(int i, int j, int k)
{
    return (i*Nelem[1] + j)*Nelem[2] + k;
}


/*
inline int ElementBlock::mem_id_elem(int i, int j, int k)
{
    return ElementBlock::id_elem(i,j,k) * Ns_elem;
}
*/

inline int ElementBlock::ids_full(int i, int j, int k, int ie, int je, int ke)
{
    return ElementBlock::ids(i,j,k) + Ns_elem * ElementBlock::id_elem(ie,je,ke);
}


inline int ElementBlock::cyclic_add(int i, int add)
{
    int ia = i + add;

    if (ia > 2) 
        ia -= 3;

    return ia;
}
#endif
