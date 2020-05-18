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
    int Ns_elem;  // Total no. of soln points in each element
    int Nf_elem;  // Total no. of flux points "   "     "
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
    real_t *rf[3]; // Flux points...

    real_t *fields; // Field data on solution points

    /* Methods */
    inline int ids(int i, int j, int k); // soln-point index map
    inline int idf(int i, int j, int k); // flux-point index map
    inline int id_elem(int i, int j, int k);
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
inline int ElementBlock::idf(int i, int j, int k)
{
    return (i*Nf[1] + j)*Nf[2] + k;
}


/* Turn 3D indices of elements in the block into a 1D index */
inline int ElementBlock::id_elem(int i, int j, int k)
{
    return (i*Nelem[1] + j)*Nelem[2] + k;
}

#endif
