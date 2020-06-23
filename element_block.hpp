#ifndef ELEMENT_BLOCK_HPP
#define ELEMENT_BLOCK_HPP

#include "common.hpp"
#include "edge.hpp"
#include "metric.hpp"
#include "tensor_field.hpp"

class ElementBlock
{
    private:
    /* These only called by the setup() method */
    void allocate_on_host();
    void set_computational_coords();
    void set_physical_coords();
    void fill_spectral_difference_matrices();

    public:
    /* Data */
    int Nelem[3];   // No. of elements in each direction
    int Ns[3];      // No. of soln points in each element, in each direction
    int Nf[3];      // No. of flux points "   "      "     "   "       "
    int Ns_tot[3];  // Total number of soln points in each direction
    int Ns_elem;    // Total no. of soln points in each element
    int Nf_elem;    // Total no. of flux points in each element
    int Nf_dir[3];  // No. of i-direction flux points, in each element
    int Nelem_block;// No. of elements in the block
    int Ns_block;   // No. of soln points in the block
    int Nf_dir_block[3];
    int Nfield;

    LengthBucket lengths; // To encapsulate all lengths needed for data indexing
    

    /* Geometrical information */
    /* For now these are degenerate with the corresponding Process data */
    real_t corners[8][3]; // 3 physical-space coordinates for each of 8 corners
    Edge edges[12]; 

    /* Computational space locations: one element's worth only */
    VectorField xs; // Soln point locations in computational space, in each direction 
    VectorField xf; // Flux point locations... 

    /* Physical/real space: need the whole 3D grid for each r^i, i in [0,1,2] */
    VectorField rs;    // Soln points, stored as a 1D block for each vector component
    VectorField rf[3]; // Flux points: need separate VectorField for each transform dir
                       // i.e. rf[transform-dir](coord component, 3D space in 1D block)

    /* Will need on device. Leave as is since it's just one pointer, or replace
     * with ScalarField (or its alias) to make clear it's headed for device? */
    real_t* fields; // Field data on solution points

    /* Matrices interpolating between soln and flux points, and taking derivatives.
     * Need one matrix for each transform direction, stored in flattened blocks.
     * i.e. index with matrix(transform-direction, flattened data). */
    VectorField soln2flux;      // Interpolate from soln points to flux points
    VectorField fluxDeriv2soln; // Take derivative of values on flux points, interpolate
                                // back to soln points

    Metric metric;


    /* Methods */
    inline int ids(int i, int j, int k); // soln-point index map
    inline int idf(int n0, int n1, int n2, int dir); // flux-point index map
    inline int id_elem(int i, int j, int k);
    inline int ids_full(int i, int j, int k, int ie, int je, int ke); 
    void setup();
    void free();
    void move_to_device();
};


/* C-order, i slowest index, k fastest */
/* The i,j,k index solution points within an element */
inline int ElementBlock::ids(int i, int j, int k)
{
    //return Ns[0]*Ns[1]*Ns[2]*var + (i*Ns[1] + j)*Ns[2] + k;
    //return (i*Ns[1] + j)*Ns[2] + k;
    return (k*Ns[1] + j)*Ns[0] + i; // Changed order: June 23
}


/* 3D --> 1D indices for flux points within an element */
/* For dir = 0: n0 -> i/0, n1 -> j/1, n2 -> k/2 etc.
 * ie n2 is "two above"; for dir = 2, n2 -> j/1 */
/* This assumes CPU ordering (transform direction is fastest) ---
 * will need to use conditional compilation for reversed GPU version */
inline int ElementBlock::idf(int n0, int n1, int n2, int dir)
{
    return (n2*Ns[dir_plus_one[dir]] + n1)*Nf[dir] + n0;
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

#endif
