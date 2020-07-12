#ifndef TENSOR_FIELD_HPP
#define TENSOR_FIELD_HPP

#ifdef __CUDACC__
#define ACCEL_DECORATOR __host__ __device__
#else
#define ACCEL_DECORATOR
#endif


#include "common.hpp"

/* Simple structures for wrapping data arrays for passing to 
 * CUDA kernels. Deliberately try to keep this as simple as possible. 
 * Note: will also use these for objects that aren't tensor 
 * fields, like the Jacobian matrix at each point. */


/* Assume this is for a spatial 2-tensor. Could replace with a template
 * to allow the dimensionality to be set, but it's probably easier to 
 * just define a SpacetimeTensorField or whatever if that's useful. */
struct TensorField
{
    real_t* data[3][3]; 

    /* Returns a reference: allows both getting and setting the value */
    ACCEL_DECORATOR
    inline real_t& operator()(const int i, const int j, const int n)
    {
        return data[i][j][n];
    };
    

    /* For indexing a const TensorField */
    ACCEL_DECORATOR
    inline real_t  operator()(const int i, const int j, const int n) const
    {
        return data[i][j][n];
    };


    ACCEL_DECORATOR
    inline real_t*& operator()(const int i, const int j)
    {
        return data[i][j];
    };
};

using MatrixField = TensorField; // Use when the dimensions are
                                 // logically disconnected


struct VectorField
{
    real_t* data[3];

    ACCEL_DECORATOR
    inline real_t& operator()(const int i, const int n)
    {
        return data[i][n];
    };


    ACCEL_DECORATOR
    inline real_t operator()(const int i, const int n) const
    {
        return data[i][n];
    };


    ACCEL_DECORATOR
    inline real_t*& operator()(const int i)
    {
        return data[i];
    };
};


/* By analogy with the above. Not necessary, might remove. */
struct ScalarField
{
    real_t* data;
    
    /* Use () instead of [] for consistency with the 2D case */
    ACCEL_DECORATOR
    inline real_t& operator()(const int n)
    {
        return data[n];
    };


    ACCEL_DECORATOR
    inline real_t  operator()(const int n) const
    {
        return data[n];
    };


    ACCEL_DECORATOR
    inline real_t*& operator()()
    {
        return data;
    };

    
    ACCEL_DECORATOR 
    inline void operator=(real_t* external_data)
    {
        data = external_data;
        return;
    };
};


/* Include what is required to index flattened data arrays. For
 * passing to kernels. */
struct LengthBucket
{
    int Nfield;

    int Nelem[3];

    int Ns[3];
    int Ns_elem;
    int Ns_block;
    
    int Nf[3];
    int Nf_dir[3];       //Total # flux points, in each dir, in each element
    int Nf_dir_block[3]; //Total # flux points, in each dir, in the block
};


/* A vector field in a different sense, in the space of solution
 * or flux vectors. The Nfield vector components are flattened into
 * a single 1D data block. */
/*
struct MultiField
{
    real_t* data;

    int Nfield;
    int Nelem[3];
    int Ns[3];
    int Nf[3];
};
*/
#endif
