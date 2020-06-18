#ifndef TENSOR_FIELD_HPP
#define TENSOR_FIELD_HPP

#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
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
    CUDA_HOSTDEV
    inline real_t& operator()(const int i, const int j, const int n)
    {
        return data[i][j][n];
    };
    

    /* For indexing a const TensorField */
    CUDA_HOSTDEV
    inline real_t  operator()(const int i, const int j, const int n) const
    {
        return data[i][j][n];
    };


    CUDA_HOSTDEV
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

    CUDA_HOSTDEV
    inline real_t& operator()(const int i, const int n)
    {
        return data[i][n];
    };


    CUDA_HOSTDEV
    inline real_t operator()(const int i, const int n) const
    {
        return data[i][n];
    };


    CUDA_HOSTDEV
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
    CUDA_HOSTDEV
    inline real_t& operator()(const int n)
    {
        return data[n];
    };


    CUDA_HOSTDEV
    inline real_t  operator()(const int n) const
    {
        return data[n];
    };


    CUDA_HOSTDEV
    inline real_t*& operator()()
    {
        return data;
    };

    
    CUDA_HOSTDEV 
    inline void operator=(real_t* external_data)
    {
        data = external_data;
        return;
    };
};
#endif
