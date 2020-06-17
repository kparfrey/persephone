#ifndef TENSOR_FIELD_HPP
#define TENSOR_FIELD_HPP

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
    real_t* d[3][3]; // Array of pointers to real data

    /* For host side only for now */
    real_t& operator()(const int i, const int j, const int n)
    {
        return d[i][j][n];
    };
    
    /* For indexing a const TensorField */
    real_t  operator()(const int i, const int j, const int n) const
    {
        return d[i][j][n];
    };
};

using MatrixField = TensorField; // Use when the dimensions are
                                 // logically disconnected


struct VectorField
{
    real_t* d[3];
};


/* By analogy with the above. Not necessary, might remove. */
struct ScalarField
{
    real_t* d;
    
    /* Use () instead of [] for consistency with the 2D case */
    real_t& operator()(const int n)
    {
        return d[n];
    };

    real_t  operator()(const int n) const
    {
        return d[n];
    };
};
#endif
