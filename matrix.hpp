#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "common.hpp"

/* Class to do basic operations on individual 3x3 matrices,
 * using the direct methods and notation of:
 *     https://en.wikipedia.org/wiki/Invertible_matrix#Analytic_solution */


class Matrix
{
    private:
    bool cofactors_found;
    bool determinant_found;
    real_t a,b,c,d,e,f,g,h,i; // Elements of original matrix
    real_t A,B,C,D,E,F,G,H,I; // Elements of the cofactor matrix
                              // ie C is the cofactor corresponding to c (0,2)

    public:
    real_t  arr[3][3]; // Original matrix stored as an array
    real_t  inv[3][3]; // Inverse matrix
    real_t  det;


    Matrix();
    Matrix(real_t (&m)[3][3], bool storeArray=true);
    void fill(real_t (&m)[3][3]);
    void find_cofactors();
    void find_determinant();
    void find_inverse();
    void test();
};


#endif
