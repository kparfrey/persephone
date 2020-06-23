#include "lagrange_polynomials.hpp"

namespace lagrange
{
    static double* barycentric_weights(const real_t* const x, const int N)
    {
        double* weights = new double [N];
        
        double p;
        for (int j = 0; j < N; ++j)
        {
            p = 1.0;
            for (int i = 0; i < N; ++i)
            {
                if (i != j)
                    p *= x[j] - x[i];
            }

            weights[j] = 1.0 / p;
        }

        return weights;
    }


    /* Single interpolation matrix: from x0 points --> x1 points */
    real_t* barycentric_interpolation_matrix(const real_t* const x0, const int N0,
                                             const real_t* const x1, const int N1)
    {
        /* Store like N1 x N0 matrix: the 0/original-point index moves fastest */
        real_t* matrix = new real_t [N1*N0];

        double* w = barycentric_weights(x0, N0);

        double s;
        for (int k = 0; k < N1; ++k)
        {
            s = 0.0;
            for (int i = 0; i < N0; ++i)
                s += w[i] / (x1[k] - x0[i]);

            for (int j = 0; j < N0; ++j)
                matrix[k*N0 + j] = w[j] / (s * (x1[k] - x0[j]));
        }

        delete[] w;

        return matrix;
    }


    real_t* differentiation_matrix(const real_t* const x, const int N)
    {
        real_t* D = new double [N*N]();

        double* w = barycentric_weights(x, N);
        
        /* i != j entries */
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i)
                if (i != j)
                    D[i*N + j] = (w[j]/w[i]) / (x[i] - x[j]);

        /* i == j entries */
        double s;
        for (int j = 0; j < N; ++j)
        {
            s = 0.0;
            for (int i = 0; i < N; ++i)
                if (i != j)
                    s += D[j*N + i];
        
            D[j*N + j] = - s;
        }
        
        delete[] w;

        return D;
    }


    /* A = B.C --- A:(Nrows_B x Ncols_C)  B: (Nrows_B x Nsum)  C: (Nsum x Ncols_C) */
    real_t* matrix_matrix_product(const real_t* const B, const real_t* const C,
                                  const int Nrows_B, const int Ncols_C, const int Nsum)
    {
        real_t* A = new double [Nrows_B * Ncols_C];
        real_t sum;

        for (int i = 0; i < Nrows_B; ++i)
            for (int j = 0; j < Ncols_C; ++j)
            {
                sum = 0.0;
                for (int k = 0; k < Nsum; ++k)
                    sum += B[i*Nsum + k] * C[k*Ncols_C + j];
                A[i*Ncols_C + j] = sum;
            }

        return A;
    }
}