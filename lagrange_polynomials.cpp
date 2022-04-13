#include "lagrange_polynomials.hpp"
#include <cmath>

namespace lagrange
{
    double* barycentric_weights(const real_t* const x, const int N)
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


    /* Interpolate f to location s using grid x and weights w */
    real_t interpolate(const real_t* const __restrict__ x,
                       const real_t* const __restrict__ w,
                       const real_t* const __restrict__ f,
                       const real_t s, const int N)
    {
        real_t factor;
        real_t above = 0.0;
        real_t below = 0.0;

        for (int i = 0; i < N; ++i)
        {
            factor = w[i] / (s - x[i]);
            above += f[i] * factor;
            below +=        factor;
        }
        
        return above / below;
    }


    /* Interpolate V to location s using grid x and weights w, store in Vs */
    void interpolate_vector(const real_t* const __restrict__ x,
                            const real_t* const __restrict__ w,
                            const VectorField V,
                            const real_t s, const int N,
                                  real_t* const __restrict__ Vs)
    {
        real_t factor;
        real_t above[3] = {0.0};
        real_t below    =  0.0;

        for (int i = 0; i < N; ++i)
        {
            factor = w[i] / (s - x[i]);
            for (int d: dirs)
                above[d] += V(d,i) * factor;
            below += factor;
        }
        
        for (int d: dirs)
            Vs[d] = above[d] / below;
        
        return ;
    }


    /* Find df/dx at x = s, for s not at a node; from Kopriva */
    real_t differentiate(const real_t* const __restrict__ x,
                         const real_t* const __restrict__ w,
                         const real_t* const __restrict__ f,
                         const real_t s, const int N)
    {
        real_t factor;
        real_t above = 0.0;
        real_t below = 0.0;

        real_t fs = interpolate(x, w, f, s, N); // Interpolate f(s)

        for (int i = 0; i < N; ++i)
        {
            factor = w[i] / (s - x[i]);
            above += factor * (fs - f[i])/(s - x[i]);
            below += factor;
        }
        
        return above / below;
    }


    /* Find df/dx at x = x[i], when x[i] is an interpolation node */
    real_t diff_at_node(const real_t* const __restrict__ x,
                        const real_t* const __restrict__ w,
                        const real_t* const __restrict__ f,
                        const int i, const int N)
    {
        real_t lsum = 0.0;

        for (int j = 0; j < N; ++j)
            if (i != j)
                lsum += w[j] * (f[i] - f[j])/(x[i] - x[j]);
        
        return -lsum/w[i];
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


    /* Assumes that the data will be located at the Chebyshev-Gauss nodes */
    real_t* chebyshev_filtering_matrix(const int N)
    {
        double alpha = 0.1;
        double s     = 32.0;
        double eta;

        real_t* F = new real_t [N * N](); // Final filtering matrix

        double* M    = new double [N * N]; // Forward: nodes to modes
        double* Minv = new double [N * N]; // Reverse: modes to nodes
        double* L    = new double [N];     // Diagonal filtering matrix
        double* LM   = new double [N * N]; // L * M

        /* NB: V[i,j] --> V[i*N + j] */

        /* Forward transform: type-2 DCT */
        for (int k = 0; k < N; ++k)
            for (int j = 0; j < N; ++j)
                M[k*N + j] = 2.0 * std::cos(0.5 * pi * k * (2.0*j + 1.0) / N);

        /* Reverse transform: type-3 DCT with normalisation */
        double norm = 1.0 / (2.0 * N);
        for (int j = 0; j < N; ++j)
        {
            Minv[j*N + 0] = norm;
            for (int k = 1; k < N; ++k)
                Minv[j*N + k] = norm * 2.0 * std::cos(0.5 * pi * k * (2.0*j + 1.0) / N);
        }

        /* Diagonal components of filtering matrix */
        for (int i = 0; i < N; ++i)
        {
            eta  = (double)i/(N - 1.0);
            L[i] = std::exp(-alpha * std::pow(eta, s)); 
        }

        /* Multiply filtering and forward-transform matrices */
        for (int k = 0; k < N; ++k)
            for (int j = 0; j < N; ++j)
                LM[k*N + j] = L[k] * M[k*N + j];

        /* Multiply reverse-transform and filtered-forward matrices */
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                for (int m = 0; m < N; ++m)
                    F[i*N + j] += (real_t)(Minv[i*N + m] * LM[m*N + j]);

        delete[] M;
        delete[] Minv;
        delete[] L;
        delete[] LM;

        return F;
    }
}
