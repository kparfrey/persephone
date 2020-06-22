// Strictly speaking this file should import the absolute minimum
// --- it probably shouldn't need anything except the definition of real_t ?
#include "kernels.hpp"

namespace kernels
{

    /* Zero-initialized by default */
    real_t* alloc(int n)
    {
        real_t *ptr = new real_t [n]();
        return ptr;
    }


    /* If for some reason a non-initialized array is preferred for speed */
    real_t* alloc_raw(int n)
    {
        real_t *ptr = new real_t [n];
        return ptr;
    }


    void free(real_t* a)
    {
        delete[] a;
        return;
    }


    /* for 0 <= i < N: result[i] = a1*v1[i] + a2*v2[i] */
    /* Think I can get away with the __restrict__ keywords even when
     * result aliases v1 or v2, because the function is so simple. */
    void add_2_vectors(real_t* __restrict__ v1, real_t* __restrict__ v2, 
                       real_t               a1, real_t               a2, 
                       real_t* __restrict__ result, const int N)
    {
        for (int i = 0; i < N; ++i)
            result[i] = a1 * v1[i] + a2 * v2[i];

        return;
    }


    void soln_to_flux(const VectorField soln2flux, 
                      const real_t* const __restrict__ Q, 
                            VectorField Qf, 
                      const LengthBucket lb)
    {
        int id_elem;
        int mem_f, mem_offset_f;
        int mem_s, mem_offset_s;
        int mem_s2f;
        real_t lsum;

        for (int ie = 0; ie < lb.Nelem[0]; ++ie)
        for (int je = 0; je < lb.Nelem[1]; ++je)
        for (int ke = 0; ke < lb.Nelem[2]; ++ke)
        {
            id_elem = (ie*lb.Nelem[1] + je)*lb.Nelem[2] + ke;
            mem_offset_s = id_elem * lb.Ns_elem;

            /* Direction 0 */
            mem_offset_f = id_elem * lb.Nf_dir[0];

            for (int j   = 0; j   < lb.Ns[1]; ++j)
            for (int k   = 0; k   < lb.Ns[2]; ++k)
            for (int i_f = 0; i_f < lb.Nf[0]; ++i_f)
            {
                lsum = 0.0;
                for (int i_s = 0; i_s < lb.Ns[0]; ++i_s)
                {
                    mem_s   = mem_offset_s + (i_s * lb.Ns[1]  + j) * lb.Ns[2] + k;
                    mem_s2f = i_f * lb.Ns[0] + i_s;
                    lsum   += soln2flux(0, mem_s2f) * Q[mem_s];
                }

                mem_f = mem_offset_f + (k   * lb.Nsf[0] + j) * lb.Nf[0] + i_f;
                Qf(0, mem_f) = lsum;
            }
        }


        return;
    }
}
