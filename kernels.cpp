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


    void soln_to_flux(const real_t* const __restrict__ soln2flux, 
                      const real_t* const __restrict__ Q, 
                            real_t* const __restrict__ Qf, 
                      const LengthBucket lb, const int dir)
    {
        int id_elem;
        int mem_f, mem_offset_f;
        int mem_s, mem_offset_s;
        int mem_s2f;
        real_t lsum;

        int n0_s, n0_f, n1, n2; // cyclic indices -- n0 is transform dir, n1 is n0+1 etc
        int *i, *j, *k;         // true or fixed indices; i always points in 0 direction etc
        
        int dir1 = dir_plus_one[dir]; // array stored in common.hpp for convenience
        int dir2 = dir_plus_two[dir]; // dir is the transform direction

        /* Transform-relative lengths for clarity */
        int Nf0 = lb.Nf[dir];
        int Ns0 = lb.Ns[dir]; 
        int Ns1 = lb.Ns[dir1];

        /* The i,j,k indices are used to index solution points, which use a single fixed 
         * (direction-independent) layout. Need to set i -> 0-dir, j -> 1-dir, k-> 2-dir
         * indices independent of the transform direction. There's probably a neater
         * way to do this.... */
        switch(dir)
        {
            case 0:
                i = &n0_s;
                j = &n1;
                k = &n2;
                break;
            case 1:
                i = &n2;
                j = &n0_s;
                k = &n1;
                break;
            case 2:
                i = &n1;
                j = &n2;
                k = &n0_s;
                break;
        }

        for (int ie = 0; ie < lb.Nelem[0]; ++ie)
        for (int je = 0; je < lb.Nelem[1]; ++je)
        for (int ke = 0; ke < lb.Nelem[2]; ++ke)
        {
            id_elem = (ie*lb.Nelem[1] + je)*lb.Nelem[2] + ke;
            mem_offset_s = id_elem * lb.Ns_elem;

            mem_offset_f = id_elem * lb.Nf_dir[dir];

            /* Do this main loop using transform-direction-relative indices, since
             * the flux-point and transform-matrix arrays use these. Need something
             * like the i,j,k pointers for the fixed solution-point indices */
            for (n2   = 0; n2   < lb.Ns[dir2]; ++n2)
            for (n1   = 0; n1   < lb.Ns[dir1]; ++n1)
            for (n0_f = 0; n0_f < lb.Nf[dir];  ++n0_f)
            {
                lsum = 0.0;
                for (n0_s = 0; n0_s < lb.Ns[dir]; ++n0_s)
                {
                    mem_s   = mem_offset_s + (*i * lb.Ns[1]  + *j) * lb.Ns[2] + *k;
                    mem_s2f = n0_f * Ns0 + n0_s; 
                    lsum   += soln2flux[mem_s2f] * Q[mem_s];
                }

                mem_f = mem_offset_f + (n2 * Ns1 + n1) * Nf0 + n0_f;
                Qf[mem_f] = lsum;
            }
        }

        return;
    }


    void generate_fluxes(const real_t* const __restrict__ Qf,
                               real_t* const __restrict__ F ,
                         const LengthBucket lb, const int dir)
    {
        int id_elem;
        int mem_offset;
        int mem;
        int dir1 = dir_plus_one[dir]; 
        int dir2 = dir_plus_two[dir];

        int Nf0 = lb.Nf[dir];
        int Ns1 = lb.Ns[dir1];

        real_t wave_speed[3] = {0.25, 0.0, 0.0};

        for (int ie = 0; ie < lb.Nelem[0]; ++ie)
        for (int je = 0; je < lb.Nelem[1]; ++je)
        for (int ke = 0; ke < lb.Nelem[2]; ++ke)
        {
            id_elem = (ie*lb.Nelem[1] + je)*lb.Nelem[2] + ke;
            mem_offset = id_elem * lb.Nf_dir[dir];

            for (int n2 = 0; n2 < lb.Ns[dir2]; ++n2)
            for (int n1 = 0; n1 < lb.Ns[dir1]; ++n1)
            for (int n0 = 0; n0 < lb.Nf[dir];  ++n0)
            {
                mem = mem_offset + (n2 * Ns1 + n1) * Nf0 + n0;

                F[mem] = wave_speed[dir] * Qf[mem];
            }
        }

        return;
    }
}
