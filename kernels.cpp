// Strictly speaking this file should import the absolute minimum
// --- it probably shouldn't need anything except the definition of real_t ?
#include "kernels.hpp"
#include "physics_includes.hpp"
#include "numerical_flux.hpp"
#include "boundary_conditions.hpp"

//#include<stdio.h>

namespace kernels
{
    /* Solution points always use a fixed order (i along dir-0 etc) while 
     * flux-point data is stored relative to the transform direction. Use this
     * to create fixed-axis pointers to the correct transform-relative indices */
    inline static void relative_to_fixed_indices(int&  n0, int&  n1, int&  n2,
                                                 int*& i0, int*& i1, int*& i2,
                                                 const int dir)
    {
        switch(dir)
        {
            case 0:
                i0  = &n0;
                i1  = &n1;
                i2  = &n2;
                break;
            case 1:
                i0  = &n2;
                i1  = &n0;
                i2  = &n1;
                break;
            case 2:
                i0  = &n1;
                i1  = &n2;
                i2  = &n0;
                break;
        }

        return;
    }

   
    inline static void fluxes_phys_to_ref(const real_t      (*Fphys)[3],
                                                real_t* __restrict__ F,
                                          const VectorField          S,
                                          const int memory_loc,
                                          const int Nfield,
                                          const int Nf_tot)
    {
        real_t lsum;

        for (int field = 0; field < Nfield; ++field)
        {
            lsum = 0.0;
            for (int d: dirs) 
                lsum += S(d,memory_loc) * Fphys[field][d];

            F[memory_loc + field * Nf_tot] = lsum; // Save into reference fluxes
        }

        return;
    }


    /*** End of static functions ***/



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
    /* Removed __restrict__ from v1 and result, since one of the
     * calls has these being the same array. */
    void add_2_vectors(real_t* v1,     real_t* __restrict__ v2, 
                       real_t  a1,     real_t               a2, 
                       real_t* result, const int N)
    {
        for (int i = 0; i < N; ++i)
            result[i] = a1 * v1[i] + a2 * v2[i];

        return;
    }


    /* Can have the same pointer in v1 or v2 and result */
    void add_3_vectors(const real_t* const v1, const real_t* const v2, 
                       const real_t* const __restrict__ v3,
                       const real_t a1, const real_t a2, const real_t a3, 
                       real_t* const result, const int N)
    {
        for (int i = 0; i < N; ++i)
            result[i] = a1 * v1[i] + a2 * v2[i] + a3 * v3[i];

        return;
    }


    void add_vectors_inPlace(       real_t* const __restrict__ v,
                              const real_t* const __restrict__ v_add,
                              const int N)
    {
        for (int i = 0; i < N; ++i)
            v[i] += v_add[i];

        return;
    }


    void add_scaled_vectors_inPlace(      real_t* const __restrict__ v,
                                    const real_t* const __restrict__ v_add,
                                    const real_t                     scalar,
                                    const int N)
    {
        for (int i = 0; i < N; ++i)
            v[i] += scalar * v_add[i];

        return;
    }


    /* Element-by-element product along two arrays, giving a third
     * array of the same size. */
#if 0
    void product_2_vectors(const real_t* const __restrict__ v1,
                           const real_t* const __restrict__ v2,
                                 real_t* const __restrict__ result,
                           const int N)
    {
        for (int i = 0; i < N; ++i)
            result[i] = v1[i] * v2[i];

        return;
    }
#endif


    void multiply_by_scalar(const real_t* const __restrict__ v, 
                            const real_t                     scalar,
                                  real_t* const __restrict__ result,
                            const int                        N)
    {
        for (int i = 0; i < N; ++i)
            result[i] = scalar * v[i];

        return;
    }


    void multiply_by_scalar_inPlace(      real_t* const __restrict__ v, 
                                    const real_t                     scalar,
                                    const int                        N)
    {
        for (int i = 0; i < N; ++i)
            v[i] *= scalar;

        return;
    }


    void soln_to_flux(const real_t* const __restrict__ matrix, 
                      const real_t* const __restrict__ U, 
                            real_t* const __restrict__ Uf, 
                      const LengthBucket lb, const int dir)
    {
        int id_elem_s, id_elem_f;
        int mem_f, mem_offset_f;
        int mem_s, mem_offset_s;
        int field_offset_f, field_offset_s;
        int mem_matrix;
        real_t lsum;

        int n0_s, n0_f, n1, n2; // cyclic indices -- n0 is transform dir, n1 is n0+1 etc
        int* i = nullptr;       // fixed indices; i always points in 0 direction etc
        int* j = nullptr;
        int* k = nullptr;

        int ne0, ne1, ne2; // Transform-direction-relative indices for elements
        int* ie = nullptr; // fixed indices for element order
        int* je = nullptr;
        int* ke = nullptr;
        
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
        relative_to_fixed_indices(n0_s, n1, n2,  i,  j,  k,  dir);
        relative_to_fixed_indices(ne0, ne1, ne2, ie, je, ke, dir);

        for(int field = 0; field < lb.Nfield; ++field)
        {
            field_offset_s = field * lb.Ns_block;
            field_offset_f = field * lb.Nf_dir_block[dir];

            for(ne2 = 0; ne2 < lb.Nelem[dir2]; ++ne2)
            for(ne1 = 0; ne1 < lb.Nelem[dir1]; ++ne1)
            for(ne0 = 0; ne0 < lb.Nelem[dir];  ++ne0)
            {
                id_elem_s    = ((*ke)*lb.Nelem[1] + *je)*lb.Nelem[0] + *ie;
                mem_offset_s = field_offset_s + id_elem_s * lb.Ns_elem;

                id_elem_f    = (ne2*lb.Nelem[dir1] + ne1)*lb.Nelem[dir] + ne0;
                mem_offset_f = field_offset_f + id_elem_f * lb.Nf_dir[dir];

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
                        mem_s      = mem_offset_s + (*k * lb.Ns[1]  + *j) * lb.Ns[0] + *i;
                        mem_matrix = n0_f * Ns0 + n0_s; 
                        lsum      += matrix[mem_matrix] * U[mem_s];
                    }

                    mem_f = mem_offset_f + (n2 * Ns1 + n1) * Nf0 + n0_f;
                    Uf[mem_f] = lsum;
                }
            }
        }

        return;
    }


    void fluxDeriv_to_soln(const real_t* const __restrict__ matrix, 
                           const real_t* const __restrict__ F, 
                                 real_t* const __restrict__ dF, 
                           const LengthBucket lb, const int dir)
    {
        int id_elem_s, id_elem_f;
        int mem_f, mem_offset_f;
        int mem_s, mem_offset_s;
        int field_offset_f, field_offset_s;
        int mem_matrix;
        real_t lsum;

        int n0_s, n0_f, n1, n2; // cyclic indices -- n0 is transform dir, n1 is n0+1 etc
        int* i = nullptr;       // fixed indices; i always points in 0 direction etc
        int* j = nullptr;
        int* k = nullptr;

        int ne0, ne1, ne2; // Transform-direction-relative indices for elements
        int* ie = nullptr; // fixed indices for element order
        int* je = nullptr;
        int* ke = nullptr;
      
        int dir1 = dir_plus_one[dir]; // array stored in common.hpp for convenience
        int dir2 = dir_plus_two[dir]; // dir is the transform direction

        /* Transform-relative lengths for clarity */
        int Nf0 = lb.Nf[dir];
        int Ns1 = lb.Ns[dir1];

        relative_to_fixed_indices(n0_s, n1, n2,  i,  j,  k,  dir);
        relative_to_fixed_indices(ne0, ne1, ne2, ie, je, ke, dir);

        for(int field = 0; field < lb.Nfield; ++field)
        {
            field_offset_s = field * lb.Ns_block;
            field_offset_f = field * lb.Nf_dir_block[dir];

            for(ne2 = 0; ne2 < lb.Nelem[dir2]; ++ne2)
            for(ne1 = 0; ne1 < lb.Nelem[dir1]; ++ne1)
            for(ne0 = 0; ne0 < lb.Nelem[dir];  ++ne0)
            {
                id_elem_s    = ((*ke)*lb.Nelem[1] + *je)*lb.Nelem[0] + *ie;
                mem_offset_s = field_offset_s + id_elem_s * lb.Ns_elem;

                id_elem_f    = (ne2*lb.Nelem[dir1] + ne1)*lb.Nelem[dir] + ne0;
                mem_offset_f = field_offset_f + id_elem_f * lb.Nf_dir[dir];

                /* Do this main loop using transform-direction-relative indices, since
                 * the flux-point and transform-matrix arrays use these. Need something
                 * like the i,j,k pointers for the fixed solution-point indices */
                for (n2   = 0; n2   < lb.Ns[dir2]; ++n2)
                for (n1   = 0; n1   < lb.Ns[dir1]; ++n1)
                for (n0_s = 0; n0_s < lb.Ns[dir];  ++n0_s)
                {
                    lsum = 0.0;
                    for (n0_f = 0; n0_f < lb.Nf[dir]; ++n0_f)
                    {
                        mem_f      = mem_offset_f + (n2 * Ns1 + n1) * Nf0 + n0_f;
                        mem_matrix = n0_s * Nf0 + n0_f; 
                        lsum      += matrix[mem_matrix] * F[mem_f];
                    }

                    mem_s = mem_offset_s + (*k * lb.Ns[1]  + *j) * lb.Ns[0] + *i;
                    dF[mem_s] = lsum; 
                }
            }
        }

        return;
    }


    /* Think I can replace with a single for loop */
    void bulk_fluxes(const real_t* const __restrict__ Uf,
                           real_t* const __restrict__ F ,
                     const VectorField                S ,
                     const Physics* const __restrict__ physics,
                     //const ConservedToPrimitive*  U_to_P,
                     //const FluxesFromPrimitive* F_from_P,
                     const LengthBucket lb, const int dir)
    {
        int id_elem;
        int mem_offset;
        int mem;
        int dir1 = dir_plus_one[dir]; 
        int dir2 = dir_plus_two[dir];

        int Nf0 = lb.Nf[dir];
        int Ns1 = lb.Ns[dir1];
        int Nf_tot = lb.Nf_dir_block[dir]; //Total # of this-dir flux points in the block

        /* Conserved vars, primitive vars, and physical-direction fluxes
         * at a single point */
        real_t* Up      = new real_t [lb.Nfield];
        real_t* Pp      = new real_t [lb.Nfield];
        real_t (*Fp)[3] = new real_t [lb.Nfield][3]; // pointer to an array


        for (int ne2 = 0; ne2 < lb.Nelem[dir2]; ++ne2)
        for (int ne1 = 0; ne1 < lb.Nelem[dir1]; ++ne1)
        for (int ne0 = 0; ne0 < lb.Nelem[dir];  ++ne0)
        {
            id_elem = (ne2*lb.Nelem[dir1] + ne1)*lb.Nelem[dir] + ne0;
            mem_offset = id_elem * lb.Nf_dir[dir];

            for (int n2 = 0; n2 < lb.Ns[dir2]; ++n2)
            for (int n1 = 0; n1 < lb.Ns[dir1]; ++n1)
            for (int n0 = 0; n0 < lb.Nf[dir];  ++n0)
            {
                /* Memory location for the index-0 field */
                mem = mem_offset + (n2 * Ns1 + n1) * Nf0 + n0;

                /* Load all field variables at this location into
                 * the Up array. */
                for (int field = 0; field < lb.Nfield; ++field)
                    Up[field] = Uf[mem + field * Nf_tot];

                //(*U_to_P)(Up, Pp); // conserved -> primitive variables
                physics->ConservedToPrimitive(Up, Pp);

                /* Uf is the physical solution at flux points
                 * Calculate physical fluxes in all three dirs */
                //(*F_from_P)(Pp, Fp);
                physics->Fluxes(Pp, Fp);

                /* Transform from physical to reference-space fluxes
                 * for all fields */
                fluxes_phys_to_ref(Fp, F, S, mem, lb.Nfield, Nf_tot);
            }
        }

        delete[] Up;
        delete[] Pp;
        delete[] Fp;

        return;
    }


    void flux_divergence(const VectorField                dF,
                         const real_t* const __restrict__ Jrdetg,
                               real_t* const __restrict__ divF,
                         const LengthBucket lb)
    {
        int id_elem;
        int mem_offset;
        int field_offset;
        int mem, mem0;

        for(int field = 0; field < lb.Nfield; ++field)
        {
            field_offset = field * lb.Ns_block;

            for (int ke = 0; ke < lb.Nelem[2]; ++ke)
            for (int je = 0; je < lb.Nelem[1]; ++je)
            for (int ie = 0; ie < lb.Nelem[0]; ++ie)
            {
                id_elem = (ke*lb.Nelem[1] + je)*lb.Nelem[0] + ie;
                mem_offset = field_offset + id_elem * lb.Ns_elem;

                for (int k = 0; k < lb.Ns[2]; ++k)
                for (int j = 0; j < lb.Ns[1]; ++j)
                for (int i = 0; i < lb.Ns[0]; ++i)
                {
                    mem  = mem_offset + (k * lb.Ns[1]  + j) * lb.Ns[0] + i;
                    mem0 = mem - field_offset; // Jrdetg has only one field's worth

                    divF[mem] = (dF(0,mem) + dF(1,mem) + dF(2,mem)) / Jrdetg[mem0];
                }
            }
        }

        return;
    }


#if 0
    /* If only adding the source, could replace with add_2_vectors() */
    void scalar_field_source(      real_t* const __restrict__ divF,
                             const real_t* const __restrict__ U, 
                             const LengthBucket lb,
                             const real_t c_h, const real_t damping_rate)
    {
        int id_elem;
        int mem_offset;
        int field_offset;
        int mem;

        const int field = 8; // psi is the 9th field
        //const real_t chsq   = c_h*c_h;

        field_offset = field * lb.Ns_block;

        /* Replace by a single loop from 0 to Ns_block */
        for (int ke = 0; ke < lb.Nelem[2]; ++ke)
        for (int je = 0; je < lb.Nelem[1]; ++je)
        for (int ie = 0; ie < lb.Nelem[0]; ++ie)
        {
            id_elem = (ke*lb.Nelem[1] + je)*lb.Nelem[0] + ie;
            mem_offset = field_offset + id_elem * lb.Ns_elem;

            for (int k = 0; k < lb.Ns[2]; ++k)
            for (int j = 0; j < lb.Ns[1]; ++j)
            for (int i = 0; i < lb.Ns[0]; ++i)
            {
                mem = mem_offset + (k * lb.Ns[1]  + j) * lb.Ns[0] + i;

                //divF[mem] = chsq * divF[mem];
                //divF[mem] = chsq * divF[mem] + damping_rate * U[mem];
                divF[mem] += damping_rate * U[mem];
            }
        }

        return;
    }


    /* Replace by multiply_by_scalar() or similar */
    void store_divB(const real_t* const __restrict__ divF,
                          real_t* const __restrict__ divB, 
                    const real_t c_h,
                    const LengthBucket lb)
    {
        int id_elem;
        int mem_offset;
        int field_offset;
        int mem, mem0;

        const int field = 8; // psi is the 9th field
        const real_t over_chsq   = 1.0 / (c_h*c_h);

        field_offset = field * lb.Ns_block;

        /* Replace with single loop --- don't need to break into elements */
        for (int ke = 0; ke < lb.Nelem[2]; ++ke)
        for (int je = 0; je < lb.Nelem[1]; ++je)
        for (int ie = 0; ie < lb.Nelem[0]; ++ie)
        {
            id_elem = (ke*lb.Nelem[1] + je)*lb.Nelem[0] + ie;
            mem_offset = field_offset + id_elem * lb.Ns_elem;

            for (int k = 0; k < lb.Ns[2]; ++k)
            for (int j = 0; j < lb.Ns[1]; ++j)
            for (int i = 0; i < lb.Ns[0]; ++i)
            {
                mem  = mem_offset + (k * lb.Ns[1]  + j) * lb.Ns[0] + i;
                mem0 = mem - field_offset; // divB has only one field's worth

                /* Haven't multiplied by c_h^2 yet: divB is just psi's stored divF */
                divB[mem0] = over_chsq * divF[mem];
            }
        }

        return;
    }
#endif


    void fill_face_data(const real_t* const __restrict__ Uf,
                              FaceCommunicator           face,
                        const LengthBucket               lb)
    {
        const int dir  = face.normal_dir;
        const int dir1 = dir_plus_one[dir]; 
        const int dir2 = dir_plus_two[dir];
        int id_elem, id_elem_face;
        int mem_offset, mem_offset_face;
        int field_offset, field_offset_face;
        int mem, mem_face;

        const int Nf0 = lb.Nf[dir];
        const int Ns1 = lb.Ns[dir1];

        
        for(int field = 0; field < lb.Nfield; ++field)
        {
            field_offset      = field * lb.Nf_dir_block[dir];
            field_offset_face = field * face.Ntot;

            for (int ne2 = 0; ne2 < face.Nelem[1]; ++ne2)
            for (int ne1 = 0; ne1 < face.Nelem[0]; ++ne1)
            {
                id_elem    = (ne2*lb.Nelem[dir1] + ne1)*lb.Nelem[dir] + face.ne0;
                mem_offset = field_offset + id_elem * lb.Nf_dir[dir];

                id_elem_face    = (ne2 * lb.Nelem[dir1]) + ne1;
                mem_offset_face = field_offset_face 
                                      + id_elem_face * lb.Ns[dir2] * lb.Ns[dir1];

                /* Iterate over flux points on the relevant face
                 * of each element. */
                for (int n2 = 0; n2 < face.N[1]; ++n2)
                for (int n1 = 0; n1 < face.N[0]; ++n1)
                {
                    mem      = mem_offset      + (n2 * Ns1 + n1) * Nf0 + face.n0;
                    mem_face = mem_offset_face +  n2 * Ns1 + n1;
                    face.my_data[mem_face] = Uf[mem];
                }
            }
        }


        /* Modify data ordering for non-trivial grid connectivity, so that the
         * array's order corresponds to that expected by the matching process */
        if (face.change_data_order)
        {
            /* For the rearranged indexing into the my_data_to_send array */
            int id_elem_send;  
            int mem_offset_send;
            int mem_send;
            
            const int Ns2 = lb.Ns[dir2]; // same as face.N[1] I think --- should tidy here and above!

            for(int field = 0; field < lb.Nfield; ++field)
            {
                field_offset_face = field * face.Ntot;

                for (int ne2 = 0; ne2 < face.Nelem[1]; ++ne2)
                for (int ne1 = 0; ne1 < face.Nelem[0]; ++ne1)
                {
                    /* Original ordering */
                    id_elem_face    = ne2 * lb.Nelem[dir1] + ne1;
                    mem_offset_face = field_offset_face + id_elem_face * lb.Ns[dir2] * lb.Ns[dir1];

                    /* The 5 types of modified ordering, for elements... */
                    if (face.swap_index_order)
                    {
                        if (face.reverse_index_direction)
                        {
                            if (face.index_to_reverse == 1)
                                /* Reverse the post-swap "normal+1" index */
                                id_elem_send = ne1 * lb.Nelem[dir2] + face.Nelem[1] - 1 - ne2;
                            else
                                /* Reverse the post-swap "normal+2" index */
                                id_elem_send = (face.Nelem[0] - 1 - ne1) * lb.Nelem[dir2] + ne2;
                        }
                        else
                            /* Swapped ordering, no reversal */
                            id_elem_send = ne1 * lb.Nelem[dir2] + ne2;
                    }
                    else // No swapping, must have reversing
                    {
                        if (face.index_to_reverse == 1)
                            /* Reverse the original "normal+1" index */
                            id_elem_send = (ne2 * lb.Nelem[dir1]) + face.Nelem[0] - 1 - ne1;
                        else
                            /* Reverse the original "normal+2" index */
                            id_elem_send = (face.Nelem[1] - 1 - ne2) * lb.Nelem[dir1] + ne1;
                    }

                    mem_offset_send = field_offset_face + id_elem_send * lb.Ns[dir2] * lb.Ns[dir1];
                    
                    for (int n2 = 0; n2 < face.N[1]; ++n2)
                    for (int n1 = 0; n1 < face.N[0]; ++n1)
                    {
                        /* Original */
                        mem_face = mem_offset_face +  n2 * Ns1 + n1;
                        
                        /* The 5 types of modified ordering, on solution points ... */
                        if (face.swap_index_order)
                        {
                            if (face.reverse_index_direction)
                            {
                                if (face.index_to_reverse == 1)
                                    /* Reverse the post-swap "normal+1" index */
                                    mem_send = n1 * Ns2 + face.N[1] - 1 - n2;
                                else
                                    /* Reverse the post-swap "normal+2" index */
                                    mem_send = (face.N[0] - 1 - n1) * Ns2 + n2;
                            }
                            else
                                /* Swapped ordering, no reversal */
                                mem_send = n1 * Ns2 + n2;
                        }
                        else // No swapping, must have reversing
                        {
                            if (face.index_to_reverse == 1)
                                /* Reverse the original "normal+1" index */
                                mem_send = n2 * Ns1 + face.N[0] - 1 - n1;
                            else
                                /* Reverse the original "normal+2" index */
                                mem_send = (face.N[1] - 1 - n2) * Ns1 + n1;
                        }

                        /* Repack this field's data into face.my_data_to_send */
                        face.my_data_to_send[mem_offset_send + mem_send] = face.my_data[mem_face];
                    } // loop over solution points
                } // loop over elements
            } // end of loop over fields
        } // end of data rearrangement
        else
        {
            /* No rearrangement needed, just copy to the send array */
            for (int n = 0; n < face.Ntot_all; ++n)
                face.my_data_to_send[n] = face.my_data[n];
        }
        
        return;
    }


    void external_numerical_flux(const FaceCommunicator           face,
                                       real_t* const __restrict__ F,
                                 const NumericalFlux*             F_numerical,
                                 const VectorField                S,
                                 const LengthBucket               lb)
    {
        const int dir  = face.normal_dir;
        const int dir1 = dir_plus_one[dir]; 
        const int dir2 = dir_plus_two[dir];
        int id_elem, id_elem_face;
        int mem_offset, mem_offset_face;
        int mem, mem_face;

        const int Nf0 = lb.Nf[dir];
        const int Ns1 = lb.Ns[dir1];
        const int Nf_tot = lb.Nf_dir_block[dir]; //Total # of this-dir flux points in the block

        const real_t* __restrict__ UL_data;
        const real_t* __restrict__ UR_data;

        /* To apply external boundary conditions, set the incoming neighbour data on
         * external faces such that the fluxes are as desired. */
        if (face.external_face)
            for (int field = 0; field < lb.Nfield; ++field)
                for (int i = 0; i < face.Ntot; ++i)
                    face.neighbour_data[field*face.Ntot + i] = (*face.BC)(field, i, face.my_data);

        if (face.orientation > 0)
        {
            UL_data = face.my_data;
            UR_data = face.neighbour_data;
        }
        else
        {
            UL_data = face.neighbour_data;
            UR_data = face.my_data;
        }

        /* Conserved variables and physical fluxes at one point */
        real_t* UL = new real_t [lb.Nfield];
        real_t* UR = new real_t [lb.Nfield];
        real_t (*F_num_phys)[3] = new real_t [lb.Nfield][3]; // pointer to an array

        real_t np[3]; // The face's normal vector at a single point
        
        for (int ne2 = 0; ne2 < face.Nelem[1]; ++ne2)
        for (int ne1 = 0; ne1 < face.Nelem[0]; ++ne1)
        {
            id_elem_face    = (ne2 * lb.Nelem[dir1]) + ne1;
            mem_offset_face = id_elem_face * lb.Ns[dir2] * lb.Ns[dir1];

            id_elem    = (ne2*lb.Nelem[dir1] + ne1)*lb.Nelem[dir] + face.ne0;
            mem_offset = id_elem * lb.Nf_dir[dir];

            for (int n2 = 0; n2 < face.N[1]; ++n2)
            for (int n1 = 0; n1 < face.N[0]; ++n1)
            {
                /* Memory location for the zeroth field */
                mem_face = mem_offset_face +  n2 * Ns1 + n1;
                
                for (int field = 0; field < lb.Nfield; ++field)
                {
                    UL[field] = UL_data[mem_face + field * face.Ntot];
                    UR[field] = UR_data[mem_face + field * face.Ntot];
                }

                for (int i: dirs)
                    np[i] = face.normal(i,mem_face);

                (*F_numerical)(UL, UR, np, F_num_phys, dir);

                /* Transform from physical to reference-space fluxes.
                 * Required memory location in the full 3D flux array */
                mem = mem_offset + (n2 * Ns1 + n1) * Nf0 + face.n0;

                fluxes_phys_to_ref(F_num_phys, F, S, mem, lb.Nfield, Nf_tot);
            }
        }
        
        delete[] UL;
        delete[] UR;
        delete[] F_num_phys;

        return;
    }


    void internal_numerical_flux(const real_t* const __restrict__ Uf,
                                       real_t* const __restrict__ F,
                                 const NumericalFlux*             F_numerical,
                                 const VectorField                S,
                                 const VectorField                normal,
                                 const LengthBucket               lb,
                                 const int                        dir)
    {
        const int Nif = lb.Nelem[dir] - 1; // No. of surfaces of internal faces
                                           // normal to direction dir
        const int dir1 = dir_plus_one[dir];
        const int dir2 = dir_plus_two[dir];

        const int Nf0 = lb.Nf[dir];
        const int Ns1 = lb.Ns[dir1];

        /* Total # of this-direction flux points in the block */
        const int Nf_tot = lb.Nf_dir_block[dir]; 

        //int ne0L, ne0R; //Element indices in the normal dir on either side of face
        int id_elem_L, id_elem_R, id_elem_n;
        int mem_offset_L, mem_offset_R, mem_offset_n;
        int mem_L, mem_R, mem_n;

        /* Conserved variables and physical fluxes at one point */
        real_t* UL = new real_t [lb.Nfield];
        real_t* UR = new real_t [lb.Nfield];
        real_t (*F_num_phys)[3] = new real_t [lb.Nfield][3]; // pointer to an array

        real_t np[3]; // The face's normal vector at a single point

        for (int ne2 = 0; ne2 < lb.Nelem[dir2]; ++ne2)
        for (int ne1 = 0; ne1 < lb.Nelem[dir1]; ++ne1)
        for (int ne0 = 0; ne0 < Nif; ++ne0)
        {
            /* For each element, do the face between this one and the one to its
             * right, ie at greater dir-coord and element index in this direction */
            id_elem_L = (ne2*lb.Nelem[dir1] + ne1)*lb.Nelem[dir] + ne0;
            id_elem_R = id_elem_L + 1;
            id_elem_n = (ne2*lb.Nelem[dir1] + ne1)*(lb.Nelem[dir]+1) + ne0; // For normals

            mem_offset_L = id_elem_L * lb.Nf_dir[dir];
            mem_offset_R = id_elem_R * lb.Nf_dir[dir];
            mem_offset_n = id_elem_n * lb.Ns[dir2] * lb.Ns[dir1];
        
            for (int n2 = 0; n2 < lb.Ns[dir2]; ++n2)
            for (int n1 = 0; n1 < lb.Ns[dir1]; ++n1)
            {
                /* Memory locations for the 0th field variable */
                mem_L = mem_offset_L + (n2 * Ns1 + n1) * Nf0 + Nf0 - 1;
                mem_R = mem_offset_R + (n2 * Ns1 + n1) * Nf0 + 0;
                mem_n = mem_offset_n +  n2 * Ns1 + n1;

                for (int field = 0; field < lb.Nfield; ++field)
                {
                    UL[field] = Uf[mem_L + field * Nf_tot];
                    UR[field] = Uf[mem_R + field * Nf_tot];
                }

                for (int i: dirs)
                    np[i] = normal(i,mem_n);

                (*F_numerical)(UL, UR, np, F_num_phys, dir);

                /* Calculate reference-space flux and save into the left element */
                fluxes_phys_to_ref(F_num_phys, F, S, mem_L, lb.Nfield, Nf_tot);

                /* Copy into the right element */
                for (int field = 0; field < lb.Nfield; ++field)
                    F[mem_R + field * Nf_tot] = F[mem_L + field * Nf_tot];
            }
        }

        delete[] UL;
        delete[] UR;
        delete[] F_num_phys;

        return;
    }


    /* Average the conserved variables or their derivatives on process-external interfaces.
     * Only required when have diffusive terms. */
    void external_interface_average(const FaceCommunicator           face,
                                          real_t* const __restrict__ Uf,
                                    const LengthBucket               lb)
    {
        const int dir  = face.normal_dir;
        const int dir1 = dir_plus_one[dir]; 
        const int dir2 = dir_plus_two[dir];
        int id_elem, id_elem_face;
        int mem_offset, mem_offset_face;
        int mem, mem_face, mem_face_field;

        const int Nf0 = lb.Nf[dir];
        const int Ns1 = lb.Ns[dir1];
        const int Nf_tot = lb.Nf_dir_block[dir]; //Total # of this-dir flux points in the block

        /* To apply external boundary conditions, set the incoming neighbour data on
         * external faces such that the fluxes are as desired. */
        /** Will need to take account of this ... 
        if (face.external_face)
            for (int field = 0; field < lb.Nfield; ++field)
                for (int i = 0; i < face.Ntot; ++i)
                    face.neighbour_data[field*face.Ntot + i] = (*face.BC)(field, i, face.my_data);
         **/

        
        for (int ne2 = 0; ne2 < face.Nelem[1]; ++ne2)
        for (int ne1 = 0; ne1 < face.Nelem[0]; ++ne1)
        {
            id_elem_face    = (ne2 * lb.Nelem[dir1]) + ne1;
            mem_offset_face = id_elem_face * lb.Ns[dir2] * lb.Ns[dir1];

            id_elem    = (ne2*lb.Nelem[dir1] + ne1)*lb.Nelem[dir] + face.ne0;
            mem_offset = id_elem * lb.Nf_dir[dir];

            for (int n2 = 0; n2 < face.N[1]; ++n2)
            for (int n1 = 0; n1 < face.N[0]; ++n1)
            {
                /* Memory location on the face, for the zeroth field */
                mem_face = mem_offset_face +  n2 * Ns1 + n1;
                /* Memory location in the full 3D flux array, zeroth field */
                mem      = mem_offset + (n2 * Ns1 + n1) * Nf0 + face.n0;
                
                for (int field = 0; field < lb.Nfield; ++field)
                {
                    mem_face_field = mem_face + field * face.Ntot;
                    Uf[mem + field * Nf_tot] = 0.5 * (face.my_data[mem_face_field] +
                                                      face.neighbour_data[mem_face_field]);
                }
            }
        }
        
        return;
    }


    /* Average the conserved variables or their derivatives on process-internal interfaces.
     * Only required when have diffusive terms. */
    void internal_interface_average(      real_t* const __restrict__ Uf,
                                    const LengthBucket               lb,
                                    const int                        dir)
    {
        const int Nif = lb.Nelem[dir] - 1; // No. of surfaces of internal faces
                                           // normal to direction dir
        const int dir1 = dir_plus_one[dir];
        const int dir2 = dir_plus_two[dir];

        const int Nf0 = lb.Nf[dir];
        const int Ns1 = lb.Ns[dir1];

        /* Total # of this-direction flux points in the block */
        const int Nf_tot = lb.Nf_dir_block[dir]; 

        //int ne0L, ne0R; //Element indices in the normal dir on either side of face
        int id_elem_L, id_elem_R;
        int mem_offset_L, mem_offset_R;
        int mem_L, mem_R, mem_L0, mem_R0;

        for (int ne2 = 0; ne2 < lb.Nelem[dir2]; ++ne2)
        for (int ne1 = 0; ne1 < lb.Nelem[dir1]; ++ne1)
        for (int ne0 = 0; ne0 < Nif; ++ne0)
        {
            /* For each element, do the face between this one and the one to its
             * right, ie at greater dir-coord and element index in this direction */
            id_elem_L = (ne2*lb.Nelem[dir1] + ne1)*lb.Nelem[dir] + ne0;
            id_elem_R = id_elem_L + 1;

            mem_offset_L = id_elem_L * lb.Nf_dir[dir];
            mem_offset_R = id_elem_R * lb.Nf_dir[dir];
        
            for (int n2 = 0; n2 < lb.Ns[dir2]; ++n2)
            for (int n1 = 0; n1 < lb.Ns[dir1]; ++n1)
            {
                /* Memory locations for the 0th field variable */
                mem_L0 = mem_offset_L + (n2 * Ns1 + n1) * Nf0 + Nf0 - 1;
                mem_R0 = mem_offset_R + (n2 * Ns1 + n1) * Nf0 + 0;

                for (int field = 0; field < lb.Nfield; ++field)
                {
                    mem_L = mem_L0 + field * Nf_tot;
                    mem_R = mem_R0 + field * Nf_tot;

                    Uf[mem_L] = 0.5 * (Uf[mem_L] + Uf[mem_R]);
                    Uf[mem_R] = Uf[mem_L];
                }
            }
        }

        return;
    }

    
    void gradient_ref_to_phys(const VectorField dU_ref,
                                    VectorField dU_phys,
                              const VectorField dxdr[3],
                              const LengthBucket lb)
    {
        int field_offset;
        int mem;
        real_t lsum;

        for (int field = 0; field < lb.Nfield; ++field)
        {
            field_offset = field * lb.Ns_block;

            for (int i = 0; i < lb.Ns_block; ++i)
            {
                mem = i + field_offset;

                for (int phys: dirs)
                {
                    lsum = 0.0;
                    for (int ref: dirs)
                        lsum += dxdr[ref](phys, i) * dU_ref(ref, mem);

                    dU_phys(phys, mem) = lsum;
                }
            }
        }

        return;
    }

    
    void diffusive_flux(const real_t* const __restrict__ Uf,
                        const VectorField                dUf,
                              real_t* const __restrict__ F,
                        //const DiffusiveFluxes*           F_diff,
                        const Physics* const __restrict__ physics,
                        const real_t* const __restrict__ args,
                        const VectorField                S,
                        const LengthBucket               lb,
                        const int                        dir)
    {
        int N = lb.Nf_dir_block[dir]; //Total # of this-dir flux points in the block

        /* Conserved vars, primitive vars, and physical-direction fluxes
         * at a single point */
        real_t* Up       = new real_t [lb.Nfield];
        real_t (*dUp)[3] = new real_t [lb.Nfield][3]; 
        real_t (*F_diff_p)[3] = new real_t [lb.Nfield][3]; // pointer to an array


        for (int i = 0; i < N; ++i) // i = mem loc for zeroth field
        {
            /* Load all field variables at this location into
             * the Up array. */
            for (int field = 0; field < lb.Nfield; ++field)
            {
                Up[field] = Uf[i + field * N];
                
                for (int d: dirs)
                    dUp[field][d] = dUf(d, i + field*N);
            }

            //(*F_diff)(Up, dUp, F_diff_p, args);
            physics->DiffusiveFluxes(Up, dUp, F_diff_p, args);

            /* Transform from physical to reference-space fluxes
             * for all fields */
            fluxes_phys_to_ref(F_diff_p, F, S, i, lb.Nfield, N);
        }

        delete[] Up;
        delete[] dUp;
        delete[] F_diff_p;

        return;
    }



#if 0
    void fill_velocity_vector(      VectorField Vf,
                              const VectorField Uf,
                              const LengthBucket lb)
    {
        int N; // No. flux points in the block for a given dir
        real_t density;

        for (int d: dirs) // Fill for all 3 flux-point blocks
        {
           N = lb.Nf_dir_block[d];
           for (int i=0; i<N; ++i) // For each flux point
           {
               density = Uf(d,i); // Density is the zeroth field
               for (int vcomp=0; vcomp < 3; ++vcomp) // for each velocity component...
                   Vf(d,i + vcomp*N) = Uf(d,i + (vcomp+1)*N) / density; // mom0 is 1st field
           }
        }

        return;
    }
#endif


    /* Temporary standalone method. Should reintegrate into other operations. */
    /* Returns the timestep restriction along this reference direction.       */
    real_t local_timestep(const real_t* const __restrict__ Uf,
                          const VectorField timestep_transform,
                          const Physics* const __restrict__ physics,
                          //const ConservedToPrimitive*  U_to_P,
                          //const WaveSpeedsFromPrimitive* c_from_P,
                          const LengthBucket lb, const int dir)
    { 
        int id_elem;
        int mem_offset;
        int mem;
        int dir1 = dir_plus_one[dir]; 
        int dir2 = dir_plus_two[dir];

        int Nf0 = lb.Nf[dir];
        int Ns1 = lb.Ns[dir1];
        int Nf_tot = lb.Nf_dir_block[dir]; //Total # of this-dir flux points in the block

        /* Conserved vars, primitive vars, and max wavespeeds
         * at a single point */
        real_t* Up      = new real_t [lb.Nfield];
        real_t* Pp      = new real_t [lb.Nfield];
        real_t cp[3][2]; 
        real_t cm[3];
        real_t chi_v;
        real_t chi_v_max = 0.0;


        for (int ne2 = 0; ne2 < lb.Nelem[dir2]; ++ne2)
        for (int ne1 = 0; ne1 < lb.Nelem[dir1]; ++ne1)
        for (int ne0 = 0; ne0 < lb.Nelem[dir];  ++ne0)
        {
            id_elem = (ne2*lb.Nelem[dir1] + ne1)*lb.Nelem[dir] + ne0;
            mem_offset = id_elem * lb.Nf_dir[dir];

            for (int n2 = 0; n2 < lb.Ns[dir2]; ++n2)
            for (int n1 = 0; n1 < lb.Ns[dir1]; ++n1)
            for (int n0 = 0; n0 < lb.Nf[dir];  ++n0)
            {
                /* Memory location for the index-0 field */
                mem = mem_offset + (n2 * Ns1 + n1) * Nf0 + n0;

                /* Load all field variables at this location into
                 * the Up array. */
                for (int field = 0; field < lb.Nfield; ++field)
                    Up[field] = Uf[mem + field * Nf_tot];

                //(*U_to_P)(Up, Pp); // conserved -> primitive variables
                physics->ConservedToPrimitive(Up, Pp);

                /* Uf is the physical solution at flux points
                 * Calculate wavespeeds in all three physical dirs */
                for (int d: dirs) // iterate over phys-space directions
                {
                    //(*c_from_P)(Pp, cp[d], d);
                    physics->WaveSpeeds(Pp, cp[d], d);
                    cm[d] = MAX(cp[d][0], std::abs(cp[d][1]));
                }
                
                chi_v = 0.0;
                for (int d: dirs) // physical dir
                    chi_v += timestep_transform(d,mem) * cm[d];
                

                if (chi_v > chi_v_max)
                    chi_v_max = chi_v;
            }
        }

        delete[] Up;
        delete[] Pp;

        return(1.0 / chi_v_max);
    }

}
