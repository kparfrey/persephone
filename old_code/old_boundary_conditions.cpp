   /* Maybe a WallBC function belongs in each Physics derived class? */
    static void torus_BC(      real_t* const __restrict__ UL, 
                               real_t* const __restrict__ UR,
                         const real_t* const __restrict__ nl,
                         const Physics* const __restrict__ physics,
                         const int mem, const int orientation)
    {
        enum conserved {Density, mom0, mom1, mom2, tot_energy, B0, B1, B2, psi};
        enum primitive {density, v0  , v1  , v2,   pressure};

        /* Note: initially both UL and UR hold this element's data, so only need to
         * modify Ubc for those variables where different values are needed. */
        //real_t* Umy; // Points to the real data from this element
        //real_t* Ubc; // Points to the data for the imaginary boundary element
        
        real_t P[9];  // Primitives from incoming conserved variables

        real_t nu[3]; // contravariant components of normal vector
        physics->metric->raise(nl, nu, mem);

        /* UL and UR should be identical, so just choose L arbitrarily */
        physics->ConservedToPrimitive(UL, P, mem);

        /***
         * Remove the orientation parameter if keep this commented out.
        if (orientation > 0)
        {
            Umy = UL;
            Ubc = UR;
        }
        else
        {
            Umy = UR;
            Ubc = UL;
        }
         ***/

        //if (physics->apply_floors)
        //    UL[Density] = UR[Density] = P[density];
        
        //const real_t vdotn = nl[0]*P[v0] + nl[1]*P[v1] + nl[2]*P[v2];
        //real_t vm[3];
        const real_t Bdotn = nl[0]*P[B0] + nl[1]*P[B1] + nl[2]*P[B2];
        real_t Bm[3];

        /* Remove the normal velocity and magnetic field */
        for (int d: dirs)
        {
            //vm[d] = P[v0+d] - vdotn * nu[d]; // Impenetrable -- seems more stable than no-slip?
            //vm[d] = 0.0; // Impermeability + no slip
            Bm[d] = P[B0+d] - Bdotn * nu[d]; // Zero normal flux
        }

        /* The density, tot_energy, and psi conserved vars are unchanged.
         * Only need to overwrite the momentum and magnetic field */
        //real_t vl[3];
        //physics->metric->lower(vm, vl, mem);

        //UL[Density] = UR[Density] = P[Density] = 1.0;

        for (int d: dirs)
        {
            //UL[mom0+d] = UR[mom0+d] = P[Density] * vm[d]; // For impermeable only; need to lower vm
            UL[mom0+d] = UR[mom0+d] = 0.0; // No slip
            UL[  B0+d] = UR[  B0+d] = Bm[d];

            //Ubc[mom0+d] = - Umy[mom0+d]; // Only set non-slip for the imaginary adjoining element
                               // Should do this or v_bc = - v_my
            //Ubc[B0+d]   = Bm[d] - Bdotn * nu[d]; // B_bc = B_tangential - B_normal

            /* Right face only */
            //UR[mom0+d] = 0.0; // No slip
            //UR[  B0+d] = Bm[d];
        }

        
        // UL[psi] = UR[psi] = 0.0;

        /****
        P[pressure] = 10.0 * P[density]; // Isothermal-ish wall at set temperature
        const real_t mag_density = 0.5 * physics->metric->square(Bm, mem);
        const real_t KE_density  = 0.0; // 0.5 * P[density] * physics->metric->square(vm, mem);
        const real_t psi_density = 0.5 * P[psi] * P[psi];
        const real_t thermal_density = P[pressure]/((MHD*)physics)->gm1;
        UL[tot_energy] = UR[tot_energy] = KE_density + mag_density + thermal_density + psi_density;
        ****/

        /* Reset total energy, since you've removed some kinetic and magnetic energy by 
         * chopping off the normal velocity and magnetic field. Not clear if this is 
         * necessary and/or a good idea ... */
        /***
        const real_t mag_density = 0.5 * physics->metric->square(Bm, mem);
        const real_t KE_density  = 0.5 * P[density] * physics->metric->square(vm, mem);
        UL[tot_energy] = UR[tot_energy] = KE_density + mag_density + P[pressure]/((MHD*)physics)->gm1
                                          + 0.5 * P[psi] * P[psi];
        ***/

        /* This seems to make things worse... */
        //if (physics->apply_floors)
        if (false)
        {
            /* Since pressure might have been increased inside cons-to-prim
             * Might want to use kinetic and magnetic energy from before the BC adjustments were applied?
             * Or is it better for the total energy to be consistent? */
            const real_t mag_density = 0.5 * physics->metric->square(Bm, mem);
            //const real_t KE_density  = 0.5 * P[density] * physics->metric->square(vm, mem);
            const real_t KE_density  = 0.0; // Zero for no-slip
            UL[tot_energy] = UR[tot_energy] = KE_density + mag_density + P[pressure]/((MHD*)physics)->gm1;
        }
         

        /* Set a minimum gas pressure floor on the external boundary. 
         * Seems to help delay p < 0 a short while. */
        /***
        const real_t p_floor = 5.0;
        if (P[pressure] < p_floor)
        {
            const real_t add_internal_e = (p_floor - P[pressure])/((MHD*)physics)->gm1;
            // UR is the neighbour face's data for the torus's external faces,
            // which are face 5 having positive orientation.
            UR[tot_energy] += add_internal_e;
        }
         ***/

        /* Characteristic-variable-based outflow boundary condition for psi/
         * Not yet clear if this is better than just leaving psi alone at
         * its extrapolated value.
         * Update for DESC equilibria: seems better to not impose a BC on psi at all... */
        //const real_t c_h = std::sqrt(physics->ch_sq);
        //UL[psi] = UR[psi] = P[psi] + c_h * Bdotn;
        /* Seems better to just apply as the external data? */
        //UR[psi] = P[psi] + c_h * Bdotn;
        //UL[psi] = UR[psi] = 0.0;

        return;
    }


    /* Neumann BC is inside the interface averaging */
    /* Average the primitive variables or their derivatives on process-external interfaces.
     * Only required when have diffusive terms. */
    void external_interface_average(const FaceCommunicator           face,
                                          real_t* const __restrict__ Pf,
                                    const LengthBucket               lb,
                                    const bool                       averaging_derivs)
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


        /* The BC for derivatives -- should have already applied BCs for the solution */
        if (face.domain_external_face && averaging_derivs)
        {
            /* Don't need to set the neighbour_data to my_data here when averaging the solution itself, 
             * since this was already done in external_numerical_flux() immediately after the Dirichlet
             * boundary conditions were applied.  */
            for (int i = 0; i < face.Ntot_all; ++i)
                face.neighbour_data[i] = face.my_data[i];

            /*** For an adiabatic wall ***/
            /***/
            int mem;
            for (int i = 0; i < face.Ntot; ++i)
            {
                mem = 4 * face.Ntot + i; // index for pressure slot, holding temperature
                face.neighbour_data[mem] = face.my_data[mem] = 0.0;
            }
            /***/

            /*** Set magnetic field derivatives to zero -- equivalent to setting eta = 0 there */
            /***
            int mem;
            for (int i = 0; i < face.Ntot; ++i)
            {
                for (int comp = 0; comp < 3; comp++)
                {
                    mem = (5+comp) * face.Ntot + i; // holds d_i B^comp
                    face.neighbour_data[mem] = face.my_data[mem] = 0.0;
                }
            }
             ***/

            /* Set all derivatives to zero on the boundary ?
             * Sending zero velocity & B gradients into the diffusive flux function
             * is equivalent to setting viscosity = resistivity = 0 at the boundary */
            //for (int i = 0; i < face.Ntot_all; ++i)
                //face.my_data[i] = face.neighbour_data[i] = 0.0;

            /***
            for (int j = 0; j < face.Ntot; ++j)
            {
                int i = j + 8*face.Ntot;
                face.my_data[i] = face.neighbour_data[i] = 0.0;
            }
            ***/

            /*
            for (int field = 5; field < 8; field++) // B only
                for (int j = 0; j < face.Ntot; ++j)
                {
                    int i = j + field*face.Ntot;
                    face.my_data[i] = face.neighbour_data[i] = 0.0;
                }
             */
        }

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
                    Pf[mem + field * Nf_tot] = 0.5 * (face.my_data[mem_face_field] +
                                                      face.neighbour_data[mem_face_field]);
                }
            }
        }
        
        return;
    }


