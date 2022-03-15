From main.cpp:

....

    proc.fill_external_boundary_data();
    
    proc.divB_subsystem_iterations(500, true);
    
    write_data(proc); // Data generally lives on device

    write::message("\nFinished setup, starting time advancement \n");

...

   while(proc.time < proc.end_time)
    {
        if ((proc.time - proc.time_last_write) > proc.dt_write)
            proc.is_output_step = true;

        proc.time_advance();
        proc.divB_subsystem_iterations(6, false);

        if (proc.is_output_step)
            write_data(proc);
    }
...


from process.cpp:

void Process::divB_subsystem_iterations(const int Niter, const bool initial_cleaning)
{
    ElementBlock& eb = elements;

    /* ch is the only wavespeed in the divB subsystem, so can choose it to be anything.
     * Choose ch = 1.0 */
    if (initial_cleaning)
    {
        Physics::ch    = 1.0;
        Physics::ch_sq = 1.0;
        Physics::psi_damping_rate = 1.0 / Physics::psi_damping_const; //p_d_c = c_r from Dedner
        dt = 0.9 / tt_max_global; 
    }

    const int Nstart  = 5 * eb.Ns_block; // Beginning of B0 data
    const int Nsubsys = 4 * eb.Ns_block; // 3 B componenents + Psi
    const int Ntot    = Nfield * eb.Ns_block; 
    const int psi     = 8 * eb.Ns_block; // mem location at which the psi field begins
    const real_t one_third  = 1.0/3.0;
    const real_t two_thirds = 2.0/3.0;

    real_t* fields_inter = kernels::alloc_raw(Ntot);
    real_t* divF         = kernels::alloc(Ntot);

    /* Set Psi to 0 at beginning and end --- decouples timescales between here and main loop */
        kernels::multiply_by_scalar_inPlace(&eb.fields[psi], 0.0, eb.Ns_block);

    for (int mem = 0; mem < Ntot; ++mem)
        fields_inter[mem] = eb.fields[mem];

    if (initial_cleaning)
        write::message("Starting initial-field divergence cleaning, " + std::to_string(Niter) + " iterations");

    for (int iter = 0; iter < Niter; ++iter)
    {
        find_divF_divB_subsystem(eb.fields, divF);
        kernels::add_2_vectors(&eb.fields[Nstart], &divF[Nstart], 
                               1.0      , -dt, 
                               &fields_inter[Nstart], Nsubsys);

        find_divF_divB_subsystem(fields_inter, divF);
        kernels::add_3_vectors(&eb.fields[Nstart], &fields_inter[Nstart], &divF[Nstart], 
                               0.75     , 0.25        , -0.25*dt,  
                               &fields_inter[Nstart], Nsubsys);

        find_divF_divB_subsystem(fields_inter, divF);
        kernels::add_3_vectors(&eb.fields[Nstart], &fields_inter[Nstart], &divF[Nstart], 
                               one_third, two_thirds  , -two_thirds*dt,  
                               &eb.fields[Nstart], Nsubsys);
    }


    if (is_output_step || initial_cleaning)
    {
        find_divF_divB_subsystem(eb.fields, divF);
        kernels::multiply_by_scalar(&divF[psi], 1.0/Physics::ch_sq, eb.divB, eb.Ns_block);
    }

    if (initial_cleaning)
        kernels::multiply_by_scalar_inPlace(&eb.fields[psi], 0.0, eb.Ns_block);

    kernels::free(fields_inter);
    kernels::free(divF);

    if (initial_cleaning)
        write::message("Finished initial-field divergence cleaning");

    return;
}


/* For iterating the divergence-cleaning subsystem on its own */
void Process::find_divF_divB_subsystem(const real_t* const U, real_t* const divF)
{
    ElementBlock& eb = elements;
    const int psi = 8 * eb.Ns_block; // mem location at which the psi field begins

    /* These vectors are in "transform direction" space
     * dF is the only one that "needs" to be a VectorField, in
     * that the whole object is passed to a kernel */
    VectorField Uf; // Solution interpolated to flux points
    VectorField  F; // Fluxes
    VectorField dF; // Store d_j ( root_det_g * F(i)^j )
    
    for (int i: dirs)
    {
        Uf(i) = kernels::alloc_raw(Nfield * eb.Nf_dir_block[i]);
         F(i) = kernels::alloc_raw(Nfield * eb.Nf_dir_block[i]);
        dF(i) = kernels::alloc_raw(Nfield * eb.Ns_block);
    }

    
    for (int i: dirs)
        kernels::soln_to_flux(eb.soln2flux(i), U, Uf(i), eb.lengths, i);

    for (int i: dirs)
        divClean::bulk_fluxes(Uf(i), F(i), eb.geometry.S[i], eb.physics[i], eb.lengths, i);

    for (int i: ifaces)
        kernels::fill_face_data(Uf(faces[i].normal_dir), faces[i], eb.lengths);

    exchange_boundary_data();

    for (int i: ifaces)
    {
        int dir = faces[i].normal_dir;
        divClean::external_numerical_flux(faces[i], F(dir), F_numerical_divB_subsystem[dir],
                                             eb.geometry.S[dir], eb.lengths);
    }
    
    for (int i: dirs)
    {
        divClean::internal_numerical_flux(Uf(i), F(i), F_numerical_divB_subsystem[i],
                                         eb.geometry.S[i], eb.geometry.normal[i], eb.lengths, i);
    }

    for (int i: dirs)
        divClean::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), F(i), dF(i), eb.lengths, i);

    divClean::flux_divergence(dF, eb.geometry.Jrdetg(), divF, eb.lengths);

    /* Add the damping source term for the div-cleaning scalar field */
    kernels::add_scaled_vectors_inPlace(&divF[psi], &U[psi], 
                                        Physics::psi_damping_rate, eb.Ns_block);

    /* Add geometric source terms if not using Cartesian physical coordinates */
    if (eb.physics_soln->metric->physical_coords != cartesian)
        divClean::add_geometric_sources(divF, U, eb.physics_soln, eb.Ns_block);


    for (int i: dirs)
    {
        kernels::free(Uf(i));
        kernels::free( F(i));
        kernels::free(dF(i));
    }

    return;
}


