#include "process.hpp"

#include <string>
#include <mpi.h>
#include <cmath>
#include "params.hpp"
#include "write_screen.hpp"
#include "basic_time_integrator.hpp"
#include "boundary_conditions.hpp"
#include "physics_includes.hpp"
#include "kernels.hpp"
#include "divclean.hpp"

//#include <stdio.h>

Process::Process(Params &params)
: params(params)
{
    // The constructor just sets the params reference to the active_parameters
    // object in active_params.hpp, passed in from main.cpp.
}


void Process::write_startup_info()
{
    if (rank == 0)
        params.write_param_info();

    return;
}


void Process::setup()
{
    write::message("Setting up process --- generic setup");
    params.setup_process_generic(*this); // in params.cpp
    
    write::message("Setting up process --- parameter-type-specific setup");
    params.setup_process(*this); // in params_cartesian.cpp etc.
    
    /* Calculates the global maximum of the timestep_transform array, for use
     * in calculating the maximum stable div-cleaning speed c_h in MHD */
    MPI_Allreduce(&elements.timestep_transform_max, &tt_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    write::message("Finished setting up process\n");
    return;
}


/* Move all of the process's device-side data there, to be called at
 * the end of the setup phase. */
void Process::move_to_device()
{
    write::message("End of host-side initialization phase --- copying all required data to device");
    
    elements.move_to_device();

    write::message("Finished copying data to device and freeing memory on host");

    return;
}


void Process::time_advance()
{
    /********************************************/
    /* Calculate timestep --- temporary */
    /* Quasi-1D approach, seems to be approximately right but should
     * be able to refine. Should check how much runtime this takes.
     * Can probably restrict to each element's boundary. */
    ElementBlock& eb = elements;
    real_t dtmin_dir[3];
    real_t dtmin;
    real_t dtmin_advect;
    real_t dtmin_diff, dt_ratio;

    for (int d: dirs)
    {
        real_t* Uf = kernels::alloc_raw(Nfield*eb.Nf_dir_block[d]);
        kernels::soln_to_flux(eb.soln2flux(d), eb.fields, Uf, eb.lengths, d);
        dtmin_dir[d] = kernels::local_timestep(Uf, eb.geometry.timestep_transform[d],
                                               eb.physics[d], eb.lengths, d);                          
                                               //U_to_P, c_from_P, eb.lengths, d);                          
        delete Uf;
    }
    dtmin = MIN(dtmin_dir[0], MIN(dtmin_dir[1], dtmin_dir[2]));
    MPI_Allreduce(MPI_IN_PLACE, &dtmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dtmin_advect = cfl * dtmin;

    if (Physics::diffusive)
    {
        // For constant diffusive coefficients this can be moved to the setup phase 
        // --- same for all time steps
        // d_t_c = 1/3 stable for WaveRect at high viscosity, but 0.4 breaks sometimes...?
        // Seems to be excessively restrictive when have an inhomogeneous grid?
        // Very problem dependent: can use d_t_c = 0.9 for the MHD Alfven wave problem
        const real_t diffusion = MAX(Physics::viscosity, Physics::resistivity);
        dtmin_diff   = (Physics::diffusive_timestep_const/diffusion) * (1./(tt_max_global*tt_max_global));

        dt = MIN(dtmin_advect, dtmin_diff); 
        dt_ratio = dtmin_diff/dtmin_advect;
    }
    else
        dt = dtmin_advect;


    /* Calculate maximum stable div-cleaning wavespeed */
    if (system == mhd)
    {
        const real_t ch_divClean = 1.0 /(dt * tt_max_global); //Should be stable with 1/()
        //std::cout << ch_divClean << std::endl;
        //const real_t ch_divClean = 2.5;

        Physics::ch    = ch_divClean;
        Physics::ch_sq = ch_divClean * ch_divClean;
        //F_from_P->ch_sq = physics->c_h * physics->c_h;

        Physics::psi_damping_rate = ch_divClean / Physics::psi_damping_const; //p_d_c = c_r from Dedner
        //Physics::psi_damping_rate = Physics::psi_damping_const / dt; 
        //physics->psi_damping_exp = std::exp(-cfl * physics->psi_damping_const);
    }
    /********************************************/

    if (Physics::diffusive)
        write::message("Starting time step " + std::to_string(step) + " --- t = " + std::to_string(time)
                                             + " --- dt = " + std::to_string(dt)
                                             + " --- dt ratio: " + std::to_string(dt_ratio)
                                             + " --- last output no.: " + std::to_string(data_output_counter-1));
    else
        write::message("Starting time step " + std::to_string(step) + " --- t = " + std::to_string(time)
                                             + " --- dt = " + std::to_string(dt)
                                             + " --- last output no.: " + std::to_string(data_output_counter-1));

    /* Call the fundamental time step method we're using */
    (*time_integrator)(*this); //Since storing a pointer to a BasicTimeIntegrator

    time += dt;
    step += 1;

    return;
}


void Process::divB_subsystem_iterations(const int Niter, const bool initial_cleaning)
{
    ElementBlock& eb = elements;

    /* ch is the only wavespeed in the divB subsystem, so can choose it to be anything.
     * Choose ch = 1.0 */
    //Physics::ch    = 1.0;
    //Physics::ch_sq = 1.0;
    //Physics::psi_damping_rate = 1.0 / Physics::psi_damping_const; //p_d_c = c_r from Dedner

    //const real_t dt_divB = 0.9 / tt_max_global; 

    const int Nstart  = 5 * eb.Ns_block; // Beginning of B0 data
    const int Nsubsys = 4 * eb.Ns_block; // 3 B componenents + Psi
    const int Ntot    = Nfield * eb.Ns_block; 
    const int psi     = 8 * eb.Ns_block; // mem location at which the psi field begins
    const real_t one_third  = 1.0/3.0;
    const real_t two_thirds = 2.0/3.0;

    real_t* fields_inter = kernels::alloc(Ntot);
    real_t* divF         = kernels::alloc(Ntot);

    /* Set Psi to 0 at beginning and end --- decouples timescales between here and main loop */
    //kernels::multiply_by_scalar_inPlace(&eb.fields[psi], 0.0, eb.Ns_block);

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

        if (initial_cleaning && (iter % 10 == 0))
            std::cout << iter << "... ";
    }


    if (is_output_step || initial_cleaning)
    {
        find_divF_divB_subsystem(eb.fields, divF);
        kernels::multiply_by_scalar(&divF[psi], 1.0/Physics::ch_sq, eb.divB, eb.Ns_block);
    }

    //kernels::multiply_by_scalar_inPlace(&eb.fields[psi], 0.0, eb.Ns_block);

    kernels::free(fields_inter);
    kernels::free(divF);

    if (initial_cleaning)
        write::message("\nFinished initial-field divergence cleaning");

    return;
}


/* The fundamental function of the spectral difference method */
/* U : vector of conserved solution fields (physical space) 
 * F : vector of fluxes (reference-element space) */
void Process::find_divF(const real_t* const U, const real_t t, real_t* const divF)
{
    ElementBlock& eb = elements;

    /* These vectors are in "transform direction" space
     * dF is the only one that "needs" to be a VectorField, in
     * that the whole object is passed to a kernel */
    VectorField Uf; // Solution interpolated to flux points
    VectorField  F; // Fluxes
    VectorField dF; // Store d_j ( root_det_g * F(i)^j )
    
    /* Only used by diffusive terms. Needs to be here, rather than in add_diffusive_flux(), 
     * so can be passed to geometric sources */
    VectorField dP; // Physical-space derivatives of primitives at the solution points
    
    for (int i: dirs)
    {
        Uf(i) = kernels::alloc_raw(Nfield * eb.Nf_dir_block[i]);
         F(i) = kernels::alloc_raw(Nfield * eb.Nf_dir_block[i]);
        dF(i) = kernels::alloc_raw(Nfield * eb.Ns_block);
        
        if (Physics::diffusive)
            dP(i) = kernels::alloc_raw(Nfield * eb.Ns_block);
    }

    
    for (int i: dirs)
        kernels::soln_to_flux(eb.soln2flux(i), U, Uf(i), eb.lengths, i);

    for (int i: dirs)
        kernels::bulk_fluxes(Uf(i), F(i), eb.geometry.S[i], eb.physics[i], eb.lengths, i);
                             //U_to_P, F_from_P, eb.lengths, i);

    for (int i: ifaces)
        kernels::fill_face_data(Uf(faces[i].normal_dir), faces[i], eb.lengths);

    exchange_boundary_data();

    for (int i: ifaces)
    {
        int dir = faces[i].normal_dir;
        kernels::external_numerical_flux(faces[i], F(dir), F_numerical[dir],
                                         eb.geometry.S[dir], eb.lengths);
    }
    
    for (int i: dirs)
    {
        kernels::internal_numerical_flux(Uf(i), F(i), F_numerical[i],
                                         eb.geometry.S[i], eb.geometry.normal[i], eb.lengths, i);
    }

    /* For explicit diffusive terms, calculate the diffusive flux and add
     * to the advective fluxes before taking the flux deriv: F += F_diffusive */
    if (Physics::diffusive)
        add_diffusive_flux(Uf, dP, F);

    for (int i: dirs)
        kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), F(i), dF(i), eb.lengths, i);

    kernels::flux_divergence(dF, eb.geometry.Jrdetg(), divF, eb.lengths);

    /* Add div-cleaning scalar field's source directly here. Could alternatively do operator
     * splitting, and directly reduce psi in a subsequent substep. */
    if (system == mhd)
    {
        const int psi = 8 * eb.Ns_block; // mem location at which the psi field begins

        if (is_output_step && (substep == 1))
        {
            /* Find divB from divF[psi] and save into elements.divB */
            const real_t over_chsq = 1.0/Physics::ch_sq;
            kernels::multiply_by_scalar(&divF[psi], over_chsq, eb.divB, eb.Ns_block);
        }

        /* Add the damping source term for the div-cleaning scalar field */
        kernels::add_scaled_vectors_inPlace(&divF[psi], &U[psi], 
                                            Physics::psi_damping_rate, eb.Ns_block);
    }


    /* Add geometric source terms if not using Cartesian physical coordinates */
    if (eb.physics_soln->metric->physical_coords != cartesian)
        kernels::add_geometric_sources(divF, U, dP, eb.physics_soln, Nfield, eb.Ns_block);


    for (int i: dirs)
    {
        kernels::free(Uf(i));
        kernels::free( F(i));
        kernels::free(dF(i));
        
        if (Physics::diffusive)
            kernels::free(dP(i));
    }

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
        divClean::soln_to_flux(eb.soln2flux(i), U, Uf(i), eb.lengths, i);

    for (int i: dirs)
        divClean::bulk_fluxes(Uf(i), F(i), eb.geometry.S[i], eb.physics[i], eb.lengths, i);

    for (int i: ifaces)
        divClean::fill_face_data(Uf(faces[i].normal_dir), faces[i], eb.lengths);

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


/* Finds the diffusive flux, adds in place to the advective flux in F */
void Process::add_diffusive_flux(VectorField Uf, VectorField dP, VectorField F)
{
    VectorField Fd;     // Diffusive fluxes for all fields
    VectorField dP_ref; // Ref-space derivatives of P at the solution points
    VectorField dPf[3]; // Physical-space derivatives of P at the flux points
                        // Note: dPf[flux-point dir](deriv dir, ...)


    ElementBlock& eb = elements;

    for (int i: dirs)
    {
        /* Initialise Fd to zero - otherwise have problems if forget to set 
         * a component in the DiffusiveFluxes functor. Doesn't seem to work?? */
        Fd(i) = kernels::alloc(Nfield * eb.Nf_dir_block[i]);
        //dU(i) = kernels::alloc_raw(Nfield * eb.Ns_block);
        dP_ref(i) = kernels::alloc_raw(Nfield * eb.Ns_block);

        for (int j: dirs)
            dPf[i](j) = kernels::alloc_raw(Nfield * eb.Nf_dir_block[i]);
    }

    /* Convert all conserved variables to primitive variables, and save back in the same memory */
    for (int i: dirs)
        kernels::conserved_to_primitive_fluxpoints(Uf(i), eb.physics[i], eb.lengths, i);
    VectorField& Pf = Uf; // Reference, so don't forget arrays now contain primitives

    /* conserved -> primitive for data already exchanged & stored in faces */
    for (int i: ifaces)
        kernels::conserved_to_primitive_faces(faces[i], eb.physics[faces[i].normal_dir], eb.lengths);

    /* Average Pf on process-external faces. */
    for (int i: ifaces)
        kernels::external_interface_average(faces[i], Pf(faces[i].normal_dir), eb.lengths, false);

    /* Average Uf on process-internal faces. */
    for (int i: dirs)
        kernels::internal_interface_average(Pf(i), eb.lengths, i);

    for (int i: dirs) // i = derivative direction
    {
        // Pf(i): the P at the i-direction's flux points
        // dP_ref(i): gradients of the P in the i direction, at the solution points
        kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), Pf(i), dP_ref(i), eb.lengths, i);
    }

    /* Transform to physical-space gradient */
    kernels::gradient_ref_to_phys(dP_ref, dP, eb.geometry.dxdr, eb.lengths);

    for (int dderiv: dirs)
        for (int dflux: dirs) // for each flux-point direction...
            kernels::soln_to_flux(eb.soln2flux(dflux), dP(dderiv), dPf[dflux](dderiv), 
                                                                     eb.lengths, dflux);
    
    /* For now, do the interface averaging separately for each physical derivative direction */
    for (int dderiv: dirs)
    {
        for (int i: ifaces)
            kernels::fill_face_data(dPf[faces[i].normal_dir](dderiv), faces[i], eb.lengths);

        exchange_boundary_data(); // So have to do 3 separate MPI calls...

        for (int i: ifaces)
            kernels::external_interface_average(faces[i], dPf[faces[i].normal_dir](dderiv), eb.lengths, true);
    }
        
    for (int dflux: dirs)
        for (int dderiv: dirs)
            kernels::internal_interface_average(dPf[dflux](dderiv), eb.lengths, dflux);

    //const real_t coeffs[2] = {physics->viscosity, physics->resistivity};
    for (int i: dirs) // flux-point direction
        kernels::diffusive_flux(Pf(i), dPf[i], Fd(i), eb.physics[i], eb.geometry.S[i], eb.lengths, i);
        //kernels::diffusive_flux(Uf(i), dPf[i], Fd(i), eb.physics[i], eb.geometry.S[i], eb.lengths, i);
        //kernels::diffusive_flux(Uf(i), dUf[i], Fd(i), eb.physics[i], coeffs, eb.geometry.S[i], eb.lengths, i);
        //kernels::diffusive_flux(Uf(i), dUf[i], Fd(i), F_diff, coeffs, eb.geometry.S[i], eb.lengths, i);

    for (int i: dirs)
        kernels::add_vectors_inPlace(F(i), Fd(i), Nfield*eb.Nf_dir_block[i]);

    for (int i: dirs)
    {
        kernels::free(Fd(i));
        kernels::free(dP_ref(i));

        for (int j: dirs)
            kernels::free(dPf[i](j));
    }

    return;
}


/* Fill the stored_data array of external faces with the initial conditions 
 * extrapolated to the face. Needs to use soln_to_flux(), so everything 
 * should live on the device. */
void Process::fill_external_boundary_data()
{
    write::message("Filling initial boundary data on external faces");

    ElementBlock& eb = elements;
    VectorField Uf; // Solution interpolated to flux points

    for (int i: dirs)
    {
        Uf(i) = kernels::alloc_raw(Nfield * eb.Nf_dir_block[i]);
        kernels::soln_to_flux(eb.soln2flux(i), eb.fields, Uf(i), eb.lengths, i);
    }

    for (int i: ifaces)
        if (faces[i].external_face == true)
        {
            kernels::fill_face_data(Uf(faces[i].normal_dir), faces[i], eb.lengths);
            faces[i].BC->setup(faces[i].my_data);
        }

    for (int i: dirs)
        kernels::free(Uf(i));

    return;
}


void Process::exchange_boundary_data()
{
    MPI_Request requests[12];

    for (int i: ifaces)
        requests[i]   = faces[i].send_data();

    for (int i: ifaces)
        requests[6+i] = faces[i].receive_data();

    /* I guess I only really need to wait until the receive requests
     * are fulfilled --- could make sure all the sends are completed
     * much later? */
    MPI_Waitall(12, requests, MPI_STATUSES_IGNORE);

    return;
}
