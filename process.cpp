#include "process.hpp"

#include <string>
#include <mpi.h>
#include <cmath>
#include "params.hpp"
#include "write_screen.hpp"
#include "basic_time_integrator.hpp"
#include "boundary_conditions.hpp"
#include "physics_includes.hpp"

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
    MPI_Allreduce(&elements.l_min, &l_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    
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
        dtmin_dir[d] = kernels::local_timestep(Uf, eb.metric.timestep_transform[d],
                                               U_to_P, c_from_P, eb.lengths, d);                          
        delete Uf;
    }
    dtmin = MIN(dtmin_dir[0], MIN(dtmin_dir[1], dtmin_dir[2]));
    //dtmin = 1.0/(1/dtmin_dir[0] + 1/dtmin_dir[1] + 1/dtmin_dir[2]);
    MPI_Allreduce(MPI_IN_PLACE, &dtmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    // This factor of 0.25 seems to be the stability limit -- 0.27 breaks for the WaveRect
    // Can use up to 0.33 for a rectilinear grid
    // Seems to be excessively restrictive when have an inhomogeneous grid?
    dtmin_diff = 0.15 * l_min_global*l_min_global / system_data->viscosity; 
    //dtmin_diff   = 0.25 * (1/system_data->viscosity) * (1./(tt_max_global*tt_max_global));

    dtmin_advect = cfl * dtmin;

    dt = MIN(dtmin_advect, dtmin_diff); 
    dt_ratio = dtmin_diff/dtmin_advect;

    /* Calculate maximum stable div-cleaning wavespeed */
    if (system == mhd)
    {
        system_data->c_h = cfl /(dt * tt_max_global);
        system_data->psi_damping_rate = cfl * system_data->psi_damping_const / dt;
        //system_data->psi_damping_exp = std::exp(-cfl * system_data->psi_damping_const);
    }
    /********************************************/

    write::message("Starting time step " + std::to_string(step) + " --- t = " + std::to_string(time)
                                         + " --- dt = " + std::to_string(dt)
                                         + " --- dt ratio: " + std::to_string(dt_ratio));

    /* Call the fundamental time step method we're using */
    (*time_integrator)(*this); //Since storing a pointer to a BasicTimeIntegrator

    time += dt;
    step += 1;

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
    
    for (int i: dirs)
    {
        Uf(i) = kernels::alloc_raw(Nfield * eb.Nf_dir_block[i]);
         F(i) = kernels::alloc_raw(Nfield * eb.Nf_dir_block[i]);
        dF(i) = kernels::alloc_raw(Nfield * eb.Ns_block);
    }
    
    for (int i: dirs)
        kernels::soln_to_flux(eb.soln2flux(i), U, Uf(i), eb.lengths, i);

    for (int i: dirs)
        kernels::bulk_fluxes(Uf(i), F(i), eb.metric.S[i], 
                             U_to_P, F_from_P, eb.lengths, i);

    for (int i: ifaces)
        kernels::fill_face_data(Uf(faces[i].normal_dir), faces[i], eb.lengths);

    exchange_boundary_data();

    for (int i: ifaces)
    {
        int dir = faces[i].normal_dir;
        kernels::external_numerical_flux(faces[i], F(dir), F_numerical,
                                         eb.metric.S[dir], eb.lengths);
    }
    
    for (int i: dirs)
        kernels::internal_numerical_flux(Uf(i), F(i), F_numerical,
                                         eb.metric.S[i], eb.metric.normal[i], eb.lengths, i);



    /* For explicit diffusive terms, include calculation of the diffusive flux here
     * and add to the advective fluxes before taking the flux deriv: F += F_diffusive */
    if ((system == navier_stokes) && system_data->diffusive)
    {
        VectorField Fd;     // Diffusive fluxes for all fields
        VectorField dU_ref; // Ref-space derivatives of U at the solution points
        VectorField dU;     // Physical-space derivatives of U at the solution points
        VectorField dUf[3]; // Physical-space derivatives of U at the flux points
                            // Note: dUf[flux-point dir](deriv dir, ...)

        for (int i: dirs)
        {
            Fd(i) = kernels::alloc_raw(Nfield * eb.Nf_dir_block[i]);
            dU(i) = kernels::alloc_raw(Nfield * eb.Ns_block);
            dU_ref(i) = kernels::alloc_raw(Nfield * eb.Ns_block);

            for (int j: dirs)
                dUf[i](j) = kernels::alloc_raw(Nfield * eb.Nf_dir_block[i]);
        }

        /* Average Uf on process-external faces. 
         * Data already in FaceCommunicators from the earlier exchange */
        for (int i: ifaces)
            kernels::external_interface_average(faces[i], Uf(faces[i].normal_dir), eb.lengths);

        /* Average Uf on process-internal faces. */
        for (int i: dirs)
            kernels::internal_interface_average(Uf(i), eb.lengths, i);

        for (int i: dirs) // i = derivative direction
        {
            // Uf(i): the U at the i-direction's flux points
            // dU_ref(i): gradients of the U in the i direction, at the solution points
            kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), Uf(i), dU_ref(i), eb.lengths, i);
        }

        /* Transform to physical-space gradient */
        kernels::gradient_ref_to_phys(dU_ref, dU, eb.metric.dxdr, eb.lengths);

        for (int dderiv: dirs)
            for (int dflux: dirs) // for each flux-point direction...
                kernels::soln_to_flux(eb.soln2flux(dflux), dU(dderiv), dUf[dflux](dderiv), 
                                                                         eb.lengths, dflux);
        
        /* For now, do the interface averaging separately for each physical derivative direction */
        for (int dderiv: dirs)
        {
            for (int i: ifaces)
                kernels::fill_face_data(dUf[faces[i].normal_dir](dderiv), faces[i], eb.lengths);

            exchange_boundary_data(); // So have to do 3 separate MPI calls...

            for (int i: ifaces)
                kernels::external_interface_average(faces[i], dUf[faces[i].normal_dir](dderiv), eb.lengths);
        }
            
        for (int dflux: dirs)
            for (int dderiv: dirs)
                kernels::internal_interface_average(dUf[dflux](dderiv), eb.lengths, dflux);

        const real_t coeffs[1] = {system_data->viscosity};
        for (int i: dirs) // flux-point direction
            kernels::diffusive_flux(Uf(i), dUf[i], Fd(i), F_diff, coeffs, eb.metric.S[i], eb.lengths, i);

        for (int i: dirs)
            kernels::add_vectors_in_place(F(i), Fd(i), Nfield*eb.Nf_dir_block[i]);

        for (int i: dirs)
        {
            kernels::free(Fd(i));
            kernels::free(dU(i));
            kernels::free(dU_ref(i));

            for (int j: dirs)
                kernels::free(dUf[i](j));
        }
    } // end of diffusive section




    for (int i: dirs)
        kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), F(i), dF(i), eb.lengths, i);

    kernels::flux_divergence(dF, eb.metric.Jrdetg(), divF, eb.lengths);

    /* Add div-cleaning scalar field's source directly here. Could alternatively do operator
     * splitting, and directly reduce psi in a subsequent substep. */
    if (system == mhd)
    {
        if (is_output_step && (substep == 1))
            kernels::store_divB(divF, elements.divB, eb.lengths);

        kernels::scalar_field_source(divF, U, eb.lengths, system_data->c_h, system_data->psi_damping_rate);
    }

    for (int i: dirs)
    {
        kernels::free(Uf(i));
        kernels::free( F(i));
        kernels::free(dF(i));
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
