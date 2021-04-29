#include "process.hpp"

#include <string>
#include <mpi.h>
#include <cmath>
#include "params.hpp"
#include "write_screen.hpp"
#include "basic_time_integrator.hpp"
#include "boundary_conditions.hpp"

//#include <stdlib.h>

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
    dt = cfl * dtmin;
    /********************************************/

    write::message("Starting time step " + std::to_string(step) + " --- t = " + std::to_string(time)
                                         + " --- dt = " + std::to_string(dt));

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

    for (int i: dirs)
        kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), F(i), dF(i), eb.lengths, i);

    kernels::flux_divergence(dF, eb.metric.Jrdetg(), divF, eb.lengths);

    /* Add div-cleaning scalar field's source directly here. Could alternatively do operator
     * splitting, and directly reduce psi in a subsequent substep. */
    if (system == mhd)
    {
        const real_t c_h = 0.9;        
        const real_t c_r = 1.0; // Dedner recommends 0.18
        const real_t c_p = std::sqrt(c_r * c_h);
        kernels::store_divB(divF, elements.divB, eb.lengths); // Only do this at output stage...
        kernels::scalar_field_source(divF, U, eb.lengths, c_h, c_p);
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
