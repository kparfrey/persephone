#include "process.hpp"

#include <string>
#include <mpi.h>
#include "params.hpp"
#include "write_screen.hpp"
#include "basic_time_integrator.hpp"

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
    write::message("Starting time step " + std::to_string(step) + " --- t = " + std::to_string(time));

    /* Call the fundamental time step method we're using */
    (*time_integrator)(*this); //Since storing a pointer to a BasicTimeIntegrator

    time += dt;
    step += 1;

    return;
}




/* The fundamental function of the spectral difference method */
/* U : vector of conserved solution fields (physical space) 
 * F : vector of fluxes (reference-element space) */
void Process::find_RHS(real_t* U, real_t t, real_t* result)
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
    
    /*  Simple: assume Nfield = 1 etc */
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
        kernels::external_face_numerical_flux(faces[i], F(dir), F_numerical,
                                              eb.metric.S[dir], eb.lengths);
    }
    
    for (int i: dirs)
        kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), F(i), dF(i), eb.lengths, i);

    kernels::flux_divergence(dF, eb.metric.Jrdetg(), result, eb.lengths);

    for (int i: dirs)
    {
        kernels::free(Uf(i));
        kernels::free( F(i));
        kernels::free(dF(i));
    }

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
