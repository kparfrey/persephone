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


/*** TEMPORARY ***/
#if 0
static void singleElement_periodicBC(real_t* F, LengthBucket lb, int dir)
{
    int dir1 = dir_plus_one[dir]; 
    int dir2 = dir_plus_two[dir];

    int Nf0 = lb.Nf[dir];
    int Ns1 = lb.Ns[dir1];

    int L, R;
    real_t fnum;

    for(int n2 = 0; n2 < lb.Ns[dir2]; ++n2)
    for(int n1 = 0; n1 < lb.Ns[dir1]; ++n1)
    {
        L = (n2 * Ns1 + n1) * Nf0 + 0;
        R = (n2 * Ns1 + n1) * Nf0 + Nf0-1;

        /* Central flux */
        fnum = 0.5 * (F[L] + F[R]);

        /* Upwind */
        //fnum = F[R]; // For wave velocity in +ve 0-direction
        //             // Does seem more stable/dissipative

        F[L] = fnum;
        F[R] = fnum;
        //F[R] = fnum + (rand() / (RAND_MAX + 1.))*1e-2;
    }
    
    return;
}
#endif




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

#if 1
    for (int i: ifaces)
    {
        int dir = faces[i].normal_dir;
        kernels::external_face_numerical_flux(faces[i], F(dir), 
                                              eb.metric.S[dir], eb.lengths);
    }
#endif
        
    //for (int i: dirs)
    //    singleElement_periodicBC(F(i), eb.lengths, i);
    
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
#if 1
    MPI_Request requests[12];

    for (int i: ifaces)
        requests[i]   = faces[i].send_data();

    for (int i: ifaces)
        requests[6+i] = faces[i].receive_data();

    /* I guess I only really need to wait until the receive requests
     * are fulfilled --- could make sure all the sends are completed
     * much later? */
    MPI_Waitall(12, requests, MPI_STATUSES_IGNORE);
#endif

    /**
    for (int i: ifaces)
    {
        FaceCommunicator& face = faces[i];
        MPI_Sendrecv(face.my_data, face.N_tot, MPI_real_t, face.neighbour_rank, 
                     i,
                     face.neighbour_data, face.N_tot, MPI_real_t, face.neighbour_rank,
                     i, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
    }
    **/


    return;
}
