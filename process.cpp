#include "process.hpp"

#include <string>
#include <mpi.h>
#include <cmath>
#include <filesystem>
#include "params.hpp"
#include "write_screen.hpp"
#include "basic_time_integrator.hpp"
#include "boundary_conditions.hpp"
#include "physics_includes.hpp"
#include "kernels.hpp"

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

    write::message("Creating data directory");
    if (rank == 0)
        if (std::filesystem::exists("data") == false)
            if (std::filesystem::create_directory("data") == false)
                write::error("Failed to create data directory", destroy);
    
    write::message("Finished setting up process\n");
    
    MPI_Barrier(MPI_COMM_WORLD);

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

    Physics::ch = 0.0; // Find timestep without cleaning wave

    /* Should be able to convert this to a single search over data on the solution points? */
    real_t vmax_dir[3]; 
    real_t vmax;
    for (int d: dirs)
    {
        real_t* Uf = kernels::alloc_raw(Nfield*eb.Nf_dir_block[d]);
        kernels::soln_to_flux(eb.soln2flux(d), eb.fields, Uf, eb.lengths, d);
        dtmin_dir[d] = kernels::local_timestep(Uf, vmax_dir[d], eb.geometry.timestep_transform[d],
                                               eb.physics[d], eb.lengths, d);                          
        delete Uf;
    }
    dtmin = MIN(dtmin_dir[0], MIN(dtmin_dir[1], dtmin_dir[2]));
    vmax  = MAX(vmax_dir[0],  MAX(vmax_dir[1],  vmax_dir[2]));

    MPI_Allreduce(MPI_IN_PLACE, &dtmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dtmin_advect = cfl * dtmin;

    if (Physics::diffusive)
    {
        // For constant diffusive coefficients this can be moved to the setup phase 
        // --- same for all time steps
        // d_t_c = 1/3 stable for WaveRect at high viscosity, but 0.4 breaks sometimes...?
        // Seems to be excessively restrictive when have an inhomogeneous grid?
        // Very problem dependent: can use d_t_c = 0.9 for the MHD Alfven wave problem
        real_t diffusion = MAX(MAX(Physics::viscosity, Physics::resistivity), Physics::conductivity);
        diffusion = MAX(diffusion, 1e-15);
        dtmin_diff   = (Physics::diffusive_timestep_const/diffusion) * (1./(tt_max_global*tt_max_global));

        dt = MIN(dtmin_advect, dtmin_diff); 
        dt_ratio = dtmin_diff/dtmin_advect;
    }
    else
        dt = dtmin_advect;


    /* Calculate maximum stable div-cleaning wavespeed */
    if (system == mhd)
    {
        real_t lambda_max = 1.0 /(dt * tt_max_global); // Max allowable e'value on the grid

        real_t safety = 0.8; // 0.95;
        if (cfl < 0.4)
            safety = 0.75; // Seem to run into trouble at very high ch for small cfl

        /* For using the non-Galilean-invariant method of Derigs+ 2018 */
        Physics::ch = safety * std::sqrt(lambda_max * (lambda_max - vmax));
        //std::cout << "ch = " << Physics::ch << std::endl;

        Physics::psi_damping_rate = Physics::ch / Physics::psi_damping_const; //p_d_c = c_r from Dedner
        //std::cout << Physics::psi_damping_rate << std::endl;
    }
    /********************************************/

    if (rank == 0)
    {
        using namespace std;
        cout << scientific;
        const int defaultprec = cout.precision();
        cout.precision(5);

        if (Physics::diffusive)
            cout << "Starting time step " << step << " --- t = " << time
                 << " --- dt = " << dt 
                 << " --- dt ratio: " << dt_ratio
                 << " --- last output no.: " << data_output_counter-1 << endl;
        else
            cout << "Starting time step " << step << " --- t = " << time
                 << " --- dt = " << dt << " --- last output no.: " << data_output_counter-1 << endl;

        cout << defaultfloat;
        cout.precision(defaultprec);
    }

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
    
    /* Only used by diffusive terms. Needs to be here, rather than in add_diffusive_flux(), 
     * so can be passed to geometric sources. Also used by MHD for GLM sources -- grad(psi). */
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

    for (int i: ifaces)
        kernels::fill_face_data(Uf(faces[i].normal_dir), faces[i], eb.lengths);

    exchange_boundary_data();

    for (int i: ifaces)
    {
        if (faces[i].domain_external_face)
            kernels::dirichlet_boundary_conditions(faces[i]);

        int dir = faces[i].normal_dir;
        kernels::external_numerical_flux(faces[i], F(dir), F_numerical[dir],
                                         eb.geometry.S[dir], eb.lengths);
    }
    
    for (int i: dirs)
        kernels::internal_numerical_flux(Uf(i), F(i), F_numerical[i],
                                         eb.geometry.S[i], eb.geometry.normal[i], eb.lengths, i);

    /* For explicit diffusive terms, calculate the diffusive flux and add
     * to the advective fluxes before taking the flux deriv: F += F_diffusive */
    if (Physics::diffusive)
        add_diffusive_flux(Uf, dP, F);

    /* For testing that dP's are being calculated correctly -- put a known function *
     * in, say, pressure and dump out its derivative.                               */
    //for (int i=0; i < eb.Ns_block; i++)
    //    eb.divB[i] = dP(0,4*eb.Ns_block + i);
    //return;

    for (int i: dirs)
        kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), F(i), dF(i), eb.lengths, i);
    

    kernels::flux_divergence(dF, eb.geometry.Jrdetg(), divF, eb.lengths);

    /**** Add full GLM sources here, including Powell terms. Should move and reorganize.
     **** See Eqns 3.17 and 4.75 of Derigs+ 2018.                                    ***/
    if (system == mhd)
    {
        const int pstart = 8 * eb.Ns_block; // mem location at which the psi field begins

        /* Find divB from divF[psi] and save into elements.divB */
        //kernels::multiply_by_scalar(&divF[pstart], 1.0/Physics::ch, eb.divB, eb.Ns_block);
        
        /* Find divB independently, by averaging B components at the element interfaces.
         * This seems to work better. */
        find_divB(&U[5 * eb.Ns_block], eb.divB);

        real_t rho, vl[3], vu[3], Bu[3], Bl[3], vdotB, divB;
        real_t gradpsi_dot_v;
        const int N = eb.Ns_block;
        for (int i = 0; i < N; ++i)
        {
            divB = eb.divB[i];
            rho = U[i]; // 0th field
            for (int d: dirs)
                vl[d] = U[(1+d)*N + i] / rho;

            for (int d: dirs)
                Bu[d] = U[(5+d)*N + i];

            eb.physics_soln->metric->raise(vl, vu, i);
            eb.physics_soln->metric->lower(Bu, Bl, i);
            vdotB = eb.physics_soln->metric->dot(vu, Bu, i);

            /* Using Eulerian method of Derigs+ 2018, Sec 3.9 
               Lagrangian method seems unstable for discontinuous field loop test?
               Add the grad(psi).v terms to use the Lagrangian method */
            //gradpsi_dot_v = dP(0,pstart+i) * vu[0] + dP(1,pstart+i) * vu[1] + dP(2,pstart+i) * vu[2];

            // Momentum sources
            for (int d: dirs)
                divF[(1+d)*N + i] += divB * Bl[d];

            // Energy source
            divF[4*N + i] += divB * vdotB; // + gradpsi_dot_v * U[pstart + i];

            // Induction equation sources
            for (int d: dirs)
                divF[(5+d)*N +i] += divB * vu[d];

            // Psi equation sources
            divF[pstart + i] += Physics::psi_damping_rate * U[pstart + i]; // + gradpsi_dot_v;
        }
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


/* Separate function for finding divB --- can repurpose to find derivative of
 * any vector field. Element interfaces are joined using simple averaging
 * rather than e.g. the HLL flux.                                               */
void Process::find_divB(const real_t* const B, real_t* const divB)
{
    ElementBlock& eb = elements;

    const int hold_Nfield = eb.lengths.Nfield; // True total no. of fields

    /* These vectors are in "transform direction" space */
    VectorField Bf; // Physical B components interpolated to flux points
    VectorField B_ref; // Reference-space components in each transform direction
    VectorField dB; // Store d_j ( root_det_g * B(i)^j )
                    
    for (int i: dirs)
    {
        Bf(i)    = kernels::alloc_raw(3 * eb.Nf_dir_block[i]);
        B_ref(i) = kernels::alloc_raw(    eb.Nf_dir_block[i]);
        dB(i)    = kernels::alloc_raw(    eb.Ns_block);
    }

    eb.lengths.Nfield = 3; 

    for (int i: dirs)
        kernels::soln_to_flux(eb.soln2flux(i), B, Bf(i), eb.lengths, i);
    
    for (int i: ifaces)
    {
        kernels::fill_face_data(Bf(faces[i].normal_dir), faces[i], eb.lengths);

        if (faces[i].domain_external_face)
        {
            FaceCommunicator& f = faces[i];

            /* Just set neighbour data to existing data. This is better
             * than doing nothing */
            //for (int mem = 0; mem < f.Ntot * 3; ++mem)
            //    f.neighbour_data[mem] = f.my_data[mem];

            /* HACK: set normal flux to zero on the boundary - move to kernels */
            /***/
            real_t nl[3], nu[3];
            real_t B[3],  Bm[3];
            real_t Bdotn;
            for (int j = 0; j < f.Ntot; ++j) // for each location on the face...
            {
                for (int d: dirs)
                {
                    nl[d] = f.normal(d,j);
                    B[d]  = f.my_data[j + d * f.Ntot];
                }

                Bdotn = nl[0]*B[0] + nl[1]*B[1] + nl[2]*B[2];

                f.physics->metric->raise(nl, nu, j);
            
                for (int d: dirs)
                    Bm[d] = B[d] - Bdotn * nu[d]; // Zero normal flux
                
                for (int d: dirs)
                    f.neighbour_data[j + d * f.Ntot] = f.my_data[j + d * f.Ntot] = Bm[d];
            }
            /***/

            /* DOUBLE HACK: set normal B to 1 for the Hartmann flow test .... */
            /*** 
            for (int j = 0; j < f.Ntot; ++j) 
            {
                f.neighbour_data[j + 0 * f.Ntot] = f.my_data[j + 0 * f.Ntot] = 0.0;
                f.neighbour_data[j + 1 * f.Ntot] = f.my_data[j + 1 * f.Ntot] = 1.0;
                f.neighbour_data[j + 2 * f.Ntot] = f.my_data[j + 2 * f.Ntot] = 0.0;
            }
             ***/
        }
    }
    
    /* Exchanges more data than necessary since only have 3 "fields" here */
    exchange_boundary_data();

    /* Average Bf on process-external faces. */
    for (int i: ifaces)
        kernels::external_interface_average(faces[i], Bf(faces[i].normal_dir), eb.lengths);

    /* Average Bf on process-internal faces. */
    for (int i: dirs)
        kernels::internal_interface_average(Bf(i), eb.lengths, i);

    /* On each flux-point-block, this converts the three physical vector components into a 
     * single component of the vector, in reference space, along that flux-transform direction. */
    for (int i: dirs)
        kernels::phys_vector_to_ref_density(Bf(i), B_ref(i), eb.geometry.S[i], eb.lengths, i);
    
    eb.lengths.Nfield = 1; // Previous operation means we're effectively dealing with
                           // the flux of a single scalar field now.

    for (int i: dirs)
        kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), B_ref(i), dB(i), eb.lengths, i);
    
    kernels::flux_divergence(dB, eb.geometry.Jrdetg(), divB, eb.lengths);
                
    for (int i: dirs)
    {
        kernels::free(Bf(i));
        kernels::free(B_ref(i));
        kernels::free(dB(i));
    }

    eb.lengths.Nfield = hold_Nfield;

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
        kernels::external_interface_average(faces[i], Pf(faces[i].normal_dir), eb.lengths);

    /* Average Pf on process-internal faces. */
    for (int i: dirs)
        kernels::internal_interface_average(Pf(i), eb.lengths, i);

    for (int i: dirs) // i = derivative direction
    {
        // Pf(i): the P at the i-direction's flux points
        // dP_ref(i): gradients of the P in the i direction, at the solution points
        kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), Pf(i), dP_ref(i), eb.lengths, i);
    }

    /* Transform to physical-space gradient on the solution points */
    kernels::gradient_ref_to_phys(dP_ref, dP, eb.geometry.dxdr, eb.lengths);

    for (int dderiv: dirs)
        for (int dflux: dirs) // for each flux-point direction...
            kernels::soln_to_flux(eb.soln2flux(dflux), dP(dderiv), dPf[dflux](dderiv), 
                                                                     eb.lengths, dflux);
    
    /* For now, do the interface averaging separately for each *physical* derivative direction */
    for (int dderiv: dirs)
    {
        for (int i: ifaces)
            kernels::fill_face_data(dPf[faces[i].normal_dir](dderiv), faces[i], eb.lengths);

        exchange_boundary_data(); // Do 3 separate MPI calls, one for each physical derivative direction
                                  // This means we never have all three derivatives at the same time...

        for (int i: ifaces)
        {
            if (faces[i].domain_external_face)
                kernels::neumann_boundary_conditions(faces[i]);

            kernels::external_interface_average(faces[i], dPf[faces[i].normal_dir](dderiv), eb.lengths);
        }
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
        if (faces[i].domain_external_face)
        {
            kernels::fill_face_data(Uf(faces[i].normal_dir), faces[i], eb.lengths);
            //faces[i].BC->setup(faces[i].my_data);
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
