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



/* The fundamental function of the spectral difference method */
/* U : vector of conserved solution fields (physical space) 
 * F : vector of fluxes (reference-element space) 
 * Note: here divF is just the negative of the "right hand side", for those
 *       fields not evolved by a flux divergence. */
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

    /************************************************** 
     **** Find B from vector potential and save into U.
     **** This is very hacky for now!                 */
#if 0
    eb.lengths.operation_mode = vecpot_mode;
    const int A0 = 8;
    const int A0_mem = A0 * eb.Ns_block;

    VectorField Af; // Vec. pot. interpolated to flux points
    for (int i: dirs)
        Af(i) = kernels::alloc_raw(3 * eb.Nf_dir_block[i]);

    for (int i: dirs)
        kernels::soln_to_flux(eb.soln2flux(i), &U[A0_mem], Af(i), eb.lengths, i);

    // Problem: will also need hydro data on faces to do the upwinding in 
    // the Riemann flux....
    for (int i: ifaces)
        kernels::fill_face_data(Af(faces[i].normal_dir), faces[i], eb.lengths);

    for (int i: dirs)
        kernels::free(Af(i));

    eb.lengths.operation_mode = normal_mode;
#endif
    /**************************************************/
    
    for (int i: dirs)
        kernels::soln_to_flux(eb.soln2flux(i), U, Uf(i), eb.lengths, i);

    for (int i: dirs)
        kernels::bulk_fluxes(Uf(i), F(i), eb.geometry.S[i], eb.physics[i], eb.lengths, i);

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
        kernels::internal_numerical_flux(Uf(i), F(i), F_numerical[i],
                                         eb.geometry.S[i], eb.geometry.normal[i], eb.lengths, i);

    /* For explicit diffusive terms, calculate the diffusive flux and add
     * to the advective fluxes before taking the flux deriv: F += F_diffusive */
    if (Physics::diffusive)
        add_diffusive_flux(Uf, dP, F);

    for (int i: dirs)
        kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), F(i), dF(i), eb.lengths, i);

    kernels::flux_divergence(dF, eb.geometry.Jrdetg(), divF, eb.lengths);

    /*
     * Dump out the initial divB
    if (rank == 0)
    {
        real_t Bsq, Bu[3];
        for (int i = 0; i < eb.Ns_block; ++i)
        {
            for (int d: dirs)
                Bu[d] = U[(5+d)*eb.Ns_block + i];
            Bsq = eb.physics_soln->metric->square(Bu, i);

            std::cout << std::abs(eb.divB[i])/std::sqrt(Bsq) << std::endl;
        }

        //exit(33);
    }
    */

    /* Add geometric source terms if not using Cartesian physical coordinates */
    if (eb.physics_soln->metric->physical_coords != cartesian)
        kernels::add_geometric_sources(divF, U, dP, eb.physics_soln, Nfield, eb.Ns_block);


    /* VECTOR POTENTIAL STUFF */
    if(true)
    {
        eb.lengths.Ncons = 3;
        int mem;

        /* Here the F and dF both just store the electric field --- no derivative being taken */
        for (int i: dirs)
            kernels::fluxDeriv_to_soln(eb.flux2soln(i), &F(i,8*eb.Nf_dir_block[i]), &dF(i,8*eb.Ns_block), eb.lengths, i);
            
        for (int f = 8; f < 11; ++f)
        for (int i = 0; i < eb.Ns_block; ++i)
        {
            mem = f*eb.Ns_block + i;
            divF[mem] = (dF(0,mem) + dF(1,mem) + dF(2,mem))/3.0;
        }
        
        eb.lengths.Ncons = 8;
    }

    /* For now just try transporting A using the B already at solution points */
    if (false)
    {
        real_t rdetg;
        const int Ns = eb.Ns_block;
        real_t* Up   = new real_t [Nfield];
        real_t* Pp   = new real_t [Nfield];
        real_t* v;
        real_t* B;
        
        for (int mem = 0; mem < Ns; ++mem)
        {
            /* Load all field variables at this location into the Up array. */
            for (int field = 0; field < Nfield; ++field)
                Up[field] = U[mem + field * Ns];

            eb.physics_soln->ConservedToPrimitive(Up, Pp, mem);
            v = &Pp[1];
            B = &Pp[5];

            rdetg = eb.physics_soln->metric->rdetg[mem];
    
            divF[mem + 8*Ns] = rdetg * (v[2]*B[1] - v[1]*B[2]);
            divF[mem + 9*Ns] = rdetg * (v[0]*B[2] - v[2]*B[0]);
            divF[mem +10*Ns] = rdetg * (v[1]*B[0] - v[0]*B[1]);

            if (Physics::diffusive)
            {
                // Just implemented for Cartesian geometry for now
                divF[mem + 8*Ns] += Physics::resistivity * (dP(1,mem+7*Ns) - dP(2,mem+6*Ns));
                divF[mem + 9*Ns] += Physics::resistivity * (dP(2,mem+5*Ns) - dP(0,mem+7*Ns));
                divF[mem +10*Ns] += Physics::resistivity * (dP(0,mem+6*Ns) - dP(1,mem+5*Ns));
            }
        }
        
        delete[] Up;
        delete[] Pp;
    }


    /* Temp: Powell source terms */
    if (false)
    {
        real_t rho, vl[3], vu[3], Bu[3], Bl[3], vdotB, divB;
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

            // Momentum sources
            for (int d: dirs)
                divF[(1+d)*N + i] += divB * Bl[d];

            // Energy source
            divF[4*N + i] += divB * vdotB;

            // Induction equation sources
            //for (int d: dirs)
            //    divF[(5+d)*N +i] += divB * vu[d];
        }
    }


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

    /* Average Pf on process-internal faces. */
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

    for (int i: dirs) // flux-point direction
        kernels::diffusive_flux(Pf(i), dPf[i], Fd(i), eb.physics[i], eb.geometry.S[i], eb.lengths, i);

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


void Process::replace_B(real_t* const U)
{
    ElementBlock& eb = elements;

    eb.lengths.Nfield = 3;
    eb.lengths.Ncons = 3;
    for (int i: ifaces)
        faces[i].Ntot_all = 3 * faces[i].Ntot;

    const int Ns = eb.Ns_block;
    const int A0 = 8;
    const int A0_mem = A0 * Ns;
    
    VectorField Af;
    VectorField dA;
    VectorField dA_ref;

    for (int d: dirs)
    {
        Af(d) = kernels::alloc(3 * eb.Nf_dir_block[d]);
        dA(d) = kernels::alloc(3 * eb.Ns_block);
        dA_ref(d) = kernels::alloc(3 * eb.Ns_block);
    }

    for (int i: dirs)
        kernels::soln_to_flux(eb.soln2flux(i), &U[A0_mem], Af(i), eb.lengths, i);

    for (int i: ifaces)
        kernels::fill_face_data(Af(faces[i].normal_dir), faces[i], eb.lengths);

    exchange_boundary_data();
    
    for (int i: ifaces)
        kernels::external_interface_average(faces[i], Af(faces[i].normal_dir), eb.lengths, false);
    
    for (int i: dirs)
        kernels::internal_interface_average(Af(i), eb.lengths, i);
    
    for (int i: dirs) // i = derivative direction
    {
        // Af(i): A at the i-direction's flux points
        // dA_ref(i): gradients of A in the i direction, at the solution points
        kernels::fluxDeriv_to_soln(eb.fluxDeriv2soln(i), Af(i), dA_ref(i), eb.lengths, i);
    }

    /* Transform to physical-space gradient */
    kernels::gradient_ref_to_phys(dA_ref, dA, eb.geometry.dxdr, eb.lengths);

    /* Use physical dA to find curl(A), save into U: B = Binit + curl(A) */
    real_t rdetg;
    int c0, c1, c2; // indices for vector components
    for (int mem = 0; mem < Ns; ++mem)
    {
        rdetg = eb.physics_soln->metric->rdetg[mem];
        c0 = mem;
        c1 = mem + Ns;
        c2 = mem + 2*Ns;

        U[mem + 5*Ns] = eb.Binit[c0] + (dA(1,c2) - dA(2,c1)) / rdetg;
        U[mem + 6*Ns] = eb.Binit[c1] + (dA(2,c0) - dA(0,c2)) / rdetg;
        U[mem + 7*Ns] = eb.Binit[c2] + (dA(0,c1) - dA(1,c0)) / rdetg;
    }


    for (int d: dirs)
    {
        kernels::free(dA(d));
        kernels::free(dA_ref(d));
        kernels::free(Af(d));
    }
    
    eb.lengths.Nfield = 11;
    eb.lengths.Ncons = 8;
    for (int i: ifaces)
        faces[i].Ntot_all = 11 * faces[i].Ntot;

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
