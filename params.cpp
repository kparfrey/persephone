#include "params.hpp"

#include <string>
#include "process.hpp"
#include "write_screen.hpp"
#include "basic_time_integrator.hpp"
#include "physics_includes.hpp"
#include "numerical_flux.hpp"


/* That part of Process setup which is the same for Cartesian,
 * Spherical, toroidal etc. configurations. */
void Params::setup_process_generic(Process &proc)
{
    proc.time     = 0.0;
    proc.end_time = end_time;
    proc.dt_write = dt_write;
    proc.step     = 0;
    proc.cfl      = cfl;
    proc.data_output_counter = 0;
    proc.is_output_step = false;

    switch (time_method)
    {
        case rk2_midpoint:
            proc.time_integrator = new RK2_midpoint; 
            break;
        case rk3_ssp:
            proc.time_integrator = new RK3_SSP; 
            break;
        default:
            write::error("Time integration method not recognized.");
    }

    switch (equations)
    {
        case scalar_advection:
            for (int d: dirs)
                proc.elements.physics[d] = new ScalarAdvection;
            proc.elements.physics_soln = new ScalarAdvection;
            //proc.system_data = new SystemData_scalar_advection;
            //proc.U_to_P   = new UtoP_scalar_advection;
            //proc.c_from_P = new WaveSpeeds_scalar_advection;
            //proc.F_from_P = new Fluxes_scalar_advection;
            break;
        case navier_stokes:
            for (int d: dirs)
                proc.elements.physics[d] = new NavierStokes;
            proc.elements.physics_soln = new NavierStokes;
            //proc.system_data = new SystemData_navstokes;
            //proc.U_to_P   = new UtoP_navstokes;
            //proc.c_from_P = new WaveSpeeds_navstokes;
            //proc.F_from_P = new Fluxes_navstokes;
            //proc.F_diff   = new DiffusiveFluxes_navstokes;
            break;
        case mhd:
            for (int d: dirs)
                proc.elements.physics[d] = new MHD;
            proc.elements.physics_soln = new MHD;
            //proc.system_data = new SystemData_mhd;
            //proc.U_to_P   = new UtoP_mhd;
            //proc.c_from_P = new WaveSpeeds_mhd;
            //proc.F_from_P = new Fluxes_mhd;
            //proc.F_diff   = new DiffusiveFluxes_mhd;
            break;
        default:
            write::error("Equation system not recognised.");
    }

    proc.system = equations;
    proc.Nfield = proc.elements.physics_soln->Nfield;

    for (int d: dirs)
    {
        /* Move inside a switch once more flux choices are defined */
        proc.F_numerical[d] = new HLL;

        proc.F_numerical[d]->Nfield  = proc.Nfield;
        proc.F_numerical[d]->physics = proc.elements.physics[d]; // Convenience pointer

        /* The above convenience pointer to proc.physics[d] is the only 
         * difference between the three copies of NumericalFlux */

        /* For running the divergence-cleaning subsystem on its own */
        proc.F_numerical_divB_subsystem[d] = new HLL_divB_subsystem;
        proc.F_numerical_divB_subsystem[d]->physics = proc.elements.physics[d]; // Convenience pointer
    }

    //proc.F_numerical->U_to_P   = proc.U_to_P;   
    //proc.F_numerical->c_from_P = proc.c_from_P;
    //proc.F_numerical->F_from_P = proc.F_from_P;

    return;
}
