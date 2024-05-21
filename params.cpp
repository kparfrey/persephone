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
            proc.elements.physics_soln   = new ScalarAdvection;
            for (int f: ifaces)
                proc.faces[f].physics    = new ScalarAdvection;
            break;
        case navier_stokes:
            for (int d: dirs)
                proc.elements.physics[d] = new NavierStokes;
            proc.elements.physics_soln   = new NavierStokes;
            for (int f: ifaces)
                proc.faces[f].physics    = new NavierStokes;
            break;
        case mhd:
            for (int d: dirs)
                proc.elements.physics[d] = new MHD;
            proc.elements.physics_soln   = new MHD;
            for (int f: ifaces)
                proc.faces[f].physics    = new MHD;
            break;
        default:
            write::error("Equation system not recognised.");
    }

    proc.system = equations;
    proc.Nfield = proc.elements.physics_soln->Nfield;

    for (int d: dirs)
    {
        /* Move inside a switch once more flux choices are defined */
        //proc.F_numerical[d] = new HLL; // This is now taken care of by the templates

        proc.F_numerical[d].Nfield  = proc.Nfield;
        proc.F_numerical[d].physics = proc.elements.physics[d]; // Convenience pointer

        /* The above convenience pointer to proc.physics[d] is the only 
         * difference between the three copies of NumericalFlux */

        /* For running the divergence-cleaning subsystem on its own */
        //proc.F_numerical_divB_subsystem[d] = new HLL_divB_subsystem;
        //proc.F_numerical_divB_subsystem[d]->physics = proc.elements.physics[d]; // Convenience pointer
    }

    return;
}
