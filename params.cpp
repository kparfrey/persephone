#include <string>
#include "process.hpp"
#include "write_screen.hpp"
#include "basic_time_integrator.hpp"
#include "physics_includes.hpp"
#include "numerical_flux.hpp"

/* Note: this file is #included by params.hpp because Params is a template class. 
 * No indepedent object file is compiled. */

/* That part of Process setup which is the same for Cartesian,
 * Spherical, toroidal etc. configurations. */
template <class ParamsType>
template <class ProcType>
void Params<ParamsType>::setup_process_generic(ProcType& proc)
{
    proc.time     = 0.0;
    proc.end_time = end_time;
    proc.dt_write = dt_write;
    proc.step     = 0;
    proc.cfl      = cfl;
    proc.data_output_counter = 0;
    proc.is_output_step = false;

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
