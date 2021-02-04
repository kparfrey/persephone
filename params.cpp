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

    switch (time_method)
    {
        case rk2_midpoint:
            proc.time_integrator = new RK2_midpoint; 
            break;
        default:
            write::error("Time integration method not recognized.");
    }

    switch (equations)
    {
        case scalar_advection:
            proc.system_data = new SystemData_scalar_advection;
            proc.U_to_P   = new UtoP_scalar_advection;
            proc.c_from_P = new WaveSpeeds_scalar_advection;
            proc.F_from_P = new Fluxes_scalar_advection;
            break;
        case euler:
            proc.system_data = new SystemData_euler;
            proc.U_to_P   = new UtoP_euler;
            proc.c_from_P = new WaveSpeeds_euler;
            proc.F_from_P = new Fluxes_euler;
            break;
        default:
            write::error("Equation system not recognised.");
    }

    proc.Nfield   = proc.system_data->Nfield;

    /* Move inside a switch once more flux choices are defined */
    if (geometry == simple_geometry) // Need to change, since broadening definition of simple_geometry
        proc.F_numerical = new HLL_straight;
    else
        proc.F_numerical = new HLL;

    proc.F_numerical->Nfield   = proc.Nfield;
    proc.F_numerical->U_to_P   = proc.U_to_P;   // Convenience pointers
    proc.F_numerical->c_from_P = proc.c_from_P;
    proc.F_numerical->F_from_P = proc.F_from_P;

    return;
}
