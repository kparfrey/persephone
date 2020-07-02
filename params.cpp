#include "params.hpp"
#include "process.hpp"
#include "write_screen.hpp"
#include "basic_time_integrator.hpp"


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

    switch(time_method)
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
            proc.Nfield = 1;
            break;
        case scalar_wave:
            proc.Nfield = 1;
            break;
        default:
            write::error("Equation system not recognised.");
    }

    //write::variable<real_t>("CFL", cfl);
    //write::variable<real_t>("End time", end_time);
    //write::variable<real_t>("dt", proc.dt);
    //write::variable<int>("No. of time steps", int(end_time/proc.dt));

    return;
}
