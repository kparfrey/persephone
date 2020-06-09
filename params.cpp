#include "params.hpp"
#include "process.hpp"
#include "write_screen.hpp"

Params::Params(EqnSystem equations, TimeIntegrator time_integrator,
                                        real_t cfl, real_t end_time) 
       : equations(equations), time_integrator(time_integrator), 
                                   cfl(cfl), end_time(end_time)
{ 
    switch (equations)
    {
        case scalar_advection:
            Nfield = 1;
            break;
        case scalar_wave:
            Nfield = 1;
            break;
        default:
            write::error("Equation system not recognised.");
    }
}


/* That part of Process setup which is the same for Cartesian,
 * Spherica, toroidal etc. configurations. */
void Params::setup_process_generic(Process &proc)
{
    proc.time     = 0.0;
    proc.end_time = end_time;
    proc.cfl      = cfl;
    proc.data_output_counter = 0;

    // Should set dt here with a general-purpose method

    write::variable<real_t>("CFL", cfl);
    write::variable<real_t>("End time", end_time);
    write::variable<real_t>("dt", proc.dt);
    write::variable<int>("No. of time steps", int(end_time/proc.dt));

    proc.time_integrator = time_integrator;

    return;
}
