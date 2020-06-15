#include "params.hpp"
#include "process.hpp"
#include "write_screen.hpp"
#include "basic_time_integrator.hpp"


Params::Params(EqnSystem equations, BasicTimeMethod time_method,
                                        real_t cfl, real_t end_time) 
       : equations(equations), time_method(time_method), 
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
 * Spherical, toroidal etc. configurations. */
void Params::setup_process_generic(Process &proc)
{
    proc.time     = 0.0;
    proc.end_time = end_time;
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

    // Should set dt here with a general-purpose method

    write::variable<real_t>("CFL", cfl);
    write::variable<real_t>("End time", end_time);
    write::variable<real_t>("dt", proc.dt);
    write::variable<int>("No. of time steps", int(end_time/proc.dt));

    return;
}
