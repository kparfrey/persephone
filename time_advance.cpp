#include "time_advance.hpp"

#include "process.hpp"
#include "write_screen.hpp"

static void step_rk2_midpoint(Process& proc)
{
    return;
}



void time_advance(Process& proc)
{
    switch(proc.time_integrator)
    {
        case rk2_midpoint:
            step_rk2_midpoint(proc);
            break;
        default:
            write::error("Time integrator not recognized.");
    }

    proc.time += proc.dt;

    return;
}
