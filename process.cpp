#include "process.hpp"
#include "params.hpp"
#include "time_integrators.hpp"
#include "write_screen.hpp"

Process::Process(Params &params)
: params(params)
{
    // The constructor just sets the params reference to the active_parameters
    // object in active_params.hpp, passed in from main.cpp.
}


void Process::write_startup_info()
{
    if (rank) return; // Only write from the root proc

    params.write_param_info();

    return;
}


void Process::time_advance()
{
    //basic_timestep_method(*this);
    
    switch(time_integrator)
    {
        case rk2_midpoint:
            step_rk2_midpoint(*this);
            break;
        default:
            write::error("Time integrator not recognized.");
    }

    time += dt;

    return;
}
