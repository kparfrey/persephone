#include "process.hpp"

#include <string>
#include "params.hpp"
#include "write_screen.hpp"
#include "basic_time_integrator.hpp"

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
    write::message("Starting time step " + std::to_string(step) + " --- t = " + std::to_string(time));

    /* Call the fundamental time step method we're using */
    (*time_integrator)(*this); //Since storing a pointer to a BasicTimeIntegrator

    time += dt;
    step += 1;

    return;
}


/* The fundamental function of the spectral difference method */
/* Q : vector of primitive solution fields 
 * F : vector of fluxes */
void Process::find_RHS(real_t* Q, real_t t, real_t* result)
{
    for (int i = 0; i < elements.Ns_block; ++i)
        result[i] = 1.0;
    return;
}
