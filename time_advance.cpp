#include "time_advance.hpp"

#include "process.hpp"
#include "write_screen.hpp"
#include "kernels.hpp"


/* Dummy function to compile */
static void find_RHS(real_t* fields, real_t t, real_t* result, Process& proc)
{
    for (int i = 0; i < proc.elements.Ns_block; ++i)
        result[i] = 1.0;

    return;
}



/* Assume one ElementBlock per process */
static void step_rk2_midpoint(Process& proc)
{
    ElementBlock& eb = proc.elements;

    real_t* fields_mid = kernels::alloc_raw(eb.Ns_block);
    real_t* RHS        = kernels::alloc_raw(eb.Ns_block);

    find_RHS(eb.fields, proc.time, RHS, proc);
    kernels::add_2_vectors(eb.fields, RHS, 1.0, 0.5*proc.dt, fields_mid, eb.Ns_block);

    find_RHS(fields_mid, proc.time + 0.5*proc.dt, RHS, proc);
    kernels::add_2_vectors(eb.fields, RHS, 1.0,     proc.dt,  eb.fields, eb.Ns_block);

    kernels::free(fields_mid);
    kernels::free(RHS);

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
