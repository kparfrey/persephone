#include "process.hpp"
#include "params.hpp"

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
