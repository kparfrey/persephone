#include "process.hpp"
#include "params.hpp"

Process::Process(Params &params)
: params(params)
{
    // The constructor just sets the params reference to the active_parameters
    // object in active_params.hpp, passed in from main.cpp.
}


void Process::write_error(string error, bool destroy)
{
    cout << "Error --- rank: " << rank << ": " << error << endl;

    if (destroy == true)
        exit(1);

    return;
}


void Process::write_message(string message)
{
    if (rank == 0)
        cout << "Progress: " << message << endl;

    return;
}


void Process::write_startup_info()
{
    if (rank) return; // Only write from the root proc

    params.write_param_info();

    return;
}
