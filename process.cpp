#include "process.hpp"

Process::Process(Params &params)
: params(params)
{
    // The constructor just sets the params reference to the active_parameters
    // object in active_params.hpp, passed in from main.cpp.
}


void Process::writeError(string error, bool destroy)
{
    cout << "Error --- rank: " << rank << ": " << error << endl;

    if (destroy == true)
        exit(1);

    return;
}


void Process::writeMessage(string message)
{
    if (rank == 0)
        cout << "Progress: " << message << endl;

    return;
}
