#ifndef PROCESS_HPP
#define PROCESS_HPP

#include <iostream>
#include <string>
#include "common.hpp"
#include "params.hpp"
#include "active_params.hpp"

using std::cout;
using std::endl;
using std::string;


class Process
{
    public:

    Params params = active_parameters;

      
    /* Local data */
    int rank;
    int Nbricks; // No. of element bricks in this process

    /* Global data */
    int Nproc;


    /* Methods */
    void writeError(string error, bool destroy);
    void writeMessage(string message);
    template<typename type> void writeVariable(string message, type variable);
};


void Process::writeError(string error, bool destroy = true)
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

template<typename type>
void Process::writeVariable(string message, type variable)
{
    if (rank == 0)
        cout << "Data: " << message << ": " << variable << endl;

    return;
}


#endif
