#ifndef PROCESS_HPP
#define PROCESS_HPP

#include <iostream>
#include <string>
#include "common.hpp"

class Params;

using std::cout;
using std::endl;
using std::string;


class Process
{
    public:

    // NB: use a reference so the virtual functions work correctly...
    Params &params;
      
    /* Local data */
    int rank;

    /* Global data */
    int Nproc;


    /* Methods */
    Process(Params &params);
    void writeError(string error, bool destroy = true);
    void writeMessage(string message);
    template<typename type> void writeVariable(string message, type variable);
};


/*** Templated member functions ***/

template<typename type>
void Process::writeVariable(string message, type variable)
{
    if (rank == 0)
        cout << "Data: " << message << ": " << variable << endl;

    return;
}
#endif
