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
    int proc_idx[3];    // Set this proc's indices in the global 3D proc array
    real_t extent[3];   // Spatial extent of proc's region in each dir
    real_t edges[3][2]; 

    /* Global data */
    int Nproc;


    /* Methods */
    Process(Params &params);
    void write_error(string error, bool destroy = true);
    void write_message(string message);
    void write_startup_info();
    template<typename type> void write_variable(string message, type variable);
};


/*** Templated member functions ***/

template<typename type>
void Process::write_variable(string message, type variable)
{
    if (rank == 0)
        cout << "Data: " << message << ": " << variable << endl;

    return;
}
#endif
