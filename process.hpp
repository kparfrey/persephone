#ifndef PROCESS_HPP
#define PROCESS_HPP

#include <iostream>
#include <string>
#include "common.hpp"
#include "face_communicator.hpp"

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
    FaceCommunicator faces[6];
    real_t corners[8][3]; // 3 spatial coordinates for each of 8 corners
    real_t edges[12];     // Lengths of edges in physical units 
    
    // Don't want to build connectivity assumptions into the Process class
    // --- all of this should be set at initialization and unnecessary after
    //int proc_idx[3];    // Set this proc's indices in the global 3D proc array
    //real_t extent[3];   // Spatial extent of proc's region in each dir
    

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
