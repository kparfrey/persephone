#ifndef PROCESS_HPP
#define PROCESS_HPP

#include <iostream>
#include <string>
#include "common.hpp"
#include "face_communicator.hpp"
#include "element_block.hpp"

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
    real_t corners[8][3]; // 3 physical-space coordinates for each of 8 corners
    //real_t edges[12]; // Should become some kind of general curve object 

    ElementBlock elements; // Start with a single ElementBlock per process...
    

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
