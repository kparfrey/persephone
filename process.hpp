#ifndef PROCESS_HPP
#define PROCESS_HPP

#include "common.hpp"
#include "face_communicator.hpp"
#include "element_block.hpp"
#include "edge.hpp"

class Params;
class BasicTimeIntegrator;


class Process
{
    public:

    // NB: use a reference so the virtual functions work correctly...
    Params &params;
      
    /* Local data */
    int rank;
    FaceCommunicator faces[6];
    real_t corners[8][3]; // 3 physical-space coordinates for each of 8 corners
    Edge edges[12]; // Should become some kind of general curve object 

    ElementBlock elements; // Start with a single ElementBlock per process...


    /* Group data 
     * A general dataset can be divided into a number of logically
     * Cartesian groups. A Cartesian dataset has only one group. */
    int group;
    int group_idx[3]; // 3D indices of this process within its group
    int Nproc_group[3]; // No. of procs in each direction in this group
    

    /* Global data */
    int Nproc;
    int data_output_counter;
    real_t time;
    real_t end_time;
    int step;
    real_t cfl;
    real_t dt;
    BasicTimeIntegrator* time_integrator;


    /* Methods */
    Process(Params &params);
    void write_startup_info();
    void time_advance();
    void find_RHS(real_t* fields, real_t t, real_t* result);
    void move_to_device();
};
#endif
