#ifndef PROCESS_HPP
#define PROCESS_HPP

#include "common.hpp"
#include "face_communicator.hpp"
#include "element_block.hpp"
#include "edge.hpp"

class Params;
class BasicTimeIntegrator;
class ConservedToPrimitive;
class WaveSpeedsFromPrimitive;
class FluxesFromPrimitive;
class NumericalFlux;

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
    int Nfield;
    real_t time;
    real_t end_time;
    int step;
    real_t cfl;
    real_t dt;       // time step
    
    real_t dt_write; // time between output writes to disk
    real_t time_last_write;
    real_t step_last_write;
    int data_output_counter;
    
    BasicTimeIntegrator*     time_integrator;
    ConservedToPrimitive*    U_to_P;
    WaveSpeedsFromPrimitive* c_from_P;
    FluxesFromPrimitive*     F_from_P;
    NumericalFlux*           F_numerical;


    /* Methods */
    Process(Params &params);
    void write_startup_info();
    void time_advance();
    void find_divF(const real_t* const U, const real_t t, real_t* const divF);
    void move_to_device();
    void exchange_boundary_data();
};
#endif
