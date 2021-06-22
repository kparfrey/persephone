#ifndef PROCESS_HPP
#define PROCESS_HPP

#include <mpi.h>
#include "common.hpp"
#include "face_communicator.hpp"
#include "element_block.hpp"
//#include "edge.hpp"

class Params;
class BasicTimeIntegrator;
class Physics;
class NumericalFlux;

#if 0
class SystemData;
class ConservedToPrimitive;
class WaveSpeedsFromPrimitive;
class FluxesFromPrimitive;
class DiffusiveFluxes;
#endif

class Process
{
    public:

    // NB: use a reference so the virtual functions work correctly...
    Params &params;
      
    /* Local data */
    int rank;
    FaceCommunicator faces[6];

    ElementBlock elements; // Start with a single ElementBlock per process...


    /* Group data 
     * A general dataset can be divided into a number of logically
     * Cartesian groups. A Cartesian dataset has only one group. */
    int group;
    int group_rank;     // rank of this Process within its group
    int group_idx[3];   // 3D indices of this process within its group
    int Nproc_group[3]; // No. of procs in each direction in this group 
                        // Should be renamed to Nproc_idx[3] or similar
    int Nproc_group_tot;
    
    MPI_Comm group_comm; // A communicator for MPI calls within this group


    /* Global data */
    int Ngroup;
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
    bool is_output_step; // true when writing data at the end of this step
    int substep; // The number of the substep being calculated - e.g. for
                 // RK3 will be 1, 2, or 3
    
    real_t tt_max_global;  // Global maximum of timestep transform, for finding c_h 

    BasicTimeIntegrator* time_integrator;

    EqnSystem      system;
    Physics*       physics;
    NumericalFlux* F_numerical;

    /***
    SystemData*              system_data;
    ConservedToPrimitive*    U_to_P;
    WaveSpeedsFromPrimitive* c_from_P;
    FluxesFromPrimitive*     F_from_P;
    DiffusiveFluxes*         F_diff;
    ***/


    /* Methods */
    Process(Params &params);
    void write_startup_info();
    void setup();
    void time_advance();
    void find_divF(const real_t* const U, const real_t t, real_t* const divF);
    void add_diffusive_flux(VectorField Uf, VectorField F);
    void move_to_device();
    void fill_external_boundary_data();
    void exchange_boundary_data();
};
#endif
