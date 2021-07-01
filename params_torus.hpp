#ifndef PARAMS_TORUS
#define PARAMS_TORUS

#include "params.hpp"
#include "torus_mode_pack.hpp"
#include "cerfon_freidberg.hpp"


class ParamsTorus : public Params
{
    public:

    TorusCentralPolygon central_polygon; // Use a square for now

    /* For these, interpret as N*[3] = N*[central, outer, phi]
     * where the central group(s) have extent central x central x phi
     *   and the outer   groups   have extent central x outer   x phi */ 
    int Nproc[3];  // No. of processes in each direction, in each disc zone
                   // NB: Nproc[2] gives the no. of Processes the whole way around
    int Nelem[3];  // No. of elements in each direction, in each process
    int Ns[3];     // No. of solution points in each direction, in each element

    TorusProblemType torus_problem_type;
    TorusModePack boundary_modes;    // Don't pass to constructor
    CerfonFreidbergConfig cf_config; // Holds parameters of a Cerfon Freidberg problem

    /* Secondary or derived quantities */
    int Ngroup_central;  // 1 for square
    int Ngroup_outer;    // 4 for square
    int Nproc_central;   // equivalent to Nproc[0] --- 1D measurement
    int Nproc_outer;     // equivalent to Nproc[1]
    int Nproc_group_central; // total (3D) number of procs in the central group
    int Nproc_group_outer;   
    
    

    /* General methods */
    void write_param_info() override;
    void secondary_params() override;
    void setup_process(Process &proc) override;
    void setup_elementblock(ElementBlock& elements, Process &proc) override;
    void set_initial_state(ElementBlock& elements) override;

    /* Constructor */
    ParamsTorus(EqnSystem equations,
                BasicTimeMethod time_method,
                int (& Nproc_)[3], 
                int (& Nelem_)[3], 
                int (& Ns_)[3], 
                real_t cfl = 0.8,
                real_t end_time = 1.0,
                real_t dt_write = 0.5,
                TorusCentralPolygon central_polygon = square)
    : Params(equations, time_method, cfl, end_time, dt_write), 
      central_polygon(central_polygon)
    {
        torus_problem_type = input_config_file; // Default choice

        for (int i=0; i<3; i++)
        {
            Nproc[i] = Nproc_[i];
            Nelem[i] = Nelem_[i];
            Ns[i]    = Ns_[i];
        }
    
        secondary_params();
    }
};
#endif
