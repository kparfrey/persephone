#ifndef PARAMS_TORUS
#define PARAMS_TORUS

#include "params.hpp"
#include "torus_mode_pack.hpp"
#include "cerfon_freidberg.hpp"
#include "desc.hpp"

#include <string>


class ParamsTorus : public Params<ParamsTorus>
{
    public:

    TorusProblemType problem_type;
    TorusGridMethod grid_method;
    std::string input_file;


    /* For these, interpret as N*[3] = N*[central, outer, phi]
     * where the central group(s) have extent central x central x phi
     *   and the outer   groups   have extent central x outer   x phi */ 
    int Nproc[3];  // No. of processes in each direction, in each disc zone
                   // NB: Nproc[2] gives the no. of Processes the whole way around
    int Nelem[3];  // No. of elements in each direction, in each process
    int Ns[3];     // No. of solution points in each direction, in each element

    TorusCentralPolygon central_polygon; 

    //TorusModePack boundary_modes;     // Don't pass to constructor -- should be a pointer??
    //CerfonFreidbergConfig* cf_config; // Holds parameters of a Cerfon Freidberg problem
    //DescConfig* desc_config;          // .... or of a DESC problem
    int desc_iteration;               // The iteration number you want to start from

    TorusConfig* torus_config;

    /* Secondary or derived quantities */
    int Ngroup_central;  // 1 for square
    int Ngroup_outer;    // 4 for square
    int Nproc_central;   // equivalent to Nproc[0] --- 1D measurement
    int Nproc_outer;     // equivalent to Nproc[1]
    int Nproc_group_central; // total (3D) number of procs in the central group
    int Nproc_group_outer;   
    
    
    void secondary_params();


    /* Called from base class interface functions */
    void write_param_info_();

    template <class ProcType>
    void setup_process_(ProcType& proc);
    
    template <class ProcType>
    void setup_elementblock_(ElementBlock &elements, ProcType& proc);
    
    void set_initial_state_(ElementBlock &elements);


    /* Constructor */
    ParamsTorus(EqnSystem equations,
                BasicTimeMethod time_method,
                int (& Nproc_)[3], 
                int (& Nelem_)[3], 
                int (& Ns_)[3], 
                real_t cfl = 0.8,
                real_t end_time = 1.0,
                real_t dt_write = 0.5,
                TorusProblemType problem_type = cerfon_freidberg,
                TorusGridMethod grid_method = internal_surface_expansion,
                std::string input_file = "output.h5",
                int desc_iteration = 1)
    : Params(equations, time_method, cfl, end_time, dt_write),
      problem_type(problem_type), grid_method(grid_method), input_file(input_file),
      desc_iteration(desc_iteration)
    {
        central_polygon = square; 

        for (int i=0; i<3; i++)
        {
            Nproc[i] = Nproc_[i];
            Nelem[i] = Nelem_[i];
            Ns[i]    = Ns_[i];
        }
    
        secondary_params();
    }
};


#include "params_torus.cpp"

#endif
