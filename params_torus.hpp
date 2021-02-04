#include "params.hpp"

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

    /* Secondary or derived quantities */
    int Ngroup;  // No. of disc groups or disc zones
    
    /* Remove from torus params, since these now depend on which group the proc is in
    int Nproc_group;
    int Nf[3];   // Where is this set?
    int Ns_elem; // Total numbers of soln/flux points per element
    int Nf_elem;
    int Nelem_proc;
     */


    /* General methods */
    virtual void write_param_info();
    virtual void secondary_params();
    virtual void setup_process(Process &proc);
    virtual void setup_elementblock(ElementBlock& elements, Process &proc);
    virtual void set_initial_state(ElementBlock& elements, EqnSystem equations);
    real_t set_dt_basic(ElementBlock& elements);


    /* Constructor */
    ParamsTorus(EqnSystem equations,
                    BasicTimeMethod time_method,
                    int (& Nproc_)[3], 
                    int (& Nelem_)[3], 
                    int (& Ns_)[3], 
                    TorusCentralPolygon central_polygon,
                    GeometryClass geometry = simple_geometry,
                    real_t cfl = 0.8,
                    real_t end_time = 1.0,
                    real_t dt_write = 0.5)
    : Params(equations, time_method, geometry, cfl, end_time, dt_write), 
      central_polygon(central_polygon)
    {
        for (int i=0; i<3; i++)
        {
            Nproc[i] = Nproc_[i];
            Nelem[i] = Nelem_[i];
            Ns[i]    = Ns_[i];
        }
    
        secondary_params();
    }
};
