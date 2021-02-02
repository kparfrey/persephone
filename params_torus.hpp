#include "params.hpp"

class ParamsTorus : public Params
{
    public:

    int Nproc[3];  // No. of processes in each direction in the whole domain
    int Nelem[3];  // No. of elements in each direction, in each process
    int Ns[3];     // No. of solution points in each direction, in each element

    /* Secondary or derived quantities */
    int Nf[3];
    int Ns_elem; // Total numbers of soln/flux points per element
    int Nf_elem;
    int Nelem_proc;


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
                    real_t (& domain_limits_)[3][2],
                    real_t cfl = 0.8,
                    real_t end_time = 1.0,
                    real_t dt_write = 0.5)
    : Params(equations, time_method, cfl, end_time, dt_write)
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
