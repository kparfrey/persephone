#include "params.hpp"

/* Simple parameter holder for homogeneous Cartesian domain */
class ParamsCartesian : public Params
{
    public:

    int Nproc[3];  // No. of processes in each direction in the whole domain
    int Nelem[3];  // No. of elements in each direction, in each process
    int Ns[3];     // No. of solution points in each direction, in each element

    real_t domain_edge[3][2];
    

    /* Secondary or derived quantities */
    int Nf[3];
    int Ns_elem; // Total numbers of soln/flux points per element
    int Nf_elem;
    int Nelem_proc;


    /* General methods */
    virtual void write_param_info();
    virtual void secondary_params();
    virtual void setup_process(Process &proc);
    virtual void setup_elementblock(ElementBlock &elements, Process &proc);
    virtual void set_initial_state(ElementBlock &elements);


    /* Constructor */
    ParamsCartesian(EqnSystem equations,
                    int (& Nproc_)[3], 
                    int (& Nelem_)[3], 
                    int (& Ns_)[3], 
                    real_t (& domain_edge_)[3][2],
                    real_t cfl = 0.8,
                    real_t end_time = 1.0)
    : Params(equations, cfl, end_time)
    {
        for (int i=0; i<3; i++)
        {
            Nproc[i] = Nproc_[i];
            Nelem[i] = Nelem_[i];
            Ns[i]    = Ns_[i];
            domain_edge[i][0] = domain_edge_[i][0];
            domain_edge[i][1] = domain_edge_[i][1];
        }
    
        secondary_params();
    }
};
