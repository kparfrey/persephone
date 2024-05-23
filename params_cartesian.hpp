#include "params.hpp"


/* Simple parameter holder for homogeneous Cartesian domain */
class ParamsCartesian : public Params<ParamsCartesian>
{
    public:

    int Nproc[3];  // No. of processes in each direction in the whole domain
    int Nelem[3];  // No. of elements in each direction, in each process
    int Ns[3];     // No. of solution points in each direction, in each element

    bool periodic; // Whether we're dealing with a fully periodic domain

    /* Secondary or derived quantities */
    int Ns_elem; // Total numbers of soln/flux points per element
    int Nf_elem;
    int Nelem_proc; // Is this necessary?


    void secondary_params();
    

    /* Called from base class interface functions */
    void write_param_info_();

    template <class ProcType>
    void setup_process_(ProcType& proc);
    
    template <class ProcType>
    void setup_elementblock_(ElementBlock& elements, ProcType& proc);

    void set_initial_state_(ElementBlock& elements);


    /* Constructor */
    ParamsCartesian(EqnSystem equations,
                    int (& Nproc_)[3], 
                    int (& Nelem_)[3], 
                    int (& Ns_)[3], 
                    real_t cfl = 0.8,
                    real_t end_time = 1.0,
                    real_t dt_write = 0.5,
                    bool periodic = true)
    : Params(equations, cfl, end_time, dt_write), periodic(periodic)
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
