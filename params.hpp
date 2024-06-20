#ifndef PARAMS_HPP
#define PARAMS_HPP


class ElementBlock;
class Process;


/* Base parameter class */

template <class ParamsType>
class Params
{
    public:

    /* Total quantities referring to whole domain */
    int Ngroup; 
    int Nproc_domain;
    int Nelem_domain;
    int Ns_domain;

    real_t cfl;
    real_t end_time;
    real_t dt_write; // elapsed time between writing to disk

    /* Constructor */
    Params(real_t cfl, real_t end_time, real_t dt_write)
          : cfl(cfl), end_time(end_time), dt_write(dt_write){}
 

    /* Provided by this class, in params.cpp */
    void setup_process_generic(Process& proc);


    /* Defined in derived classes */
    void write_param_info()
    {
        static_cast<ParamsType*>(this)->write_param_info_();
        return;
    }

    void setup_process(Process& proc)
    {
        static_cast<ParamsType*>(this)->setup_process_(proc);
        return;
    }
    
    void setup_elementblock(ElementBlock& elements, Process& proc)
    {
        static_cast<ParamsType*>(this)->setup_elementblock_(elements, proc);
        return;
    }
    
    void set_initial_state(ElementBlock& elements)
    {
        static_cast<ParamsType*>(this)->set_initial_state_(elements);
        return;
    }
};

#endif
