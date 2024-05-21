#ifndef PARAMS_HPP
#define PARAMS_HPP

#include "common.hpp"

class ElementBlock;
class Physics;


/* Base parameter class */

template <class ParamsType>
class Params
{
    public:

    EqnSystem equations;
    BasicTimeMethod time_method;

    /* Total quantities referring to whole domain */
    int Ngroup; 
    int Nproc_domain;
    int Nelem_domain;
    int Ns_domain;

    real_t cfl;
    real_t end_time;
    real_t dt_write; // elapsed time between writing to disk

    /* Constructor */
    Params(EqnSystem equations, BasicTimeMethod time_method, 
           real_t cfl, real_t end_time, real_t dt_write)
          : equations(equations), time_method(time_method),
            cfl(cfl), end_time(end_time), dt_write(dt_write){}
 

    /* Provided below */
    template <class ProcType>
    void setup_process_generic(ProcType& proc);


    /* Defined in derived classes, which hide these functions */
    template <class ProcType>
    void setup_process(ProcType& proc){}
    
    template <class ProcType>
    void setup_elementblock(ElementBlock &elements, ProcType& proc){}
    
    void write_param_info(){}
    void secondary_params(){}
    void set_initial_state(ElementBlock &elements){}

};

#include "params.cpp"

#endif
