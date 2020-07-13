#ifndef PARAMS_HPP
#define PARAMS_HPP

#include "common.hpp"

class Process; // Reference to Process type will be used below
class ElementBlock;


/* Base parameter class */

class Params
{
    public:

    EqnSystem equations;
    BasicTimeMethod time_method;

    /* Total quantities referring to whole domain */
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
 

    /* Methods: all pure virtual ones must be defined in derived classes */
    virtual void write_param_info() = 0;
    virtual void secondary_params() = 0;
    virtual void setup_process(Process &proc) = 0;
    virtual void setup_elementblock(ElementBlock &elements, Process &proc) = 0;
    virtual void set_initial_state(ElementBlock &elements, EqnSystem equations) = 0;

    void setup_process_generic(Process& proc);
};
#endif
