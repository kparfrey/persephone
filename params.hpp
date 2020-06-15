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
    TimeIntegrator time_integrator;

    /* Total quantities referring to whole domain */
    int Nproc_domain;
    int Nelem_domain;
    int Ns_domain;

    real_t cfl;
    real_t end_time;
    int Nfield;

    /* Constructor */
    Params(EqnSystem equations, TimeIntegrator time_integrator,
                                   real_t cfl, real_t end_time);

    /* Methods: all pure virtual ones must be defined in derived classes */
    virtual void write_param_info() = 0;
    virtual void secondary_params() = 0;
    virtual void setup_process(Process &proc) = 0;
    virtual void setup_elementblock(ElementBlock &elements, Process &proc) = 0;
    virtual void set_initial_state(ElementBlock &elements) = 0;

    void setup_process_generic(Process& proc);
};
#endif
