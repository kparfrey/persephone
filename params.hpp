#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <iostream>
#include "common.hpp"

class Process; // Process type will be used below


/* Base parameter class */

class Params
{
    public:

    /* Total quantities referring to whole domain */
    int Nproc_domain;
    int Nelem_domain;
    int Ns_domain;

    real_t cfl;
    real_t end_time;

    /* Methods: all must be defined in derived classes */
    virtual void secondary_params(){};
    virtual void setup_process(Process &proc){};
    virtual void write_param_info(){std::cout << "Calling wrong method!" << std::endl;};

    Params(real_t cfl, real_t end_time) : cfl(cfl), end_time(end_time) { }
};

#endif
