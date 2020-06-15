#ifndef BASIC_TIME_INTEGRATOR_HPP
#define BASIC_TIME_INTEGRATOR_HPP

#include "process.hpp"
#include "write_screen.hpp"
#include "kernels.hpp"


/* Abstract base class for underlying fundamental time step functor */
class BasicTimeIntegrator
{
    public:
    virtual void operator()(Process& proc) = 0;
};


class RK2_midpoint : public BasicTimeIntegrator
{
    public:
    virtual void operator()(Process& proc)
    {
        ElementBlock& eb = proc.elements;

        real_t* fields_mid = kernels::alloc_raw(eb.Ns_block);
        real_t* RHS        = kernels::alloc_raw(eb.Ns_block);

        proc.find_RHS(eb.fields, proc.time, RHS);
        kernels::add_2_vectors(eb.fields, RHS, 1.0, 0.5*proc.dt, fields_mid, eb.Ns_block);

        proc.find_RHS(fields_mid, proc.time + 0.5*proc.dt, RHS);
        kernels::add_2_vectors(eb.fields, RHS, 1.0,     proc.dt,  eb.fields, eb.Ns_block);

        kernels::free(fields_mid);
        kernels::free(RHS);

        return;
    }
};
#endif
