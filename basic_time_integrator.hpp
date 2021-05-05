#ifndef BASIC_TIME_INTEGRATOR_HPP
#define BASIC_TIME_INTEGRATOR_HPP

#include "process.hpp"
#include "kernels.hpp"
//#include "physics_includes.hpp"


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
        const int Ntot = eb.Nfield * eb.Ns_block;

        real_t* fields_mid = kernels::alloc_raw(Ntot);
        real_t* divF       = kernels::alloc_raw(Ntot);

        /* dU/dt = - divF */
        proc.substep = 1;
        proc.find_divF(eb.fields, proc.time, divF);
        kernels::add_2_vectors(eb.fields, divF, 1.0, -0.5*proc.dt, fields_mid, Ntot);

        proc.substep = 2;
        proc.find_divF(fields_mid, proc.time + 0.5*proc.dt, divF);
        kernels::add_2_vectors(eb.fields, divF, 1.0,     -proc.dt,  eb.fields, Ntot);

        /* psi source term as a separate substep */
        //kernels::multiply_by_scalar(&eb.fields[8*eb.Ns_block], 
        //        proc.system_data->psi_damping_exp, eb.Ns_block);

        kernels::free(fields_mid);
        kernels::free(divF);

        return;
    }
};
#endif
