#include "common.hpp"
#include "rk3_ssp.hpp"
#include "process.hpp"
#include "kernels.hpp"
#include "element_block.hpp"


void RK3_SSP::takeStep(Process& proc) const
{
    ElementBlock& eb = proc.elements;
    const int Ntot = eb.Nfield * eb.Ns_block;

    /* Intermediate fields, reused */
    real_t* fields_inter = kernels::alloc_raw(Ntot);

    real_t* divF         = kernels::alloc_raw(Ntot);

    const real_t one_third  = 1.0/3.0;
    const real_t two_thirds = 2.0/3.0;

    proc.substep = 1;
    proc.find_divF(eb.fields, proc.time, divF);
    kernels::add_2_vectors(eb.fields, divF, 
                           1.0      , -proc.dt, 
                           fields_inter, Ntot);
    kernels::floors(fields_inter, eb.physics_soln, eb.lengths);

    //kernels::filter_field(&fields_inter[8*eb.Ns_block], eb.chebyshev_filter, eb.lengths);


    proc.substep = 2;
    proc.find_divF(fields_inter, proc.time + proc.dt, divF);
    kernels::add_3_vectors(eb.fields, fields_inter, divF, 
                           0.75     , 0.25        , -0.25*proc.dt,  
                           fields_inter, Ntot);
    kernels::floors(fields_inter, eb.physics_soln, eb.lengths);

    //kernels::filter_field(&fields_inter[8*eb.Ns_block], eb.chebyshev_filter, eb.lengths);


    proc.substep = 3;
    proc.find_divF(fields_inter, proc.time + 0.5*proc.dt, divF);
    kernels::add_3_vectors(eb.fields, fields_inter, divF, 
                           one_third, two_thirds  , -two_thirds*proc.dt,  
                           eb.fields, Ntot);
    kernels::floors(eb.fields, eb.physics_soln, eb.lengths);


    //for (int f = 0; f < 5; ++f)
    //    kernels::filter_field(&eb.fields[f*eb.Ns_block], eb.chebyshev_filter, eb.lengths);
    

    kernels::free(fields_inter);
    kernels::free(divF);

    return;
}
