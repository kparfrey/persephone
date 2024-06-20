#ifndef RK3_SSP_HPP
#define RK3_SSP_HPP

#include "basic_time_integrator.hpp"

class RK3_SSP : public BasicTimeIntegrator<RK3_SSP>
{
    public:

    void takeStep(Process& proc) const;
};

#endif

