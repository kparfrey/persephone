#ifndef RK2_MIDPOINT_HPP
#define RK2_MIDPOINT_HPP

#include "basic_time_integrator.hpp"

class RK2_midpoint : public BasicTimeIntegrator<RK2_midpoint>
{
    public:
    void takeStep(Process& proc) const;
};

#endif
