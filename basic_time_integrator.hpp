#ifndef BASIC_TIME_INTEGRATOR_HPP
#define BASIC_TIME_INTEGRATOR_HPP

class Process;


template <class TimeStepType>
class BasicTimeIntegrator
{
    public:

    void operator()(Process& proc) const
    {
        static_cast<const TimeStepType*>(this)->takeStep(proc);
        return;
    }
};

#endif
