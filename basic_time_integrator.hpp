#ifndef BASIC_TIME_INTEGRATOR_HPP
#define BASIC_TIME_INTEGRATOR_HPP

class Process;


template <class TimeStepType>
class BasicTimeIntegrator
{
    public:
    void operator()(Process& proc)
    {
        static_cast<TimeStepType*>(this)->takeStep(proc);
        return;
    }
};

#endif
