#ifndef TIME_INTEGRATORS_HPP
#define TIME_INTEGRATORS_HPP

#include "common.hpp"

class Process;

void step_rk2_midpoint(Process& proc);

//void time_advance(Process&);
#endif
