#include "params.hpp"
#include "write_screen.hpp"

Params::Params(EqnSystem equations, TimeIntegrator time_integrator,
                                        real_t cfl, real_t end_time) 
       : equations(equations), time_integrator(time_integrator), 
                                   cfl(cfl), end_time(end_time)
{ 
    switch (equations)
    {
        case scalar_advection:
            Nfield = 1;
            break;
        case scalar_wave:
            Nfield = 1;
            break;
        default:
            write::error("Equation system not recognised.");
    }
}
