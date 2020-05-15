#include "common.hpp"
#include "params.hpp"

Params::Params(EqnSystem equations, real_t cfl, real_t end_time) 
       : equations(equations), cfl(cfl), end_time(end_time)
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
            std::cout << "Equation system not recognised." << std::endl;
            exit(1);
    }
}
