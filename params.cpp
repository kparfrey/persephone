#include <string>
#include "common.hpp"
#include "params.hpp"
//#include "write_screen.hpp"
//#include "basic_time_integrator.hpp"
//#include "physics_includes.hpp"
//#include "numerical_flux.hpp"
#include "process.hpp"


/* That part of Process setup which is the same for Cartesian,
 * Spherical, toroidal etc. configurations. */
template <class ParamsType>
void Params<ParamsType>::setup_process_generic(Process& proc)
{
    proc.time     = 0.0;
    proc.end_time = end_time;
    proc.dt_write = dt_write;
    proc.step     = 0;
    proc.cfl      = cfl;
    proc.data_output_counter = 0;
    proc.is_output_step = false;

    proc.system = PhysicsType::system; // Store in Process for convenience
    proc.Nfield = proc.elements.physics_soln.Nfield;

    for (int d: dirs)
    {
        proc.F_numerical[d].physics = &(proc.elements.physics[d]); // Convenience pointer
        
        proc.F_numerical[d].Nfield  = proc.Nfield;

        /* The above convenience pointer to proc.physics[d] is the only 
         * difference between the three copies of NumericalFlux */

        /* For running the divergence-cleaning subsystem on its own */
        //proc.F_numerical_divB_subsystem[d] = new HLL_divB_subsystem;
        //proc.F_numerical_divB_subsystem[d]->physics = proc.elements.physics[d]; // Convenience pointer
    }

    return;
}

/* Explicit template instantiation, so can have Params in its own compilation unit / object file. */
template class Params<ParamsCartesian>;
template class Params<ParamsTorus>;
