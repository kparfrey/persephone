#include "params_cartesian.hpp"

/* Rather than choosing a problem-specific method, may want to pass in the 
 * various parameter changes directly here in the constructor. Then each 
 * specific problem can just be stored as a separate file, and soft-linked to:
 * ln -s specific_problem_version.hpp active_params.hpp
 */

/* For discussion of non-standard features see active_params_torus_0.hpp */

static int Nproc[3] = {2,2,1};
static int Nelem[3] = {5,5,1};
static int Ns[3]    = {4,4,1};

// 0.6 stable for Euler, ~0.5 for NS? RK3_SSP more stable.
static real_t cfl      = 0.5; 
static real_t end_time = 30.0;
static real_t dt_write = 0.5;

static ParamsCartesian active_parameters(Nproc, Nelem, Ns, cfl, end_time, dt_write);
