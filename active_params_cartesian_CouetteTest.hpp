#include "params_cartesian.hpp"

/* For discussion of non-standard features see active_params_torus_0.hpp */

static int Nproc[3] = {1,2,1};
static int Nelem[3] = {1,5,1};
static int Ns[3]    = {1,5,1};

static real_t cfl      = 0.1; 
static real_t end_time = 150.0;
static real_t dt_write = 0.5;

static bool periodic = false; // Only fully periodic in the 2-direction

static ParamsCartesian active_parameters(Nproc, Nelem, Ns, cfl, end_time, dt_write, periodic);
