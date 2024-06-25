#include "params_cartesian.hpp"

/* For discussion of non-standard features see active_params_torus_0.hpp */

static int Nproc[3] = {1,4,1};
static int Nelem[3] = {1,8,1};
static int Ns[3]    = {1,6,1};

static real_t cfl      = 0.8; 
static real_t end_time = 20.0;
static real_t dt_write = 0.1;

static bool periodic = false; // Non-periodic in the y (1) direction

static ParamsCartesian active_parameters(Nproc, Nelem, Ns, cfl, end_time, dt_write, periodic);
