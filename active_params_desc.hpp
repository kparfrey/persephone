#include "common.hpp"
#include "params_torus.hpp"

/* Rather than choosing a problem-specific method, may want to pass in the 
 * various parameter changes directly here in the constructor. Then each 
 * specific problem can just be stored as a separate file, and soft-linked to:
 * ln -s specific_problem_version.hpp active_params.hpp
 */

/* Why static variables in a header? This file is only ever to be included ONCE,
 * in main.cpp via a softlink from active_params.hpp. The static variables are 
 * to make it explicit that these names aren't exposed to the linker and therefore
 * don't pollute the global namespace. 
 * In general, the code aims to have zero extern variables.
 * That this file is only to be included once is reflected in the lack of include
 * guards.
 */

static EqnSystem equations = mhd;
static BasicTimeMethod time_method = rk3_ssp;
static TorusProblemType problem_type = desc_input;

static int Nproc[3] = {1,1,1};
static int Nelem[3] = {4,4,4};
static int Ns[3]    = {6,6,6};

static real_t cfl      = 0.9; // Seems stable to 1.0 for ideal MHD
static real_t end_time = 10.0;
static real_t dt_write = 0.01;

static ParamsTorus active_parameters(equations, time_method, 
                                     Nproc, Nelem, Ns, 
                                     cfl, end_time, dt_write,
                                     problem_type,
                                     "NCSX_output.h5", 7);
                                     //"DSHAPE_output.h5", 2);
                                     //"AXISYM_output.h5", 20);
