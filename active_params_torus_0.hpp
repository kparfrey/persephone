#ifndef ACTIVE_PARAMS_HPP
#define ACTIVE_PARAMS_HPP

#include "common.hpp"
#include "params_torus.hpp"

/* Rather than choosing a problem-specific method, may want to pass in the 
 * various parameter changes directly here in the constructor. Then each 
 * specific problem can just be stored as a separate file, and soft-linked to:
 * ln -s specific_problem_version.hpp active_params.hpp
 */

static EqnSystem equations = euler;
static BasicTimeMethod time_method = rk2_midpoint;
static TorusCentralPolygon central_polygon = square;

static int Nproc[3] = {1,1,1};
static int Nelem[3] = {3,3,2};
static int Ns[3]    = {8,8,5};

static real_t cfl      = 0.8;
static real_t end_time = 0.5;
static real_t dt_write = 0.25;

ParamsTorus active_parameters(equations, time_method, 
                              Nproc, Nelem, Ns, 
                              cfl, end_time, dt_write,
                              central_polygon);

#endif
