#ifndef ACTIVE_PARAMS_HPP
#define ACTIVE_PARAMS_HPP

#include "common.hpp"
#include "params_cartesian.hpp"

/* Rather than choosing a problem-specific method, may want to pass in the 
 * various parameter changes directly here in the constructor. Then each 
 * specific problem can just be stored as a separate file, and soft-linked to:
 * ln -s specific_problem_version.hpp active_params.hpp
 */

static EqnSystem equations = scalar_advection;
static TimeIntegrator time_integrator = rk2_midpoint;

static int Nproc[3] = {4,2,1};
static int Nelem[3] = {2,2,1};
static int Ns[3]    = {16,8,1};

static real_t edges[3][2] = {{0.0,2.0}, {0.0,1.0}, {0.0, 1.0}};

static real_t cfl      = 0.8;
static real_t end_time = 1.5;

ParamsCartesian active_parameters(equations, time_integrator, 
                                  Nproc, Nelem, Ns, edges,
                                  cfl, end_time);

#endif
