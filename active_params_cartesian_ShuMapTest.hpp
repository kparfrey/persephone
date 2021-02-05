#ifndef ACTIVE_PARAMS_HPP
#define ACTIVE_PARAMS_HPP

#include "common.hpp"
#include "params_cartesian.hpp"

/* Rather than choosing a problem-specific method, may want to pass in the 
 * various parameter changes directly here in the constructor. Then each 
 * specific problem can just be stored as a separate file, and soft-linked to:
 * ln -s specific_problem_version.hpp active_params.hpp
 */

static EqnSystem equations = euler;
static BasicTimeMethod time_method = rk2_midpoint;

static int Nproc[3] = {2,2,1};
static int Nelem[3] = {3,3,1};
static int Ns[3]    = {12,8,1};

static GeometryClass geometry = full_geometry;
static real_t limits[3][2] = {{-5.,5.}, {-3.,3.}, {-0.1, 0.1}};

static real_t cfl      = 0.8;
static real_t end_time = 75.0;
static real_t dt_write = 0.25;

ParamsCartesian active_parameters(equations, time_method, 
                                  Nproc, Nelem, Ns,
                                  cfl, end_time, dt_write,
                                  geometry, limits);

#endif
