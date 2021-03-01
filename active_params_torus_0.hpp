#ifndef ACTIVE_PARAMS_HPP
#define ACTIVE_PARAMS_HPP

#include "common.hpp"
#include "params_torus.hpp"
#include "torus_mode_pack.hpp"

/* Rather than choosing a problem-specific method, may want to pass in the 
 * various parameter changes directly here in the constructor. Then each 
 * specific problem can just be stored as a separate file, and soft-linked to:
 * ln -s specific_problem_version.hpp active_params.hpp
 */

static EqnSystem equations = euler;
static BasicTimeMethod time_method = rk2_midpoint;
static TorusCentralPolygon central_polygon = square;

static int Nproc[3] = {1,1,1};
static int Nelem[3] = {4,3,32};
static int Ns[3]    = {8,6,10};

static real_t cfl      = 0.8;
static real_t end_time = 0.5;
static real_t dt_write = 0.25;

constexpr int Nm = 3;
constexpr int Nk = 3;
static real_t Rmk[Nm][Nk] = {{2.0, 0.289, 5.73e-2}, {0.8, 0.15, 6.89e-2}, {-0.2, -3.8e-3, 4.72e-2}};
static real_t Zmk[Nm][Nk] = {{0.0, 0.45, 4.8e-3}, {1.6, -0.53, -9.28e-3}, {-0.2,5.73e-3, -2.04e-2}};
//static real_t Rmk[Nm][Nk] = {{2.0}, {1.5}, {0.0}, {0.0}};
//static real_t Zmk[Nm][Nk] = {{0.0}, {0.8}, {0.0}, {0.0}};
TorusModePack TMP(Rmk, Zmk);

ParamsTorus active_parameters(equations, time_method, 
                              Nproc, Nelem, Ns, 
                              cfl, end_time, dt_write,
                              central_polygon, TMP);

#endif
