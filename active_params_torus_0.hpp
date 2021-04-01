#include "common.hpp"
#include "params_torus.hpp"
#include "torus_mode_pack.hpp"

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

static EqnSystem equations = euler;
static BasicTimeMethod time_method = rk2_midpoint;
static TorusCentralPolygon central_polygon = square;

static int Nproc[3] = {2,2,1};
static int Nelem[3] = {4,3,2};
static int Ns[3]    = {8,6,4};

static real_t cfl      = 0.8;
static real_t end_time = 0.200;
static real_t dt_write = 0.010;

//constexpr static int Nm = 3;
//constexpr static int Nk = 3;
//static real_t Rmk[Nm][Nk] = {{2.0, 0.289, 5.73e-2}, {0.8, 0.15, 6.89e-2}, {-0.2, -3.8e-3, 4.72e-2}};
//static real_t Zmk[Nm][Nk] = {{0.0, 0.45, 4.8e-3}, {1.6, -0.53, -9.28e-3}, {-0.2,5.73e-3, -2.04e-2}};
////static real_t Rmk[Nm][Nk] = {{2.0}, {1.5}, {0.0}, {0.0}};
////static real_t Zmk[Nm][Nk] = {{0.0}, {0.8}, {0.0}, {0.0}};

/* Circular unit-radius cylinder centred at the origin */
static real_t Rmk[2][1] = {{0.0}, {1.0}};
static real_t Zmk[2][1] = {{0.0}, {1.0}};

static TorusModePack boundary_modes(Rmk, Zmk);

static ParamsTorus active_parameters(equations, time_method, 
                                     Nproc, Nelem, Ns, 
                                     cfl, end_time, dt_write,
                                     central_polygon, boundary_modes);
