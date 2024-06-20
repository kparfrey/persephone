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


/*** Constructing the params object ***/
static TorusProblemType problem_type = desc_input;
static TorusGridMethod  grid_method  = internal_surface_expansion; // blend_boundary_to_axis; 

static int Nproc[3] = {1,1,1};
static int Nelem[3] = {8,8,1};
static int Ns[3]    = {5,5,5};

static real_t cfl      = 0.8; // Seems stable to 1.0 for ideal MHD
static real_t end_time = 0.01;
static real_t dt_write = 2e-4;

static ParamsTorus active_parameters(Nproc, Nelem, Ns, 
                                     cfl, end_time, dt_write,
                                     problem_type, grid_method,
                                     //"easy_eqlm_cheb1_output.h5", 2);
                                     //"W7-X_output.h5", 2);
                                     //"DESC_example_stellarator.h5",0);
                                     //"WISTELL-A_output.h5", 4);
                                     //"NCSX_quad_output.h5", 7);
                                     //"NCSX_cheb1_output.h5", 7);
                                     "DSHAPE_output.h5", 2);
                                     //"DSHAPE_cheb1_output.h5", 2);
                                     //"AXISYM_output.h5", 20);
