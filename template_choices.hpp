#ifndef TEMPLATE_CHOICES_HPP
#define TEMPLATE_CHOICES_HPP


/*** Specific choices made for the various template objects. 
 *** Will be included by common.hpp, and therefore these types are available everywhere.
 *** This reduces the need to make every class containing a templated class object a template
 *** itself in order to pass the template type down the hierarchy. ***/

#if 1

#include "physics_includes.hpp"
using PhysicsType  = MHD;

#include "numerical_flux.hpp"
using NumFluxType  = HLL;

#include "params_torus.hpp"
#include "params_cartesian.hpp"
using ParamsType   = ParamsTorus;

#include "time_integrator_includes.hpp"
using TimeStepType = RK3_SSP;
#endif

#if 0
#define ParamsType ParamsTorus
#define TimeStepType RK3_SSP
#define NumFluxType HLL
#define PhysicsType MHD
#endif

#endif
