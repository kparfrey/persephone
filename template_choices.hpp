#ifndef TEMPLATE_CHOICES_HPP
#define TEMPLATE_CHOICES_HPP


/*** Specific choices made for the various template objects. 
 *** Will be included by common.hpp, and therefore these types are available everywhere.
 *** This reduces the need to make every class containing a templated class object a template
 *** itself in order to pass the template type down the hierarchy. ***/


#include "physics_includes.hpp"
using PhysicsType  = MHD; // Choices: ScalarAdvection - NavierStokes - MHD

#include "numerical_flux.hpp"
using NumFluxType  = HLL;

#include "params_torus.hpp"
#include "params_cartesian.hpp"
using ParamsType   = ParamsTorus; // Choices: ParamsCartesian - ParamsTorus

#include "time_integrator_includes.hpp"
using TimeStepType = RK3_SSP; // Choices: RK2_midpoint - RK3_SSP

#include "boundary_conditions.hpp"
using BCType       = WallBC_NoSlip_FixedNormalB; // WallBC_NoSlip_FixedNormalBC
                                                 // PeriodicPressureBC
                                                 // CouettePlateBC
                                                 // HartmannPlateBC

#endif
