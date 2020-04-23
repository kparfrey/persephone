#ifndef ACTIVE_PARAMS_HPP
#define ACTIVE_PARAMS_HPP

#include "common.hpp"
#include "params.hpp"

/* Rather than choosing a problem-specific method, may want to pass in the 
 * various parameter changes directly here in the constructor. Then each 
 * specific problem can just be stored as a separate file, and soft-linked to:
 * ln -s specific_problem_version.hpp active_params.hpp
 */

ParamsCartesian active_parameters;

#endif
