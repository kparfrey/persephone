#ifndef COMMON_HPP
#define COMMON_HPP

/*** Problem specifications, definitions, and conditional compilation ***
 *** that is needed everywhere. Will be eventually split into several ***
 *** files for convenience.                                           ***/


#define FLOAT_POINT_PRECISION DOUBLE

constexpr int Nvar = 1; // Total number of field variables


/***************************************************/
/**** Options above this line, do not change below */


#define SINGLE 1
#define DOUBLE 2

#if (FLOAT_POINT_PRECISION == SINGLE)
  #define real_t    float
  #define realMPI_t MPI_FLOAT
#else
  #define real_t    double
  #define realMPI_t MPI_DOUBLE
#endif


#include <iostream>
#include <iomanip>
#include <string>
#include "element_block.hpp"
#include "process.hpp"

#endif
