#ifndef COMMON_HPP
#define COMMON_HPP

#define FLOAT_POINT_PRECISION DOUBLE


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


#include "element_brick.hpp"
#include "process.hpp"

#endif
