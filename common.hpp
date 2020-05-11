#ifndef COMMON_HPP
#define COMMON_HPP

/*** Problem specifications, definitions, and conditional compilation ***
 *** that is needed everywhere. Will be eventually split into several ***
 *** files for convenience.                                           ***/

/*** Avoid includes here without a good reason....                    ***/


#define FLOAT_POINT_PRECISION DOUBLE

constexpr int Nvar = 1; // Total number of field variables

constexpr int dirs[3] = {0, 1, 2}; // The three spatial directions

constexpr int ifaces[6] = {0,1,2,3,4,5}; // Indices for faces (proc, elem, etc.)


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

#endif
