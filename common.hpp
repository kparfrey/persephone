#ifndef COMMON_HPP
#define COMMON_HPP

#include <math.h>


/*** Problem specifications, definitions, and conditional compilation ***
 *** that is needed everywhere. Will be eventually split into several ***
 *** files for convenience.                                           ***/

/*** Avoid includes here without a good reason....                    ***/


#define FLOAT_POINT_PRECISION DOUBLE



/***************************************************/
/**** Options above this line, do not change below */

enum EqnSystem {scalar_advection, scalar_wave};

constexpr int dirs[3] = {0, 1, 2}; // The three spatial directions

constexpr int ifaces[6] = {0,1,2,3,4,5}; // Indices for faces (proc, elem, etc.)

constexpr int icorners[8] = {0,1,2,3,4,5,6,7}; 

constexpr int iedges[12] = {0,1,2,3,4,5,6,7,8,9,10,11}; 

/* Coordinate directions spanning each of the six faces 
 * In each of the 3 pairs, the first/lower-index one is at lower values 
 * of the missing/normal coordinate */
constexpr int face_coords[6][2] = {{0,1},{0,1},{0,2},{0,2},{1,2},{1,2}};

/* Relative locations of the 8 corners for a unit cube */
constexpr int corner_coords[8][3] = {{0, 0, 0},  // 0
                                     {1, 0, 0},  // 1
                                     {1, 1, 0},  // 2
                                     {0, 1, 0},  // 3
                                     {0, 0, 1},  // 4
                                     {1, 0, 1},  // 5
                                     {1, 1, 1},  // 6
                                     {0, 1, 1}}; // 7

#define SINGLE 1
#define DOUBLE 2

#if (FLOAT_POINT_PRECISION == SINGLE)
  #define real_t    float
  #define realMPI_t MPI_FLOAT
#else
  #define real_t    double
  #define realMPI_t MPI_DOUBLE
#endif


#define pi M_PI

#endif
