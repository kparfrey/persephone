#ifndef COMMON_HPP
#define COMMON_HPP


/*** Problem specifications, definitions, and conditional compilation ***
 *** that is needed everywhere. Will be eventually split into several ***
 *** files for convenience.                                           ***/

/*** Avoid includes here without a good reason....                    ***/


#define FLOAT_POINT_PRECISION DOUBLE



/***************************************************/
/**** Options above this line, do not change below */

enum EqnSystem {scalar_advection, scalar_wave};
enum BasicTimeMethod {rk2_midpoint, rk3_classic};

enum Components {covariant, contravariant};
enum CoordTransDir {phys2ref, ref2phys};
enum WhoWrites {root_only, all_ranks, one_rank};
enum Prognosis {destroy, survive};
//enum Operation {soln_to_flux, fluxDeriv_to_soln};

constexpr int dirs[3] = {0, 1, 2}; // The three spatial directions

constexpr int dir_plus_one[3] = {1, 2, 0}; // Cyclic addition for direction-
constexpr int dir_plus_two[3] = {2, 0, 1}; // dependent indexing

constexpr int ifaces[6] = {0,1,2,3,4,5}; // Indices for faces (proc, elem, etc.)

constexpr int icorners[8] = {0,1,2,3,4,5,6,7}; 

constexpr int iedges[12] = {0,1,2,3,4,5,6,7,8,9,10,11}; 

/* Coordinate directions spanning each of the six faces 
 * In each of the 3 pairs, the first/lower-index one is at lower values 
 * of the missing/normal coordinate */
//constexpr int face_coords[6][2] = {{0,1},{0,1},{0,2},{0,2},{1,2},{1,2}};

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
  #define MPI_real_t MPI_FLOAT
#else
  #define real_t    double
  #define MPI_real_t MPI_DOUBLE
#endif

/* Should all be double --- will need to redefine for float */
#define pi   M_PI
#define pi_2 M_PI_2
#define pi_4 M_PI_4
#define one_pi M_1_PI

#endif
