#ifndef COMMON_HPP
#define COMMON_HPP


/*** Problem specifications, definitions, and conditional compilation ***
 *** that is needed everywhere. Will be eventually split into several ***
 *** files for convenience.                                           ***/

/*** Avoid includes here without a good reason....                    ***/


#define FLOAT_POINT_PRECISION DOUBLE_PRECISION_MY
//enum Precision {single_precision, double_precision}
//const Precision precision = double_precision

/***************************************************/
/**** Options above this line, do not change below */

enum EqnSystem {scalar_advection, navier_stokes, mhd};
enum BasicTimeMethod {rk2_midpoint, rk3_ssp};
enum PhysicalCoords {cartesian, cylindrical}; // Minkowski space unless noted

enum WhoWrites {root_only, all_ranks, one_rank};
enum Prognosis {destroy, survive};
enum TorusCentralPolygon {square, hexagon, octagon};
enum TorusProblemType {cerfon_freidberg, explicit_modes, desc_input}; 

enum OpMode {normal_mode, vecpot_mode}; // Whether currently in "normal" mode or
                                        // dealing with the vector potential

//enum Components {covariant, contravariant};
//enum CoordTransDir {phys2ref, ref2phys};
//enum Operation {soln_to_flux, fluxDeriv_to_soln};

constexpr int dirs[3] = {0, 1, 2}; // The three spatial directions

constexpr int dir_plus_one[3] = {1, 2, 0}; // Cyclic addition for direction-
constexpr int dir_plus_two[3] = {2, 0, 1}; // dependent indexing

constexpr int ifaces[6] = {0,1,2,3,4,5}; // Indices for faces (proc, elem, etc.)

constexpr int icorners[8] = {0,1,2,3,4,5,6,7}; 

constexpr int iedges[12] = {0,1,2,3,4,5,6,7,8,9,10,11}; 


#ifdef __CUDACC__
#define ACCEL_DECORATOR __host__ __device__
#else
#define ACCEL_DECORATOR
#endif


#define SINGLE_PRECISION_MY 1
#define DOUBLE_PRECISION_MY 2

#if (FLOAT_POINT_PRECISION == SINGLE_PRECISION_MY)
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


#define TINY 1e-15
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define SIGN(A)  ((A) >  0  ?  1  : -1 )

#endif
