#ifndef COMMON_HPP
#define COMMON_HPP


/*** Problem specifications, definitions, and conditional compilation ***
 *** that is needed everywhere. Will be eventually split into several ***
 *** files for convenience.                                           ***/

/*** Avoid includes here without a good reason....                    ***/


#define FLOAT_POINT_PRECISION DOUBLE



/***************************************************/
/**** Options above this line, do not change below */

enum EqnSystem {scalar_advection, euler};
enum BasicTimeMethod {rk2_midpoint, rk3_classic};


enum WhoWrites {root_only, all_ranks, one_rank};
enum Prognosis {destroy, survive};
enum GeometryClass {simple_geometry, full_geometry};
enum TorusCentralPolygon {square, hexagon, octagon};

//enum Components {covariant, contravariant};
//enum CoordTransDir {phys2ref, ref2phys};
//enum Operation {soln_to_flux, fluxDeriv_to_soln};

constexpr int dirs[3] = {0, 1, 2}; // The three spatial directions

constexpr int dir_plus_one[3] = {1, 2, 0}; // Cyclic addition for direction-
constexpr int dir_plus_two[3] = {2, 0, 1}; // dependent indexing

constexpr int ifaces[6] = {0,1,2,3,4,5}; // Indices for faces (proc, elem, etc.)

constexpr int icorners[8] = {0,1,2,3,4,5,6,7}; 

constexpr int iedges[12] = {0,1,2,3,4,5,6,7,8,9,10,11}; 

/* A lot (all?) of the following are only used in element_block.cpp, when 
 * setting the physical coords. Break into a separate hpp file? Something 
 * like geometry_labels.hpp maybe. Might be able to use some of this in
 * FaceCommunicator too, for consistency. */ 

/* Coordinate directions spanning each of the six faces. Order is by cyclic
 * permutation: {n plus 1, n plus 2} where n is the dir normal to the face */
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

/* Corner-edge relationship: is corner n at the x = 0 or x = 1 point of edge n,
 * where x is that edge's own internal 1D coordinate? */
constexpr int edge_to_corner[8] = {0, 0, 1, 1, 0, 0, 1, 1};

/* Reference coordinate direction along which each edge is aligned / anti-aligned */
constexpr int edge_dir[12] = {0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2};

/* Reference coord direction normal to each of the 6 faces */
constexpr int face_normal[6] = {2, 2, 0, 0, 1, 1};

#ifdef __CUDACC__
#define ACCEL_DECORATOR __host__ __device__
#else
#define ACCEL_DECORATOR
#endif


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


#define TINY 1e-15
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define SIGN(A)  ((A) >  0  ?  1  : -1 )

#endif
