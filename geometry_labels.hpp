#ifndef GEOMETRY_LABELS_HPP
#define GEOMETRY_LABELS_HPP

/* This file stores useful integer arrays for organising the various 
 * relationships between corners, edges, faces, and directions. Only
 * used in the setup phase, so will only be included by a few .cpp files.
 */


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


/*** corner_map and edge_map only used in transfinite_map.cpp ***/

/* 3D label of each of the 6 faces' 4 corners. Ordering is via cyclic permutation from
 * the face-normal direction; ie face 0 has normal-dir = 2, so its face coords are 
 * {0, 1} with the 0 direction moving first. */
constexpr int corner_map[6][4] = {{0, 1, 2, 3},  // face 0
                                  {4, 5, 6, 7},  //      1
                                  {0, 3, 7, 4},  //      2
                                  {1, 2, 6, 5},  //      3
                                  {0, 4, 5, 1},  //      4
                                  {3, 7, 6, 2}}; //      5


/* 3D label of each of the 6 faces' 4 edges. Ordering is as for corner_map ---
 * these two much have identical ordering systems so that the edges and corners
 * create the correct 2D edge-corner diagram. */
constexpr int edge_map[6][4] = {{ 0,  1,  2,  3},  // face 0
                                { 4,  5,  6,  7},  //      1
                                { 3, 11,  7,  8},  //      2
                                { 1, 10,  5,  9},  //      3
                                { 8,  4,  9,  0},  //      4
                                {11,  6, 10,  2}}; //      5

#endif
