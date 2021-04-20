#ifndef TRANSFINITE_MAP_HPP
#define TRANSFINITE_MAP_HPP
 
#include "common.hpp"

class DomainMap;
class Edge;


void analytic_transfinite_map_2D(const real_t x[2], 
                                 DomainMap* const map,
                                 const real_t corners[4][3], 
                                       real_t r[3]);

void analytic_transfinite_map_3D(const real_t x[3], 
                                 DomainMap* const map,
                                 const real_t corners[8][3], 
                                       real_t r[3]);

void polynomial_transfinite_map_2D(const real_t x[2],
                                   const Edge edges[4], 
                                   const real_t corners[4][3], 
                                         real_t r[3]);

void polynomial_transfinite_map_3D(const real_t x[3],
                                   const Edge edges[12], 
                                   const real_t corners[8][3], 
                                         real_t r[3]);

void drdx_transfinite_map_2D(const int dir, 
                             const real_t x[2],
                             const Edge edges[4], 
                             const real_t corners[4][3], 
                                   real_t dr[3]);

void drdx_transfinite_map_3D(const real_t x[3], const Edge edges[12],
                             const real_t corners[8][3], 
                                   real_t dr[3][3]);

#endif
