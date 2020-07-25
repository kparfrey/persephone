#ifndef TRANSFINITE_MAP_HPP
#define TRANSFINITE_MAP_HPP
 
#include "common.hpp"

class DomainMap;
class Edge;


void analytic_transfinite_map_2D(const real_t x[2], DomainMap* const map,
                                 const real_t corners[4][3], real_t r[3]);

void polynomial_transfinite_map_2D(const real_t x[2], const int point_idx[2],
                                   Edge edges[4], const real_t corners[4][3], real_t r[3]);

#endif
