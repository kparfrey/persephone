#ifndef DOMAIN_MAP_TORUS_HPP
#define DOMAIN_MAP_TORUS_HPP

#include "common.hpp"
#include "domain_map.hpp"

#include <cmath>



class BasicSquareTorusMap : public DomainMap
{
    public:

    virtual void operator()(const int n, const real_t x, real_t r[3])
    {
    }
};

#endif
