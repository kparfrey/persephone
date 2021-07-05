#ifndef CERFON_FREIDBERG_HPP
#define CERFON_FREIDBERG_HPP

#include "common.hpp"
#include "domain_map.hpp"

class CerfonFreidbergConfig
{
    public:

    enum Machine {tftr, jet, iter, mast, nstx, spheromak, frc};

    const Machine machine = frc;

    const real_t A = 0.0;
    
    real_t epsilon, kappa, delta, alpha;

    real_t c[13] = {}; // Coeffs in the potential expansion, initialised to zero


    CerfonFreidbergConfig();
};


/**
class CerfonFreidbergMap : public DomainMap
{
    public:

};
**/

#endif
