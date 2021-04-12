#ifndef BOUNDARY_CONDITIONS_HPP 
#define BOUNDARY_CONDITIONS_HPP

#include "common.hpp"

class BoundaryConditions
{
    public:
    const int Ntot;   // Total number of points on the face
    const int Nfield;

    real_t* stored_data; // Usually the initial conditions

    BoundaryConditions(const int Ntot, const int Nfield)
        : Ntot(Ntot), Nfield(Nfield) 
    {
        stored_data = new real_t [Nfield * Ntot];
    }

    /* Make sure to call after the face's data has been
     * filled with the initial conditions and moved to device */
    void setup(const real_t* const data)
    {
        for(int i = 0; i < Nfield*Ntot; ++i)
            stored_data[i] = data[i];
        return;
    }

    virtual void apply(real_t* const __restrict__ data) = 0;
};


class ImplosionTestBC : public BoundaryConditions
{
    public:

    virtual void apply(real_t* const __restrict__ data)
    {
    }
};

#endif
