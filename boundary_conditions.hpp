#ifndef BOUNDARY_CONDITIONS_HPP 
#define BOUNDARY_CONDITIONS_HPP

#include "common.hpp"
#include "kernels.hpp"

class BoundaryConditions
{
    public:
    const int Ntot;   // Total number of points on the face
    const int Nfield;

    real_t* stored_data; // Usually the initial conditions

    BoundaryConditions(const int Ntot, const int Nfield)
        : Ntot(Ntot), Nfield(Nfield) {}

    void setup(const real_t* const data);

    /* Use a functor to set the external BCs for each field and face location.
     * Can this be made more efficient and/or cleaner? */
    virtual real_t operator()(const int field, const int i,
                              const real_t* const __restrict__ my_data) = 0;

    /* Original method, with loop in this file. Now use the functor approach 
     * above, with the loop in kernels */
    //virtual void apply(const real_t* const __restrict__ my_data,
    //                         real_t* const __restrict__ neighbour_data) = 0;
};


/* Make sure to call after the face's data has been
 * filled with the initial conditions and moved to device */
inline void BoundaryConditions::setup(const real_t* const data)
{
    stored_data = kernels::alloc_raw(Nfield * Ntot);

    for(int i = 0; i < Nfield*Ntot; ++i)
        stored_data[i] = data[i];

    return;
}


class ImplosionTestBC : public BoundaryConditions
{
    enum conserved {density, mom0, mom1, mom2, tot_energy};

    public:

    ImplosionTestBC(const int Ntot, const int Nfield)
        : BoundaryConditions(Ntot, Nfield) {}

    
    virtual real_t operator()(const int field, const int i,
                              const real_t* const __restrict__ my_data)
    {
        switch(field)
        {
            case density:
                return stored_data[density*Ntot + i];
            case mom0:
                return - my_data[mom0*Ntot + i]; // Radial wall
            case mom1:
            case mom2:
                return 0.0;
            case tot_energy:
                return stored_data[tot_energy*Ntot + i];
        }

        return 0.0;
    }

    /***
    virtual void apply(const real_t* const __restrict__ my_data,
                             real_t* const __restrict__ neighbour_data)
    {
        for (int i = 0; i < Ntot; ++i)
        {
            neighbour_data[density*Ntot + i]    = stored_data[density*Ntot + i];
            neighbour_data[mom0*Ntot + i]       = - my_data[mom0*Ntot + i]; 
            neighbour_data[mom1*Ntot + i]       = 0.0;
            neighbour_data[mom2*Ntot + i]       = 0.0;
            neighbour_data[tot_energy*Ntot + i] = stored_data[tot_energy*Ntot + i];
        }

        return;
    }
     ***/
};

#endif
