#ifndef BOUNDARY_CONDITIONS_HPP 
#define BOUNDARY_CONDITIONS_HPP

#include "common.hpp"
#include "kernels.hpp"

class BoundaryConditions
{
    public:
    
    BoundaryConditions(){}

    virtual void dirichlet(      real_t* const __restrict__ U, 
                           const real_t* const __restrict__ nl,
                           const Physics* const __restrict__ physics,
                           const int mem) const = 0;
    
    virtual void neumann(real_t* const __restrict__ dP) const = 0;
};



class WallBC_NoSlip_ZeroNormalB : public BoundaryConditions
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy, B0, B1, B2, psi};
    enum primitive {density, v0  , v1  , v2,   pressure};

    public:

    WallBC_NoSlip_ZeroNormalB(){}

    void dirichlet(      real_t* const __restrict__ U, 
                   const real_t* const __restrict__ nl,
                   const Physics* const __restrict__ physics,
                   const int mem) const override
    {
        real_t P[9];  // Primitives from incoming conserved variables
        real_t nu[3]; // contravariant components of normal vector
                      
        physics->metric->raise(nl, nu, mem);

        physics->ConservedToPrimitive(U, P, mem);
        
        //const real_t vdotn = nl[0]*P[v0] + nl[1]*P[v1] + nl[2]*P[v2];
        //real_t vm[3];
        const real_t Bdotn = nl[0]*P[B0] + nl[1]*P[B1] + nl[2]*P[B2];
        real_t Bm[3];

        /* Remove the normal velocity and magnetic field */
        for (int d: dirs)
        {
            //vm[d] = P[v0+d] - vdotn * nu[d]; // Impenetrable -- seems more stable than no-slip?
            //vm[d] = 0.0; // Impermeability + no slip
            Bm[d] = P[B0+d] - Bdotn * nu[d]; // Zero normal magnetic flux
        }
                
        for (int d: dirs)
        {
            //U[mom0+d] = P[Density] * vm[d]; // For impermeable only; need to lower vm first
            U[mom0+d] = 0.0; // No slip
            U[  B0+d] = Bm[d];
        }

        return;
    }

    
    /*** For an adiabatic wall ***/
    void neumann(real_t* const __restrict__ dP) const override
    {
        /* Using the fact that the pressure slot holds temperature (and the Psi slot holds B^2)
         * from kernels::conserved_to_primitive_fluxpoints() */

        dP[pressure] = 0.0; // Zero total (not just normal) temperature gradient
        
        return;
    }
};


class PeriodicPressureBC : public BoundaryConditions
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure  };

    public:
    const int orientation;
    real_t dPdx = 5e-2;

    PeriodicPressureBC(const int orientation): orientation(orientation) {}

    
    void dirichlet(      real_t* const __restrict__ U, 
                   const real_t* const __restrict__ ,
                   const Physics* const __restrict__ physics,
                   const int mem) const override
    {
        real_t* P = new real_t [physics->Nfield];

        physics->ConservedToPrimitive(U, P, mem);

        /* Have only 1 grid point in x-dir, so the value extrapolated
         * to the faces is equal to the central solution value */
        P[pressure] += orientation * 0.5 * dPdx;

        physics->PrimitiveToConserved(P, U, mem);

        delete[] P;
        return;
    }


    void neumann(real_t* const __restrict__ ) const override
    {
        return;
    }
};


class CouettePlateBC : public BoundaryConditions
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure  };

    public:
    const int orientation;
    real_t v_top = 1.0;

    CouettePlateBC(const int orientation): orientation(orientation) {}
    
    void dirichlet(      real_t* const __restrict__ U, 
                   const real_t* const __restrict__ ,
                   const Physics* const __restrict__ physics,
                   const int mem) const override
    {
        real_t* P = new real_t [physics->Nfield];

        physics->ConservedToPrimitive(U, P, mem);

        /* orientation =  1: v0 = v_top
         * orientation = -1: v0 = 0     */
        P[v0] = 0.5 * (1.0 + orientation) * v_top;
        P[v1] = 0.0;
        P[v2] = 0.0;

        physics->PrimitiveToConserved(P, U, mem);

        delete[] P;
        return;
    }


    void neumann(real_t* const __restrict__ ) const override
    {
        return;
    }
};


/* To match analytic solution need to make incompressible: set fluxes of
 * density and y-momentum (mom1) to zero. */
class HartmannPlateBC : public BoundaryConditions
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy, B0, B1, B2, psi};
    enum primitive {density, v0  , v1  , v2,   pressure  };

    public:
    const int orientation;
    real_t v_top    = 1.0;

    HartmannPlateBC(const int orientation): orientation(orientation) {}
    
    void dirichlet(      real_t* const __restrict__ U, 
                   const real_t* const __restrict__ ,
                   const Physics* const __restrict__ physics,
                   const int mem) const override
    {
        real_t* P = new real_t [physics->Nfield];

        physics->ConservedToPrimitive(U, P, mem);

        /* orientation =  1: v0 = v_top
         * orientation = -1: v0 = - v_top     */
        P[v0] = orientation * v_top;
        P[v1] = 0.0;
        P[v2] = 0.0;

        P[B0] = 0.0;
        P[B1] = 1.0;
        P[B2] = 0.0;

        physics->PrimitiveToConserved(P, U, mem);

        delete[] P;
        return;
    }


    void neumann(real_t* const __restrict__ ) const override
    {
        return;
    }
};

#endif
