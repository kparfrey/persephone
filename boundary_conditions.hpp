#ifndef BOUNDARY_CONDITIONS_HPP 
#define BOUNDARY_CONDITIONS_HPP


template <class T>
class BoundaryConditions
{
    public:
    
    bool stores_data; // Whether or not this object needs memory allocated
    int  orientation; // Not used in all BC types

    BoundaryConditions(){stores_data = false;}

    void dirichlet(      real_t* const __restrict__ U, 
                   const real_t* const __restrict__ nl,
                   const PhysicsType                physics,
                   const int mem) const
    {
        static_cast<T const*>(this)->dirichlet_(U, nl, physics, mem);
        return;
    }
    

    void neumann(real_t* const __restrict__ dP) const
    {
        static_cast<T const*>(this)->neumann_(dP);
        return;
    }
};



class WallBC_NoSlip_FixedNormalB : public BoundaryConditions<WallBC_NoSlip_FixedNormalB>
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy, B0, B1, B2, psi};
    enum primitive {density, v0  , v1  , v2,   pressure};

    public:

    real_t* Bdotn_initial; // Store n.B from the initial conditions

    WallBC_NoSlip_FixedNormalB() {stores_data = true;}


    void allocate(const int N)
    {
        Bdotn_initial = new real_t [N](); // Set to zero by default
    }


    void dirichlet_(      real_t* const __restrict__ U, 
                    const real_t* const __restrict__ nl,
                    const PhysicsType                physics,
                    const int mem) const
    {
        real_t P[9];  // Primitives from incoming conserved variables
        real_t nu[3]; // contravariant components of normal vector
                      
        physics.metric.raise(nl, nu, mem);

        physics.ConservedToPrimitive(U, P, mem);
        
        //const real_t vdotn = nl[0]*P[v0] + nl[1]*P[v1] + nl[2]*P[v2];
        const real_t Bdotn = nl[0]*P[B0] + nl[1]*P[B1] + nl[2]*P[B2];
                
        for (int d: dirs)
        {
            //P[v0+d] -= - vdotn * nu[d]; // Impenetrable only
            P[v0+d]  = 0.0; // No slip
            P[B0+d] += (Bdotn_initial[mem] - Bdotn) * nu[d];
        }

        //P[psi] = 0.0;

        /* Want to set the primitives and then reconstruct the conserved variables in case
         * the pressure or density was reset by the flooring procedure in the earlier call
         * to MHD::ConservedToPrimitive() */
        physics.PrimitiveToConserved(P, U, mem);

        return;
    }

    
    /*** For an adiabatic wall ***/
    void neumann_(real_t* const __restrict__ dP) const
    {
        /* Using the fact that the pressure slot holds temperature (and the Psi slot holds B^2)
         * from kernels::conserved_to_primitive_fluxpoints() */

        dP[pressure] = 0.0; // Zero total (not just normal) temperature gradient
        
        return;
    }
};


class PeriodicPressureBC : public BoundaryConditions<PeriodicPressureBC>
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure  };

    public:
    real_t dPdx = 5e-2;

    PeriodicPressureBC(){stores_data = false;}
    //PeriodicPressureBC(const int orientation): orientation(orientation) {}

    
    void dirichlet_(      real_t* const __restrict__ U, 
                    const real_t* const __restrict__ ,
                    const PhysicsType                physics,
                    const int mem) const
    {
        real_t* P = new real_t [physics.Nfield];

        physics.ConservedToPrimitive(U, P, mem);

        /* Have only 1 grid point in x-dir, so the value extrapolated
         * to the faces is equal to the central solution value */
        P[pressure] += orientation * 0.5 * dPdx;

        physics.PrimitiveToConserved(P, U, mem);

        delete[] P;
        return;
    }


    void neumann_(real_t* const __restrict__ ) const
    {
        return;
    }
};


class CouettePlateBC : public BoundaryConditions<CouettePlateBC>
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure  };

    public:
    real_t v_top = 1.0;

    CouettePlateBC(){stores_data = false;}
    //CouettePlateBC(const int orientation): orientation(orientation) {}
    
    void dirichlet_(      real_t* const __restrict__ U, 
                    const real_t* const __restrict__ ,
                    const PhysicsType                physics,
                    const int mem) const
    {
        real_t* P = new real_t [physics.Nfield];

        physics.ConservedToPrimitive(U, P, mem);

        /* orientation =  1: v0 = v_top
         * orientation = -1: v0 = 0     */
        P[v0] = 0.5 * (1.0 + orientation) * v_top;
        P[v1] = 0.0;
        P[v2] = 0.0;

        physics.PrimitiveToConserved(P, U, mem);

        delete[] P;
        return;
    }


    void neumann_(real_t* const __restrict__ ) const
    {
        return;
    }
};


/* To match analytic solution need to make incompressible: set fluxes of
 * density and y-momentum (mom1) to zero. */
class HartmannPlateBC : public BoundaryConditions<HartmannPlateBC>
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy, B0, B1, B2, psi};
    enum primitive {density, v0  , v1  , v2,   pressure  };

    public:
    real_t v_top    = 1.0;

    HartmannPlateBC(){stores_data = false;}
    //HartmannPlateBC(const int orientation): orientation(orientation) {}
    
    void dirichlet_(      real_t* const __restrict__ U, 
                    const real_t* const __restrict__ ,
                    const PhysicsType                physics,
                    const int mem) const
    {
        real_t* P = new real_t [physics.Nfield];

        physics.ConservedToPrimitive(U, P, mem);

        /* orientation =  1: v0 = v_top
         * orientation = -1: v0 = - v_top     */
        P[v0] = orientation * v_top;
        P[v1] = 0.0;
        P[v2] = 0.0;

        P[B0] = 0.0;
        P[B1] = 1.0;
        P[B2] = 0.0;

        physics.PrimitiveToConserved(P, U, mem);

        delete[] P;
        return;
    }


    void neumann_(real_t* const __restrict__ ) const
    {
        return;
    }
};

#endif
