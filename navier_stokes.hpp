#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include <cmath>
#include "physics.hpp"

/* U --- conserved
 * 0 rho
 * 1 rho v_0
 * 2 rho v_1
 * 3 rho v_2
 * 4 E  ---  total energy density
 *
 * P --- primitive
 * 0 rho
 * 1 v^0
 * 2 v^1
 * 3 v^2
 * 4 p  ---  pressure
 */

class NavierStokes : public Physics<NavierStokes>
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure  };

    public:
    const real_t gamma = 5.0/3.0;
    const real_t gm1   = gamma - 1.0;

    NavierStokes()
    {
        system = navier_stokes;
        Nfield = 5;
        variables = new string [5];
        variables[0] = "rho";
        variables[1] = "v0";
        variables[2] = "v1";
        variables[3] = "v2";
        variables[4] = "p";

        diffusive = true;
        viscosity = 1e-2;
        
        diffusive_timestep_const = 0.3; 
    }


    void ConservedToPrimitive_(const real_t* const __restrict__ U, 
                                     real_t* const __restrict__ P,
                               const int mem) const ;

    void PrimitiveToConserved_(const real_t* const __restrict__ P, 
                                     real_t* const __restrict__ U,
                               const int mem) const ;

    void WaveSpeeds_(const real_t* const __restrict__ P, 
                           real_t* const __restrict__ c,
                     const int dir,
                     const int mem) const ;

    void Fluxes_(const real_t* const __restrict__ P, 
                       real_t (*__restrict__ F)[3],
                 const int mem) const ;

    void DiffusiveFluxes_(const real_t* const __restrict__ U, 
                          const real_t (* const __restrict__ dU)[3],
                                real_t (*__restrict__ F)[3],
                          const int mem) const ;

    void OrthonormaliseVectors_(real_t* const P, const int mem) const ;
    
    void Floors_(real_t* const __restrict__ U, const int mem) const 
    {
        return;
    }
};


inline void NavierStokes::ConservedToPrimitive_(const real_t* const __restrict__ U, 
                                                      real_t* const __restrict__ P,
                                                const int mem) const
{
    P[density] = U[density]; 

    real_t vl[3]; // Covariant velocity components

    vl[0] = U[mom0] / U[density];
    vl[1] = U[mom1] / U[density];
    vl[2] = U[mom2] / U[density];

    /* Store the contravariant components into P */
    metric.raise(vl, &P[v0], mem);
    
    const real_t KE_density = 0.5 * P[density] * 
                              (vl[0]*P[v0] + vl[1]*P[v1] + vl[2]*P[v2]);

    P[pressure] = gm1 * (U[tot_energy] - KE_density);

    return;
}



inline void NavierStokes::PrimitiveToConserved_(const real_t* const __restrict__ P, 
                                                      real_t* const __restrict__ U,
                                                const int mem) const
{
    const real_t rho = P[density];

    U[density] = rho;

    real_t vl[3];
    metric.lower(&P[v0], vl, mem);

    U[mom0] = rho * vl[0];
    U[mom1] = rho * vl[1];
    U[mom2] = rho * vl[2];
    
    const real_t KE_density = 0.5 * rho * (vl[0]*P[v0] + vl[1]*P[v1] + vl[2]*P[v2]);

    U[tot_energy] = KE_density + P[pressure] / gm1;

    return;
}


/* Returns the contravariant component of the max wave speed */
inline void NavierStokes::WaveSpeeds_(const real_t* const __restrict__ P, 
                                            real_t* const __restrict__ c,
                                      const int dir,
                                      const int mem) const 
{
    const real_t sound_speed_sq = gamma * P[pressure] / P[density];
    
    /* Assume diagonal metric for now... Tag:DIAGONAL */
    /* This is the contravariant comp of the sound speed */
    const real_t sound_speed = std::sqrt(sound_speed_sq / metric.g[dir][mem]);

    c[0] = MAX(0.0, P[v0+dir] + sound_speed);
    c[1] = MIN(0.0, P[v0+dir] - sound_speed);

    return;
}


inline void NavierStokes::Fluxes_(const real_t* const __restrict__ P, 
                                        real_t (*__restrict__ F)[3],
                                  const int mem) const
{
    /* Pointer directly to contravariant velocity components */
    const real_t* const vu = &P[v0];

    /* Covariant component of velocity */
    real_t vl[3];
    metric.lower(vu, vl, mem);

    const real_t KE_density = 0.5 * P[density] * (vl[0]*vu[0] + vl[1]*vu[1] + vl[2]*vu[2]); 
    const real_t E          = KE_density + P[pressure] / gm1;

    real_t v;

    for (int d: dirs)
    {
        v = vu[d]; // velocity in this "flux direction"

        F[density][d]    = v * P[density];
        F[tot_energy][d] = v * (E + P[pressure]);

        F[mom0][d] = P[density] * v * vl[0];
        F[mom1][d] = P[density] * v * vl[1];
        F[mom2][d] = P[density] * v * vl[2];

        F[mom0+d][d] += P[pressure];
    }

    return;
}


/* Still needs to be finished for general coords --- could build in from MHD object */
inline void NavierStokes::DiffusiveFluxes_(const real_t* const __restrict__    P, 
                                           const real_t (* const __restrict__ dP)[3],
                                                 real_t (*__restrict__ F)[3],
                                           const int mem) const
{
    const real_t mu = P[density] * viscosity; // mu = dynamic viscosity
    //const real_t mu = viscosity; // For Couette test, which is defined for const. *dynamic* viscosity

    const real_t lambda = - (2.0/3.0) * mu; // from Stokes hypothesis: zero bulk viscosity (ζ)
    
    real_t v[3];
    real_t tau[3][3];
    real_t dv[3][3]; // dv[component][deriv dir]
    real_t divv;
    real_t rdetg_deriv[3]; // (1/rdetg)*d_j(rdetg) at a single point

    for (int d: dirs)
        v[d] = P[v0+d]; // v[d] = v^d -- P stores the contravariant components

    for (int comp: dirs)
        for (int deriv: dirs)
            dv[comp][deriv] = dP[v0+comp][deriv]; // dv[i][j] = d_j v^i
            
    for (int d: dirs)
        rdetg_deriv[d] = metric.rdetg_deriv[d][mem];

    /* div(v) = d_i v^i + v^j d_j(rdetg) / rdetg */
    divv = dv[0][0] + dv[1][1] + dv[2][2] + 
           v[0] * rdetg_deriv[0] + v[1] * rdetg_deriv[1] + v[2] * rdetg_deriv[2];

    for (int d: dirs)
        tau[d][d] = 2*mu*dv[d][d] + lambda*divv;

    tau[0][1] = tau[1][0] = mu * (dv[0][1] + dv[1][0]);
    tau[0][2] = tau[2][0] = mu * (dv[0][2] + dv[2][0]);
    tau[1][2] = tau[2][1] = mu * (dv[1][2] + dv[2][1]);

    for (int d: dirs) // flux direction
    {
        F[Density][d] = 0.0;

        for (int i: dirs)
            F[mom0+i][d] = - tau[d][i];

        /* F(tot E)^d = - v^i tau_i^d
         * i.e. dE/dt = ... - div( v.tau )
         * Contracting on tau's lower index, so don't need any metric terms */
        F[tot_energy][d] = - (v[0]*tau[0][d] + v[1]*tau[1][d] + v[2]*tau[2][d]); 
    }

    return;
}


inline void NavierStokes::OrthonormaliseVectors_(real_t* const P, const int mem) const
{
    metric.orthonormals(&P[v0], mem);

    return;
}
#endif
