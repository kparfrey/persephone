#ifndef MHD_HPP
#define MHD_HPP

#include <cmath>
#include "physics.hpp"
#include <iostream>

/* U --- conserved
 * 0 rho
 * 1 rho v_0
 * 2 rho v_1
 * 3 rho v_2
 * 4 E  ---  total energy density
 * 5 B^0
 * 6 B^1
 * 7 B^2
 * 8 psi --- div-cleaning scalar field
 *
 * P --- primitive
 * 0 rho
 * 1 v^0
 * 2 v^1
 * 3 v^2
 * 4 p  ---  pressure
 * 5 B^0
 * 6 B^1
 * 7 B^2
 * 8 psi
 */

class MHD : public Physics
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy, B0, B1, B2, psi};
    enum primitive {density, v0  , v1  , v2,   pressure};

    public:
    const real_t gamma = 5.0/3.0;
    const real_t gm1 = gamma - 1.0;

    /* The code uses Heaviside-Lorentz units internally: E_tot = 0.5 * B^2 + p/(gamma - 1) */
    const real_t sqrt_mu0 = std::sqrt(4*pi*1e-7); // B input & output in Teslas
    //const real_t sqrt_mu0 = 1.0; // Input & output in Heaviside-Lorentz: E_mag = 0.5 * B^2
    
    //const real_t p_floor = 1e-3; // pressure floor, for recovering from inversion errors...

    MHD()
    {
        system = mhd;
        Nfield = 9;
        variables = new string [9];
        variables[0] = "rho";
        variables[1] = "v0";
        variables[2] = "v1";
        variables[3] = "v2";
        variables[4] = "p";
        variables[5] = "B0";
        variables[6] = "B1";
        variables[7] = "B2";
        variables[8] = "psi";

        diffusive   = true;
        viscosity   = 1e-6;
        resistivity = 1e-6;
        diffusive_timestep_const = 1.0; // Default: 1/3, but larger can be more stable?!

        /* Divergence-cleaning parameters */
        psi_damping_const = 1.0; // c_r --- Dedner suggests 0.18
        //psi_damping_const = 0.03; // 0 < p_d_const < 1

        apply_floors = true;
    }


    void ConservedToPrimitive(const real_t* const __restrict__ U, 
                                    real_t* const __restrict__ P,
                              const int mem) const override;

    void WaveSpeeds(const real_t* const __restrict__ P, 
                          real_t* const __restrict__ c,
                    const int dir,
                    const int mem) const override; 

    void Fluxes(const real_t* const __restrict__ P, 
                      real_t (*__restrict__ F)[3],
                const int mem) const override;

    void Fluxes_divB_subsystem(const real_t* const __restrict__ P, 
                                     real_t (*__restrict__ F)[3],
                               const int mem) const override;

    void DiffusiveFluxes(const real_t* const __restrict__ U, 
                         const real_t (* const __restrict__ dU)[3],
                               real_t (*__restrict__ F)[3],
                         const int mem) const override;

    void OrthonormaliseVectors(real_t* const P, const int mem) const override;
    
    void Floors(real_t* const __restrict__ U, const int mem) const override;
};


inline void MHD::ConservedToPrimitive(const real_t* const __restrict__ U, 
                                            real_t* const __restrict__ P,
                                      const int mem) const
{
    P[density] = U[density]; 

    if (P[density] < 1e-3)
    {
        if (apply_floors)
            P[density] = 1e-3;
        else if (P[density] < 0)
            exit(99);
    }

    real_t vl[3]; // Covariant velocity components

    vl[0] = U[mom0] / P[density]; // Momenta are covariant components
    vl[1] = U[mom1] / P[density];
    vl[2] = U[mom2] / P[density];

    /* Store the contravariant components into P */
    metric->raise(vl, &P[v0], mem);
    
    /* Magnetic fields are contravariant components */
    P[B0] = U[B0];
    P[B1] = U[B1];
    P[B2] = U[B2];

    P[psi] = U[psi];

    const real_t KE_density = 0.5 * P[density] * 
                              (vl[0]*P[v0] + vl[1]*P[v1] + vl[2]*P[v2]);

    const real_t mag_density = 0.5 * metric->square(&P[B0], mem);
    
    P[pressure] = gm1 * (U[tot_energy] - KE_density - mag_density);

    if (apply_floors)
    {
        const real_t beta = P[pressure] / mag_density;
        const real_t beta_min = 1e-4;

        if (beta < beta_min)
            P[pressure] = beta_min * mag_density;
    }

    /**
    if (KE_density + mag_density > U[tot_energy])
    {
        std::cout << "Inversion error:  R = " << metric->rdetg[mem] << "  Mag = " << mag_density << "  Energy = " << U[tot_energy] << "  KE = " << KE_density << "   p = " << U[tot_energy] - KE_density - mag_density << std::endl;
        //P[pressure] = p_floor;
        exit(45);
    }
    else
        P[pressure] = gm1 * (U[tot_energy] - KE_density - mag_density);
     **/

    return;
}


/* Floors just calls cons-to-prim -- put the actual flooring mechanism there,
 * since you'll want to use it also on flux points, at the boundary etc. */
inline void MHD::Floors(real_t* const __restrict__ U, const int mem) const
{
    real_t* P = new real_t [9];

    MHD::ConservedToPrimitive(U, P, mem);

    U[Density] = P[density];

    /* Having to recalculate vsq and Bsq here -- can probably do better */
    const real_t KE_density  = 0.5 * P[density] * metric->square(&P[v0], mem);
    const real_t mag_density = 0.5 * metric->square(&P[B0], mem);

    U[tot_energy] = mag_density + KE_density + P[pressure]/gm1;

    delete[] P;
    return;
}


inline void MHD::WaveSpeeds(const real_t* const __restrict__ P, 
                                  real_t* const __restrict__ c,
                            const int dir,
                            const int mem) const 
{
    const real_t p = P[pressure];

    /* cs_sq = (sound speed)^2 */
    const real_t cs_sq = gamma * p / P[density];

    /* va_sq = (alfven speed)^2 */
    const real_t* const Bu = &P[B0];
    const real_t Bsq   = metric->square(Bu, mem);
    const real_t va_sq =  Bsq / P[density]; 
    
    /* This version from e.g. Stone+ 2008, Appendix A, Eqn A10 */
#if 0
    const real_t va_dir_sq = P[B0+dir]*P[B0+dir]/P[density]; 
    const real_t factor = std::sqrt(std::pow(cs_sq+va_sq,2.0) - 4.0*cs_sq*va_dir_sq);
    const real_t fast_speed = std::sqrt(0.5*(cs_sq + va_sq + factor));
#endif

    /* From Spuit, Essential MHD for Astrophysics, https://arxiv.org/pdf/1301.5572.pdf */
    real_t fast_speed;
    if (va_sq > 1e-12)
    {
        real_t ku[3] = {0.0,0.0,0.0}; // Contravariant components of wavevector
        real_t kl[3];
        ku[dir] = 1.0; // Note the normalisation here is taken care of at the end
        metric->lower(ku, kl, mem);
        const real_t kdotB = kl[0]*Bu[0] + kl[1]*Bu[1] + kl[2]*Bu[2];
        const real_t ksq   = kl[0]*ku[0] + kl[1]*ku[1] + kl[2]*ku[2];
        const real_t cost  = kdotB / std::sqrt(ksq * Bsq);

        const real_t b = std::sqrt(cs_sq/va_sq) + std::sqrt(va_sq/cs_sq);

        /* The square of the fast speed */
        const real_t usq = 0.5 * (cs_sq + va_sq)*(1 + std::sqrt(1 - 4*cost*cost/(b*b)));

        /* This is the contra component. Assume diagonal metric for now... Tag:DIAGONAL */
        fast_speed = std::sqrt(usq / ((DiagonalSpatialMetric*)metric)->g[dir][mem]);
    }
    else
        fast_speed = std::sqrt(cs_sq);


    c[0] = MAX(0.0, P[v0+dir] + fast_speed); // P stores the contra velocity components
    c[1] = MIN(0.0, P[v0+dir] - fast_speed);

    return;
}


inline void MHD::Fluxes(const real_t* const __restrict__ P, 
                              real_t (*__restrict__ F)[3],
                        const int mem) const
{
    /* Pointers directly to contravariant components */
    const real_t* const vu = &P[v0];
    const real_t* const Bu = &P[B0];

    /* Covariant components */
    real_t vl[3];
    real_t Bl[3];
    metric->lower(vu, vl, mem);
    metric->lower(Bu, Bl, mem);

    const real_t KE_density = 0.5 * P[density] * (vl[0]*vu[0] + vl[1]*vu[1] + vl[2]*vu[2]); 
    const real_t Bsq   = Bl[0]*Bu[0] + Bl[1]*Bu[1] + Bl[2]*Bu[2];
    const real_t E     = KE_density + 0.5*Bsq + P[pressure] / gm1;
    const real_t Ptot  = P[pressure] + 0.5*Bsq;
    const real_t Bdotv = Bl[0]*vu[0] + Bl[1]*vu[1] + Bl[2]*vu[2];

    real_t v;
    real_t B;

    for (int d: dirs)
    {
        v = vu[d]; // velocity in this direction
        B = Bu[d]; // B in this direction

        F[density][d]    = v * P[density];
        F[tot_energy][d] = v * (E + Ptot) - Bdotv*B;

        F[mom0][d] = P[density] * v * vl[0] - B*Bl[0];
        F[mom1][d] = P[density] * v * vl[1] - B*Bl[1];
        F[mom2][d] = P[density] * v * vl[2] - B*Bl[2];

        F[mom0+d][d] += Ptot;

        /* This can be tidied; the mom'm one too */
        F[B0][d] = Bu[0]*v - B*vu[0]; 
        F[B1][d] = Bu[1]*v - B*vu[1]; 
        F[B2][d] = Bu[2]*v - B*vu[2]; 

        F[B0+d][d] = P[psi] * ((DiagonalSpatialMetric*)metric)->ginv[d][mem]; // Overwrite the above

        F[psi][d] = ch_sq * B; 
    }

    return;
}


inline void MHD::Fluxes_divB_subsystem(const real_t* const __restrict__ P, 
                                             real_t (*__restrict__ F)[3],
                                       const int mem) const
{
    for (int field = 0; field < Nfield; ++field)
        for (int d: dirs)
            F[field][d] = 0.0;

    for (int d: dirs)
    {
        F[B0+d][d] = P[psi] * ((DiagonalSpatialMetric*)metric)->ginv[d][mem]; 
        F[psi][d]  = ch_sq * P[B0+d]; // P holds contravariant components
    }

    return;
}


inline void MHD::DiffusiveFluxes(const real_t* const __restrict__   P, 
                                 const real_t (* const __restrict__ dP)[3],
                                       real_t (*__restrict__ F)[3],
                                 const int mem) const
{
    const real_t mu = P[density] * viscosity; // mu = dynamic viscosity
    const real_t lambda = - (2.0/3.0) * mu;   // from Stokes hypothesis
    
    real_t v[3];
    real_t tau[3][3];
    real_t dv[3][3]; // dv[component][deriv dir] = d_deriv v^component
    real_t rdetg_deriv[3]; // (1/rdetg)*d_j(rdetg) at a single point
    real_t divv;
    real_t duB[3][3]; // duB[j][i] = d^i B^j -- 1st index raised, for Etot flux
    real_t Bl[3];     // need B_i for Etot flux
    real_t JxB;       // for Etot flux

    int dp1, dp2; // cyclic addition d+1, d+2 -- used in Etot flux
    

    for (int d: dirs)
        v[d] = P[v0+d];

    for (int comp: dirs)
        for (int deriv: dirs)
            dv[comp][deriv] = dP[v0+comp][deriv]; // dv[i][j] = d_j v^i
            
    for (int d: dirs)
        rdetg_deriv[d] = metric->rdetg_deriv[d][mem];

    /* div(v) = d_i v^i + v^j d_j(rdetg) / rdetg */
    divv = dv[0][0] + dv[1][1] + dv[2][2] + 
           v[0] * rdetg_deriv[0] + v[1] * rdetg_deriv[1] + v[2] * rdetg_deriv[2];

    for (int d: dirs)
        tau[d][d] = 2*mu*dv[d][d] + lambda*divv;

    tau[0][1] = tau[1][0] = mu * (dv[0][1] + dv[1][0]);
    tau[0][2] = tau[2][0] = mu * (dv[0][2] + dv[2][0]);
    tau[1][2] = tau[2][1] = mu * (dv[1][2] + dv[2][1]);

    /* Raise index on the derivative operator to form duB[j][i] = d^i B^j */
    for (int d: dirs)
        metric->raise(dP[B0+d], duB[d], mem);
    
    metric->lower(&P[B0], Bl, mem); // P[Bi] = B^i 

    for (int d: dirs) // flux direction
    {
        /* Hydro part */
        F[Density][d] = 0.0;

        for (int i: dirs)
            F[mom0+i][d] = - tau[d][i];

        /* Shouldn't this have a contribution from resistivity too? 
         * Yes, should have an eta*JxB contribution from the resistive part 
         * of the electric field */
        /* F(tot E)^d = - v^i tau_i^d + eta (J x B)^d
         * Contracting on tau's lower index, so don't need any metric terms */
        dp1 = dir_plus_one[d];
        dp2 = dir_plus_two[d];

        /* This should be (j x B)^d */
        JxB = Bl[dp2]*(duB[d][dp2] - duB[dp2][d]) - Bl[dp1]*(duB[dp1][d] - duB[d][dp1]);

        F[tot_energy][d] = - (v[0]*tau[0][d] + v[1]*tau[1][d] + v[2]*tau[2][d]) 
                           + resistivity * JxB;

        /* Magnetic part - resistivity is really a magnetic diffusivity */
        F[B0][d] = - resistivity * (dP[B0][d] - dP[B0+d][0]); 
        F[B1][d] = - resistivity * (dP[B1][d] - dP[B0+d][1]); 
        F[B2][d] = - resistivity * (dP[B2][d] - dP[B0+d][2]); 

        F[B0+d][d] = 0.0; // Overwrite the above

        F[psi][d] = 0.0;
        //F[psi][d] = - 0.1 * dP[psi][d];
    }

    return;
}


inline void MHD::OrthonormaliseVectors(real_t* const P, const int mem) const
{
    metric->orthonormals(&P[v0], mem);
    metric->orthonormals(&P[B0], mem);

    return;
}

/****************************************************************/


#if 0
class SystemData_mhd : public SystemData
{
    public:

    SystemData_mhd()
    {
        Nfield = 9;
        variables = new string [9];
        variables[0] = "rho";
        variables[1] = "v0";
        variables[2] = "v1";
        variables[3] = "v2";
        variables[4] = "p";
        variables[5] = "B0";
        variables[6] = "B1";
        variables[7] = "B2";
        variables[8] = "psi";

        diffusive = true;
        viscosity = 3e-3;
        resistivity = 1e-3;

        /* Divergence-cleaning parameters */
        psi_damping_const = 0.01; // 0 < p_d_const < 1
    }
};

class UtoP_mhd : public ConservedToPrimitive
{
    public:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure, B0, B1, B2, psi};

    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ U, 
                                         real_t* const __restrict__ P) const
    {
        P[density] = U[density]; 

        P[v0] = U[mom0] / U[density];
        P[v1] = U[mom1] / U[density];
        P[v2] = U[mom2] / U[density];

        P[B0] = U[B0];
        P[B1] = U[B1];
        P[B2] = U[B2];

        P[psi] = U[psi];

        const real_t KE_density  = 0.5 * P[density] 
                                   * (P[v0]*P[v0] + P[v1]*P[v1] + P[v2]*P[v2]);
        const real_t mag_density = 0.5 * (P[B0]*P[B0] + P[B1]*P[B1] + P[B2]*P[B2]);

        P[pressure] = gm1_mhd * (U[tot_energy] - KE_density - mag_density);

        return;
    }
};


class WaveSpeeds_mhd : public WaveSpeedsFromPrimitive
{
    public:
    enum primitive {density, v0, v1, v2, pressure, B0, B1, B2, psi};

    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ P, 
                                         real_t* const __restrict__ c,
                                   const int dir) const 
    {
        /* asq = (sound speed)^2 */
        const real_t asq = gamma_mhd * P[pressure] / P[density];

        const real_t bsq    = (P[B0]*P[B0]+P[B1]*P[B1]+P[B2]*P[B2])/P[density]; 
        const real_t bdirsq = P[B0+dir]*P[B0+dir]/P[density]; 

        const real_t factor = std::sqrt(std::pow(asq+bsq,2.0) - 4.0*asq*bdirsq);
        const real_t fast_speed = std::sqrt(0.5*(asq + bsq + factor));

        c[0] = MAX(0.0, P[v0+dir] + fast_speed);
        c[1] = MIN(0.0, P[v0+dir] - fast_speed);

        return;
    }
};


class Fluxes_mhd : public FluxesFromPrimitive
{
    public:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure, B0, B1, B2, psi};

    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ P, 
                                         real_t (*__restrict__ F)[3]) const
    {
        const real_t KE_density = 0.5 * P[density] 
                                   * (P[v0]*P[v0] + P[v1]*P[v1] + P[v2]*P[v2]);
        const real_t Bsq = P[B0]*P[B0] + P[B1]*P[B1] + P[B2]*P[B2];
        const real_t E = KE_density + 0.5*Bsq + P[pressure] / gm1_mhd;
        const real_t Ptot = P[pressure] + 0.5*Bsq;
        const real_t Bdotv = P[B0]*P[v0] + P[B1]*P[v1] + P[B2]*P[v2];

        real_t v;
        real_t B;

        for (int d: dirs)
        {
            v = P[v0+d]; // velocity in this direction
            B = P[B0+d]; // B in this direction

            F[density][d]    = v * P[density];
            F[tot_energy][d] = v * (E + Ptot) - Bdotv*B;

            F[mom0][d] = P[density] * v * P[v0] - B*P[B0];
            F[mom1][d] = P[density] * v * P[v1] - B*P[B1];
            F[mom2][d] = P[density] * v * P[v2] - B*P[B2];

            F[mom0+d][d] += Ptot;

            /* This can be tidied; the mom'm one too */
            F[B0][d] = P[B0]*v - B*P[v0]; 
            F[B1][d] = P[B1]*v - B*P[v1]; 
            F[B2][d] = P[B2]*v - B*P[v2]; 

            F[B0+d][d] = P[psi]; // Overwrite the above

            F[psi][d] = ch_sq * B; 
        }

        return;
    }
};


class DiffusiveFluxes_mhd : public DiffusiveFluxes
{
    public:
    enum conserved {Density, mom0, mom1, mom2, tot_energy, B0, B1, B2, psi};

    ACCEL_DECORATOR
    inline virtual void operator()(const real_t* const __restrict__ U, 
                                   const real_t (* const __restrict__ dU)[3],
                                         real_t (*__restrict__ F)[3],
                                   const real_t* const __restrict__ args) const 
    {
        const real_t mu = U[Density] * args[0]; // dynamic viscosity
        const real_t lambda = - (2.0/3.0) * mu; // from Stokes hypothesis
        const real_t eta = args[1]; // magnetic diffusivity
        
        real_t v[3];
        real_t tau[3][3];
        real_t dv[3][3]; // dv[component][deriv dir]
        real_t divv;

        for (int d: dirs)
            v[d] = U[mom0+d] / U[Density];

        for (int comp: dirs)
            for (int deriv: dirs)
                dv[comp][deriv] = (dU[mom0+comp][deriv] - v[comp]*dU[Density][deriv])
                                                                             / U[Density];
        divv = dv[0][0] + dv[1][1] + dv[2][2]; // Assuming Cartesian...

        for (int d: dirs)
            tau[d][d] = 2*mu*dv[d][d] + lambda*divv;

        tau[0][1] = tau[1][0] = mu * (dv[0][1] + dv[1][0]);
        tau[0][2] = tau[2][0] = mu * (dv[0][2] + dv[2][0]);
        tau[1][2] = tau[2][1] = mu * (dv[1][2] + dv[2][1]);

        for (int d: dirs) // flux direction
        {
            /* Hydro part */
            F[Density][d] = 0.0;

            for (int i: dirs)
                F[mom0+i][d] = - tau[d][i];

            F[tot_energy][d] = - (v[0]*tau[0][d] + v[1]*tau[1][d] + v[2]*tau[2][d]);

            /* Magnetic part */
            F[B0][d] = - eta * (dU[B0][d] - dU[B0+d][0]); 
            F[B1][d] = - eta * (dU[B1][d] - dU[B0+d][1]); 
            F[B2][d] = - eta * (dU[B2][d] - dU[B0+d][2]); 

            F[B0+d][d] = 0.0; // Overwrite the above

            F[psi][d] = 0.0;
        }

        return;
    }
};
#endif
#endif
