#ifndef MHD_HPP
#define MHD_HPP

#include <cmath>
#include "physics.hpp"

/* U --- conserved
 * 0 rho
 * 1 rho v_0
 * 2 rho v_1
 * 3 rho v_2
 * 4 E  ---  total energy density
 * 5 B_0
 * 6 B_1
 * 7 B_2
 * 8 psi --- div-cleaning scalar field
 *
 * P --- primitive
 * 0 rho
 * 1 v_0
 * 2 v_1
 * 3 v_2
 * 4 p  ---  pressure
 * 5 B_0
 * 6 B_1
 * 7 B_2
 * 8 psi
 */

class MHD : public Physics
{
    private:
    enum conserved {Density, mom0, mom1, mom2, tot_energy};
    enum primitive {density, v0  , v1  , v2,   pressure, B0, B1, B2, psi};

    public:
    const real_t gamma = 5.0/3.0;
    const real_t gm1 = gamma - 1.0;


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
        viscosity   = 3e-3;
        resistivity = 1e-3;

        /* Divergence-cleaning parameters */
        psi_damping_const = 0.01; // 0 < p_d_const < 1
    }


    void ConservedToPrimitive(const real_t* const __restrict__ U, 
                                    real_t* const __restrict__ P) const override;

    void WaveSpeeds(const real_t* const __restrict__ P, 
                          real_t* const __restrict__ c,
                    const int dir) const override;

    void Fluxes(const real_t* const __restrict__ P, 
                      real_t (*__restrict__ F)[3]) const override;

    void DiffusiveFluxes(const real_t* const __restrict__ U, 
                         const real_t (* const __restrict__ dU)[3],
                               real_t (*__restrict__ F)[3]) const override;
                         //const real_t* const __restrict__ args) const override;
};


inline void MHD::ConservedToPrimitive(const real_t* const __restrict__ U, 
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

    P[pressure] = gm1 * (U[tot_energy] - KE_density - mag_density);

    return;
}


inline void MHD::WaveSpeeds(const real_t* const __restrict__ P, 
                                  real_t* const __restrict__ c,
                            const int dir) const 
{
    /* asq = (sound speed)^2 */
    const real_t asq = gamma * P[pressure] / P[density];

    const real_t bsq    = (P[B0]*P[B0]+P[B1]*P[B1]+P[B2]*P[B2])/P[density]; 
    const real_t bdirsq = P[B0+dir]*P[B0+dir]/P[density]; 

    const real_t factor = std::sqrt(std::pow(asq+bsq,2.0) - 4.0*asq*bdirsq);
    const real_t fast_speed = std::sqrt(0.5*(asq + bsq + factor));

    c[0] = MAX(0.0, P[v0+dir] + fast_speed);
    c[1] = MIN(0.0, P[v0+dir] - fast_speed);

    return;
}


inline void MHD::Fluxes(const real_t* const __restrict__ P, 
                              real_t (*__restrict__ F)[3]) const
{
    const real_t KE_density = 0.5 * P[density] 
                               * (P[v0]*P[v0] + P[v1]*P[v1] + P[v2]*P[v2]);
    const real_t Bsq = P[B0]*P[B0] + P[B1]*P[B1] + P[B2]*P[B2];
    const real_t E = KE_density + 0.5*Bsq + P[pressure] / gm1;
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


inline void MHD::DiffusiveFluxes(const real_t* const __restrict__ U, 
                                 const real_t (* const __restrict__ dU)[3],
                                       real_t (*__restrict__ F)[3]) const
                                 //const real_t* const __restrict__ args) const 
{
    const real_t mu = U[Density] * viscosity; // mu = dynamic viscosity
    const real_t lambda = - (2.0/3.0) * mu; // from Stokes hypothesis
    //const real_t eta = args[1]; // magnetic diffusivity
    
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

        /* Magnetic part - resistivity is really a magnetic diffusivity */
        F[B0][d] = - resistivity * (dU[B0][d] - dU[B0+d][0]); 
        F[B1][d] = - resistivity * (dU[B1][d] - dU[B0+d][1]); 
        F[B2][d] = - resistivity * (dU[B2][d] - dU[B0+d][2]); 

        F[B0+d][d] = 0.0; // Overwrite the above

        F[psi][d] = 0.0;
    }

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
