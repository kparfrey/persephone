#include "cerfon_freidberg.hpp"

#include <cmath>
#include "write_screen.hpp"

static real_t psi_fn(const real_t x, const real_t y, const real_t A, const real_t* const c)
{
    const real_t x2 = x * x;
    const real_t x4 = x2 * x2;
    const real_t x6 = x2 * x4;

    const real_t y2 = y * y;
    const real_t y4 = y2 * y2;
    const real_t y6 = y2 * y4;

    const real_t lnx = std::log(x);

    const real_t psi_p = 0.5 * A * x2 * lnx  + 0.125 * (1-A) * x4;

    real_t psi[8] = {}; // the psi_i

    psi[1] = 1.0;
    psi[2] = x2;
    psi[3] = y2 - x2 * lnx; // Antoine's code has the negative of this: might need to reverse coeff sign
    psi[4] = x4 - 4 * x2 * y2;
    psi[5] = 2 * y4 - 9 * y2 * x2 + 3 * x4 * lnx - 12 * x2 * y2 * lnx;
    psi[6] = x6 - 12 * x4 * y2 + 8 * x2 * y4;
    psi[7] = 8 * y6 - 140 * y4 * x2 + 75 * y2 * x4 - 15 * x6 * lnx  
             + 180 * x4 * y2 * lnx - 120 * x2 * y4 * lnx;

    real_t psi_h = 0.0;
    for (int i = 1; i < 8; ++i)
        psi_h += c[i] * psi[i];

    real_t psi_tot = psi_p + psi_h;

    return psi_tot;
}


static real_t dpsi_dx(const real_t x, const real_t y, const real_t A, const real_t* const c)
{
    const real_t x2 = x * x;
    const real_t x3 = x * x2;
    const real_t x5 = x2 * x3;

    const real_t y2 = y * y;
    const real_t y4 = y2 * y2;

    const real_t lnx = std::log(x);

    /* psi_p = d(psi_p)/dx */
    const real_t psi_p = A * x * lnx + 0.5 * A * x  + 0.5 * (1-A) * x3;

    /* psi[i] = d(psi_i)/dx */
    real_t psi[8] = {}; 

    psi[2] = 2 * x;
    psi[3] = - 2 * x * lnx - x; 
    psi[4] = 3 * x3 - 8 * x * y2;
    psi[5] = - 18 * y2 * x + 12 * x3 * lnx + 3 * x3 - 24 * x * y2 * lnx - 12 * x * y2;
    psi[6] = 6 * x5 - 48 * x3 * y2 + 16 * x * y4;
    psi[7] = - 280 * y4 * x + 300 * y2 * x3 - 90 * x5 * lnx - 15 * x5 
             + 720 * x3 * y2 * lnx + 180 * x3 * y2 - 240 * x * y4 * lnx - 120 * x * y4;

    real_t psi_h = 0.0;
    for (int i = 2; i < 8; ++i)
        psi_h += c[i] * psi[i];

    real_t psi_tot = psi_p + psi_h;

    return psi_tot;
}


static real_t dpsi_dy(const real_t x, const real_t y, const real_t A, const real_t* const c)
{
    const real_t x2 = x * x;
    const real_t x4 = x2 * x2;

    const real_t y2 = y * y;
    const real_t y3 = y * y2;
    const real_t y5 = y2 * y3;

    const real_t lnx = std::log(x);

    /* psi_p = d(psi_p)/dy */
    const real_t psi_p = 0.0;

    /* psi[i] = d(psi_i)/dy */
    real_t psi[8] = {}; 

    psi[3] = 2 * y;
    psi[4] = - 8 * x2 * y;
    psi[5] = 8 * y3 - 18 * y * x2 - 24 * x2 * y * lnx;
    psi[6] = - 24 * x4 * y + 32 * x2 * y3;
    psi[7] = 48 * y5 - 560 * y3 * x2 + 150 * y * x4  
             + 360 * x4 * y * lnx - 480 * x2 * y3 * lnx;

    real_t psi_h = 0.0;
    for (int i = 3; i < 8; ++i)
        psi_h += c[i] * psi[i];

    real_t psi_tot = psi_p + psi_h;

    return psi_tot;
}


void CerfonFreidbergConfig::construct_equilibrium(const real_t r[3],
                                                        real_t U[9],
                                                  const real_t gamma)
{
    enum conserved {Density, mom0, mom1, mom2, tot_energy, B0, B1, B2, div_scalar};
    
    const real_t R = r[0]; // equiv to x
    const real_t Z = r[1]; // equiv to y

    const real_t psi = psi_fn(R, Z, A, c);

    const real_t BR   = - dpsi_dy(R, Z, A, c) / R; // orthonormal components
    const real_t BZ   =   dpsi_dx(R, Z, A, c) / R;
    const real_t Bphi = (B0*B0 - 2 * A * psi) / (R*R);

    const real_t p = (A - 1) * psi;

    const real_t Bsq = BR*BR + BZ*BZ + Bphi*Bphi;

    U[Density] = 1.0; // Uniform density seems reasonable?
    U[mom0]    = 0.0; // Since it's an equilibrium
    U[mom1]    = 0.0;
    U[mom2]    = 0.0;

    U[tot_energy] = 0.5 * Bsq + p/(gamma - 1.0);

    U[B0] = BR;
    U[B1] = BZ;
    U[B2] = Bphi / R; // U[B2] is the contravariant component

    U[div_scalar] = 0.0;

    return;
}


CerfonFreidbergConfig::CerfonFreidbergConfig()
{
    /* Set all the coefficients, found using Antoine's matlab script
     * Using c[1] for c_1 etc. */

    switch(machine)
    {
        case tftr:
            epsilon = 0.87 / 2.5; // R0 = 2.5
            kappa   = 1.0;
            delta   = 0.0;
            c[1] =  3.363214222958483e-02;
            c[2] = -1.288496594601578e-01;
            c[3] = -5.920909679287922e-02;
            c[4] = -5.823137092324818e-02;
            c[5] =  6.696122694686874e-03;
            c[6] = -1.479830335643033e-03;
            c[7] = -4.613291591508187e-05;
            break;
        case nstx:
            epsilon = 0.78;
            kappa   = 2.0;
            delta   = 0.35;
            break;
        case spheromak:
            epsilon = 0.95;
            kappa   = 1.0;
            delta   = 0.2;
            break;
        case frc:
            epsilon = 0.99;
            kappa   = 10.0;
            delta   = 0.7;
            break;
        default:
            write::error("Chosen machine not defined yet...", destroy);
            break;
    }

    alpha = std::asin(delta);

    return;
}
