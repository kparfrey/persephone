#include "cerfon_freidberg.hpp"

#include <cmath>
#include "write_screen.hpp"
#include <iostream>

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
    psi[4] = 4 * x3 - 8 * x * y2;
    psi[5] = - 30 * y2 * x + 12 * x3 * lnx + 3 * x3 - 24 * x * y2 * lnx;
    psi[6] = 6 * x5 - 48 * x3 * y2 + 16 * x * y4;
    psi[7] = - 400 * y4 * x + 480 * y2 * x3 - 90 * x5 * lnx - 15 * x5 
             + 720 * x3 * y2 * lnx - 240 * x * y4 * lnx;

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
    /* Note using b0 etc. so don't hide the B0 parameter in this object */
    enum conserved {density, mom0, mom1, mom2, tot_energy, b0, b1, b2, div_scalar};
    
    const real_t R = r[0]; 
    const real_t Z = r[1];
    const real_t R0 = 1.0; // Always --- maybe remove all R0 factors?
    const real_t x = R/R0;
    const real_t y = Z/R0;

    const real_t R0sq = R0 * R0;
    const real_t R04  = R0sq * R0sq;

    const real_t psi = psi_fn(x, y, A, c);

    const real_t BR   = - dpsi_dy(x, y, A, c) / (R0*R); // orthonormal components
    const real_t BZ   =   dpsi_dx(x, y, A, c) / (R0*R);
    const real_t Bphi = (R0/R) * std::sqrt(B0*B0 - 2 * A * psi/R04);

    const real_t p = (A - 1) * psi / R04 + 1e-3;

    const real_t Bsq = BR*BR + BZ*BZ + Bphi*Bphi;

    U[density] = 1.0; // Uniform density seems reasonable?
    U[mom0]    = 0.0; // Since it's an equilibrium
    U[mom1]    = 0.0;
    U[mom2]    = 0.0;

    U[tot_energy] = 0.5 * Bsq + p/(gamma - 1.0);

    U[b0] = BR;
    U[b1] = BZ;
    U[b2] = Bphi / R; // U[B2] is the contravariant component

    U[div_scalar] = 0.0;

    return;
}


CerfonFreidbergConfig::CerfonFreidbergConfig()
{
    /* Set all the coefficients, found using Antoine's matlab script
     * Using c[1] for c_1 etc. 
     * Need to take the negative of his script's c[3], since in the script
     * he defined psi_3 with the opposite sign to the paper. */

    switch(machine)
    {
        case tftr:
            epsilon = 0.87 / 2.5; // R0 = 2.5
            kappa   = 1.0;
            delta   = 0.0;
            c[1] =  3.363214222958483e-02;
            c[2] = -1.288496594601578e-01;
            c[3] =  5.920909679287922e-02;
            c[4] = -5.823137092324818e-02;
            c[5] =  6.696122694686874e-03;
            c[6] = -1.479830335643033e-03;
            c[7] = -4.613291591508187e-05;
            break;
        case nstx:
            epsilon = 0.67/0.86;
            kappa   = 2.0;
            delta   = 0.35;
            c[1] =  1.515626015086087e-02;
            c[2] = -3.219112995316729e-01;
            c[3] =  4.304434935073897e-03;
            c[4] = -2.337438273912475e-02;
            c[5] =  2.790357839469228e-04;
            c[6] = -3.833470165175256e-04;
            c[7] = -3.029439496605592e-06;
            break;
        case spheromak:
            epsilon = 0.95;
            kappa   = 1.0;
            delta   = 0.2;
            c[1] =  5.877704044228160e-04;
            c[2] = -2.560362608900638e-01;
            c[3] =  6.979227356204132e-03;
            c[4] = -6.078983140122448e-02;
            c[5] =  6.243471441125277e-03;
            c[6] = -3.066243397601454e-03;
            c[7] = -9.081518541524376e-05;
            break;
        case frc:
            epsilon = 0.99;
            kappa   = 10.0;
            delta   = 0.7;
            c[1] =  4.874020982193714e-05;
            c[2] = -4.874223413173798e-01;
            c[3] =  1.722686095823331e-06;
            c[4] = -1.904450361046603e-03;
            c[5] = -2.002800367444697e-07;
            c[6] = -3.681071366340895e-06;
            c[7] =  6.291815809141939e-10;
            break;
        case iter:
            epsilon = 2.0/6.2;
            kappa   = 1.7;
            delta   = 0.33;
            c[1] =  8.235805781421381e-02;
            c[2] = -1.932149717154881e-01;
            c[3] = -4.791944197513237e-02;
            c[4] = -4.666382427999496e-02;
            c[5] =  4.741712757718168e-03;
            c[6] = -4.119698821529545e-03;
            c[7] = -1.056698610272738e-04;
            break;
        case jet:
            epsilon = 1.0/3.0;
            kappa   = 1.7;
            delta   = 0.25;
            c[1] =  7.650174912175325e-02;
            c[2] = -2.027143233944073e-01;
            c[3] = -1.847588959046809e-02;
            c[4] = -3.658675468981453e-02;
            c[5] =  1.522777623705320e-03;
            c[6] = -1.782489256943540e-03;
            c[7] = -3.555779071815003e-05;
            break;
        case mast:
            epsilon = 0.65/0.85;
            kappa   = 2.45;
            delta   = 0.5;
            c[1] =  1.805778616813774e-02;
            c[2] = -3.309464564355566e-01;
            c[3] = -7.082532551454974e-04;
            c[4] = -1.972741308898579e-02;
            c[5] =  1.017045415679966e-04;
            c[6] = -3.903514551270571e-04;
            c[7] = -1.853186423352871e-06;
            break;
        default:
            write::error("Chosen machine not defined yet...", destroy);
            break;
    }

    alpha = std::asin(delta);

    return;
}
