#include "cerfon_freidberg.hpp"

#include <cmath>
#include "write_screen.hpp"

static real_t psi(const real_t x, const real_t y, const real_t A, const real_t* const c)
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
    psi[3] = y2 - x2 * lnx;
    psi[4] = x4 - 4 * x2 * y2;
    psi[5] = 2 * y4 - 9 * y2 * x2 + 3 * x4 * lnx - 12 * x2 * y2 * lnx;
    psi[6] = x6 - 12 * x4 * y2 + 8 * x2 * y4;
    psi[7] = 8 * y6 - 140 * y4 * x2 + 75 * y2 * x4 - 15 * x6 * lnx + 
             180 * x4 * y2 * lnx - 120 * x2 * y4 * lnx;

    real_t psi_h = 0.0;
    for (int i = 1; i < 8; ++i)
        psi_h += c[i] * psi[i];

    real_t psi_tot = psi_p + psi_h;

    return psi_tot;
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
