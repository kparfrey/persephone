#include "cerfon_freidberg.hpp"

#include <cmath>
#include "write_screen.hpp"


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
        default:
            write::error("Chosen machine not defined yet...", destroy);
            break;
    }

    return;
}
