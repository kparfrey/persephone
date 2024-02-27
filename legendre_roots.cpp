#include "legendre_roots.hpp"
#include <cmath>

namespace legendre
{
    /* This n is the number of internal points needed, and the order of the Legendre
     * polynomial whose roots are found; n = N_soln - 1. 
     * x[] will hold all the flux points, so x[0] = 0, x[-1] = 1
     * See Numerical Recipes (3rd Ed) Sec 4.6.1.         */
    void find_roots(real_t x[], const int n)
    {
        real_t z, z1, p1, p2, p3, pp;

        const real_t eps = 1e-15;

        /* Lower and upper bounds of the desired interval - the roots returned
         * will be entirely internal to this interval, i.e. x_root != x1, x2 */
        const real_t x1 = 0.0; 
        const real_t x2 = 1.0; 

        const real_t xm = 0.5 * (x2 + x1);
        const real_t xl = 0.5 * (x2 - x1);

        const int m = (n + 1) / 2;

        for (int i = 0; i < m; i++)
        {
            z = std::cos(pi * (i + 0.75)/(n + 0.5));

            do
            {
                p1 = 1.0;
                p2 = 0.0;
                for (int j = 0; j < n; j++)
                {
                    p3 = p2;
                    p2 = p1;
                    p1 = ((2.0*j + 1.0) * z * p2 - j * p3) / (j + 1.0);
                }
                pp = n * (z * p1 - p2)/(z*z - 1.0);
                z1 = z;
                z = z1 - p1/pp;
            } while (std::abs(z - z1) > eps);

            x[i + 1] = xm - xl * z; // leave x[0] for left-edge fluxpoint
            x[n - i] = xm + xl * z; // and shift this up one memory slot too
        }

        return;
    }
}
