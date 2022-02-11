#include "desc.hpp"

#include <highfive/H5Easy.hpp>
#include <cmath>
#include <iostream>
//#include <boost/math/special_functions/jacobi.hpp>
#include "jacobi.hpp"
#include "write_screen.hpp"

/*
inline static long int factorial(const int n)
{
    long int fac = 1;

    if (n > 1)
        for (int i = 2; i <=n ; ++i)
            fac *= i;

    return fac;
}
*/


#if 0
static long double zernike_DIY(const real_t rho, const real_t theta, const int l, const int m)
{
    using std::pow;

    long double Z; // Zernike polynomial Z_l^m(rho,theta) --- don't confuse with cylindrical Z
    long double R = 0.0L; // radial function R_l^m(rho)    --- don't confuse with cylindrical radius

    const int mabs = std::abs(m);
    const int s_max = (l - mabs)/2;

    if ( (l-mabs)%2 != 0 )
        exit(31);

    /* Works with this, but blows up whenever l is allowed to exceed 20, even for M < 20 */
    if (l > 20)
        return 0.0;

    int sign;
    long int top, bot;

    /* Direct method */
    /*
    for (int s = 0; s <= s_max; ++s)
    {
        if (s % 2)
            sign = -1;
        else
            sign = 1;

        top = sign * factorial(l-s);
        bot = factorial(s) * factorial((l+mabs)/2 - s) * factorial((l-mabs)/2 - s);
        R += std::pow((long double)rho, (long double)(l-2*s)) * (long double)top/((long double)bot);
    }
    */

    /* Horner's method */
    long double rho_sq = (long double)rho * (long double)rho;
    for (int s = 0; s <= s_max; ++s)
    {
        if (s % 2)
            sign = -1;
        else
            sign = 1;

        top = sign * factorial(l-s);
        bot = factorial(s) * factorial((l+mabs)/2 - s) * factorial((l-mabs)/2 - s);

        R = rho_sq * ((long double)top/((long double)bot) + R);
    }

    R *= std::pow((long double)rho, (long double)(l - (2*s_max+2)));


    if (m >= 0)
        Z = R * std::cos( m * theta);
    else
        Z = R * std::sin(-m * theta);

    return Z;
}
#endif


/* Use the Jacobi polynomial from Boost to find the radial function:
 * https://mathworld.wolfram.com/ZernikePolynomial.html
 * https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/sf_poly/jacobi.html */
double DescConfig::zernike(const real_t rho, const real_t theta, const int l, const int m, deriv_dir deriv) const
{
    double R; // Radial function for Zernike polynomial --- don't confuse with cylindrical R
    double Z; // Zernike polynomial Z_l^m(rho,theta) --- don't confuse with cylindrical Z

    /* Radial function */
    //using boost::math::jacobi;
    const unsigned int mabs = std::abs(m);
    const unsigned int lm2  = (l - mabs) / 2;
    
    double sign = 1.0;
    if (lm2 % 2)
        sign = -1.0;

    const double J = jacobi(lm2, (double)mabs, 0.0, 1.0 - 2.0*rho*rho);

    if (deriv == rho_d) 
    {
        const double Jp = jacobi_prime(lm2, (double)mabs, 0.0, 1.0 - 2.0*rho*rho);
        R = sign * std::pow(rho, mabs) * (mabs * J / rho - 4.0 * rho * Jp);
    }
    else // don't differentiate the radial function - standard Z polynomial
        R = sign * std::pow(rho, mabs) * J;

    /* Add poloidal function */
    if (deriv == theta_d)
    {
        if (m >= 0)
            Z = - R * mabs * std::sin( mabs * theta);
        else
            Z =   R * mabs * std::cos( mabs * theta);
    }
    else // don't differentiate - standard Zernike
    {
        if (m >= 0)
            Z = R * std::cos( mabs * theta);
        else
            Z = R * std::sin( mabs * theta);
    }

    return Z;
}


void DescConfig::surface_polynomial_expansion(real_t& f, const real_t r[3], const std::vector<double>& f_lmn, 
                                              const std::vector<std::vector<int>>& f_modes,
                                              deriv_dir deriv) const
{
    real_t rho   = r[0]; // i.e. r_uds = sqrt(toroidal flux function)
    real_t theta = r[1]; // The "curly theta" from https://desc-docs.readthedocs.io/en/latest/output.html
    real_t zeta  = r[2]; // toroidal angle

    int l, m, n; // 3D modes indices for each separate mode
    real_t F_toroidal;

    int N = f_lmn.size(); // Total number of modes

    double lsum = 0.0;

    for (int i = 0; i < N; ++i)
    {
        l = f_modes[i][0];
        m = f_modes[i][1];
        n = f_modes[i][2];

        if (deriv == zeta_d)
        {
            if (n >= 0)
                F_toroidal = - n * Nfp * std::sin( n * Nfp * zeta);
            else
                F_toroidal = - n * Nfp * std::cos(-n * Nfp * zeta);
        }
        else // don't differentiate in zeta
        {
            if (n >= 0)
                F_toroidal = std::cos( n * Nfp * zeta);
            else
                F_toroidal = std::sin(-n * Nfp * zeta);
        }

        lsum += f_lmn[i] * zernike(rho, theta, l, m, deriv) * F_toroidal;
    }

    f = (real_t) lsum;

    return;
}


/* Using the flux surface shape for all r, even r=1. In other words, the outer boundary
 * will be the shape of the last flux surface, not the input target boundary shape. 
 * Assuming that r_uds = rho = sqrt(psi/psi_a) --- the UDS radial coord is identified
 * with the square root of the toroidal flux function */
void DescConfig::unit_disc_to_physical_space(real_t r[3]) const
{
    real_t R = 0.0; // Cylindrical coords in physical space
    real_t Z = 0.0;

    deriv_dir deriv = none;

    surface_polynomial_expansion(R, r, R_lmn, R_modes, deriv);
    surface_polynomial_expansion(Z, r, Z_lmn, Z_modes, deriv);

    r[0] = R;
    r[1] = Z;

    return;
}


/* Need to pass in "unit disc space" coords. These are stored as (rho, theta, phi). */
void DescConfig::construct_equilibrium(const real_t r[3], real_t U[9]) const
{
    enum conserved {density, mom0, mom1, mom2, tot_energy, B0, B1, B2, div_scalar};

    const real_t rho = r[0]; // UDS radial coord, rho = sqrt(psi/psi_a) -- sqrt of flux

    real_t gamma = 5.0/3.0; // Want to pipe in from MHD object

    real_t p = 0.0;
    for (int i = 0; i < N_pressure; ++i)
        p += pressure[i] * std::pow(rho, i);

    real_t rot_trans = 0.0;
    for (int i = 0; i < N_iota; ++i)
        rot_trans += iota[i] * std::pow(rho, i);

    U[density] = 1.0; // Uniform density seems reasonable?
    U[mom0]    = 0.0; // Since it's an equilibrium
    U[mom1]    = 0.0;
    U[mom2]    = 0.0;

    /* Construct the magnetic field */
    real_t Btheta, Bzeta; // B in straight-field-line coords, contra. components; B_rho = 0
    real_t BR, BZ; // B in cylindrical coords, contravariant & orthonormal are identical
    real_t R;                            // R surface function --- recalculate
    real_t dR_drho, dR_dtheta, dR_dzeta; // derivatives of surface functions wrt SFL coords
    real_t dZ_drho, dZ_dtheta, dZ_dzeta;

    surface_polynomial_expansion(R, r, R_lmn, R_modes, none);
    surface_polynomial_expansion(dR_drho,   r, R_lmn, R_modes, rho_d);
    surface_polynomial_expansion(dR_dtheta, r, R_lmn, R_modes, theta_d);
    surface_polynomial_expansion(dR_dzeta,  r, R_lmn, R_modes, zeta_d);
    surface_polynomial_expansion(dZ_drho,   r, Z_lmn, Z_modes, rho_d);
    surface_polynomial_expansion(dZ_dtheta, r, Z_lmn, Z_modes, theta_d);
    surface_polynomial_expansion(dZ_dzeta,  r, Z_lmn, Z_modes, zeta_d);

    const real_t rdetg = (dZ_drho * dR_dtheta - dZ_dtheta * dR_drho) * R;
    const real_t prefactor = psi_a * rho / (pi * rdetg);

    Btheta = prefactor * rot_trans;
    Bzeta  = prefactor;

    BR = Btheta * dR_dtheta + Bzeta * dR_dzeta;
    BZ = Btheta * dZ_dtheta + Bzeta * dZ_dzeta;

    U[B0] = BR;
    U[B1] = BZ;
    U[B2] = Bzeta; // U[B2] is the contravariant component = Bphi/R

    const real_t Bsq = BR*BR + BZ*BZ + R*R*Bzeta*Bzeta;

    U[tot_energy] = 0.5 * Bsq + p/(gamma - 1.0);

    U[div_scalar] = 0.0;

    return;
}


DescConfig::DescConfig(std::string input_file, const int iteration)
{
    using std::vector;

    write::message("\nReading DESC data from " + input_file + ", iteration " + std::to_string(iteration));

    H5Easy::File data(input_file, H5Easy::File::ReadOnly);

    std::string base = "/_equilibria/" + std::to_string(iteration) + "/";

    L = H5Easy::load<int>(data, base + "_L");
    M = H5Easy::load<int>(data, base + "_M");
    N = H5Easy::load<int>(data, base + "_N");

    L_lmn = H5Easy::load<vector<double>>(data, base + "_L_lmn");
    R_lmn = H5Easy::load<vector<double>>(data, base + "_R_lmn");
    Z_lmn = H5Easy::load<vector<double>>(data, base + "_Z_lmn");

    L_modes = H5Easy::load<vector<vector<int>>>(data, base + "_L_basis/_modes");
    R_modes = H5Easy::load<vector<vector<int>>>(data, base + "_R_basis/_modes");
    Z_modes = H5Easy::load<vector<vector<int>>>(data, base + "_Z_basis/_modes");

    iota     = H5Easy::load<vector<double>>(data, base + "_pressure/_params");
    pressure = H5Easy::load<vector<double>>(data, base + "_iota/_params");

    Nfp   = H5Easy::load<double>(data, base + "_NFP");
    psi_a = H5Easy::load<double>(data, base + "_Psi");

    N_L = L_lmn.size();
    N_R = R_lmn.size();
    N_Z = Z_lmn.size();

    N_iota     = iota.size();
    N_pressure = pressure.size();

    write::variable<int>("Radial (lambda) polynomial resolution    - L ", L);
    write::variable<int>("Surface polynomials, poloidal resolution - M ", M);
    write::variable<int>("Surface polynomials, toroidal resolution - N ", N);
    write::variable<int>("No. of modes in lambda polynomial ", N_L);
    write::variable<int>("No. of modes in R polynomial      ", N_R);
    write::variable<int>("No. of modes in Z polynomial      ", N_Z);
    write::variable<int>("No. of modes in iota polynomial   ", N_iota);
    write::variable<int>("No. of modes in pressure poly.    ", N_pressure);
    write::variable<double>("No. of field periods ", Nfp);

    write::message("Finished loading input data\n");

    return;
}
