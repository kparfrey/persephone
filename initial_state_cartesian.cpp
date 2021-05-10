#include "initial_state_cartesian.hpp"

#include <cmath>
#include "element_block.hpp"
#include "physics_includes.hpp"


static real_t scalar_function(const real_t r[3])
{
    real_t sigma[3] = {0.2, 0.2, 0.2};
    real_t r0[3]    = {0.0, 0.0, 0.0};

    real_t rf;
    real_t f;
    real_t arg = 0.0;
    for (int i: dirs)
    {
        rf   = (r[i] - r0[i]) / sigma[i];
        arg += 0.5 * rf * rf ;
    }
    
    f = std::exp(-arg);

    return f;
}


/* From Sun, Wang, Liu (2007), CiCP */
static void shu_vortex(const real_t r[3],
                             real_t* const __restrict__ fields,
                       const int loc0,
                       const int Ns_block)
{
    enum conserved {Density, mom0, mom1, mom2, energy};
    real_t density, v0, v1, v2, pressure; // primitive variables

    const real_t eps = 5.0;
    const real_t rsq = r[0]*r[0] + r[1]*r[1];

    /* Velocities */
    const real_t dv = 0.5*eps*one_pi*std::exp(0.5*(1.0-rsq));
    v0 = 1.0 - dv * r[1];
    v1 = 1.0 + dv * r[0];
    v2 = 0.0;

    const real_t specific_KE = 0.5 * (v0*v0 + v1*v1 + v2*v2);

    /* Density and pressure */
    const real_t dT = - gm1_navstokes * eps*eps * std::exp(1.0-rsq) 
                                        / (8 * gamma_navstokes * pi*pi);
    density  = std::pow(1.0 + dT, 1.0/gm1_navstokes);
    pressure = std::pow(density, gamma_navstokes); // Since entropy = 1

    /* Convert to conserved variables */
    fields[loc0 + Density*Ns_block] = density;
    fields[loc0 +    mom0*Ns_block] = density * v0;
    fields[loc0 +    mom1*Ns_block] = density * v1;
    fields[loc0 +    mom2*Ns_block] = density * v2;
    fields[loc0 +  energy*Ns_block] = density * specific_KE + pressure / gm1_navstokes;

    return;
}


/* From Toth (2000), section 6.3.1 */
static void alfven_wave(const real_t r[3],
                              real_t* const __restrict__ fields,
                        const int loc0,
                        const int Ns_block)
{
    enum conserved {Density, mom0, mom1, mom2, energy, B0, B1, B2, psi};
    real_t density, v0, v1, v2, pressure, b0, b1, b2; // primitive variables

    real_t x, y, alpha, cosa, sina;
    real_t bpar, bperp, vpar, vperp;
    real_t kinetic_energy, magnetic_energy;

    density  = 1.0;
    pressure = 0.1;

    alpha = pi / 6.0;
    cosa  = std::cos(alpha);
    sina  = std::sin(alpha);
    x = r[0];
    y = r[1];

    vpar = 0.0; // vpar = 1 gives standing wave
    bpar = 1.0;

    vperp = 0.1 * std::sin(2*pi*(x*cosa + y*sina));
    bperp = vperp;
    v2    = 0.1 * std::cos(2*pi*(x*cosa + y*sina));
    b2    = v2;

    v0 = vpar*cosa - vperp*sina;
    v1 = vpar*sina + vperp*cosa;
    b0 = bpar*cosa - bperp*sina;
    b1 = bpar*sina + bperp*cosa;

    kinetic_energy  = 0.5 * density * (v0*v0 + v1*v1 + v2*v2);
    magnetic_energy = 0.5 * (b0*b0 + b1*b1 + b2*b2);

    /* Convert to conserved variables */
    fields[loc0 + Density*Ns_block] = density;
    fields[loc0 +    mom0*Ns_block] = density * v0;
    fields[loc0 +    mom1*Ns_block] = density * v1;
    fields[loc0 +    mom2*Ns_block] = density * v2;
    fields[loc0 +  energy*Ns_block] = kinetic_energy + magnetic_energy + pressure / gm1_mhd;
    fields[loc0 +      B0*Ns_block] = b0;
    fields[loc0 +      B1*Ns_block] = b1;
    fields[loc0 +      B2*Ns_block] = b2;
    fields[loc0 +     psi*Ns_block] = 0.0;

    return;
}


void set_initial_state_cartesian(ElementBlock& eb, EqnSystem system)
{
    int mem_offset;
    int loc0; // Memory location for the 0th field
    real_t r[3];

    for (int ke = 0; ke < eb.Nelem[2]; ++ke)
    for (int je = 0; je < eb.Nelem[1]; ++je)
    for (int ie = 0; ie < eb.Nelem[0]; ++ie)
    {
        mem_offset = eb.id_elem(ie, je, ke) * eb.Ns_elem;

        for (int k = 0; k < eb.Ns[2]; ++k)
        for (int j = 0; j < eb.Ns[1]; ++j)
        for (int i = 0; i < eb.Ns[0]; ++i)
        {
            loc0 = eb.ids(i,j,k) + mem_offset;
            for (int d: dirs)
                r[d] = eb.rs(d,loc0);

            switch(system)
            {
                case scalar_advection:
                    eb.fields[loc0] = scalar_function(r);
                    break;
                case navier_stokes:
                    shu_vortex(r, eb.fields, loc0, eb.Ns_block);
                    break;
                case mhd:
                    alfven_wave(r, eb.fields, loc0, eb.Ns_block);
                    break;
            }
        }
    }

    return;
}
