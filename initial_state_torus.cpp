#include "initial_state_torus.hpp"

#include <cmath>
#include "element_block.hpp"
#include "physics_includes.hpp"
#include "cerfon_freidberg.hpp"


static void implosion(const real_t r[3],
                            real_t* const __restrict__ fields,
                      const int loc0,
                      const int Ns_block,
                      const NavierStokes* const __restrict__ physics)
{
    enum conserved {Density, mom0, mom1, mom2, energy};
    real_t density, v0, v1, v2, pressure; // primitive variables

    const real_t radial = std::sqrt(r[0]*r[0] + r[1]*r[1]);

    /* Velocities */
    v0 = 0.0; 
    v1 = 0.0;
    v2 = 0.0;
    
    const real_t specific_KE = 0.5 * (v0*v0 + v1*v1 + v2*v2);

    /* Density and pressure */
    const real_t rramp  = 0.6;
    const real_t density_floor = 0.01;
    real_t reff = 1.0;
    if (radial < rramp)
        reff = radial / rramp;

    density  = density_floor + 0.5*(1.0 - std::cos(reff * pi));
    pressure = std::pow(density, physics->gamma); // Since entropy = 1

    /* Convert to conserved variables */
    fields[loc0 + Density*Ns_block] = density;
    fields[loc0 +    mom0*Ns_block] = density * v0;
    fields[loc0 +    mom1*Ns_block] = density * v1;
    fields[loc0 +    mom2*Ns_block] = density * v2;
    fields[loc0 +  energy*Ns_block] = density * specific_KE + pressure / physics->gm1;

    return;
}


void set_euler_torus(ElementBlock& eb)
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

            implosion(r, eb.fields, loc0, eb.Ns_block, (NavierStokes*)eb.physics_soln);
        }
    }

    return;
}

#if 0
void set_CerfonFreidberg(ElementBlock& eb, CerfonFreidbergConfig& cf_config)
{
    real_t r[3];
    real_t U[9]; // Conserved variables at a point

    cf_config.gamma = ((MHD*)eb.physics_soln)->gamma;

    /* loc0 is memory location for the 0th field */
    for (int loc0 = 0; loc0 < eb.Ns_block; ++loc0)
    {
        for (int d: dirs)
            r[d] = eb.rs(d,loc0);

        cf_config.construct_equilibrium(r, U);
    
        for (int i = 0; i < eb.Nfield; ++i)
            eb.fields[loc0 + i*eb.Ns_block] = U[i];
    }

    return;
}
#endif


void set_torus_initial_state(ElementBlock& eb, TorusConfig& config)
{
    real_t r_uds[3];
    real_t r_phys[3];
    real_t U[9]; // Conserved variables at a point

    config.gamma = ((MHD*)eb.physics_soln)->gamma;

    /* loc0 is memory location for the 0th field */
    for (int loc0 = 0; loc0 < eb.Ns_block; ++loc0)
    {
        for (int d: dirs)
        {
            r_uds[d]  = eb.rs_pre_transform(d,loc0);
            r_phys[d] = eb.rs(d,loc0);
        }

        config.construct_equilibrium(r_uds, r_phys, U);
    
        for (int i = 0; i < eb.Nfield; ++i)
            eb.fields[loc0 + i*eb.Ns_block] = U[i];
    }

    return;
}


/* Just a placeholder, used when testing grid construction */
void set_uniform(ElementBlock& eb)
{
    for (int loc0 = 0; loc0 < eb.Ns_block; ++loc0)
    {
        for (int i = 0; i < eb.Nfield; ++i)
            eb.fields[loc0 + i*eb.Ns_block] = 1.0;
    }

    return;
}
