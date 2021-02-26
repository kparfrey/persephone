#ifndef TORUS_MODE_PACK_HPP
#define TORUS_MODE_PACK_HPP

#include "common.hpp"

/* Stores the R and Z mode coefficients for doing the  
 * transformation fromunit-disc space to physical space */    
class BoundaryModePack
{
    public:
    int Nm; // No. of poloidal modes
    int Nk; // No. of toroidal modes
    
    real_t **Rmk; // R modes: R ~ r Sum Rmk[m][k] Cos(m theta) Cos(k phi)
    real_t **Zmk; // Z modes: Z ~ r Sum Zmk[m][k] Sin(m theta) Cos(k phi)

    /* Use a templated constructor to read off Nm and Nk from the arrays */
    template <int Nm, int Nk>
    BoundaryModePack(real_t (&Rmk)[Nm][Nk], real_t (&Zmk)[Nm][Nk])
    {
        this->Nm = Nm;
        this->Nk = Nk;

        setup();

        for (int m = 0; m < Nm; ++m)
        for (int k = 0; k < Nk; ++k)
        {
           this->Rmk[m][k] = Rmk[m][k]; 
           this->Zmk[m][k] = Zmk[m][k]; 
        }
    }

    void setup();
};


void BoundaryModePack::setup()
{
    Rmk = new real_t* [Nm];
    Zmk = new real_t* [Nm];

    for (int i = 0; i < Nm; ++i)
    {
        Rmk[i] = new real_t [Nk];
        Zmk[i] = new real_t [Nk];
    }

    return;
}
#endif
