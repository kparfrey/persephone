#ifndef ELEMENT_BRICK_HPP
#define ELEMENT_BRICK_HPP

#include "common.hpp"




class ElementBrick
{
    public:

    /* Data */
    int Ns[3];  // No. soln points in each element, in each direction
    int Nf[3];  // No. flux points "   "      "     "   "       "
    int Ns_tot; // Total no. of soln points in each element
    int Nf_tot; // Total no. of flux points "   "     "
    int Nv;     // No. of field variables being evolved


    /* Methods */
    inline int sidx(int var, int i, int j, int k); // soln-point index map
    void allocate();
    void free();
};


//ElementBrick::ElementBrick(LocalParams &lp)
//{
//}


/*** C-order, i slowest index, k fastest ***/
inline int ElementBrick::sidx(int var, int i, int j, int k)
{
    return Ns[0]*Ns[1]*Ns[2]*var + (i*Ns[1] + j)*Ns[2] + k;
}


// Should this be "setup", or even run as the constructor?
void ElementBrick::allocate()
{
    return;
}


void ElementBrick::free()
{
    return;
}
#endif
