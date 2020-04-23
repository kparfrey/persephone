#ifndef PARAMS_HPP
#define PARAMS_HPP

#include "common.hpp"

class Process; // Process type will be used below


/* Base parameter class. Will eventually likely have this, ParamsCartesian, 
 * ParamsToroidal, ParamsSpherical etc., all in separate files.
 */
class Params
{
    public:

    /* Total quantities referring to whole domain */
    int Nproc_tot;
    int Nbricks_tot;

    int Nelems_tot_domain;
    int Ns_tot_domain;

    real_t cfl      = 0.8;
    real_t end_time = 1.0;

    /* Methods: all must be defined in derived classes */
    virtual void setup_params(){};
    virtual void secondary_params(){};
    virtual void setup_process(Process &proc){};
};




/* Simple parameter holder for homogeneous Cartesian domain */
/* Declare with some default values                         */
class ParamsCartesian : public Params
{
    public:

    int Nproc[3]   = {1, 1, 1}; // No. of processes in each direction in the whole domain
    int Nbricks[3] = {1, 1, 1}; // No. of element bricks in each direction, in each process
    int Nelems[3]  = {1, 1, 1}; // No. of elements in each direction, in each brick
    int Ns[3]      = {8, 8, 8}; // No. of solution points in each direction, in each element

    real_t domain_edge[3][2] = {{0.0,1.0}, {0.0,1.0}, {0.0, 1.0}};
    

    /* Secondary or derived quantities */
    int Nf[3];
    int Ns_tot; // Total numbers of soln/flux points per element
    int Nf_tot;

    /* General methods */
    virtual void setup_params(); 
    virtual void secondary_params();
    virtual void setup_process(Process &proc);

    /* List defined problems here */
    void test1();
};


void ParamsCartesian::secondary_params()
{
    Nproc_tot = Nproc[0]*Nproc[1]*Nproc[2];
    Ns_tot    = Ns[0]*Ns[1]*Ns[2];
    Nf[0]  = Ns[0]+1;
    Nf[1]  = Ns[1]+1;
    Nf[2]  = Ns[2]+1;
    Nf_tot = Ns[0]*Ns[1]*Nf[2] + Ns[0]*Ns[2]*Nf[1] + Ns[1]*Ns[2]*Nf[0]; // 3 Nf Ns^2 if all equal
}


void ParamsCartesian::setup_params()
{
    test1();

    secondary_params();
};


void ParamsCartesian::setup_process(Process &proc)
{
    return;
}



/***************************************************************************/
/*** Defined problems ******************************************************/

void ParamsCartesian::test1()
{
    Nproc[0] = 4; 
}
#endif
