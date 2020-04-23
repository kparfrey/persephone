#ifndef PARAMS_HPP
#define PARAMS_HPP

#include "common.hpp"

using std::cout;
using std::endl;
    
class Process; // Process type will be used below


/* Base parameter class. Will eventually likely have this, ParamsCartesian, 
 * ParamsToroidal, ParamsSpherical etc., all in separate files.
 */
class Params
{
    public:

    /* Total quantities referring to whole domain */
    int Nproc_tot;
    int Nbrick_tot_domain;
    int Nelem_tot_domain;
    int Ns_tot_domain;

    real_t cfl      = 0.8;
    real_t end_time = 1.0;

    /* Methods: all must be defined in derived classes */
    virtual void secondary_params(){};
    virtual void setup_process(Process &proc){};
    virtual void write_param_info(){std::cout << "Calling wrong method!" << std::endl;};
};




/* Simple parameter holder for homogeneous Cartesian domain */
/* Declare with some default values                         */
class ParamsCartesian : public Params
{
    public:

    int Nproc[3];  // No. of processes in each direction in the whole domain
    int Nbrick[3]; // No. of element bricks in each direction, in each process
    int Nelem[3];  // No. of elements in each direction, in each brick
    int Ns[3];     // No. of solution points in each direction, in each element

    real_t domain_edge[3][2];
    

    /* Secondary or derived quantities */
    int Nf[3];
    int Ns_tot; // Total numbers of soln/flux points per element
    int Nf_tot;
    int Nelem_tot;

    /* General methods */
    virtual void secondary_params();
    virtual void setup_process(Process &proc);
    virtual void write_param_info();


    /* Constructor */
    ParamsCartesian(int (& Nproc_)[3] , 
                    int (& Nbrick_)[3],
                    int (& Nelem_)[3] , 
                    int (& Ns_)[3]    , 
                    real_t (& domain_edge_)[3][2])
    {
        for (int i=0; i<3; i++)
        {
            Nproc[i]  = Nproc_[i];
            Nbrick[i] = Nbrick_[i];
            Nelem[i]  = Nelem_[i];
            Ns[i]     = Ns_[i];
            domain_edge[i][0] = domain_edge_[i][0];
            domain_edge[i][1] = domain_edge_[i][1];
        }
    
        secondary_params();
    }
};


void ParamsCartesian::secondary_params()
{
    Nproc_tot         = Nproc[0] * Nproc[1]  * Nproc[2];
    Nbrick_tot_domain = Nbrick[0]* Nbrick[1] * Nbrick[2] * Nproc_tot;
    Nelem_tot         = Nelem[0] * Nelem[1]  * Nelem[2];
    Nelem_tot_domain  = Nelem_tot * Nbrick_tot_domain;
    Ns_tot            = Ns[0]*Ns[1]*Ns[2];
    Ns_tot_domain     = Ns_tot * Nelem_tot_domain;

    Nf[0]  = Ns[0]+1;
    Nf[1]  = Ns[1]+1;
    Nf[2]  = Ns[2]+1;
    Nf_tot = Ns[0]*Ns[1]*Nf[2] + Ns[0]*Ns[2]*Nf[1] + Ns[1]*Ns[2]*Nf[0]; // 3 Nf Ns^2 if all equal
}


void ParamsCartesian::write_param_info()
{
    cout << endl;
    cout << "***** Parameters ******************************" << endl;
    cout << "Procs:    ";
    for (int i: {0,1,2}) cout << Nproc[i] << "  ";
    cout << "---  Total: " << Nproc_tot << endl;

    cout << "Bricks:   ";
    for (int i: {0,1,2}) cout << Nbrick[i] << "  ";
    cout << " per proc  ---  Total: " << Nbrick_tot_domain << endl;

    cout << "Elements: ";
    for (int i: {0,1,2}) cout << Nelem[i] << "  ";
    cout << " per brick ---  Total: " << Nelem_tot_domain << endl;

    cout << "Solution: ";
    for (int i: {0,1,2}) cout << Ns[i] << "  ";
    cout << " per elem  ---  Total: " << Ns_tot_domain << endl;

    cout << endl;

    cout << "Domain:   ";
    for (int i: {0,1,2}) cout << domain_edge[i][0] << " --> " << domain_edge[i][1] << "   ";
    cout << endl;

    cout << "***********************************************" << endl;
    cout << endl;

    return;
}


void ParamsCartesian::setup_process(Process &proc)
{
    return;
}


#endif
