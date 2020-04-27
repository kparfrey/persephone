#include "params_cartesian.hpp"

using std::cout;
using std::endl;
    
void ParamsCartesian::secondary_params()
{
    Nproc_tot        = Nproc[0] * Nproc[1] * Nproc[2];
    Nbloc_tot_domain = Nbloc[0] * Nbloc[1] * Nbloc[2] * Nproc_tot;
    Nelem_tot        = Nelem[0] * Nelem[1] * Nelem[2];
    Nelem_tot_domain = Nelem_tot * Nbloc_tot_domain;
    Ns_tot           = Ns[0] * Ns[1] * Ns[2];
    Ns_tot_domain    = Ns_tot * Nelem_tot_domain;

    Nf[0]  = Ns[0]+1;
    Nf[1]  = Ns[1]+1;
    Nf[2]  = Ns[2]+1;
    Nf_tot = Ns[0]*Ns[1]*Nf[2] + Ns[0]*Ns[2]*Nf[1] + Ns[1]*Ns[2]*Nf[0]; // 3 Nf Ns^2 if all equal
}


void ParamsCartesian::write_param_info()
{
    cout << endl;
    cout << "***** Parameters ******************************" << endl;
    cout << "Processes: ";
    for (int i: {0,1,2}) cout << Nproc[i] << "  ";
    cout << "           ---  Total: " << Nproc_tot << endl;

    cout << "Blocks:    ";
    for (int i: {0,1,2}) cout << Nbloc[i] << "  ";
    cout << " per proc  ---  Total: " << Nbloc_tot_domain << endl;

    cout << "Elements:  ";
    for (int i: {0,1,2}) cout << Nelem[i] << "  ";
    cout << " per block ---  Total: " << Nelem_tot_domain << endl;

    cout << "Soln pnts: ";
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
