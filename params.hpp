#ifndef PARAMS_HPP
#define PARAMS_HPP


/* Simple parameter holder for homogeneous Cartesian domain */
/* Declare with some default values                         */
struct Params
{
    int Nproc[3]   = {1, 1, 1}; // No. of processes in each direction in the whole domain
    int Nbricks[3] = {1, 1, 1}; // No. of element bricks in each direction, in each process
    int Nelems[3]  = {1, 1, 1}; // No. of elements in each direction, in each brick
    int Ns[3]      = {8, 8, 8}; // No. of solution points in each direction, in each element

    real_t domain_edge[3][2] = {{0.0,1.0}, {0.0,1.0}, {0.0, 1.0}};
    
    real_t cfl      = 0.8;
    real_t end_time = 1.0;


    /* Secondary or derived quantities */
    int Nproc_tot;
    int Ns_tot;
    int Nf[3];
    int Nf_tot;

    /* Generic methods */
    void setup_params(); // Modify to choose defined problem from the list below
    void secondary_params();

    /* List defined problems here */
    void test1();
};


void Params::secondary_params()
{
    Nproc_tot = Nproc[0]*Nproc[1]*Nproc[2];
    Ns_tot    = Ns[0]*Ns[1]*Ns[2];
    Nf[0] = Ns[0]+1;
    Nf[1] = Ns[1]+1;
    Nf[2] = Ns[2]+1;
    Nf_tot = Ns[0]*Ns[1]*Nf[2] + Ns[0]*Ns[2]*Nf[1] + Ns[1]*Ns[2]*Nf[0]; // 3 Nf Ns^2 if all equal
}


void Params::setup_params()
{
    test1();

    secondary_params();
};


/***************************************************************************/
/*** Defined problems ******************************************************/

void Params::test1()
{
    Nproc[0] = 4; 
}


#endif
