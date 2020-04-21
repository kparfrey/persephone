#ifndef PARAMS_HPP
#define PARAMS_HPP


/* Simple parameter holder for homogeneous Cartesian domain */
struct Params
{
    int Nproc[3];   // No. of processes in each direction in the whole domain
    int Nbricks[3]; // No. of element bricks in each direction, in each process
    int Nelems[3];  // No. of elements in each direction, in each brick
    int Ns[3];      // No. of solution points in each direction, in each element

    real_t domain_edge[3][2];
    
    real_t cfl;
    real_t end_time;
};
#endif
