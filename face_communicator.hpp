#ifndef FACE_COMMUNICATOR_HPP
#define FACE_COMMUNICATOR_HPP

#include "common.hpp"
class Process;

class FaceCommunicator
{
    public:

    int my_rank;        // This process's rank
    int neighbour_rank; // Rank of MPI process this face is paired with 
    int my_idx;         // Index/label from 0 to 5 for this face
    int neighbour_idx;
    int normal_dir; // Direction of the normal to the face: 0,1,2 in ref space
    int orientation;// +1 if the elem's OUTWARD normal points along the 0/1/2 axis
                    // -1 if it points along the negative of the 0/1/2 axis

    bool external_face; // True if this face needs a boundary condition rather
                        // than communication from another process

    int N[2]; // Total no. of (solution) points in each dir on the face
              // Indices assigned cyclically: N[0] -> normal_dir + 1 cyclic
              //                              N[1] -> normal_dir + 2 cyclic
    int N_tot; // Shorthand for N[0] * N[1]

    real_t *my_data;
    real_t *neighbour_data;

    real_t *my_data_host;
    real_t *neighbour_data_host;


    /* Member functions */
    void setup(Process& proc, int face_idx);
};

#endif
