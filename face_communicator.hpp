#ifndef FACE_COMMUNICATOR_HPP
#define FACE_COMMUNICATOR_HPP

#include "common.hpp"

class FaceCommunicator
{
    public:

    int N[2]; // Total no. of (solution) points in each dir on the face
              // Lower directional index will always be in "0" slot

    int my_rank;        // This process's rank
    int neighbour_rank; // Rank of MPI process this face is paired with 
    int my_idx;         // Index/label from 0 to 5 for this face
    int neighbour_idx;

    bool external_face; // True if this face needs a boundary condition rather
                        // than communication from another process

    real_t *my_data;
    real_t *neighbour_data;
};

#endif
