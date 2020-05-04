#ifndef FACE_COMMUNICATOR_HPP
#define FACE_COMMUNICATOR_HPP

#include "common.hpp"

class FaceCommunicator
{
    public:

    int N[2]; // No. of (solution) points in each dir on the face
              // Lower directional index will always be in "0" slot

    int my_rank;        // This process's rank
    int neighbour_rank; // Rank of MPI process this face is paired with 

    bool external_face; // True if this face needs a boundary condition rather
                        // than communication from another process

    real_t *my_data;
    real_t *neighbour_data;
};

#endif
