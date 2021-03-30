#ifndef FACE_COMMUNICATOR_HPP
#define FACE_COMMUNICATOR_HPP

#include <mpi.h>
#include "common.hpp"
#include "tensor_field.hpp"

class Process;

class FaceCommunicator
{
    private:
    void allocate();

    public:

    int my_rank;          // This process's rank
    int my_group_rank;    // This process's rank within its group
    int neighbour_rank;   // Rank of MPI process this face is paired with 
    int my_id;            // Index/label from 0 to 5 for this face
    int neighbour_id;     // 0-5 index label for this face's neighbour
    int neighbour_idx[3]; // 3D groupwise index of this face's neighbours
    int neighbour_group;  // The group of this face's neighbour
    int normal_dir;       // Direction of the normal to the face: 0,1,2 in ref space
    int orientation;      // +1 if the elem's OUTWARD normal points along the 0/1/2 axis
                          // -1 if it points along the negative of the 0/1/2 axis

    bool external_face; // True if this face needs a boundary condition rather
                        // than communication from another process

    int n0;  // Normal-direction flux-point index of this face
    int ne0; // Normal-direction element index of this face

    int Nelem[2]; // No. of elems in each direction on the face
    int N[2]; // No. of (solution) points in each dir on the face, in each elem
              // Indices assigned cyclically: N[0] -> normal_dir + 1 cyclic
              //                              N[1] -> normal_dir + 2 cyclic
    int Ntot; // Total no. of points on face: N[0] * N[1] * Nelem[0] * Nelem[1]
    int Ntot_all; // Total no. of data values: Nfield * N_tot

    VectorField normal;

    real_t* my_data;
    real_t* neighbour_data;

    real_t* my_data_host;
    real_t* neighbour_data_host;


    /* For handling non-trivial connectivity */
    bool change_data_order;  // Set to true to use any of the following
    bool swap_index_order;
    bool reverse_index_direction;
    int index_to_reverse;


    /* Member functions */
    void setup(Process& proc, int face_id);
    MPI_Request send_data();
    MPI_Request receive_data();
};

#endif
