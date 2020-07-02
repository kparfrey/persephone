#include "face_communicator.hpp"
#include "process.hpp"

void FaceCommunicator::setup(Process& proc, int face_idx)
{
    my_idx  = face_idx;
    my_rank = proc.rank;

    /* neighbour_idx corresponds to a simple Cartesian connectivity.
     * Need to modify later if not true for this face. */
    switch (my_idx)
    {
        case 0:
            normal_dir  =  2;
            orientation = -1; 
            neighbour_idx = 1;
            break;
        case 1:
            normal_dir  =  2;
            orientation =  1;
            neighbour_idx = 0;
            break;
        case 2:
            normal_dir  =  0;
            orientation = -1;
            neighbour_idx = 3;
            break;
        case 3:
            normal_dir  =  0;
            orientation =  1;
            neighbour_idx = 2;
            break;
        case 4:
            normal_dir  =  1;
            orientation = -1;
            neighbour_idx = 5;
            break;
        case 5:
            normal_dir  =  1;
            orientation =  1;
            neighbour_idx = 4;
            break;
    }

    int normal_p1 = dir_plus_one[normal_dir];
    int normal_p2 = dir_plus_two[normal_dir];

    N[0] = proc.elements.Ns[normal_p1];
    N[1] = proc.elements.Ns[normal_p2];

    N_tot = N[0] * N[1];

    external_face = false; // By default
    
    return;
}
