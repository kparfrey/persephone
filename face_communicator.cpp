#include "face_communicator.hpp"

#include "process.hpp"
#include "kernels.hpp"
#include "write_screen.hpp"

void FaceCommunicator::setup(Process& proc, int face_id)
{
    my_id   = face_id;
    my_rank = proc.rank;

    /* neighbour_id corresponds to a simple Cartesian connectivity.
     * Need to modify later if not true for this face. */
    switch (my_id)
    {
        case 0:
            normal_dir  =  2;
            orientation = -1; 
            neighbour_id = 1;
            break;
        case 1:
            normal_dir  =  2;
            orientation =  1;
            neighbour_id = 0;
            break;
        case 2:
            normal_dir  =  0;
            orientation = -1;
            neighbour_id = 3;
            break;
        case 3:
            normal_dir  =  0;
            orientation =  1;
            neighbour_id = 2;
            break;
        case 4:
            normal_dir  =  1;
            orientation = -1;
            neighbour_id = 5;
            break;
        case 5:
            normal_dir  =  1;
            orientation =  1;
            neighbour_id = 4;
            break;
    }

    /* Normal-direction indices of this face: flux-point and element */
    if (orientation == -1)
    {
        n0  = 0; 
        ne0 = 0;
    }
    else
    {
        n0  = proc.elements.Nf[normal_dir] - 1;
        ne0 = proc.elements.Nelem[normal_dir] -1;
    }
    

    int normal_p1 = dir_plus_one[normal_dir];
    int normal_p2 = dir_plus_two[normal_dir];

    Nelem[0] = proc.elements.Nelem[normal_p1];
    Nelem[1] = proc.elements.Nelem[normal_p2];
    N[0] = proc.elements.Ns[normal_p1];
    N[1] = proc.elements.Ns[normal_p2];

    N_tot     = N[0] * N[1] * Nelem[0] * Nelem[1];
    N_tot_all = proc.Nfield * N_tot;

    external_face = false; // By default
   
    allocate();

    return;
}


void FaceCommunicator::allocate()
{
    my_data        = kernels::alloc(N_tot_all);
    neighbour_data = kernels::alloc(N_tot_all);

#if USING_ACCEL
    my_data_host        = new real_t [N_tot_all];
    neighbour_data_host = new real_t [N_tot_all];
#else
    /* If not using a GPU etc just point directly
     * to the device data arrays */
    my_data_host        = my_data;
    neighbour_data_host = neighbour_data;
#endif

    return;
}


MPI_Request FaceCommunicator::send_data()
{
    /* Will need to add up/down data movement when using accelerators */

    MPI_Request request;
    int send_tag = my_id;

    MPI_Isend(my_data_host, N_tot_all, MPI_real_t, neighbour_rank, send_tag, 
                                               MPI_COMM_WORLD, &request);

    return request;
}



MPI_Request FaceCommunicator::receive_data()
{
    /* Will need to add up/down data movement when using accelerators */

    MPI_Request request;
    int recv_tag = neighbour_id;
    
    MPI_Irecv(neighbour_data_host, N_tot_all, MPI_real_t, neighbour_rank,
                                  recv_tag, MPI_COMM_WORLD, &request);

    return request;
}

