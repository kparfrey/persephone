#include "face_communicator.hpp"

#include "process.hpp"
#include "kernels.hpp"
#include "write_screen.hpp"

void FaceCommunicator::setup(Process& proc, int face_id)
{
    my_id   = face_id;
    my_rank = proc.rank;
    my_group_rank = proc.group_rank;

    /* neighbour_id corresponds to a simple Cartesian connectivity.
     * Need to modify later if not true for this face (e.g. the face
     * is on a group boundary. */
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

    Ntot     = N[0] * N[1] * Nelem[0] * Nelem[1];
    Ntot_all = proc.Nfield * Ntot;

    external_face = false; // By default
    change_data_order = false;
    swap_index_order  = false;
    reverse_index_direction = false;

   
    allocate();

    /* Fill normal array from ElementBlock */
    int id_elem_face, id_elem_normals;
    int mem_offset_face, mem_offset_normals;
    int mem_face, mem_normals;

    for (int ne2 = 0; ne2 < Nelem[1]; ++ne2)
    for (int ne1 = 0; ne1 < Nelem[0]; ++ne1)
    {
        /* Offsets for this face */
        id_elem_face    = (ne2 * Nelem[0]) + ne1;
        mem_offset_face = id_elem_face * N[0] * N[1];

        /* Offsets for the normals array */
        if (orientation == -1)
            id_elem_normals = (ne2*Nelem[0] + ne1)*(proc.elements.Nelem[normal_dir]+1) + 0;
        else /* Want rightmost face: ne0 + 1 */
            id_elem_normals = (ne2*Nelem[0] + ne1)*(proc.elements.Nelem[normal_dir]+1) + ne0+1;
        mem_offset_normals = id_elem_normals * N[0] * N[1];

        for (int n2 = 0; n2 < N[1]; ++n2)
        for (int n1 = 0; n1 < N[0]; ++n1)
        {
            mem_face    = mem_offset_face    + n2 * N[0] + n1;
            mem_normals = mem_offset_normals + n2 * N[0] + n1;  

            for (int i: dirs)
                normal(i,mem_face) = proc.elements.metric.normal[normal_dir](i,mem_normals);
        }
    }

    return;
}


void FaceCommunicator::allocate()
{
    my_data        = kernels::alloc(Ntot_all);
    neighbour_data = kernels::alloc(Ntot_all);

#if USING_ACCEL
    my_data_host        = new real_t [Ntot_all];
    neighbour_data_host = new real_t [Ntot_all];
#else
    /* If not using a GPU etc just point directly
     * to the device data arrays */
    my_data_host        = my_data;
    neighbour_data_host = neighbour_data;
#endif

    for (int i: dirs)
        normal(i) = new real_t [Ntot];

    return;
}


MPI_Request FaceCommunicator::send_data()
{
    /* Will need to add up/down data movement when using accelerators */

    MPI_Request request = MPI_REQUEST_NULL;
    int send_tag = my_id;

    if (!external_face)
        MPI_Isend(my_data_host, Ntot_all, MPI_real_t, neighbour_rank, send_tag, 
                                                   MPI_COMM_WORLD, &request);

    return request;
}



MPI_Request FaceCommunicator::receive_data()
{
    /* Will need to add up/down data movement when using accelerators */

    MPI_Request request = MPI_REQUEST_NULL;
    int recv_tag = neighbour_id;
    
    if (!external_face)
        MPI_Irecv(neighbour_data_host, Ntot_all, MPI_real_t, neighbour_rank,
                                      recv_tag, MPI_COMM_WORLD, &request);

    return request;
}

