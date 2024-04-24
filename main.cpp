#include <mpi.h>
#include "common.hpp"
#include "process.hpp"
#include "active_params.hpp"
#include "write_mesh.hpp"
#include "write_data.hpp"
#include "write_screen.hpp"

#include "element_block.hpp"
#include "physics.hpp"
#include "mhd.hpp"


static void startMPI(int argc, char *argv[], Process &proc)
{
    int error = MPI_Init(&argc, &argv);
    
    if (error != MPI_SUCCESS)
    {
        write::basic("Error starting MPI");
        exit(1);
    }
 
    MPI_Comm_rank(MPI_COMM_WORLD, &(proc.rank));   // Get rank of this process    
    MPI_Comm_size(MPI_COMM_WORLD, &(proc.Nproc));  // Number of running processes

    if (proc.rank == 0)
        if (proc.Nproc != proc.params.Nproc_domain)
        {
            write::basic("Error: number of running processes different from number expected from params.");
            write::basic("Expected: " + write::str(proc.params.Nproc_domain) + " processes");
            exit(1);
        }

    write::store_rank(proc.rank);
    write::message("MPI started");
    write::variable<int>("Active processes", proc.Nproc);

    MPI_Barrier(MPI_COMM_WORLD);
 
    return;
}


int main(int argc, char *argv[])
{
    Process proc(active_parameters);

    startMPI(argc, &(*argv), proc);

    proc.write_startup_info();
    
    proc.setup();

    write_mesh(proc); // Mesh data still on host

    proc.elements.free_setup_memory(); // Keep edges until mesh has been written

    proc.move_to_device();

    proc.fill_external_boundary_data();
    
    write_data(proc); // Data generally lives on device

    write::message("\nFinished setup, starting time advancement \n");

    /*** Run divB cleaner a few times first ***/
    int Nclean = 1000;
    int clean_output_freq = 100;
    ((MHD*)proc.elements.physics_soln)->divB_subsystem_only = true;
    for (int d: dirs)
        ((MHD*)proc.elements.physics[d])->divB_subsystem_only = true;
    write::message("Starting divB cleaning.....");
    for (int i = 0; i < Nclean; i++)
    {
        proc.time_advance();

        if (proc.step % clean_output_freq == 0)
            write_data(proc);
    }
    ((MHD*)proc.elements.physics_soln)->divB_subsystem_only = false;
    for (int d: dirs)
        ((MHD*)proc.elements.physics[d])->divB_subsystem_only = false;

    write::message("Starting standard MHD evolution.....");
    /******************************************/


    while(proc.time < proc.end_time)
    {
        if ((proc.time - proc.time_last_write) > proc.dt_write)
            proc.is_output_step = true;

        proc.time_advance();

        if (proc.is_output_step)
            write_data(proc);
    }

    if (proc.step > proc.step_last_write + 1)
        write_data(proc);

    MPI_Finalize();
    exit(0);
}
