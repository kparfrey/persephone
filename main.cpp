#include <mpi.h>
#include "common.hpp"
#include "process.hpp"
#include "active_params.hpp"
#include "write_mesh.hpp"


static void startMPI(int argc, char *argv[], Process &proc)
{
    int error = MPI_Init(&argc, &argv);
    
    if (error != MPI_SUCCESS)
        proc.write_error("Error starting MPI");
 
    MPI_Comm_rank(MPI_COMM_WORLD, &(proc.rank));   // Get rank of this process    
    MPI_Comm_size(MPI_COMM_WORLD, &(proc.Nproc));  // Number of running processes

    proc.write_message("MPI started");
    proc.write_variable<int>("number of active processes", proc.Nproc);

    MPI_Barrier(MPI_COMM_WORLD);
 
    return;
}


int main(int argc, char *argv[])
{
    Process proc(active_parameters);

    startMPI(argc, &(*argv), proc);

    proc.write_startup_info();
    
    proc.params.setup_process(proc); // Need to do something about this...
    
    write_mesh(proc);


    MPI_Finalize();
    exit(0);
}
