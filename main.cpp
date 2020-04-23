#include <mpi.h>
#include "common.hpp"


static void startMPI(int argc, char *argv[], Process &proc)
{
    int error = MPI_Init(&argc, &argv);
    
    if (error != MPI_SUCCESS)
        proc.writeError("Error starting MPI");
 
    MPI_Comm_rank(MPI_COMM_WORLD, &(proc.rank));   // Get rank of this process    
    MPI_Comm_size(MPI_COMM_WORLD, &(proc.Nproc));  // Number of running processes

    proc.writeMessage("MPI started");
    proc.writeVariable<int>("number of active processes", proc.Nproc);
 
    return;
}


int main(int argc, char *argv[])
{
    Process proc;

    startMPI(argc, &(*argv), proc);
    
    if (proc.rank == 0) proc.params.write_param_info();

    MPI_Finalize();
    exit(0);
}
