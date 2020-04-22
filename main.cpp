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
 
    return;
}


int main(int argc, char *argv[])
{
    Process proc;
    proc.params.setup_params();

    startMPI(argc, &(*argv), proc);

    MPI_Finalize();
    exit(0);
}
