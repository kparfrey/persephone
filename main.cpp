#include <mpi.h>
#include "common.hpp"
#include "process.hpp"
#include "active_params.hpp"
#include "write_mesh.hpp"

#include "matrix.hpp"

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

    if (proc.rank==0)
    {
        //real_t m[3][3] = {{1.0,4.5,3.1},{-0.1,9.1,-7.7},{5.0,-3.0,1.8}};
        //real_t m[3][3] = {{8.0,-4.5,2.1},{-0.1,5.1,-1.7},{5.0,-3.0,3.8}};
        //real_t m[3][3] = {{0.8,-1.5,2.1},{-0.8,1.1,-1.7},{0.9,-2.0,1.6}};
        real_t m[3][3] = {{8000.,-15000.,21000.},{-8000.,11000,-17000},{9000,-20000,16000.}};
        Matrix matrix(m);
        matrix.test();
    }


    MPI_Finalize();
    exit(0);
}
