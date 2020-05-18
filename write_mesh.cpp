#include "write_mesh.hpp"

using std::vector;
using HighFive::File;


static real_t*** alloc_3D_array(vector<size_t> dims)
{
    real_t ***array = new real_t** [dims[0]];

    for (size_t i = 0; i < dims[0]; ++i)
    {
        array[i] = new real_t* [dims[1]];
        for (size_t j = 0; j < dims[1]; ++j)
            array[i][j] = new real_t [dims[2]];
    }

    return array;
}


static void threeDify(vector<size_t> dims, real_t *oneD, real_t ***threeD)
{
    for (size_t i = 0; i < dims[0]; ++i)
    for (size_t j = 0; j < dims[1]; ++j)
    for (size_t k = 0; k < dims[2]; ++k)
    {
        /* i,j,k iterate straightforwardly <<across elements>> 
         * in a single Process/ElementBlock */
        threeD[i][j][k] = i + 2.0*j + 3.0*k;
    }

    return;
}


void write_mesh(Process &proc)
{
    vector<size_t> local_dims(3);
    vector<size_t> global_dims(3);
    vector<size_t> offset(3);

    proc.write_message("Creating mesh file.");
    File meshfile("mesh.h5", File::ReadWrite | File::Create | File::Truncate,
                                HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));


    /* Solution point locations */
    /* Basic structured Cartesian grid: assume every process is identical etc. */
    for (int i: dirs)
    {
        local_dims[i]  = std::size_t(proc.elements.Ns_tot[i]);
        global_dims[i] = std::size_t(local_dims[i] * proc.Nproc_group[i]);
        offset[i] = std::size_t(proc.group_idx[i] * local_dims[i]);
    }

    for (int i: dirs)
        proc.write_variable<size_t>("Global dims", global_dims[i]);

    proc.write_message("Creating dataset.");
    HighFive::DataSet dataset = meshfile.createDataSet<real_t>("test_dataset", 
                                       HighFive::DataSpace(global_dims));

    proc.write_message("Creating 3D array.");
    real_t ***array3D = alloc_3D_array(local_dims);

    proc.write_message("Converting from 1D to 3D array.");
    threeDify(local_dims, proc.elements.rs[0], array3D);

    proc.write_message("Writing data to disk.");
    dataset.select(offset, local_dims).write(array3D);

    proc.write_message("Finished writing to disk.");
    return;
}
