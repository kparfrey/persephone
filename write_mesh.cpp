#include "write_mesh.hpp"

#include <vector>
#include <mpi.h>
#include <filesystem>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include "process.hpp"
#include "element_block.hpp"
#include "write_screen.hpp"
#include "write_file_utils.hpp"
#include "tensor_field.hpp"


using std::string;


/* Very basic mesh output for simplest Cartesian topology */
void write_mesh(Process &proc)
{
    ElementBlock& eb = proc.elements;
    std::vector<size_t> local_dims(3);
    std::vector<size_t> global_dims(3);
    std::vector<size_t> offset(3);

    string filename = "data/mesh.h5";

    if (proc.rank == 0)
    {
        if (std::filesystem::exists(filename))
            if (std::filesystem::remove(filename))
                write::message("Deleted existing mesh file.");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    write::message("Creating mesh file");
    HighFive::File meshfile(filename, HighFive::File::Create,
                                HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));


    /* Solution point locations */
    /* Basic structured Cartesian grid: assume every process is identical etc. */
    for (int i: dirs)
    {
        local_dims[i]  = size_t(eb.Ns_tot[i]);
        global_dims[i] = size_t(local_dims[i] * proc.Nproc_group[i]);
        offset[i]      = size_t(proc.group_idx[i] * local_dims[i]);
    }

    /* For now, just write all data as separate scalar variables */
    constexpr int Nvec   = 1;
    constexpr int Nscal  = 0;
    constexpr int Nwrite = Nscal + 3*Nvec;
    int count = 0; // For keeping track of the variable lists

    /* Automating the organization of vector components --- will want
     * to break this into a function for use with data output too. */
    string names[Nwrite];
    real_t* datalist[Nwrite];

    string group              = "/coords";
    string vecnames[Nvec]     = {"r"};
    VectorField veclist[Nvec] = {eb.rs};

    vector_organization(meshfile, Nvec, count, group, vecnames, veclist, names, datalist);

    /* Allocate data buffer for repacking each component */
    real_t* data = new real_t [local_dims[0]*local_dims[1]*local_dims[2]];

    for (int i = 0; i < Nwrite; i++)
    {
        HighFive::DataSet dataset = meshfile.createDataSet<real_t>(names[i], 
                                           HighFive::DataSpace(global_dims));

        /* Need to repack data from native ordering to elementblock-wise logical */
        repack(datalist[i], data, eb);

        /* Pass the repacked 1D array cast as a triple pointer */
        write::message("Writing " + names[i]);
        dataset.select(offset, local_dims).write((real_t***)data);
    }

    delete[] data;

    write::message("Finished writing mesh file to disk.");

    return;
}
