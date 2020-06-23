#include "write_data.hpp"

#include <vector>
#include <mpi.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include "process.hpp"
#include "element_block.hpp"
#include "write_file_utils.hpp"
#include "write_screen.hpp"


using std::string;


/* Very basic data output for simple Cartesian mesh */
void write_data(Process &proc)
{
    ElementBlock& eb = proc.elements;
    std::vector<size_t> local_dims(3);
    std::vector<size_t> global_dims(3);
    std::vector<size_t> offset(3);


    /* Create filename */
    std::stringstream filenum;
    filenum << std::setw(4) << std::setfill('0') << proc.data_output_counter;
    string filename = "data/data" + filenum.str() + ".h5";

    if (proc.rank == 0)
    {
        if (std::filesystem::exists(filename))
            if (std::filesystem::remove(filename))
                write::message("Deleted existing data file " + filename);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    write::message("Creating data file " + std::to_string(proc.data_output_counter));
    HighFive::File datafile(filename, HighFive::File::Create,
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
    constexpr int Nvec   = 0;
    constexpr int Nscal  = 1;
    constexpr int Nwrite = Nscal + 3*Nvec;
    //int count = 0; // For keeping track of the variable lists

    string names[Nwrite];
    real_t* datalist[Nwrite];

    //string group           = "/coords";
    //string vecnames[Nvec]  = {"r"};
    //real_t** veclist[Nvec] = {eb.rs};
    //vector_organization(meshfile, Nvec, count, group, vecnames, veclist, names, datalist);

    names[0] = "phi";
    datalist[0] = eb.fields;

    /* Allocate data buffer for repacking each component */
    real_t* data = new real_t [local_dims[0]*local_dims[1]*local_dims[2]];

    for (int i = 0; i < Nwrite; i++)
    {
        HighFive::DataSet dataset = datafile.createDataSet<real_t>(names[i], 
                                             HighFive::DataSpace(global_dims));

        /* Need to repack data from native ordering to elementblock-wise logical */
        repack(datalist[i], data, eb);

        /* Pass the repacked 1D array cast as a triple pointer */
        write::message("Writing " + names[i]);
        dataset.select(offset, local_dims).write((real_t***)data);
    }

    delete[] data;

    proc.data_output_counter++;
    proc.time_last_write = proc.time;
    proc.step_last_write = proc.step;

    write::message("Finished writing data file to disk.");
    return;
}
