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
#include <highfive/H5Easy.hpp>

#include "process.hpp"
#include "element_block.hpp"
#include "physics_includes.hpp"
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

    string groupstring = std::to_string(proc.group);


    /* Create filename */
    std::stringstream filenum;
    filenum << std::setw(4) << std::setfill('0') << proc.data_output_counter;
    string filename = "data/data" + filenum.str() + ".h5";

    if (proc.rank == 0)
    {
        if (std::filesystem::exists(filename))
            if (std::filesystem::remove(filename))
                write::message("Deleted existing data file " + filename);
        
        H5Easy::File datafile(filename, H5Easy::File::Create);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* MPI collective part */
    for (int ig = 0; ig < proc.Ngroup; ++ig)
    {
        if (ig == proc.group)
        { 
            write::message("Creating data file " + std::to_string(proc.data_output_counter));
            HighFive::File datafile(filename, HighFive::File::ReadWrite,
                                        HighFive::MPIOFileDriver(proc.group_comm, MPI_INFO_NULL));

            /* Create HDF5 groups for our Groups of processes */
            datafile.createGroup(groupstring);

            /* Data at solution points */
            /* Assume every Process in a given Group is identical */
            for (int i: dirs)
            {
                local_dims[i]  = size_t(eb.Ns_tot[i]);
                global_dims[i] = size_t(local_dims[i] * proc.Nproc_group[i]);
                offset[i]      = size_t(proc.group_idx[i] * local_dims[i]);
            }

            /* Transform from conserved to primitive variables */
            real_t* primitives = new real_t [proc.Nfield * eb.Ns_block];
            
            real_t* Up = new real_t [proc.Nfield];
            real_t* Pp = new real_t [proc.Nfield];
            for (int n = 0; n < eb.Ns_block; ++n)
            {
                for (int field = 0; field < proc.Nfield; ++field)
                    Up[field] = eb.fields[n + field * eb.Ns_block];

                eb.physics_soln->ConservedToPrimitive(Up, Pp, n);
                eb.physics_soln->OrthonormaliseVectors(Pp, n);

                for (int field = 0; field < proc.Nfield; ++field)
                    primitives[n + field * eb.Ns_block] = Pp[field];

                /* Convert B units, e.g. to Teslas */
                if (proc.system == mhd)
                    for (int field = 5; field <= 8; ++field)
                        primitives[n + field * eb.Ns_block] *= ((MHD*)(proc.elements.physics_soln))->sqrt_mu0;
            }
            delete[] Up;
            delete[] Pp;

            /* Allocate data buffer for repacking each component */
            real_t* data = new real_t [local_dims[0]*local_dims[1]*local_dims[2]];

            for (int i = 0; i < proc.Nfield; ++i)
            {
                string name = groupstring + "/" + eb.physics_soln->variables[i];
                HighFive::DataSet dataset = datafile.createDataSet<real_t>(name, 
                                                     HighFive::DataSpace(global_dims));

                /* Need to repack data from native ordering to elementblock-wise logical */
                repack(&primitives[i*eb.Ns_block], data, eb);

                /* Pass the repacked 1D array cast as a triple pointer */
                // write::message("Writing " + name);
                dataset.select(offset, local_dims).write((real_t***)data);
            }
            
            delete[] primitives;

            /* Diagnostic field output: just divB for now */
            if (proc.system == mhd)
            {
                string name = groupstring + "/divB";
                HighFive::DataSet dataset = datafile.createDataSet<real_t>(name, 
                                                     HighFive::DataSpace(global_dims));

                /* Need to repack data from native ordering to elementblock-wise logical */
                repack(eb.divB, data, eb);

                /* Pass the repacked 1D array cast as a triple pointer */
                // write::message("Writing " + name);
                dataset.select(offset, local_dims).write((real_t***)data);
            }

            delete[] data;
        }// closes: if (ig == proc.Ngroup) 

        MPI_Barrier(MPI_COMM_WORLD);
    }// End of MPI collective part -- closes for loop over groups



    /* Write time etc. from the root process */
    /* Move to top... */
    if (proc.rank == 0)
    {
        //write::message("Writing time etc.");

        /* Reopen the mesh file in single-processor mode */
        H5Easy::File datafile(filename, H5Easy::File::ReadWrite);
        H5Easy::dump(datafile, "time", proc.time);
        H5Easy::dump(datafile, "step", proc.step);
        H5Easy::dump(datafile, "dt",   proc.dt);

        /**
        HighFive::File datafile(filename, HighFive::File::ReadWrite);
        
        HighFive::DataSet ds0 = datafile.createDataSet<real_t>("time", HighFive::DataSpace::From(proc.time));
        ds0.write(proc.time);

        HighFive::DataSet ds1 = datafile.createDataSet<int>("step", HighFive::DataSpace::From(proc.step));
        ds1.write(proc.step);

        HighFive::DataSet ds2 = datafile.createDataSet<real_t>("dt", HighFive::DataSpace::From(proc.dt));
        ds2.write(proc.dt);
        **/
    }

    write::message("Finished writing data file to disk.");

    proc.data_output_counter++;
    proc.time_last_write = proc.time;
    proc.step_last_write = proc.step;
    proc.is_output_step  = false; 

    return;
}
