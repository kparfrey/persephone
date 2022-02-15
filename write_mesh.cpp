#include "write_mesh.hpp"

#include <vector>
#include <mpi.h>
#include <filesystem>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

#include "process.hpp"
#include "element_block.hpp"
#include "write_screen.hpp"
#include "write_file_utils.hpp"
#include "tensor_field.hpp"
#include "edge.hpp"
#include "params.hpp"

using std::string;


/* Very basic mesh output for simplest Cartesian topology */
void write_mesh(Process &proc)
{
    ElementBlock& eb = proc.elements;
    std::vector<size_t> local_dims(3);
    std::vector<size_t> global_dims(3);
    std::vector<size_t> offset(3);

    string filename = "data/mesh.h5";
    string groupstring = std::to_string(proc.group);


    if (proc.rank == 0)
    {
        if (std::filesystem::exists(filename))
            if (std::filesystem::remove(filename))
                write::message("Deleted existing mesh file.");
    
        /* Create file and write global data 
         * Use H5Easy to reduce amount of code for writing simple datasets */
        write::message("Creating mesh file");
        
        H5Easy::File meshfile(filename, H5Easy::File::Create);
        H5Easy::dump(meshfile, "/Ngroup", proc.Ngroup);
        H5Easy::dump(meshfile, "/Nproc",  proc.Nproc);
        H5Easy::dump(meshfile, "/Nfield", proc.Nfield);
    }


    /* Start MPI output phase */
    /* Do groups one at a time. Need a full MPI barrier after each group. */
    for (int ig = 0; ig < proc.Ngroup; ++ig)
    {
        if (ig == proc.group)
        {
            HighFive::File meshfile(filename, HighFive::File::ReadWrite,
                                        HighFive::MPIOFileDriver(proc.group_comm, MPI_INFO_NULL));
            //HighFive::File meshfile(filename, HighFive::File::Create,
            //                            HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));

            /* Create HDF5 groups for our Groups of processes */
            meshfile.createGroup(groupstring);

            /* Solution point locations */
            /* Assume every Process in a given Group is identical */
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

            string group              = groupstring; // + "/coords";
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
                // write::message("Writing " + names[i]);
                dataset.select(offset, local_dims).write((real_t***)data);
            }

            delete[] data;


            /* Write edges */
            meshfile.createGroup(groupstring + "/edges");

            //for (int i = 0; i < 4; ++i) // for every edge index...
            for (int i: iedges) // for every edge index...
            {
                string name = groupstring + "/edges/" + std::to_string(i);

                int N = eb.edges[0][i].N; // N along edge; use 0th elem as prototype

                local_dims[0] = size_t(eb.Nelem_block); // 1D element index
                local_dims[1] = 3;                      // coordinate components
                local_dims[2] = size_t(N); 

                global_dims[0] = size_t(local_dims[0] * proc.Nproc_group_tot);
                global_dims[1] = local_dims[1]; // i.e. only "stack" new data along the 0 axis
                global_dims[2] = local_dims[2];

                offset[0] = size_t(proc.group_rank * local_dims[0]); // Again assume just 1 group
                offset[1] = 0;
                offset[2] = 0;

                real_t* edge_data = new real_t [local_dims[0]*local_dims[1]*local_dims[2]];

                /* Store edge data as
                 *    file[group][edge][elem_id, coord, point] 
                 * where edge is 0 .. 11
                 *       elem_id is from 0 to (no. of elems in the group) - 1 
                 *       coord is 0 .. 2
                 *       point indexes the points along the edge, 0 .. Nf-1  */
                int mem;
                for (int nelem = 0; nelem < eb.Nelem_block; ++nelem)
                    for (int d: dirs) // coordinate component
                        for (int j = 0; j < N; ++j) // edge node 
                        {
                            mem = (nelem * 3 + d) * N + j;
                            edge_data[mem] = eb.edges[nelem][i].r(d,j);
                        }

                // write::message("Writing " + name);

                HighFive::DataSet dataset = meshfile.createDataSet<real_t>(name, 
                                                   HighFive::DataSpace(global_dims));
                dataset.select(offset, local_dims).write((real_t***)edge_data);

                delete[] edge_data;
            } // closes: loop over this element's 12 edges
        } // closes: if (ig == proc.Ngroup)

        MPI_Barrier(MPI_COMM_WORLD);
    } // End of MPI output phase - closes for loop over groups


    /* Write group properties --- separate from collective MPI part
     * Assume all Processes in the Group are identical, so write from one proc
     * Could probably replace with the H5Easy interface at some point */
    for (int ig = 0; ig < proc.Ngroup; ++ig)
    {
        if ((proc.group_rank == 0) && (ig == proc.group))
        {
            // write::message("Writing properties for group " + groupstring);

            /* Reopen the mesh file in single-processor mode */
            HighFive::File meshfile(filename, HighFive::File::ReadWrite);
            
            int scalar = 0;

            /* Nelem_tot: total number of elements in the group */
            HighFive::DataSet dset = meshfile.createDataSet<int>(groupstring + "/Nelem_tot", HighFive::DataSpace::From(scalar));
            dset.write(eb.Nelem_block * proc.Nproc_group_tot);

            std::vector<int> Nelem(3); // No. of elements in the group, in each direction
            std::vector<int> Ns(3);    // No. of solution points in each elem, in each dir
            std::vector<int> Nf(3); 

            for (int i: dirs)
            {
                Nelem[i] = eb.Nelem[i] * proc.Nproc_group[i];
                Ns[i]    = eb.Ns[i];
                Nf[i]    = eb.Nf[i];
            }

            HighFive::DataSet ds0 = meshfile.createDataSet<int>(groupstring + "/Nelem", HighFive::DataSpace::From(Nelem));
            ds0.write(Nelem);

            HighFive::DataSet ds1 = meshfile.createDataSet<int>(groupstring + "/Ns", HighFive::DataSpace::From(Ns));
            ds1.write(Ns);

            HighFive::DataSet ds2 = meshfile.createDataSet<int>(groupstring + "/Nf", HighFive::DataSpace::From(Nf));
            ds2.write(Nf);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }


    write::message("Finished writing mesh file to disk.");

    return;
}
