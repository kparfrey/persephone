#include "write_mesh.hpp"

using std::string;



/* Repack a 1D data array from native (elementwise etc.) ordering to simple
 * order which is logically Cartesian across the elementblock for a single variable */
/* This works but needs to be made more efficient for use in data output */
static void repack(real_t* original, real_t* eb_logical, ElementBlock& eb)
{
    int il, jl, kl;          // elementblock-logical indices
    int Nl1 = eb.Ns_tot[1];  // Logical dimensions are the same as Ns_tot[3]
    int Nl2 = eb.Ns_tot[2];

    for (int ie = 0; ie < eb.Nelem[0]; ie++)
    for (int je = 0; je < eb.Nelem[1]; je++)
    for (int ke = 0; ke < eb.Nelem[2]; ke++)
    for (int i  = 0; i  < eb.Ns[0];    i++)
    for (int j  = 0; j  < eb.Ns[1];    j++)
    for (int k  = 0; k  < eb.Ns[2];    k++)
    {
        // Assume going from solution points for now...
        il = ie * eb.Ns[0] + i;
        jl = je * eb.Ns[1] + j;
        kl = ke * eb.Ns[2] + k;
        eb_logical[(il*Nl1 + jl)*Nl2 + kl] = original[eb.ids_full(i,j,k,ie,je,ke)];
    }

    return;
}


/* Create groups in the file, and fill in names & datalist from vecnames & veclist.
 * Pass in ref to count, the index in names/datalist at which the vector components begin. */
static void vector_organization(HighFive::File file,
                                const int Nvec, int& count, string group,
                                string vecnames[], real_t** veclist[],
                                string names[], real_t* datalist[])
{
    string subgroup; // The name of the vector within its parent group

    for (int i = 0; i < Nvec; i++)
    {
        subgroup = group + "/" + vecnames[i];
        file.createGroup(subgroup);
        for (int d: dirs)
        {
            names[count]    = subgroup + "/" + std::to_string(d);
            datalist[count] = veclist[i][d];
            count++;
        }
    }

    return;
}


/* Very basic mesh output for simplest Cartesian topology */
void write_mesh(Process &proc)
{
    ElementBlock& eb = proc.elements;
    std::vector<size_t> local_dims(3);
    std::vector<size_t> global_dims(3);
    std::vector<size_t> offset(3);

    proc.write_message("Creating mesh file...");
    HighFive::File meshfile("mesh.h5", HighFive::File::Overwrite,
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

    string group           = "/coords";
    string vecnames[Nvec]  = {"r"};
    real_t** veclist[Nvec] = {eb.rs};

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
        proc.write_message("Writing " + names[i] + "...");
        dataset.select(offset, local_dims).write((real_t***)data);
    }

    delete[] data;

    proc.write_message("Finished writing mesh file to disk.");
    return;
}
