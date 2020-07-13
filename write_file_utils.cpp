#include "write_file_utils.hpp"

/* Repack a 1D data array from native (elementwise etc.) ordering to simple
 * order which is logically Cartesian across the elementblock for a single variable */
/* This works but needs to be made more efficient for use in data output */
void repack(real_t* original, real_t* eb_logical, ElementBlock& eb)
{
    int il, jl, kl;          // elementblock-logical indices
    int Nl1 = eb.Ns_tot[1];  // Logical dimensions are the same as Ns_tot[3]
    int Nl2 = eb.Ns_tot[2];

    for (int ke = 0; ke < eb.Nelem[2]; ke++)
    for (int je = 0; je < eb.Nelem[1]; je++)
    for (int ie = 0; ie < eb.Nelem[0]; ie++)
    for (int k  = 0; k  < eb.Ns[2];    k++)
    for (int j  = 0; j  < eb.Ns[1];    j++)
    for (int i  = 0; i  < eb.Ns[0];    i++)
    {
        // Assume going from solution points for now...
        il = ie * eb.Ns[0] + i;
        jl = je * eb.Ns[1] + j;
        kl = ke * eb.Ns[2] + k;

        /* When converted to 3D array, will be indexable as array[dir0][dir1][dir2]
         * Swap the i and k below (+ Nl2->Nl0) if prefer array[dir2][dir1][dir0] */
        eb_logical[(il*Nl1 + jl)*Nl2 + kl] = original[eb.ids_full(i,j,k,ie,je,ke)];
    }

    return;
}


/* Create groups in the file, and fill in names & datalist from vecnames & veclist.
 * Pass in ref to count, the index in names/datalist at which the vector components begin. */
void vector_organization(HighFive::File file,
                         const int Nvec, int& count, string group,
                         string vecnames[], VectorField veclist[],
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
            datalist[count] = veclist[i](d);
            count++;
        }
    }

    return;
}

