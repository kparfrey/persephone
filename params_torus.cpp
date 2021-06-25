#include "params_torus.hpp"

#include <mpi.h>
#include <cmath>
#include "process.hpp"
#include "element_block.hpp"
#include "initial_state_torus.hpp"
#include "write_screen.hpp"
#include "domain_map_torus.hpp"
#include "boundary_conditions.hpp"
#include "spatial_metric.hpp"
#include "physics_includes.hpp"

using std::cout;
using std::endl;

#include <string>
using std::string;
using std::to_string;
    
/* Called by the ParamsTorus constructor */
void ParamsTorus::secondary_params()
{
    switch (central_polygon)
    {
        case square:
            Ngroup = 5;
            Ngroup_central = 1;
            Ngroup_outer   = 4;
            break;
        default:
            write::error("Torus central polygon not supported", destroy);
    }

    //Nelem_proc   = Nelem[0] * Nelem[1] * Nelem[2];
    //Ns_elem      = Ns[0] * Ns[1] * Ns[2];
    
    Nproc_central = Nproc[0];
    Nproc_outer   = Nproc[1];

    Nproc_group_central = Nproc[0] * Nproc[0] * Nproc[2];
    Nproc_group_outer   = Nproc[0] * Nproc[1] * Nproc[2];

    Nproc_domain = Ngroup_central * Nproc_group_central + Ngroup_outer * Nproc_group_outer;
    Nelem_domain = (Ngroup_central*Nelem[0]*Nelem[0] + Ngroup_outer*Nelem[0]*Nelem[1]) * Nelem[2];
    Ns_domain    = (Ngroup_central*Ns[0]*Ns[0]       + Ngroup_outer*Ns[0]*Ns[1])       * Ns[2];
}


void ParamsTorus::write_param_info()
{
    cout << endl;
    cout << "***** Parameters ******************************" << endl;
    cout << "***** Torus domain ****************************" << endl;
    cout << "***** Need to update properly ... *************" << endl;
    cout << "Processes: ";
    for (int i: dirs) cout << Nproc[i] << "  ";
    cout << "          ---  Total: " << Nproc_domain << endl;

    cout << "Elements:  ";
    for (int i: dirs) cout << Nelem[i] << "  ";
    cout << " per proc ---  Total: " << Nelem_domain << endl;

    cout << "Soln pnts: ";
    for (int i: dirs) cout << Ns[i] << "  ";
    cout << " per elem ---  Total: " << Ns_domain << endl;

    cout << endl;

    cout << "***********************************************" << endl;
    cout << endl;

    return;
}


void ParamsTorus::setup_process(Process &proc)
{
    int Nproc_group_central = Nproc[0] * Nproc[0] * Nproc[2]; // Move to params_torus variables?
    int Nproc_group_outer   = Nproc[0] * Nproc[1] * Nproc[2];

    proc.Ngroup = Ngroup;

    if (proc.rank < Nproc_group_central)
    {
        proc.group = 0;
        proc.group_rank = proc.rank;
        proc.Nproc_group[0] = Nproc[0];
        proc.Nproc_group[1] = Nproc[0];
    }
    else // outer group
    {
        /* Use floor division to separate ranks into groups */
        proc.group = 1 + (proc.rank - Nproc_group_central) / Nproc_group_outer; 
        proc.group_rank = proc.rank - Nproc_group_central - (proc.group - 1) * Nproc_group_outer;
        proc.Nproc_group[0] = Nproc[1]; // dir-0 is outward, so associated with Nproc[..,outer,..] 
        proc.Nproc_group[1] = Nproc[0]; // dir-1 is like theta, so goes with Nproc[central,...]
    }

    /* direction-2 is the same for all groups */
    proc.Nproc_group[2] = Nproc[2];

    proc.Nproc_group_tot = proc.Nproc_group[0] * proc.Nproc_group[1] * proc.Nproc_group[2];

    /* Create new MPI communicator for this group */
    MPI_Comm_split(MPI_COMM_WORLD, proc.group, proc.group_rank, &proc.group_comm);


    /* Groupwise indices of this process */
    proc.group_idx[2] = proc.group_rank % proc.Nproc_group[2];
    proc.group_idx[1] = (proc.group_rank / proc.Nproc_group[2]) % proc.Nproc_group[1];
    proc.group_idx[0] = proc.group_rank / (proc.Nproc_group[1] * proc.Nproc_group[2]);

    for (int i: dirs)
    {
        /* Copy data into ElementBlock for later convenience */
        proc.elements.group_idx[i]   = proc.group_idx[i];
        proc.elements.Nproc_group[i] = proc.Nproc_group[i];
    }

    setup_elementblock(proc.elements, proc);
    
    set_initial_state(proc.elements);

    write::variable<real_t>("CFL", proc.cfl);
    write::variable<real_t>("End time", proc.end_time);


    /*************** Setup Faces ******************************/
    for (int i: ifaces)
    {
        FaceCommunicator& face = proc.faces[i];

        face.setup(proc, i);

        /* Inter-process connectivity */
        int normal = face.normal_dir;
        face.neighbour_group = proc.group; // Modify later if this Process is on a group boundary

        /* Start by setting neighbour_idx[3] relative to the group with a straightfoward
         * Cartesian connectivity */
        for (int d: dirs)
            face.neighbour_idx[d] = proc.group_idx[d];
        
        face.neighbour_idx[normal] += face.orientation;
    }

    /* Make periodic in the 2-direction (phi) */
    if (proc.faces[0].neighbour_idx[2] < 0)
        proc.faces[0].neighbour_idx[2] = proc.Nproc_group[2] - 1;
    
    if (proc.faces[1].neighbour_idx[2] > proc.Nproc_group[2] - 1)
        proc.faces[1].neighbour_idx[2] = 0;

        
    /* Adjust neighbour indices and neighbour groups on group boundaries,
     * for faces with normals in the 0-1 (poloidal) plane */
    switch (central_polygon)
    {
        case square: // Put the following into a separate private/static function?
            if (proc.group == 0) // central square
            {
                if (proc.faces[2].neighbour_idx[0] < 0)
                {
                    FaceCommunicator& f = proc.faces[2];
                    f.neighbour_idx[0] = 0;
                    f.neighbour_idx[1] = proc.Nproc_group[1] - proc.group_idx[1] - 1; // dir reverses
                    f.neighbour_group  = 4;
                    f.neighbour_id     = 2;
                    f.change_data_order = true;
                    f.reverse_index_direction = true;
                    f.index_to_reverse = 1; // reverse normal_dir + 1
                }

                if (proc.faces[3].neighbour_idx[0] > proc.Nproc_group[0] - 1)
                {
                    FaceCommunicator& f = proc.faces[3];
                    f.neighbour_idx[0] = 0;
                    f.neighbour_group  = 2;
                    f.neighbour_id     = 2;
                }

                if (proc.faces[4].neighbour_idx[1] < 0)
                {
                    // NB these are the indices in group-1, in which the coords have exchanged connectivity
                    FaceCommunicator& f = proc.faces[4];
                    f.neighbour_idx[0] = 0;
                    f.neighbour_idx[1] = proc.group_idx[0]; // group-1:dir-1 aligned with group-0:dir-0
                    f.neighbour_group  = 1;
                    f.neighbour_id     = 2;
                    f.change_data_order = true;
                    f.swap_index_order  = true;
                }

                if (proc.faces[5].neighbour_idx[1] > proc.Nproc_group[1] - 1)
                {
                    FaceCommunicator& f = proc.faces[5];
                    f.neighbour_idx[0] = 0;
                    f.neighbour_idx[1] = proc.Nproc_group[0] - proc.group_idx[0] - 1;
                    f.neighbour_group  = 3;
                    f.neighbour_id     = 2;
                    f.change_data_order = true;
                    f.swap_index_order  = true;
                    f.reverse_index_direction = true;
                    f.index_to_reverse = 1; // reverse normal_dir + 1
                }
            }
            else // the four outer quads
            {
                /* All face-2s point into the central square */
                if (proc.faces[2].neighbour_idx[0] < 0)
                {
                    FaceCommunicator& f = proc.faces[2];
                    f.neighbour_group = 0;

                    switch (proc.group)
                    {
                        case 1:
                            f.neighbour_idx[0] = proc.group_idx[1];
                            f.neighbour_idx[1] = 0;
                            f.neighbour_id     = 4;
                            f.change_data_order = true;
                            f.swap_index_order  = true;
                            break;
                        case 2:
                            f.neighbour_idx[0] = Nproc_central - 1;
                            f.neighbour_id     = 3;
                            break;
                        case 3:
                            f.neighbour_idx[0] = Nproc_central - 1 - proc.group_idx[1];
                            f.neighbour_idx[1] = Nproc_central - 1;
                            f.neighbour_id     = 5;
                            f.change_data_order = true;
                            f.swap_index_order  = true;
                            f.reverse_index_direction = true;
                            f.index_to_reverse = 2; // reverse normal_dir + 2
                            break;
                        case 4:
                            f.neighbour_idx[0] = 0;
                            f.neighbour_idx[1] = Nproc_central - 1 - proc.group_idx[1];
                            f.neighbour_id     = 2;
                            f.change_data_order = true;
                            f.reverse_index_direction = true;
                            f.index_to_reverse = 1; // reverse normal_dir + 1
                            break;
                    }
                }
                
                if (proc.faces[3].neighbour_idx[0] > proc.Nproc_group[0] - 1)
                {
                    FaceCommunicator& f = proc.faces[3];
                    f.neighbour_idx[0] = -1;
                    f.neighbour_group  = -1;
                    f.neighbour_id     = -1;
                    f.external_face = true;

#if 0
                    /* Make periodic for testing */
                    f.neighbour_idx[0] = Nproc_outer - 1;
                    f.neighbour_idx[1] = Nproc_central - 1 - proc.group_idx[1];
                    f.neighbour_id     = 3;
                    f.change_data_order = true;
                    f.reverse_index_direction = true;
                    f.index_to_reverse = 1; // reverse normal_dir + 1
                    
                    switch(proc.group)
                    {
                        case 1:
                            f.neighbour_group = 3;
                            break;
                        case 2:
                            f.neighbour_group = 4;
                            break;
                        case 3:
                            f.neighbour_group = 1;
                            break;
                        case 4:
                            f.neighbour_group = 2;
                            break;
                    }
#endif
                }
                
                if (proc.faces[4].neighbour_idx[1] < 0)
                {
                    FaceCommunicator& f = proc.faces[4];
                    f.neighbour_idx[1] = Nproc_central - 1;
                    if (proc.group == 1)
                        f.neighbour_group = 4;
                    else
                        f.neighbour_group = proc.group - 1;
                }

                if (proc.faces[5].neighbour_idx[1] > proc.Nproc_group[1] - 1)
                {
                    FaceCommunicator& f = proc.faces[5];
                    f.neighbour_idx[1] = 0;
                    if (proc.group == 4)
                        f.neighbour_group = 1;
                    else
                        f.neighbour_group = proc.group + 1;
                }
            }
            break;
        default:
            write::error("Torus central polygon not supported", destroy);
    }

    for (int i: ifaces)
    {
        FaceCommunicator& f = proc.faces[i];


        /* Translate the groupwise neighbour_idx[3] into a global rank 
         * For the torus, Nproc_group[1] = Nproc_central for all groups, central and outer */
        int neighbour_rank_groupwise = (f.neighbour_idx[0]  * proc.Nproc_group[1]
                                      + f.neighbour_idx[1]) * proc.Nproc_group[2]
                                      + f.neighbour_idx[2];
        
        int groupwise_to_global_offset;
        int nbgrp = f.neighbour_group;
        if (nbgrp < Ngroup_central)
            groupwise_to_global_offset = nbgrp * Nproc_group_central;
        else
            groupwise_to_global_offset = Ngroup_central * Nproc_group_central
                                         + (nbgrp - Ngroup_central) * Nproc_group_outer;

        f.neighbour_rank = groupwise_to_global_offset + neighbour_rank_groupwise;


        /* Set up external boundary conditions */
        if (f.external_face == true)
            f.BC = new ImplosionTestBC(f.Ntot, proc.Nfield, equations);
    }

    /* Write out connectivity info for testing */
#if 0
    string sp = "  ";
    string grouprank = sp + "Group rank:" + to_string(proc.group_rank) + sp;
    string globalrank = sp + "Global rank:" + to_string(proc.rank) + sp;
    for (int ig = 0; ig < 5; ++ig)
    {
        if (proc.group == ig)
        {
            string group_idx = sp + to_string(proc.group_idx[0])+sp+to_string(proc.group_idx[1]) + sp;
            cout << "Group:"+to_string(proc.group) + group_idx + grouprank + globalrank << endl;
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (proc.rank == 18 || proc.rank == 19)
    {
        cout << globalrank + to_string(proc.faces[2].neighbour_rank) << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
#endif

    return;
}


void ParamsTorus::setup_elementblock(ElementBlock &elements, Process &proc)
{
    if (proc.group < Ngroup_central)
    {
        /* In a central group */
        elements.Nelem[0] = elements.Nelem[1] = Nelem[0];
        elements.Ns[0] = elements.Ns[1] = Ns[0];
    }
    else 
    {
        /* In an outer group */
        /* NB: ParamsTorus.Nelem[3] = Nelem[central, outer, phi], while
         * ElementBlock.Nelem[3] refer to the three local coordinate directions */
        elements.Nelem[0] = Nelem[1];
        elements.Nelem[1] = Nelem[0];
        elements.Ns[0] = Ns[1];
        elements.Ns[1] = Ns[0];
    }

    elements.Nelem[2] = Nelem[2];
    elements.Ns[2]    = Ns[2];

    // Can this be moved to ElementBlock::setup()? 
    elements.Nelem_block = elements.Nelem[0] * elements.Nelem[1] * elements.Nelem[2];

    elements.Nfield   = proc.Nfield;

    elements.map = new BasicSquareTorusMap(proc.group, boundary_modes);

    elements.physics_soln->metric = new DiagonalSpatialMetric(cylindrical);
    for (int d: dirs)
        elements.physics[d]->metric = new DiagonalSpatialMetric(cylindrical);

    /* At this point all external information is present, and the internal
     * setup method can take over. */
    elements.setup();

    return;
}


/* Should this be moved to initial_state_torus? */
void ParamsTorus::set_initial_state(ElementBlock &elements)
{
    switch (elements.physics_soln->system)
    {
        case navier_stokes:
            set_euler_torus(elements);
            break;
        case mhd:
            write::error("MHD not implemented for torus", destroy);
            break;
        case scalar_advection:
            write::error("Scalar advection not implemented for torus", destroy);
            break;
    }

    return;
}
