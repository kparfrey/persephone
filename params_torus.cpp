#include "params_torus.hpp"

#include <cmath>
#include "process.hpp"
#include "element_block.hpp"
#include "initial_state_cartesian.hpp"
#include "write_screen.hpp"
#include "domain_map.hpp"

using std::cout;
using std::endl;
    
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

    //Nproc_group  = Nproc[0] * Nproc[1] * Nproc[2];
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

    if (proc.rank < Nproc_group_central)
    {
        proc.group = 0
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
    
    set_initial_state(proc.elements, equations);

    /* Should be made general and moved to a separate member function */
    proc.dt = 5e-3; //cfl * set_dt_basic(proc.elements); 

    write::variable<real_t>("CFL", proc.cfl);
    write::variable<real_t>("End time", proc.end_time);
    write::variable<real_t>("dt", proc.dt);
    write::variable<int>("No. of time steps", int(end_time/proc.dt));


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
            neighbour_idx[d] = proc.group_idx[d];
        
        neighbour_idx[normal] += face.orientation;
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
                    f.neighbour_idx[1] = f.my_group_rank[0]; // group-1:dir-1 aligned with group-0:dir-0
                    f.neighbour_group  = 1;
                    f.neighbour_id     = 2;
                }

                if (proc.faces[5].neighbour_idx[1] > proc.Nproc_group[1] - 1)
                {
                    FaceCommunicator& f = proc.faces[5];
                    f.neighbour_idx[0] = 0;
                    f.neighbour_idx[1] = proc.Nproc_group[0] - proc.group_idx[0] - 1;
                    f.neighbour_group  = 3;
                    f.neighbour_id     = 2;
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
                            f.neighbour_idx[0] = proc.group_rank[1];
                            f.neighbour_idx[1] = 0;
                            f.neighbour_id     = 4;
                            break;
                        case 2:
                            f.neighbour_idx[0] = Nproc_central - 1;
                            f.neighbour_id     = 3;
                            break;
                        case 3:
                            f.neighbour_idx[0] = Nproc_central - 1 - proc.group_rank[1];
                            f.neighbour_idx[1] = Nproc_central - 1;
                            f.neighbour_id     = 5;
                            break;
                        case 4:
                            f.neighbour_idx[0] = 0;
                            f.neighbour_idx[1] = Nproc_central - 1 - proc.group_rank[1];
                            f.neighbour_id     = 2;
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

    /* Translate the groupwise neighbour_idx[3] into a global rank */
    for (int i: ifaces)
    {
        /* For the torus, Nproc_group[1] = Nproc_central for all groups, central and outer */
        int neighbour_rank_groupwise = (neighbour_idx[0]  * proc.Nproc_group[1]
                                      + neighbour_idx[1]) * proc.Nproc_group[2]
                                      + neighbour_idx[2];
        
        int groupwise_to_global_offset;
        int nbgrp = proc.faces[i].neighbour_group;
        if (nbgrp < Ngroup_central)
            groupwise_to_global_offset = nbgroup * Nproc_group_central;
        else
            groupwise_to_global_offset = Ngroup_central * Nproc_group_central
                                         + (nbgrp - Ngroup_central) * Nproc_group_outer;

        proc.faces[i].neighbour_rank = groupwise_to_global_offset + neighbour_rank_groupwise;
    }

    return;
}


void ParamsTorus::setup_elementblock(ElementBlock &elements, Process &proc)
{
    for (int i: dirs)
    {
        elements.Nelem[i] = Nelem[i];
        elements.Ns[i] = Ns[i];
        elements.Nf[i] = Nf[i];
    }

    elements.Nelem_block = Nelem_proc; 
    elements.Nfield      = proc.Nfield;

    elements.geometry = geometry;

    if (geometry == simple_geometry)
    {
        /* Basic Cartesian shape --- 90 degree angles etc. 
         * Move more of this to ElementBlock.set_physical_coords_simple().
         * Corners not really needed in this simple case. */
        real_t domain_length;
        real_t length_per_proc[3];
        for (int i: dirs) 
        {
            domain_length = domain_limits[i][1] - domain_limits[i][0];
            length_per_proc[i] = domain_length / Nproc[i]; // Evenly tile
        }

        real_t proc_origin[3]; // Coordinates of corner 0
        for (int i: dirs) 
            proc_origin[i] = domain_limits[i][0] + proc.group_idx[i] * length_per_proc[i];

        for (int i: icorners)
            for (int j: dirs)
                elements.corners[i][j] = proc_origin[j] + corner_coords[i][j]*length_per_proc[j];
    }
    else // full_geometry
        elements.map = new WaveRect2D; //new QuarterAnnulusMap; // specify manually for now...


    /* At this point all external information is present, and the internal
     * setup method can take over. */
    elements.setup();

    return;
}


/* Should this be moved to initial_state_torus? */
void ParamsTorus::set_initial_state(ElementBlock &elements, EqnSystem equations)
{
    switch (equations)
    {
        case scalar_advection:
            set_scalar(elements);
            break;
        case euler:
            set_euler(elements);
            break;
    }

    return;
}


/* Very basic method to set the time step. Assumes every ElementBlock
 * is indentical, so only useful for trivial Cartesian grids, and assumes
 * the wave speed is a constant */
real_t ParamsTorus::set_dt_basic(ElementBlock& eb)
{
    real_t dt_max;
    real_t wave_speed = 1.0;
    real_t dmin[3];
    real_t speed_over_dx;

    /* Timestep based on minimal flux-point spacing */
    for (int d: dirs)
        dmin[d] = std::fabs(eb.rf[d](d,1) - eb.rf[d](d,0));

    speed_over_dx = 0.0;
    for (int d: dirs)
        speed_over_dx += wave_speed / dmin[d];


    dt_max = 1.0 / speed_over_dx;

    return dt_max;
}
