#include "params_cartesian.hpp"

#include <cmath>
#include "process.hpp"
#include "element_block.hpp"
#include "initial_state_cartesian.hpp"
#include "write_screen.hpp"

using std::cout;
using std::endl;
    
/* Called by the ParamsCartesian constructor */
void ParamsCartesian::secondary_params()
{
    Nproc_domain = Nproc[0] * Nproc[1] * Nproc[2];
    Nelem_proc   = Nelem[0] * Nelem[1] * Nelem[2];
    Nelem_domain = Nelem_proc * Nproc_domain;
    Ns_elem      = Ns[0] * Ns[1] * Ns[2];
    Ns_domain    = Ns_elem * Nelem_domain;
}


void ParamsCartesian::write_param_info()
{
    cout << endl;
    cout << "***** Parameters ******************************" << endl;
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

    cout << "Domain:   ";
    for (int i: dirs) cout << domain_edge[i][0] << " --> " << domain_edge[i][1] << "   ";
    cout << endl;

    cout << "***********************************************" << endl;
    cout << endl;

    return;
}


void ParamsCartesian::setup_process(Process &proc)
{
    /* Setup common to all geometries. */
    setup_process_generic(proc);


    /* This process's set of indices in the global 3D array of processes   */
    int proc_idx[3];

    proc_idx[2] = proc.rank % Nproc[2];
    proc_idx[1] = (proc.rank / Nproc[2]) % Nproc[1];
    proc_idx[0] = proc.rank / (Nproc[1]*Nproc[2]);

    /* A Cartesian dataset needs only one Cartesian group */
    proc.group = 0;
    for (int i: dirs)
    {
        proc.group_idx[i] = proc_idx[i];
        proc.Nproc_group[i] = Nproc[i];
    }



    /* For now: assume basic Cartesian shape --- 90 degree angles etc. */
    real_t domain_length;
    real_t length_per_proc[3];
    for (int i: dirs) 
    {
        domain_length = domain_edge[i][1] - domain_edge[i][0];
        length_per_proc[i] = domain_length / Nproc[i]; // Evenly tile
    }

    real_t proc_origin[3]; // Coordinates of corner 0
    for (int i: dirs) 
        proc_origin[i] = domain_edge[i][0] + proc_idx[i] * length_per_proc[i];

    for (int i: icorners)
        for (int j: dirs)
            proc.corners[i][j] = proc_origin[j] + corner_coords[i][j]*length_per_proc[j];


    setup_elementblock(proc.elements, proc);
    
    set_initial_state(proc.elements);

    /* Should be made general and moved to a separate member function */
    proc.dt = cfl * set_dt_basic(proc.elements); 

    write::variable<real_t>("CFL", proc.cfl);
    write::variable<real_t>("End time", proc.end_time);
    write::variable<real_t>("dt", proc.dt);
    write::variable<int>("No. of time steps", int(end_time/proc.dt));


    /* Setup faces */
    for (int i: ifaces)
    {
        FaceCommunicator& face = proc.faces[i];

        face.setup(proc, i);

        /* Inter-process connectivity */
        /* Make periodic for now */
        int normal = face.normal_dir;
        int neighbour_idx[3];
        for (int d: dirs)
            neighbour_idx[d] = proc.group_idx[d];
        
        neighbour_idx[normal] += face.orientation;

        /* Make periodic */
        if (neighbour_idx[normal] < 0)
            neighbour_idx[normal] = proc.Nproc_group[normal] - 1;
        
        if (neighbour_idx[normal] > proc.Nproc_group[normal] - 1)
            neighbour_idx[normal] = 0;


        face.neighbour_rank = (neighbour_idx[0]  * proc.Nproc_group[1]
                             + neighbour_idx[1]) * proc.Nproc_group[2]
                             + neighbour_idx[2];
    }


    /***
    if (proc_idx[2] == 0) proc_base = Nproc[2];
    else proc_base = proc_idx[2];
    proc.faces[0].neighbour_rank = proc_base - 1;
    
    if (proc_idx[2] == Nproc[2]-1) proc_base = -1;
    else proc_base = proc_idx[2];
    proc.faces[1].neighbour_rank = proc_base + 1;

    if (proc_idx[0] == 0) proc_base = Nproc[0];
    else proc_base = proc_idx[0];
    proc.faces[2].neighbour_rank = proc_base - 1;

    if (proc_idx[0] == Nproc[0]-1) proc_base = -1;
    else proc_base = proc_idx[0];
    proc.faces[3].neighbour_rank = proc_base + 1;

    if (proc_idx[1] == 0) proc_base = Nproc[1];
    else proc_base = proc_idx[1];
    proc.faces[4].neighbour_rank = proc_base - 1;

    if (proc_idx[1] == Nproc[1]-1) proc_base = -1;
    else proc_base = proc_idx[1];
    proc.faces[5].neighbour_rank = proc_base + 1;
    ***/


    return;
}


void ParamsCartesian::setup_elementblock(ElementBlock &elements, Process &proc)
{
    for (int i: dirs)
    {
        elements.Nelem[i] = Nelem[i];
        elements.Ns[i] = Ns[i];
        elements.Nf[i] = Nf[i];
    }

    elements.Nelem_block = Nelem_proc; 
    elements.Nfield      = proc.Nfield;


    /* Set geometrical information: trivial for now since Process and
     * ElementBlock boundaries are assumed to coincide */
    for (int i: icorners)
        for (int j: dirs)
            elements.corners[i][j] = proc.corners[i][j];

    for (int i: iedges)
        elements.edges[i] = proc.edges[i];


    /* At this point all external information is present, and the internal
     * setup method can take over. */
    elements.setup();

    return;
}


void ParamsCartesian::set_initial_state(ElementBlock &elements)
{
    set_scalar(elements);

    return;
}


/* Very basic method to set the time step. Assumes every ElementBlock
 * is indentical, so only useful for trivial Cartesian grids, and assumes
 * the wave speed is a constant */
real_t ParamsCartesian::set_dt_basic(ElementBlock& eb)
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
