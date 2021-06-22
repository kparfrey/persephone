#include "params_cartesian.hpp"

#include <cmath>
#include <mpi.h>
#include "process.hpp"
#include "element_block.hpp"
#include "initial_state_cartesian.hpp"
#include "write_screen.hpp"
#include "domain_map.hpp"
#include "geometry_labels.hpp"
#include "spatial_metric.hpp"

using std::cout;
using std::endl;
    
/* Called by the ParamsCartesian constructor */
void ParamsCartesian::secondary_params()
{
    Ngroup       = 1;
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

    cout << "***********************************************" << endl;
    cout << endl;

    return;
}


void ParamsCartesian::setup_process(Process &proc)
{
    /* This process's set of indices in the global 3D array of processes   */
    int proc_idx[3];

    proc_idx[2] = proc.rank % Nproc[2];
    proc_idx[1] = (proc.rank / Nproc[2]) % Nproc[1];
    proc_idx[0] = proc.rank / (Nproc[1]*Nproc[2]);

    /* A Cartesian dataset needs only one Cartesian group */
    proc.Ngroup = 1;
    proc.group  = 0;
    proc.group_rank = proc.rank;

    for (int i: dirs)
    {
        proc.group_idx[i]   = proc_idx[i];
        proc.Nproc_group[i] = Nproc[i];

        proc.elements.group_idx[i]   = proc_idx[i];
        proc.elements.Nproc_group[i] = Nproc[i];
    }
    
    /* Used by mesh and data output routines. */
    MPI_Comm_dup(MPI_COMM_WORLD, &proc.group_comm);

    proc.Nproc_group_tot = proc.Nproc_group[0] * proc.Nproc_group[1] * proc.Nproc_group[2];

    setup_elementblock(proc.elements, proc);
    
    set_initial_state(proc.elements, proc.physics);

    write::variable<real_t>("CFL", proc.cfl);
    write::variable<real_t>("End time", proc.end_time);


    /* Setup faces */
    for (int i: ifaces)
    {
        FaceCommunicator& face = proc.faces[i];

        write::variable<int>("Setting up face:", i);

        face.setup(proc, i);

        /* Inter-process connectivity */
        /* Make periodic for now */
        int normal = face.normal_dir;
        for (int d: dirs)
            face.neighbour_idx[d] = proc.group_idx[d];
        
        face.neighbour_idx[normal] += face.orientation;

        /* Make periodic */
        if (face.neighbour_idx[normal] < 0)
            face.neighbour_idx[normal] = proc.Nproc_group[normal] - 1;
        
        if (face.neighbour_idx[normal] > proc.Nproc_group[normal] - 1)
            face.neighbour_idx[normal] = 0;


        face.neighbour_rank = (face.neighbour_idx[0]  * proc.Nproc_group[1]
                             + face.neighbour_idx[1]) * proc.Nproc_group[2]
                             + face.neighbour_idx[2];
    }

    write::message("Finished setup_process()");

    return;
}


void ParamsCartesian::setup_elementblock(ElementBlock &elements, Process &proc)
{
    for (int i: dirs)
    {
        elements.Nelem[i] = Nelem[i];
        elements.Ns[i] = Ns[i];
    }

    elements.Nelem_block = Nelem_proc; 
    elements.Nfield      = proc.Nfield;

    //elements.map = new BasicRect2D;
    elements.map = new WaveRect2D; //new QuarterAnnulusMap; // specify manually for now...

    elements.geometry.metric_s = new DiagonalSpatialMetric(cartesian);
    for (int d: dirs)
        elements.geometry.metric_f[d] = new DiagonalSpatialMetric(cartesian);

    /* At this point all external information is present, and the internal
     * setup method can take over. */
    elements.setup();

    return;
}


void ParamsCartesian::set_initial_state(ElementBlock &elements, Physics* physics)
{
    set_initial_state_cartesian(elements, physics);

    return;
}
