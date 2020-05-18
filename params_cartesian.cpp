#include "params_cartesian.hpp"
#include "process.hpp"

using std::cout;
using std::endl;
    
void ParamsCartesian::secondary_params()
{
    Nproc_domain = Nproc[0] * Nproc[1] * Nproc[2];
    Nelem_proc   = Nelem[0] * Nelem[1] * Nelem[2];
    Nelem_domain = Nelem_proc * Nproc_domain;
    Ns_elem      = Ns[0] * Ns[1] * Ns[2];
    Ns_domain    = Ns_elem * Nelem_domain;

    Nf[0]   = Ns[0]+1;
    Nf[1]   = Ns[1]+1;
    Nf[2]   = Ns[2]+1;
    Nf_elem = Ns[0]*Ns[1]*Nf[2] + Ns[0]*Ns[2]*Nf[1] + Ns[1]*Ns[2]*Nf[0]; // 3 Nf Ns^2 if all equal
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
    /* This process's set of indices in the global 3D array of processes   *
     * Assume row-major flattening order : eg ElementBlock::sidx() for CPU */
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


    for (int i: ifaces)
    {
        proc.faces[i].my_rank = proc.rank;
        proc.faces[i].my_idx  = i;
        proc.faces[i].external_face = false; // Fully periodic for now...
        proc.faces[i].N[0] = Ns[face_coords[i][0]]; // # soln points on the face
        proc.faces[i].N[1] = Ns[face_coords[i][1]];
    }

    /* Inter-process connectivity */
    // Make periodic for now
    int proc_base;

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


    proc.faces[0].neighbour_idx = 1; // Simple globally Cartesian connectivity 
    proc.faces[1].neighbour_idx = 0; 
    proc.faces[2].neighbour_idx = 3; 
    proc.faces[3].neighbour_idx = 2; 
    proc.faces[4].neighbour_idx = 5; 
    proc.faces[5].neighbour_idx = 4; 


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
        proc_origin[i] = proc_idx[i] * length_per_proc[i];

    for (int i: icorners)
        for (int j: dirs)
            proc.corners[i][j] = proc_origin[j] + corner_coords[i][j]*length_per_proc[j];


    setup_elementblock(proc.elements, proc);

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

    elements.Ns_elem = Ns_elem; // Should these be set here rather
    elements.Nf_elem = Nf_elem; // than by params object?
    elements.Nelem_block = Nelem_proc; 
    elements.Nfield  = Nfield;


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
