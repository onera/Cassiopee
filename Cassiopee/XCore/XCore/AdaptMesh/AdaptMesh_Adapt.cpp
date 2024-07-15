#include "Mesh.h"
#include "Hexa.h"

static
void cache_prerefinement_data(Mesh *M)
{
    M->gnf_old = 0;
    
    Int ndup = 0, nfree = 0;
    for (Int i = 0; i < M->nf; i++) {
        if (M->face_to_ppatch.find(i) != M->face_to_ppatch.end()) {
            ndup++;
        } else {
            nfree++;
        }
    }

    Int gnfree, gndup;
    MPI_Allreduce(&nfree, &gnfree, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ndup, &gndup, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
    M->gnf_old = gnfree + gndup/2; 

    M->nc_old = M->nc;
    M->nf_old = M->nf;
}

static
void print_postrefinement_data(Mesh *M)
{
    Int gnc = 0;
    MPI_Allreduce(&M->nc, &gnc, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (M->pid == 0) {
        printf("    Total cells after refinement: %d\n", gnc);
    }

    Float balanced = gnc / (Float) M->npc;
    Float my_imbalance = fabs((M->nc - balanced) / (Float)balanced * 100.0);
    Float max_imbalance;
    MPI_Allreduce(&my_imbalance, &max_imbalance, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (M->pid == 0) printf("    Max cell imbalance: %.2f%%\n", max_imbalance);
}

PyObject *K_XCORE::AdaptMesh_Adapt(PyObject *self, PyObject *args)
{
    PyObject *MESH;

    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "AdaptMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(MESH, "AdaptMesh");

    if (M->pid == 0) puts("Adapting...");

    // Orient mesh
    if (M->pid == 0) puts("    Orienting mesh...");
    Mesh_set_orientation(M);

    // Isolate cells/faces/edges to be refined
    
    std::vector<Int> ref_cells, ref_faces;
    std::set<UEdge> ref_edges;

    if (M->pid == 0) puts("    Isolating ref entities...");
    Mesh_get_ref_entities(M, ref_cells, ref_faces, ref_edges);

    // Allocate for refinement

    if (M->pid == 0) puts("    Resizing for refinement...");
    Mesh_resize_for_refinement(M, ref_cells, ref_faces, ref_edges);

    // Refine lower-level cells/faces first
    if (M->pid == 0) puts("    Sorting ref entities by level...");
    Mesh_sort_ref_entities_by_level(M, ref_cells, ref_faces, ref_edges);

    // Cache the pre-refinement number of cells and faces

    cache_prerefinement_data(M);

    if (M->pid == 0) puts("    Refining...");
    Mesh_refine(M, ref_cells, ref_faces, ref_edges);

    print_postrefinement_data(M);

    if (M->pid == 0) puts("    Updating global cell ids...");
    Mesh_update_global_cell_ids(M);

    if (M->pid == 0) puts("    Updating comm patches...");
    Mesh_update_ppatches(M);

    if (M->pid == 0) puts("    Updating boundary patches...");
    Mesh_update_bpatches(M);

    if (M->pid == 0) puts("    Updating global face ids...");
    Mesh_update_global_face_ids(M);

    if (M->pid == 0) puts("    Conformizing face edges...");
    Mesh_conformize_face_edge(M);
 
    //if (M->pid == 0) puts("    Exporting CGNS array...");
    //PyObject *karray = Mesh_export_karray(M);

    if (M->pid == 0) puts("    Done.");

    return Py_None;
    //return karray;
}
