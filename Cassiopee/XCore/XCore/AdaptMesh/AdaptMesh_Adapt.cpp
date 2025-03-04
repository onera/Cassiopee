/*    
    Copyright 2013-2025 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "Mesh.h"
#include "Hexa.h"

static
void cache_prerefinement_data(Mesh *M)
{
    M->gnf_old = 0;
    
    E_Int ndup = 0, nfree = 0;
    for (E_Int i = 0; i < M->nf; i++) {
        if (M->face_to_ppatch.find(i) != M->face_to_ppatch.end()) {
            ndup++;
        } else {
            nfree++;
        }
    }

    E_Int gnfree, gndup;
    MPI_Allreduce(&nfree, &gnfree, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ndup, &gndup, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
    M->gnf_old = gnfree + gndup/2; 

    M->nc_old = M->nc;
    M->nf_old = M->nf;
}

static
void print_postrefinement_data(Mesh *M)
{
    E_Int gnc = M->nc;
    MPI_Allreduce(&M->nc, &gnc, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    E_Int gnf = Mesh_get_global_face_count(M);

    if (M->pid == 0) {
        printf("    Total cells after refinement: " SF_D_ "\n", gnc);
        printf("    Total faces after refinement: " SF_D_ "\n", gnf);
    }

    if (M->npc == 1) return;

    E_Float balanced = gnc / (E_Float) M->npc;
    E_Float my_imbalance = fabs((M->nc - balanced) / (E_Float)balanced * 100.0);
    E_Float max_imbalance;
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

    // Isolate cells/faces/edges to be refined
    
    std::vector<E_Int> ref_cells, ref_faces;
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

    if (M->pid == 0) puts("    Done.");

    return Py_None;
}
