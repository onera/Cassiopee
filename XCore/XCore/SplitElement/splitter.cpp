//#if defined(USE_STLMAP)
#include <unordered_map>
#include <unordered_set>

#include "kPython.h"
using std::unordered_map;
using std::unordered_set;
//#else
//#include "Map/robin_hood.h"
// using robin_hood::unordered_map;
//#endif

#include "Memory/vector_view.hpp"
#include "Numpy/Numpy.h"
#include "Numpy/vector.hpp"
#include "SplitElement/splitter.h"
#include "String/kstring.h"
#include "paradigma/mesh/pdm_dmesh_nodal.h"
#include "paradigma/ppart/pdm_part.h"
#include "xmpi/context.hpp"
#include "xmpi/xmpi.hpp"

#if PY_MAJOR_VERSION >= 3
#define CString2PyString PyUnicode_FromString
#else
#define CString2PyString PyString_FromString
#endif

namespace
{
    std::pair<int, int>
    section_type(const char *eltType)
    {
        if (K_STRING::cmp(eltType, "TRI") == 0 || K_STRING::cmp(eltType, "TRI*") == 0)
            return {PDM_MESH_NODAL_TRIA3, 3};
        else if (K_STRING::cmp(eltType, "QUAD") == 0 || K_STRING::cmp(eltType, "QUAD*") == 0)
            return {PDM_MESH_NODAL_QUAD4, 4};
        else if (K_STRING::cmp(eltType, "TETRA") == 0 || K_STRING::cmp(eltType, "TETRA*") == 0)
            return {PDM_MESH_NODAL_TETRA4, 4};
        else if (K_STRING::cmp(eltType, "PYRA") == 0 || K_STRING::cmp(eltType, "PYRA*") == 0)
            return {PDM_MESH_NODAL_PYRAMID5, 5};
        else if (K_STRING::cmp(eltType, "PENTA") == 0 || K_STRING::cmp(eltType, "PENTA*") == 0)
            return {PDM_MESH_NODAL_PRISM6, 6};
        else if (K_STRING::cmp(eltType, "HEXA") == 0 || K_STRING::cmp(eltType, "HEXA*") == 0)
            return {PDM_MESH_NODAL_HEXA8, 8};
        else if (K_STRING::cmp(eltType, "BAR") == 0 || K_STRING::cmp(eltType, "BAR*") == 0)
            return {PDM_MESH_NODAL_BAR2, 2};
        else
            throw(std::runtime_error("Bad type of element"));
        return {-1, -1};
    }
    // ----------------------------------------------------------------------------------------
    int
    nb_verts_for_element(int pdm_elt_type)
    {
        switch (pdm_elt_type) {
            case PDM_MESH_NODAL_BAR2:
                return 2;
                break;
            case PDM_MESH_NODAL_TRIA3:
                return 3;
                break;
            case PDM_MESH_NODAL_QUAD4:
                return 4;
                break;
            case PDM_MESH_NODAL_TETRA4:
                return 4;
                break;
            case PDM_MESH_NODAL_PYRAMID5:
                return 5;
                break;
            case PDM_MESH_NODAL_PRISM6:
                return 6;
                break;
            case PDM_MESH_NODAL_HEXA8:
                return 8;
                break;
            default:
                return -1;
        }
        return -1;
    }
}  // namespace

static unordered_map<int, std::vector<PDM_g_num_t>> elt2global_vertices;

struct type_element_connectivity {
    int                               type_element, nb_verts_per_elt;
    KCore::Numpy::vector<PDM_g_num_t> elt2vertsglob;
    type_element_connectivity(int tp_elt, int nb_v_per_e, PDM_g_num_t *e2v_glob, int nb_elt)
      : type_element(tp_elt), nb_verts_per_elt(nb_v_per_e), elt2vertsglob{nb_elt}
    {
        for (int i = 0; i < nb_elt; ++i) elt2vertsglob[i] = e2v_glob[i];
    }
};

struct zone_connectivity_t {
    E_Int                                  nb_vertices;          /* Nombre de sommets contenu dans une zone */
    E_Int                                  nb_elements;          /* Nombre total d'éléments contenus dans la zone */
    E_Int                                  first_ind_glob_vert;  // Premier index global des sommets pour cette zone
    E_Int                                  first_ind_glob_elts;  // Premier index global des éléments pour cette zone
    std::vector<type_element_connectivity> per_element_connectivities;  // Connectivité par élément
    KCore::Numpy::vector<double> vertices_coordinates; // 3 * nb_vertices en taille, stocké (x0,y0,z0,x1,y1,z1,...)
};

using domain_connectivity_t = std::vector<zone_connectivity_t>;

struct py_splitted_zone {
    PyObject_HEAD KCore::Numpy::vector<int> cellTag, cellFaceIdx, cellFace;
    KCore::Numpy::vector<int>               faceTag, faceCell, faceVertexIdx, faceVertex;
    KCore::Numpy::vector<int>               facePartBound, facePartBoundProcIdx, facePartBoundPartIdx;
    KCore::Numpy::vector<int>               vertexTag, faceGroupIdx, faceGroup;
    KCore::Numpy::vector<PDM_g_num_t>       cell_loc2glob, face_loc2glob, vertex_loc2glob, faceGroup_loc2glob;
    KCore::Numpy::vector<double>            vertex;
};
// ########################################################################################
static PyObject *
build_splitted_zone(int hdl_part, int ipart)
{
    xcore::communicator comm = xcore::context::globalCommunicator();
    int nCell, nProc, nTPart, nFace, nFacePartBound, nVertex, sCellFace, sFaceVertex, sFaceGroup, nFaceGroup;

    PDM_part_part_dim_get(hdl_part, ipart, &nCell, &nFace, &nFacePartBound, &nVertex, &nProc, &nTPart, &sCellFace,
                          &sFaceVertex, &sFaceGroup, &nFaceGroup);
    int *        cellTag, *cellFaceIdx, *cellFace;
    PDM_g_num_t *cellLNToGN, *faceLNToGN, *vertexLNToGN, *faceGroupLNToGN;
    int *        faceTag, *faceCell, *faceVertexIdx, *faceVertex;
    int *        facePartBound, *facePartBoundProcIdx, *facePartBoundPartIdx;
    int *        vertexTag, *faceGroupIdx, *faceGroup;
    double *     vertex;

    PDM_part_part_val_get(hdl_part, ipart, &cellTag, &cellFaceIdx, &cellFace, &cellLNToGN, &faceTag, &faceCell,
                          &faceVertexIdx, &faceVertex, &faceLNToGN, &facePartBoundProcIdx, &facePartBoundPartIdx,
                          &facePartBound, &vertexTag, &vertex, &vertexLNToGN, &faceGroupIdx, &faceGroup,
                          &faceGroupLNToGN);
    PyObject *splt_zone = PyDict_New();
    if (cellTag) 
    {
        KCore::Numpy::vector<int> pcellTag(cellTag, cellTag + nCell);
        PyDict_SetItem(splt_zone, CString2PyString("cellTag"), (PyObject *)KCore::Numpy::shared_with_python(pcellTag));
    }
    if (cellFaceIdx) {
        KCore::Numpy::vector<int> pCellFaceIdx(cellFaceIdx, cellFaceIdx + nCell + 1);
        PyDict_SetItem(splt_zone, CString2PyString("cellFaceIdx"),
                       (PyObject *)KCore::Numpy::shared_with_python(pCellFaceIdx));
    }
    if (cellFace) {
        KCore::Numpy::vector<int> pCellFace(cellFace, cellFace + cellFaceIdx[nCell]);
        PyDict_SetItem(splt_zone, CString2PyString("cellFace"),
                       (PyObject *)KCore::Numpy::shared_with_python(pCellFace));
    }
    if (cellLNToGN) {
        KCore::Numpy::vector<PDM_g_num_t> pCellLN2GN(cellLNToGN, cellLNToGN + nCell);
        PyDict_SetItem(splt_zone, CString2PyString("cellLN2GN"),
                       (PyObject *)KCore::Numpy::shared_with_python(pCellLN2GN));
    }
    if (faceLNToGN) {
        KCore::Numpy::vector<PDM_g_num_t> pFaceLNToGN(faceLNToGN, faceLNToGN + nFace);
        PyDict_SetItem(splt_zone, CString2PyString("faceLN2GN"),
                       (PyObject *)KCore::Numpy::shared_with_python(pFaceLNToGN));
    }
        if (vertexLNToGN) {
            KCore::Numpy::vector<PDM_g_num_t> pVertexLN2GN(vertexLNToGN, vertexLNToGN + nVertex);
            PyDict_SetItem(splt_zone, CString2PyString("vertexLN2GN"),
                           (PyObject *)KCore::Numpy::shared_with_python(pVertexLN2GN));
        }
        if (faceTag) {
            KCore::Numpy::vector<int> pFaceTag(faceTag, faceTag + nFace);
            PyDict_SetItem(splt_zone, CString2PyString("faceTag"),
                           (PyObject *)KCore::Numpy::shared_with_python(pFaceTag));
        }
        if (faceCell) {
            KCore::Numpy::vector<int> pFaceCell(faceCell, faceCell + 2 * nFace);
            PyDict_SetItem(splt_zone, CString2PyString("faceCell"),
                           (PyObject *)KCore::Numpy::shared_with_python(pFaceCell));
        }
        if (faceVertexIdx) {
            KCore::Numpy::vector<int> pFaceVertexIdx(faceVertexIdx, faceVertexIdx + nFace + 1);
            PyDict_SetItem(splt_zone, CString2PyString("faceVertexIdx"),
                           (PyObject *)KCore::Numpy::shared_with_python(pFaceVertexIdx));
        }
        if (faceVertex) {
            KCore::Numpy::vector<int> pFaceVertex(faceVertex, faceVertex + faceVertexIdx[nFace]);
            PyDict_SetItem(splt_zone, CString2PyString("faceVertex"),
                           (PyObject *)KCore::Numpy::shared_with_python(pFaceVertex));
        }
        if (facePartBoundProcIdx) {
            KCore::Numpy::vector<int> pFacePartBoundProcIdx(facePartBoundProcIdx, facePartBoundProcIdx + comm.size + 1);
            PyDict_SetItem(splt_zone, CString2PyString("facePartBoundProcIdx"),
                           (PyObject *)KCore::Numpy::shared_with_python(pFacePartBoundProcIdx));
        }
        if (facePartBoundPartIdx) {
            KCore::Numpy::vector<int> pFacePartBoundPartIdx(facePartBoundPartIdx, facePartBoundPartIdx + nTPart + 1);
            PyDict_SetItem(splt_zone, CString2PyString("facePartBoundPartIdx"),
                           (PyObject *)KCore::Numpy::shared_with_python(pFacePartBoundPartIdx));
        }
        if (facePartBound) {
            KCore::Numpy::vector<int> pFacePartBound(facePartBound, facePartBound + 4 * nFacePartBound);
            PyDict_SetItem(splt_zone, CString2PyString("facePartBound"),
                           (PyObject *)KCore::Numpy::shared_with_python(pFacePartBound));
        }
        if (vertexTag) {
            KCore::Numpy::vector<int> pVertexTag(vertexTag, vertexTag + nVertex);
            PyDict_SetItem(splt_zone, CString2PyString("vertexTag"),
                           (PyObject *)KCore::Numpy::shared_with_python(pVertexTag));
        }
        if (vertex) {
            KCore::Numpy::vector<double> pVertex({nVertex,3});
            memcpy(pVertex.data(), vertex, 3*nVertex*sizeof(double));
            PyDict_SetItem(splt_zone, CString2PyString("vertex"),
                           (PyObject *)KCore::Numpy::shared_with_python(pVertex));
        }
        if (faceGroupIdx) {
            KCore::Numpy::vector<int> pfaceGroupIdx(faceGroupIdx, faceGroupIdx + nFaceGroup + 1);
            PyDict_SetItem(splt_zone, CString2PyString("faceGroupIdx"),
                           (PyObject *)KCore::Numpy::shared_with_python(pfaceGroupIdx));
        }
        if (faceGroup) {
            KCore::Numpy::vector<int> pFaceGroup(faceGroup, faceGroup + sFaceGroup);
            PyDict_SetItem(splt_zone, CString2PyString("faceGroup"),
                           (PyObject *)KCore::Numpy::shared_with_python(pFaceGroup));
        }
        if (faceGroupLNToGN) {
            KCore::Numpy::vector<PDM_g_num_t> pFaceGroupLNToGN(faceGroupLNToGN, faceGroupLNToGN + sFaceGroup);
            PyDict_SetItem(splt_zone, CString2PyString("faceGroupLN2GN"),
                           (PyObject *)KCore::Numpy::shared_with_python(pFaceGroupLNToGN));
        }
        return Py_BuildValue("O", splt_zone);
    }
    // ----------------------------------------------------------------------------------------
    PyObject *split_elements(PyObject * self, PyObject * args)
    {
        int       nbparts = 1;
        PyObject *set_of_zones;
        // Attention, nbparts => nombre de partitions par procs !!!!
        // Découpage se base sur le nombre de processus * nombre de partitions
        if (not PYPARSETUPLEI(args, "O|l", "O|i", &set_of_zones, &nbparts)) return NULL;
        if (not PyList_Check(set_of_zones)) {
            PyErr_SetString(PyExc_RuntimeError, "Set of Zones must be stored in a python list");
            return NULL;
        }
        Py_ssize_t nb_zones = PyList_Size(set_of_zones);

        domain_connectivity_t domain;
        domain.reserve(nb_zones);
        xcore::communicator comm = xcore::context::globalCommunicator();

        // Pour chaque zone du domaine de calcul :
        for (Py_ssize_t i_zone = 0; i_zone < nb_zones; ++i_zone) {
            // Récupération d'une zone courante : la connectivité de la zone + nombre total de sommets dans la zone,
            // nombre total d'éléments dans la zone
            PyObject *pt_pyzone = PyList_GetItem(set_of_zones, i_zone);
            if (not PyTuple_Check(pt_pyzone)) {
                PyErr_SetString(PyExc_TypeError,
                                "Each element of the zone list must be a tuple : (zone, number of "
                                "vertices, number of elements)");
                return NULL;
            }

            zone_connectivity_t zconnect;

            PyObject *py_zone;
            E_Int     nb_total_vertices, nb_total_elements;
            if (!PYPARSETUPLEI(pt_pyzone, "Oll", "Oii", &py_zone, &nb_total_vertices, &nb_total_elements)) {
                PyErr_SetString(PyExc_TypeError,
                                "The tuple must contains three values : the zone connectivity, "
                                "the total number of vertices and total number of elements on the zone");
                return NULL;
            }
            if (not PyList_Check(py_zone)) {
                PyErr_SetString(PyExc_TypeError,
                                "The zone object must be a list of connectivities per type of element");
                return NULL;
            }
            zconnect.nb_vertices = nb_total_vertices;
            zconnect.nb_elements = nb_total_elements;
            zconnect.per_element_connectivities.reserve(PyList_Size(py_zone));

            Py_ssize_t nb_elt_type = PyList_Size(py_zone);  // Nombre de types d'éléments dans la zone
            // Pour chaque type d'éléments contenus dans la zone :
            for (Py_ssize_t itype = 0; itype < nb_elt_type; ++itype) {
                // Récupération d'un type d'élément : : ("Type d'élément (chaîne de caractère)", connectivité élt ->
                // sommet,  )
                PyObject *pt_pyelt_type = PyList_GetItem(py_zone, itype);
                if (not PyTuple_Check(pt_pyelt_type)) {
                    PyErr_SetString(
                        PyExc_TypeError,
                        "Each element of the elements list must be a tuple (\"element type\", connectivity");
                    return NULL;
                }
                // Array contains the kind of element...
                const char *   elt_type;
                PyArrayObject *pt_e2v, *pt_coords;
                if (!PYPARSETUPLEI(pt_pyelt_type, "sOO", "sOO", &elt_type, &pt_e2v, &pt_coords)) {
                    PyErr_SetString(PyExc_TypeError,
                                    "The tuple must contain a string given type of element, "
                                    "and the connectivity array (numpy)");
                    return NULL;
                }
                int type_elt, nb_verts_per_elt;
                std::tie(type_elt, nb_verts_per_elt) = section_type(elt_type);
                E_Int  size, nfld;
                E_Int *e2v;
                E_Int coord_size, coord_nfld;
                if (not K_NUMPY::getFromNumpyArray((PyObject *)pt_e2v, e2v, size, nfld)) {
                    PyErr_SetString(PyExc_TypeError,
                                    "The element to vertices array must be a "
                                    "numpy array containing E_Int values");
                    return NULL;
                }
                assert(nfld==1);
                PDM_g_num_t *pt_elt2vert;
                if (sizeof(E_Int) == sizeof(PDM_g_num_t))
                    pt_elt2vert = (PDM_g_num_t *)e2v;
                else {
                    //std::cout << "Copie d'indices avec representation entiere differente" << std::endl;
                    pt_elt2vert = new PDM_g_num_t[size];
                    std::copy(e2v, e2v + size, pt_elt2vert);
                }
                zconnect.per_element_connectivities.emplace_back(type_elt, nb_verts_per_elt, pt_elt2vert, int(size));
                zconnect.vertices_coordinates = std::move(KCore::Numpy::vector<double>(pt_coords));
            }  // Fin pour itype
            domain.push_back(zconnect);
        }  // Fin pour i_zone
        // On va compter le nombre total de sommets, le nombre total d'éléments et renuméroter globalement les
        // sommets et les éléments :
        int nb_tot_vertices = 0, nb_tot_elements = 0;
        int nb_loc_vertices = 0, nb_loc_elements = 0;  // dnVerts et dnCells
        for (auto &zone : domain) {
            zone.first_ind_glob_vert = nb_tot_vertices;
            zone.first_ind_glob_elts = nb_tot_elements;
            nb_tot_elements += zone.nb_elements;
            nb_tot_vertices += zone.nb_vertices;
            for (auto &elts_type : zone.per_element_connectivities)
                for (auto &index : elts_type.elt2vertsglob) index += zone.first_ind_glob_vert;
        }
        double *     dVtxCoord   = new double[3 * nb_tot_vertices];  // NULL;
        int ind_vert = 0;
        int nb_loc_stored_verts = 0;
        for (auto &zone : domain) {
            double* coords;
            E_Int coord_size, coord_nfld;
            PyArrayObject* py_coords = KCore::Numpy::shared_with_python(zone.vertices_coordinates);
            if (not K_NUMPY::getFromNumpyArray((PyObject *)py_coords, coords, coord_size, coord_nfld))
            {
                PyErr_SetString(PyExc_TypeError, 
                    R"RAW(The coordinates of vertices array must be a numpy array containing double values)RAW");
                return NULL;                        
            }
            nb_loc_stored_verts += coord_size;
            for ( E_Int i = 0; i < coord_size; ++i )
            {
                dVtxCoord[3*ind_vert+3*i+0] = coords[i + 0*coord_size];
                dVtxCoord[3*ind_vert+3*i+1] = coords[i + 1*coord_size];
                dVtxCoord[3*ind_vert+3*i+2] = coords[i + 2*coord_size];
            }
            ind_vert += zone.nb_vertices;
            for (auto &elts_type : zone.per_element_connectivities) {
                unordered_set<PDM_g_num_t> set_of_vertices;
                set_of_vertices.insert(elts_type.elt2vertsglob.begin(), elts_type.elt2vertsglob.end());
                //===========================================================================
                // int                        type_element, nb_verts_per_elt;
                // K_MEMORY::vector_view<PDM_g_num_t> elt2vertsglob;
                //===========================================================================
                nb_loc_elements += elts_type.elt2vertsglob.size() / elts_type.nb_verts_per_elt;
                nb_loc_vertices += set_of_vertices.size();
                if (elt2global_vertices.find(elts_type.type_element) == elt2global_vertices.end())
                    elt2global_vertices[elts_type.type_element] =
                        std::vector<PDM_g_num_t>(elts_type.elt2vertsglob.begin(), elts_type.elt2vertsglob.end());
                else
                    elt2global_vertices[elts_type.type_element].insert(
                        elt2global_vertices[elts_type.type_element].end(), elts_type.elt2vertsglob.begin(),
                        elts_type.elt2vertsglob.end());
            }
        }
        // Construction de la structure DMesh de Paradigma pour le découpage :
        xcore::Ext_Communicator &ext_comm = comm.get_implementation();
        PDM_MPI_Comm             pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&ext_comm);
        int hdl_dmesh = PDM_DMesh_nodal_create(pdm_comm, nb_tot_vertices, nb_tot_elements);

        for (auto &elt_type : elt2global_vertices) {
            int type_elt = elt_type.first;
            int nb_vert_per_elt = nb_verts_for_element(type_elt);
            int id_section = PDM_DMesh_nodal_section_add(hdl_dmesh, (PDM_Mesh_nodal_elt_t)type_elt);
            PDM_DMesh_nodal_section_std_set(hdl_dmesh, id_section, elt_type.second.size() / nb_vert_per_elt,
                                            elt_type.second.data());
        }
        PDM_DMesh_nodal_coord_set(hdl_dmesh, nb_loc_stored_verts, dVtxCoord);        
        comm.barrier();
        PDM_DMesh_nodal_cell_face_compute(hdl_dmesh);

        int       hdl_ppart;
        const int nb_properties_for_cells = 3;
        const int nCellPerCache = 1024, is_asynchrone = 1, is_vectorisation = 1;
        int       cell_properties[nb_properties_for_cells] = {nCellPerCache, is_asynchrone, is_vectorisation};

        int  nFaceGroupe   = 0;  // Nombre de types de conditition limite dont les raccords...
        int *dFaceGroupIdx = NULL, *dFaceGroup = NULL;

        int have_dCellPart = 0;
        // nb local d'éléments
        std::vector<int> dCellPart(nb_tot_elements);

        int nb_face_tot = PDM_DMesh_nodal_total_n_face_get(hdl_dmesh);
        PDM_g_num_t *dFaceCell;
        int nb_faces = PDM_DMesh_nodal_face_cell_get(hdl_dmesh, &dFaceCell);
        int nb_vertices = PDM_DMesh_nodal_n_vtx_get(hdl_dmesh);
        const double* d_coords = PDM_DMesh_nodal_vtx_get(hdl_dmesh); 
        int *        dFaceVtxIdx;
        PDM_g_num_t *dFaceVtx    = NULL;
        int          nb_face_loc = PDM_DMesh_nodal_face_vtx_get(hdl_dmesh, &dFaceVtxIdx, &dFaceVtx);
        PDM_part_create(&hdl_ppart, pdm_comm, PDM_PART_SPLIT_PTSCOTCH, "PDM_PART_RENUM_CELL_NONE",
                        "PDM_PART_RENUM_FACE_NONE", nb_properties_for_cells, cell_properties, 0,
                        NULL,  // Two last for face properties ( eq. to cell properties but for faces )
                        nbparts, nb_loc_elements, nb_face_loc, nb_vertices, nFaceGroupe, NULL, NULL, NULL, NULL,
                        have_dCellPart, dCellPart.data(), dFaceCell, dFaceVtxIdx, dFaceVtx, 
                        NULL, d_coords, NULL, dFaceGroupIdx, dFaceGroup);

        PyObject *py_part_list = PyList_New(nbparts);
        for (int i = 0; i < nbparts; ++i)
            PyList_SetItem(py_part_list, i, (PyObject *)build_splitted_zone(hdl_ppart, i));
        delete[] dVtxCoord;
        PDM_DMesh_nodal_free(hdl_dmesh);
        return py_part_list;
    }  // Fin fonction split_elements
