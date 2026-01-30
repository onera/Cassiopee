/*
    Copyright 2013-2026 ONERA.

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
#include <unordered_map>
#include <unordered_set>
#include <array>

#include "kPython.h"
using std::unordered_map;
using std::unordered_set;

// stub
PyObject *splitElements(PyObject *self, PyObject *args)
{
    Py_INCREF(Py_None);
    return Py_None;
}
// #include "Memory/vector_view.hpp"
// #include "Numpy/Numpy.h"
// #include "Numpy/vector.hpp"
// #include "SplitElement/splitter.h"
// #include "String/kstring.h"
// #include "paradigma/mesh/pdm_dmesh_nodal.h"
// #include "paradigma/ppart/pdm_part.h"
// #include "xmpi/context.hpp"
// #include "xmpi/xmpi.hpp"

// #if PY_MAJOR_VERSION >= 3
// #define CString2PyString PyUnicode_FromString
// #else
// #define CString2PyString PyString_FromString
// #endif

// //#   define SPLITTER_TRACE

// //namespace
// //{
//     std::pair<int, int>
//     section_type( const char *eltType )
//     {
//         if ( K_STRING::cmp( eltType, "TRI" ) == 0 || K_STRING::cmp( eltType, "TRI*" ) == 0 )
//             return {PDM_MESH_NODAL_TRIA3, 3};
//         else if ( K_STRING::cmp( eltType, "QUAD" ) == 0 || K_STRING::cmp( eltType, "QUAD*" ) == 0 )
//             return {PDM_MESH_NODAL_QUAD4, 4};
//         else if ( K_STRING::cmp( eltType, "TETRA" ) == 0 || K_STRING::cmp( eltType, "TETRA*" ) == 0 )
//             return {PDM_MESH_NODAL_TETRA4, 4};
//         else if ( K_STRING::cmp( eltType, "PYRA" ) == 0 || K_STRING::cmp( eltType, "PYRA*" ) == 0 )
//             return {PDM_MESH_NODAL_PYRAMID5, 5};
//         else if ( K_STRING::cmp( eltType, "PENTA" ) == 0 || K_STRING::cmp( eltType, "PENTA*" ) == 0 )
//             return {PDM_MESH_NODAL_PRISM6, 6};
//         else if ( K_STRING::cmp( eltType, "HEXA" ) == 0 || K_STRING::cmp( eltType, "HEXA*" ) == 0 )
//             return {PDM_MESH_NODAL_HEXA8, 8};
//         else if ( K_STRING::cmp( eltType, "BAR" ) == 0 || K_STRING::cmp( eltType, "BAR*" ) == 0 )
//             return {PDM_MESH_NODAL_BAR2, 2};
//         else
//             throw( std::runtime_error( "Bad type of element" ) );
//         return {-1, -1};
//     }
//     // ----------------------------------------------------------------------------------------
//     int
//     nb_verts_for_element( int pdm_elt_type )
//     {
//         switch ( pdm_elt_type ) {
//             case PDM_MESH_NODAL_BAR2:
//                 return 2;
//                 break;
//             case PDM_MESH_NODAL_TRIA3:
//                 return 3;
//                 break;
//             case PDM_MESH_NODAL_QUAD4:
//                 return 4;
//                 break;
//             case PDM_MESH_NODAL_TETRA4:
//                 return 4;
//                 break;
//             case PDM_MESH_NODAL_PYRAMID5:
//                 return 5;
//                 break;
//             case PDM_MESH_NODAL_PRISM6:
//                 return 6;
//                 break;
//             case PDM_MESH_NODAL_HEXA8:
//                 return 8;
//                 break;
//             default:
//                 return -1;
//         }
//         return -1;
//     }
// //}  // namespace

// static unordered_map<int, std::vector<PDM_g_num_t>> elt2global_vertices;
// static unordered_map<int, std::vector<int>> elttag_global;

// struct type_element_connectivity {
//     int type_element, nb_verts_per_elt, tag;
//     KCore::Numpy::vector<PDM_g_num_t> elt2vertsglob;
//     type_element_connectivity( int tp_elt, int nb_v_per_e, PDM_g_num_t *e2v_glob, int nb_elt, int num_zone )
//       : type_element( tp_elt ), nb_verts_per_elt( nb_v_per_e ), elt2vertsglob{(unsigned int)nb_elt}, tag{num_zone}
//     {
//         for ( int i = 0; i < nb_elt; ++i )
//             elt2vertsglob[ i ] = e2v_glob[ i ];

//     }
// };

// struct zone_connectivity_t {
//     E_Int num_zone;
//     E_Int nb_vertices;                                                  /* Nombre de sommets contenu dans une zone */
//     E_Int nb_elements;                                                  /* Nombre total d'éléments contenus dans la zone */
//     E_Int first_ind_glob_vert;                                          // Premier index global des sommets pour cette zone
//     E_Int first_ind_glob_elts;                                          // Premier index global des éléments pour cette zone
//     std::vector<type_element_connectivity> per_element_connectivities;  // Connectivité par élément
//     PyObject *vertices_coordinates;
//     //KCore::Numpy::vector<double> vertices_coordinates; // 3 * nb_vertices en taille, stocké (x0,y0,z0,x1,y1,z1,...)
// };

// using domain_connectivity_t = std::vector<zone_connectivity_t>;

// struct py_splitted_zone {
//     PyObject_HEAD KCore::Numpy::vector<int> cellTag, cellFaceIdx, cellFace;
//     KCore::Numpy::vector<int> faceTag, faceCell, faceVertexIdx, faceVertex;
//     KCore::Numpy::vector<int> facePartBound, facePartBoundProcIdx, facePartBoundPartIdx;
//     KCore::Numpy::vector<int> vertexTag, faceGroupIdx, faceGroup;
//     KCore::Numpy::vector<PDM_g_num_t> cell_loc2glob, face_loc2glob, vertex_loc2glob, faceGroup_loc2glob;
//     KCore::Numpy::vector<double> vertex;
// };
// // ########################################################################################
// static PyObject *
// build_splitted_zone( int hdl_part, int ipart )
// {
//     xcore::communicator comm = xcore::context::globalCommunicator();
//     int nCell, nProc, nTPart, nFace, nFacePartBound, nVertex, sCellFace, sFaceVertex, sFaceGroup, nFaceGroup;

//     PDM_part_part_dim_get( hdl_part, ipart, &nCell, &nFace, &nFacePartBound, &nVertex, &nProc, &nTPart, &sCellFace,
//                            &sFaceVertex, &sFaceGroup, &nFaceGroup );
// #   if defined( SPLITTER_TRACE )
//     std::cout << "Splitted data : nCell =  " << nCell << ", nFace = " << nFace << ", nVertex = " << nVertex << std::endl;
//     std::flush(std::cout);
// #   endif
//     int *cellTag, *cellFaceIdx, *cellFace;
//     PDM_g_num_t *cellLNToGN, *faceLNToGN, *vertexLNToGN, *faceGroupLNToGN;
//     int *faceTag, *faceCell, *faceVertexIdx, *faceVertex;
//     int *facePartBound, *facePartBoundProcIdx, *facePartBoundPartIdx;
//     int *vertexTag, *faceGroupIdx, *faceGroup;
//     double *vertex;

//     PDM_part_part_val_get( hdl_part, ipart, &cellTag, &cellFaceIdx, &cellFace, &cellLNToGN, &faceTag, &faceCell,
//                            &faceVertexIdx, &faceVertex, &faceLNToGN, &facePartBoundProcIdx, &facePartBoundPartIdx,
//                            &facePartBound, &vertexTag, &vertex, &vertexLNToGN, &faceGroupIdx, &faceGroup,
//                            &faceGroupLNToGN );
// #   if defined( SPLITTER_TRACE )
//     std::cout << "Splitted data => donnees recuperees..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     PyObject *splt_zone = PyDict_New();
//     if ( cellTag ) {
//         KCore::Numpy::vector<int> pcellTag( cellTag, cellTag + nCell );
//         PyDict_SetItem( splt_zone, CString2PyString( "cellTag" ), (PyObject *)KCore::Numpy::shared_with_python( pcellTag ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "cellTag recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( cellFaceIdx ) {
//         KCore::Numpy::vector<int> pCellFaceIdx( cellFaceIdx, cellFaceIdx + nCell + 1 );
//         PyDict_SetItem( splt_zone, CString2PyString( "cellFaceIdx" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pCellFaceIdx ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "cellFaceIdx recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( cellFace ) {
//         KCore::Numpy::vector<int> pCellFace( cellFace, cellFace + cellFaceIdx[ nCell ] );
//         PyDict_SetItem( splt_zone, CString2PyString( "cellFace" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pCellFace ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "cellFace recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( cellLNToGN ) {
//         KCore::Numpy::vector<PDM_g_num_t> pCellLN2GN( cellLNToGN, cellLNToGN + nCell );
//         PyDict_SetItem( splt_zone, CString2PyString( "cellLN2GN" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pCellLN2GN ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "cellLNToGN recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( faceLNToGN ) {
//         KCore::Numpy::vector<PDM_g_num_t> pFaceLNToGN( faceLNToGN, faceLNToGN + nFace );
//         PyDict_SetItem( splt_zone, CString2PyString( "faceLN2GN" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pFaceLNToGN ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "faceLNToGN recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( vertexLNToGN ) {
//         KCore::Numpy::vector<PDM_g_num_t> pVertexLN2GN( vertexLNToGN, vertexLNToGN + nVertex );
//         PyDict_SetItem( splt_zone, CString2PyString( "vertexLN2GN" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pVertexLN2GN ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "vertexLNToGN recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( faceTag ) {
//         KCore::Numpy::vector<int> pFaceTag( faceTag, faceTag + nFace );
//         PyDict_SetItem( splt_zone, CString2PyString( "faceTag" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pFaceTag ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "faceTag recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( faceCell ) {
//         KCore::Numpy::vector<int> pFaceCell( faceCell, faceCell + 2 * nFace );
//         PyDict_SetItem( splt_zone, CString2PyString( "faceCell" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pFaceCell ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "faceCell recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( faceVertexIdx ) {
//         KCore::Numpy::vector<int> pFaceVertexIdx( faceVertexIdx, faceVertexIdx + nFace + 1 );
//         PyDict_SetItem( splt_zone, CString2PyString( "faceVertexIdx" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pFaceVertexIdx ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "faceVertexIdx recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( faceVertex ) {
//         KCore::Numpy::vector<int> pFaceVertex( faceVertex, faceVertex + faceVertexIdx[ nFace ] );
//         PyDict_SetItem( splt_zone, CString2PyString( "faceVertex" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pFaceVertex ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "faceVertex recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( facePartBoundProcIdx ) {
//         KCore::Numpy::vector<int> pFacePartBoundProcIdx( facePartBoundProcIdx, facePartBoundProcIdx + comm.size + 1 );
//         PyDict_SetItem( splt_zone, CString2PyString( "facePartBoundProcIdx" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pFacePartBoundProcIdx ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "facePartBoundProcIdx recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( facePartBoundPartIdx ) {
//         KCore::Numpy::vector<int> pFacePartBoundPartIdx( facePartBoundPartIdx, facePartBoundPartIdx + nTPart + 1 );
//         PyDict_SetItem( splt_zone, CString2PyString( "facePartBoundPartIdx" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pFacePartBoundPartIdx ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "facePartBoundPartIdx recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( facePartBound ) {
//         KCore::Numpy::vector<int> pFacePartBound( facePartBound, facePartBound + 4 * nFacePartBound );
//         PyDict_SetItem( splt_zone, CString2PyString( "facePartBound" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pFacePartBound ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "facePartBound recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( vertexTag ) {
//         KCore::Numpy::vector<int> pVertexTag( vertexTag, vertexTag + nVertex );
//         PyDict_SetItem( splt_zone, CString2PyString( "vertexTag" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pVertexTag ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "vertexTag recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( vertex ) {
//       KCore::Numpy::vector<double> pVertex( {(unsigned long)nVertex, 3} );
//         memcpy( pVertex.data(), vertex, 3 * nVertex * sizeof( double ) );
//         PyDict_SetItem( splt_zone, CString2PyString( "vertex" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pVertex ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "nVertex recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( faceGroupIdx ) {
//         KCore::Numpy::vector<int> pfaceGroupIdx( faceGroupIdx, faceGroupIdx + nFaceGroup + 1 );
//         PyDict_SetItem( splt_zone, CString2PyString( "faceGroupIdx" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pfaceGroupIdx ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "faceGroudIdx recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     if ( faceGroup ) {
//         KCore::Numpy::vector<int> pFaceGroup( faceGroup, faceGroup + sFaceGroup );
//         PyDict_SetItem( splt_zone, CString2PyString( "faceGroup" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pFaceGroup ) );
//     }
//     if ( faceGroupLNToGN ) {
//         KCore::Numpy::vector<PDM_g_num_t> pFaceGroupLNToGN( faceGroupLNToGN, faceGroupLNToGN + sFaceGroup );
//         PyDict_SetItem( splt_zone, CString2PyString( "faceGroupLN2GN" ),
//                         (PyObject *)KCore::Numpy::shared_with_python( pFaceGroupLNToGN ) );
//     }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "faceGroupLNToGN recupere..." << std::endl;
//     std::flush(std::cout);
// #   endif
//     // Construction d'un cellVertex et d'un cellVertexIdx plus cellTypeIdx et d'un cellType ?
//     if ( cellFaceIdx && cellFace && faceVertexIdx && faceVertex )
//     {
//         // On calcule le nombre d'éléments de chaque type. Pour l'instant : tetra, hexa, pyra, prism, triangle, quadrangle, line
//         int decal_type = -PDM_MESH_NODAL_BAR2;
//         std::array<int,PDM_MESH_NODAL_HEXA8> nb_elts_per_type;
//         nb_elts_per_type.fill(0);
//         for ( int i = 0; i < nCell; ++i )
//         {
//             int nbFacesInCell = cellFaceIdx[i+1]-cellFaceIdx[i];
//             int nbVertsOnFistFace = faceVertexIdx[cellFaceIdx[i]+1] - faceVertexIdx[cellFaceIdx[i]];
//             if ( nbFacesInCell == 2) // Seul cas : c'est une ligne
//                 nb_elts_per_type[PDM_MESH_NODAL_BAR2+decal_type] += 1;
//             else if (nbFacesInCell == 3) // Seul cas : le triangle
//                 nb_elts_per_type[PDM_MESH_NODAL_TRIA3+decal_type] += 1;
//             else if ( nbFacesInCell == 4 ) // Quadrangle ou Tetra ?
//             {
//                 if ( nbVertsOnFistFace == 2 ) // Quadrangle !
//                     nb_elts_per_type[PDM_MESH_NODAL_QUAD4+decal_type] += 1;
//                 else
//                 {
//                     assert(nbVertsOnFistFace == 3);
//                     nb_elts_per_type[PDM_MESH_NODAL_TETRA4+decal_type] += 1;
//                 }
//             }
//             else if ( nbFacesInCell == 5) // Pyramide ou Prism ?
//             {
//                 int nb_quad = 0;
//                 int nb_trig = 0;
//                 for (int iface = cellFaceIdx[i]; iface < cellFaceIdx[i+1]; ++iface )
//                 {
//                     if (faceVertexIdx[iface+1]-faceVertexIdx[iface] == 4) nb_quad += 1;
//                     else
//                     {
//                         assert(faceVertexIdx[iface+1]-faceVertexIdx[iface] == 3);
//                         nb_trig += 1;
//                     }
//                 }
//                 assert(nb_quad + nb_trig == 5);
//                 if ( nb_quad == 1)
//                     nb_elts_per_type[PDM_MESH_NODAL_PYRAMID5+decal_type] += 1;
//                 else
//                     nb_elts_per_type[PDM_MESH_NODAL_PRISM6+decal_type] += 1;
//             }
//             else // Hexa
//             {
//                 assert(nbFacesInCell == 6);
//                 nb_elts_per_type[PDM_MESH_NODAL_HEXA8+decal_type] += 1;
//             }
//         }
//         // Contruction des listes de cellVertex par type d'éléments
//         // Construction de cellVertexIdx :
//         constexpr const std::array<int, PDM_MESH_NODAL_HEXA8> nb_verts_per_elt =
//             { 2, 3, 4, 0, 4, 5, 6, 8 }; 
//         const std::array<std::string, PDM_MESH_NODAL_HEXA8> type_name = {"BAR", "TRI", "QUAD", "POLY2D", 
//                                                                          "TETRA", "PYRA", "PENTA", "HEXA"};
//         PyObject* py_cellvertex_per_elt_type = PyDict_New();
//         for (int i_elt_type = PDM_MESH_NODAL_BAR2; i_elt_type <= PDM_MESH_NODAL_HEXA8; ++i_elt_type )
//         {
//             //assert(i_elt_type+decal_type+1 < beg_elt_per_type.size());
//             if ( nb_elts_per_type[i_elt_type+decal_type] > 0 )
//             {
//                 KCore::Numpy::vector<int> pcellVertex( nb_elts_per_type[i_elt_type+decal_type] *  nb_verts_per_elt[i_elt_type+decal_type]);
//                 KCore::Numpy::vector<int> pcellTag   ( nb_elts_per_type[i_elt_type+decal_type] );
//                 PyObject* tuple_cell = PyTuple_New(2);
//                 PyTuple_SetItem(tuple_cell, 0, (PyObject *)KCore::Numpy::shared_with_python(pcellTag));
//                 PyTuple_SetItem(tuple_cell, 1, (PyObject *)KCore::Numpy::shared_with_python(pcellVertex));
//                 PyDict_SetItem(py_cellvertex_per_elt_type, PyUnicode_FromString(type_name[i_elt_type+decal_type].c_str()), 
//                                tuple_cell);
//             }
//         }
//         // On construit un tableau de renumérotation des éléments splittés pour les ranger par type :
//         std::array<int,PDM_MESH_NODAL_HEXA8> ind_elt_per_type;
//         ind_elt_per_type.fill(0);
//         //std::vector<E_Int> renum_per_type_2_splitted_renum(nCell, -1);
//         //std::vector<E_Int> splitted_renum_2_renum_per_type(nCell, -1);
//         for ( int i = 0; i < nCell; ++i )
//         {
//             int nbFacesInCell = cellFaceIdx[i+1]-cellFaceIdx[i];
//             int nbVertsOnFistFace = faceVertexIdx[cellFaceIdx[i]+1] - faceVertexIdx[cellFaceIdx[i]];
//             int ind_type_elt;
//             if ( nbFacesInCell == 2) // Seul cas : c'est une ligne
//                 ind_type_elt = PDM_MESH_NODAL_BAR2+decal_type;
//             else if (nbFacesInCell == 3) // Seul cas : le triangle
//                 ind_type_elt = PDM_MESH_NODAL_TRIA3+decal_type;
//             else if ( nbFacesInCell == 4 ) // Quadrangle ou Tetra ?
//             {
//                 if ( nbVertsOnFistFace == 2 ) // Quadrangle !
//                     ind_type_elt = PDM_MESH_NODAL_QUAD4+decal_type;
//                 else
//                     ind_type_elt = PDM_MESH_NODAL_TETRA4+decal_type;
//             }
//             else if ( nbFacesInCell == 5) // Pyramide ou Prism ?
//             {
//                 int nb_quad = 0;
//                 int nb_trig = 0;
//                 for (int iface = cellFaceIdx[i]; iface < cellFaceIdx[i+1]; ++iface )
//                 {
//                     if (faceVertexIdx[iface+1]-faceVertexIdx[iface] == 4) nb_quad += 1;
//                     else
//                     {
//                         assert(faceVertexIdx[iface+1]-faceVertexIdx[iface] == 3);
//                         nb_trig += 1;
//                     }
//                 }
//                 assert(nb_quad + nb_trig == 5);
//                 if ( nb_quad == 1)
//                     ind_type_elt = PDM_MESH_NODAL_PYRAMID5+decal_type;
//                 else
//                     ind_type_elt = PDM_MESH_NODAL_PRISM6+decal_type;
//             }
//             else // Hexa
//             {
//                 assert(nbFacesInCell == 6);
//                 ind_type_elt = PDM_MESH_NODAL_HEXA8+decal_type;
//             }
//             int ind_renum = ind_elt_per_type[ind_type_elt];
//             assert(ind_type_elt >= 0);
//             assert(ind_type_elt < ind_elt_per_type.size());
//             ind_elt_per_type[ind_type_elt] += 1;
//             PyObject* py_cell = PyDict_GetItem(py_cellvertex_per_elt_type, PyUnicode_FromString(type_name[ind_type_elt].c_str()));
//             PyObject* py_celltag    = PyTuple_GetItem(py_cell, 0);
//             KCore::Numpy::vector<int> tag((PyArrayObject*)py_celltag);
//             tag[ind_renum] = cellTag[i];
//             PyObject* py_cellvertex = PyTuple_GetItem(py_cell, 1);
//             KCore::Numpy::vector<int> indices((PyArrayObject*)py_cellvertex);
//             if (ind_type_elt == PDM_MESH_NODAL_BAR2 + decal_type)
//             {
//                 int pt_iFace = cellFaceIdx[i];
//                 int ind_cellvertex= 2*ind_renum;
//                 int iFace1   = cellFace[pt_iFace  ] - 1;
//                 int iFace2   = cellFace[pt_iFace+1] - 1;
//                 int pt_vert1 = faceVertexIdx[iFace1];
//                 int pt_vert2 = faceVertexIdx[iFace2];
//                 indices[ind_cellvertex  ] = faceVertex[pt_vert1];
//                 indices[ind_cellvertex+1] = faceVertex[pt_vert2];
//             }
//             if (ind_type_elt == PDM_MESH_NODAL_TRIA3 + decal_type)
//             {
//                 int ind_cellvertex= 3*ind_renum;
//                 int pt_iFace = cellFaceIdx[i];
//                 int iFace1   = cellFace[pt_iFace  ] - 1;
//                 int pt_vert1 = faceVertexIdx[iFace1];
//                 int iFace2   = cellFace[pt_iFace+1] - 1;
//                 int pt_vert2 = faceVertexIdx[iFace2];
//                 int iFace3   = cellFace[pt_iFace+2] - 1;
//                 int pt_vert3 = faceVertexIdx[iFace3];
//                 indices[ind_cellvertex  ] = faceVertex[pt_vert1  ];
//                 indices[ind_cellvertex+1] = faceVertex[pt_vert1+1];
//                 indices[ind_cellvertex+2] = ( faceVertex[pt_vert2  ] + faceVertex[pt_vert2+1] +
//                                               faceVertex[pt_vert3  ] + faceVertex[pt_vert3+1] -
//                                             ( faceVertex[pt_vert1  ] + faceVertex[pt_vert1+1]))/2;
//             }
//             if (ind_type_elt == PDM_MESH_NODAL_QUAD4 + decal_type)
//             {
//                 int ind_cellvertex = 4 * ind_renum;
//                 int pt_iFace = cellFaceIdx[i];
//                 int iFace1   = cellFace[pt_iFace  ] - 1;
//                 int pt_vert1 = faceVertexIdx[iFace1];
//                 int iFace2   = cellFace[pt_iFace+1] - 1;
//                 int pt_vert2 = faceVertexIdx[iFace2];
//                 int iFace3   = cellFace[pt_iFace+2] - 1;
//                 int pt_vert3 = faceVertexIdx[iFace3];
//                 int iFace4   = cellFace[pt_iFace+3] - 1;
//                 int pt_vert4 = faceVertexIdx[iFace4];
//                 indices[ind_cellvertex  ] = faceVertex[pt_vert1  ];
//                 indices[ind_cellvertex+1] = faceVertex[pt_vert1+1];
//                 indices[ind_cellvertex+2] = (faceVertex[pt_vert1+1] == faceVertex[pt_vert2])*faceVertex[pt_vert2+1] +
//                                             (faceVertex[pt_vert1+1] == faceVertex[pt_vert3])*faceVertex[pt_vert3+1] +
//                                             (faceVertex[pt_vert1+1] == faceVertex[pt_vert4])*faceVertex[pt_vert4+1];
//                 indices[ind_cellvertex+3] = (indices[ind_cellvertex+2] == faceVertex[pt_vert2])*faceVertex[pt_vert2+1] +
//                                             (indices[ind_cellvertex+2] == faceVertex[pt_vert3])*faceVertex[pt_vert3+1] +
//                                             (indices[ind_cellvertex+2] == faceVertex[pt_vert4])*faceVertex[pt_vert4+1];
//             }
//             if (ind_type_elt == PDM_MESH_NODAL_TETRA4 + decal_type)
//             {
//                 int ind_cellvertex = 4 * ind_renum;
//                 int pt_iFace = cellFaceIdx[i];
//                 int iFace1   = cellFace[pt_iFace  ] - 1;
//                 int pt_vert1 = faceVertexIdx[iFace1];
//                 int iFace2   = cellFace[pt_iFace+1] - 1;
//                 int pt_vert2 = faceVertexIdx[iFace2];
//                 int iFace3   = cellFace[pt_iFace+2] - 1;
//                 int pt_vert3 = faceVertexIdx[iFace3];
//                 int iFace4   = cellFace[pt_iFace+3] - 1;
//                 int pt_vert4 = faceVertexIdx[iFace4];
//                 indices[ind_cellvertex  ] = faceVertex[pt_vert1  ];
//                 indices[ind_cellvertex+1] = faceVertex[pt_vert1+1];
//                 indices[ind_cellvertex+2] = faceVertex[pt_vert1+2];
//                 indices[ind_cellvertex+3] = (faceVertex[pt_vert2] + faceVertex[pt_vert2+1] + faceVertex[pt_vert2+2] +
//                                              faceVertex[pt_vert3] + faceVertex[pt_vert3+1] + faceVertex[pt_vert3+2] +
//                                              faceVertex[pt_vert4] + faceVertex[pt_vert4+1] + faceVertex[pt_vert4+2] - 2 *
//                                              (faceVertex[pt_vert1] + faceVertex[pt_vert1+1] + faceVertex[pt_vert1+2]))/3;
//             }
//             if (ind_type_elt == PDM_MESH_NODAL_PYRAMID5 + decal_type)
//             {
//                 int ind_cellvertex = 5 * ind_renum;
//                 int pt_iFace = cellFaceIdx[i];

//                 std::array<int,5> ifaces;
//                 int i_face_trig = 1;
//                 for ( int j = 0; j < 5; ++j )
//                 {
//                     int iFace = cellFace[pt_iFace+j] - 1;
//                     int nb_nodes = faceVertexIdx[iFace+1]-faceVertexIdx[iFace];
//                     if (nb_nodes == 4) ifaces[0] = iFace;
//                     else {
//                         ifaces[i_face_trig] = iFace;
//                         i_face_trig += 1;
//                     }
//                 }
//                 int iFace1 = ifaces[0];
//                 int pt_vert1 = faceVertexIdx[iFace1];
//                 indices[ind_cellvertex  ] = faceVertex[pt_vert1  ];
//                 indices[ind_cellvertex+1] = faceVertex[pt_vert1+1];
//                 indices[ind_cellvertex+2] = faceVertex[pt_vert1+2];
//                 indices[ind_cellvertex+3] = faceVertex[pt_vert1+3];
//                 int acu = 0;
//                 for ( int j = 1; j < 5; ++j )
//                 {
//                     int pt_vert = faceVertexIdx[ifaces[j]];
//                     acu += faceVertex[pt_vert]+faceVertex[pt_vert+1]+faceVertex[pt_vert+2];
//                 }
//                 acu -= 2*(indices[ind_cellvertex  ]+indices[ind_cellvertex+1]+
//                           indices[ind_cellvertex+2]+indices[ind_cellvertex+3]);
//                 indices[ind_cellvertex+4] = acu/4;
//             }
//             if (ind_type_elt == PDM_MESH_NODAL_PRISM6 + decal_type)
//             {
//                 int ind_cellvertex = 6 * ind_renum;
//                 int pt_iFace = cellFaceIdx[i];

//                 std::array<int,2> ifaces;
//                 int i_face_trig = 0;
//                 for ( int j = 0; j < 6; ++j )
//                 {
//                     int iFace = cellFace[pt_iFace+j] - 1;
//                     int nb_nodes = faceVertexIdx[iFace+1]-faceVertexIdx[iFace];
//                     if (nb_nodes == 3) {
//                         ifaces[i_face_trig] = iFace;
//                         i_face_trig += 1;
//                     }
//                 }
//                 int iFace = ifaces[0];
//                 int pt_vert = faceVertexIdx[iFace];
//                 indices[ind_cellvertex  ] = faceVertex[pt_vert  ];
//                 indices[ind_cellvertex+1] = faceVertex[pt_vert+2];// Pour être conforme CGNS ?
//                 indices[ind_cellvertex+2] = faceVertex[pt_vert+1];
//                 iFace = ifaces[1];
//                 pt_vert = faceVertexIdx[iFace];
//                 indices[ind_cellvertex+3] = faceVertex[pt_vert  ];
//                 indices[ind_cellvertex+4] = faceVertex[pt_vert+1];
//                 indices[ind_cellvertex+5] = faceVertex[pt_vert+2];
//             }
//             if (ind_type_elt == PDM_MESH_NODAL_HEXA8 + decal_type)
//             {
//                 int ind_cellvertex = 8 * ind_renum;
//                 int pt_iFace      = cellFaceIdx[i];
//                 int iFace         = cellFace[pt_iFace+0] - 1;
//                 int pt_vert       = faceVertexIdx[iFace];
//                 indices[ind_cellvertex+0] = faceVertex[pt_vert  ];
//                 indices[ind_cellvertex+1] = faceVertex[pt_vert+3];
//                 indices[ind_cellvertex+2] = faceVertex[pt_vert+2];
//                 indices[ind_cellvertex+3] = faceVertex[pt_vert+1];
//                 // Recherche de la face ayant les autres sommets : 
//                 indices[ind_cellvertex+4] = -1;
//                 indices[ind_cellvertex+5] = -1;
//                 indices[ind_cellvertex+6] = -1;
//                 indices[ind_cellvertex+7] = -1;
//                 for ( int j = 1; j < 6; ++j )
//                 {
//                     iFace = cellFace[pt_iFace+j]-1;

//                     bool has_common = false;
//                     int pt_vert_base_face   = faceVertexIdx[iFace];
//                     int nb_verts_in_face    = faceVertexIdx[iFace+1] - faceVertexIdx[iFace];
//                     for ( int k = 0; k < nb_verts_in_face; ++k )
//                     {
//                         pt_vert       = pt_vert_base_face + k;
//                         int current_vertex = faceVertex[pt_vert];
//                         if ( ( current_vertex != indices[ind_cellvertex+0]) &&
//                              ( current_vertex != indices[ind_cellvertex+1]) &&
//                              ( current_vertex != indices[ind_cellvertex+2]) && 
//                              ( current_vertex != indices[ind_cellvertex+3]) )
//                         {
//                             int prec_vert = faceVertex[pt_vert_base_face + (k + nb_verts_in_face-1)%nb_verts_in_face];
//                             int next_vert = faceVertex[pt_vert_base_face + (k + 1)%nb_verts_in_face];
//                             if ( (next_vert == indices[ind_cellvertex+0]) ||
//                                  (prec_vert == indices[ind_cellvertex+0]) )
//                                 indices[ind_cellvertex+4] = current_vertex;
//                             if ( (next_vert == indices[ind_cellvertex+1]) ||
//                                  (prec_vert == indices[ind_cellvertex+1]) )
//                                 indices[ind_cellvertex+5] = current_vertex;
//                             if ( (next_vert == indices[ind_cellvertex+2]) ||
//                                  (prec_vert == indices[ind_cellvertex+2]) )
//                                 indices[ind_cellvertex+6] = current_vertex;
//                             if ( (next_vert == indices[ind_cellvertex+3]) ||
//                                  (prec_vert == indices[ind_cellvertex+3]) )
//                                 indices[ind_cellvertex+7] = current_vertex;
//                         }
//                     }
//                 }
//                 assert(indices[ind_cellvertex+4] != -1);
//                 assert(indices[ind_cellvertex+5] != -1);
//                 assert(indices[ind_cellvertex+6] != -1);
//                 assert(indices[ind_cellvertex+7] != -1);
//             }
//        }
// #   if defined( SPLITTER_TRACE )
//     std::cout << "Fin construction cellVertex..." << std::endl;
//     std::flush(std::cout);
// #   endif
//         PyDict_SetItem( splt_zone, CString2PyString( "cellVertex" ), py_cellvertex_per_elt_type);
//     }

//     return Py_BuildValue( "O", splt_zone );
// }
// // ----------------------------------------------------------------------------------------
// PyObject *splitElements( PyObject *self, PyObject *args )
// {
// #if defined( SPLITTER_TRACE )
//     std::cout << __PRETTY_FUNCTION__ << std::endl;
//     std::flush( std::cout );
// #endif
//     int nbparts = 1;
//     PyObject *set_of_zones;
//     // Attention, nbparts => nombre de partitions par procs !!!!
//     // Découpage se base sur le nombre de processus * nombre de partitions
//     if ( not PYPARSETUPLEI( args, "O|l", "O|i", &set_of_zones, &nbparts ) ) return NULL;
//     if ( not PyList_Check( set_of_zones ) ) {
//         PyErr_SetString( PyExc_RuntimeError, "Set of Zones must be stored in a python list" );
//         return NULL;
//     }
//     Py_ssize_t nb_zones = PyList_Size( set_of_zones );

//     domain_connectivity_t domain;
//     domain.reserve( nb_zones );
//     xcore::communicator comm = xcore::context::globalCommunicator();

//     // Pour chaque zone du domaine de calcul :
//     for ( Py_ssize_t i_zone = 0; i_zone < nb_zones; ++i_zone ) {
//         // Récupération d'une zone courante : la connectivité de la zone + nombre total de sommets dans la zone,
//         // nombre total d'éléments dans la zone
//         PyObject *pt_pyzone = PyList_GetItem( set_of_zones, i_zone );
//         if ( not PyTuple_Check( pt_pyzone ) ) {
//             PyErr_SetString( PyExc_TypeError,
//                              "Each element of the zone list must be a tuple : (zone, coordinates, number of "
//                              "vertices, number of elements)" );
//             return NULL;
//         }

//         zone_connectivity_t zconnect;

//         PyObject *py_zone;
//         PyArrayObject *pt_coords;
//         PyObject *py_coords;
//         E_Int nb_total_vertices, nb_total_elements;
//         if ( !PYPARSETUPLEI( pt_pyzone, "OOll", "OOii", &py_zone, &py_coords, &nb_total_vertices, &nb_total_elements ) ) {
//             PyErr_SetString( PyExc_TypeError,
//                              "The tuple must contains four values : the zone connectivity, "
//                              "the total number of vertices and total number of elements on the zone" );
//             return NULL;
//         }
//         if ( not PyList_Check( py_zone ) ) {
//             PyErr_SetString( PyExc_TypeError,
//                              "The zone object must be a list of connectivities per type of element" );
//             return NULL;
//         }
//         zconnect.nb_vertices = nb_total_vertices;
//         zconnect.nb_elements = nb_total_elements;
//         zconnect.num_zone    = i_zone;
//         zconnect.per_element_connectivities.reserve( PyList_Size( py_zone ) );

//         Py_ssize_t nb_elt_type = PyList_Size( py_zone );  // Nombre de types d'éléments dans la zone
//         // Pour chaque type d'éléments contenus dans la zone :
//         for ( Py_ssize_t itype = 0; itype < nb_elt_type; ++itype ) {
//             // Récupération d'un type d'élément : : ("Type d'élément (chaîne de caractère)", connectivité élt ->
//             // sommet,  )
//             PyObject *pt_pyelt_type = PyList_GetItem( py_zone, itype );
//             if ( not PyTuple_Check( pt_pyelt_type ) ) {
//                 PyErr_SetString(
//                     PyExc_TypeError,
//                     "Each element of the elements list must be a tuple (\"element type\", connectivity" );
//                 return NULL;
//             }
//             // Array contains the kind of element...
//             const char *elt_type;
//             PyArrayObject *pt_e2v;
//             if ( !PYPARSETUPLEI( pt_pyelt_type, "sO", "sO", &elt_type, &pt_e2v ) ) {
//                 PyErr_SetString( PyExc_TypeError,
//                                  "The tuple must contain a string given type of element, "
//                                  "and the connectivity array (numpy)" );
//                 return NULL;
//             }
//             int type_elt, nb_verts_per_elt;
//             std::tie( type_elt, nb_verts_per_elt ) = section_type( elt_type );
//             E_Int size, nfld;
//             E_Int *e2v;
//             E_Int coord_size, coord_nfld;
//             if ( not K_NUMPY::getFromNumpyArray( (PyObject *)pt_e2v, e2v, size, nfld ) ) {
//                 PyErr_SetString( PyExc_TypeError,
//                                  "The element to vertices array must be a "
//                                  "numpy array containing E_Int values" );
//                 return NULL;
//             }
//             assert( nfld == 1 );
//             PDM_g_num_t *pt_elt2vert;
//             if ( sizeof( E_Int ) == sizeof( PDM_g_num_t ) )
//                 pt_elt2vert = (PDM_g_num_t *)e2v;
//             else {
//                 //std::cout << "Copie d'indices avec representation entiere differente" << std::endl;
//                 pt_elt2vert = new PDM_g_num_t[ size ];
//                 std::copy( e2v, e2v + size, pt_elt2vert );
//             }
//             zconnect.per_element_connectivities.emplace_back( type_elt, nb_verts_per_elt, pt_elt2vert, int( size ), i_zone );
//         }                                               // Fin pour itype
//         zconnect.vertices_coordinates = py_coords;  //std::move(KCore::Numpy::vector<double>((PyArrayObject*)py_coords));
//         domain.push_back( zconnect );
//     }  // Fin pour i_zone
// #if defined( SPLITTER_TRACE )
//     std::cout << "Fini de browser les données d'entrée." << std::endl;
// #endif
//     // On va compter le nombre total de sommets, le nombre total d'éléments et renuméroter globalement les
//     // sommets et les éléments :
//     int nb_tot_vertices = 0, nb_tot_elements = 0;
//     int nb_loc_vertices = 0, nb_loc_elements = 0;  // dnVerts et dnCells
//     for ( auto &zone : domain ) {
// #if defined( SPLITTER_TRACE )
//         std::cout << "Iteration sur les zones..." << std::endl;
//     std::flush( std::cout );
// #endif
//         zone.first_ind_glob_vert = nb_tot_vertices;
//         zone.first_ind_glob_elts = nb_tot_elements;
//         nb_tot_elements += zone.nb_elements;
//         nb_tot_vertices += zone.nb_vertices;
//         for ( auto &elts_type : zone.per_element_connectivities )
//             for ( auto &index : elts_type.elt2vertsglob ) index += zone.first_ind_glob_vert;
//     }
// #if defined( SPLITTER_TRACE )
//     std::cout << "Création des coordonnées des sommets -> " << nb_tot_vertices << std::endl;
//     std::cout << "On a trouve de plus qu'il y a " << nb_tot_elements << " elements dans le probleme" << std::endl;
//     std::flush( std::cout );
// #endif
//     double *dVtxCoord = new double[ 3 * nb_tot_vertices ];  // NULL;
//     int* dVtxTag   = new int[nb_tot_vertices];
//     int ind_vert = 0;
//     int nb_loc_stored_verts = 0;
//     int izone = 0;
//     int ind_elt_glob = 0;
//     for ( auto &zone : domain ) {
// #if defined( SPLITTER_TRACE )
//         std::cout << "Récupération des coordonnées => " << zone.vertices_coordinates << std::endl;
//     std::flush( std::cout );
// #endif
//         double *coords;
//         E_Int coord_size, coord_nfld;
//         PyArrayObject *py_coords;
//         if ( not PyTuple_Check( zone.vertices_coordinates ) ) {
//             py_coords = (PyArrayObject *)zone.vertices_coordinates;
//             if ( not K_NUMPY::getFromNumpyArray( (PyObject *)py_coords, coords, coord_size, coord_nfld ) ) {
//                 PyErr_SetString( PyExc_TypeError,
//                                  R"RAW(The coordinates of vertices array must be a numpy array containing double values)RAW" );
//                 return NULL;
//             }
//             nb_loc_stored_verts += coord_size;
//             for ( E_Int i = 0; i < coord_size; ++i ) {
//                 dVtxCoord[ 3 * ind_vert + 3 * i + 0 ] = coords[ i + 0 * coord_size ];
//                 dVtxCoord[ 3 * ind_vert + 3 * i + 1 ] = coords[ i + 1 * coord_size ];
//                 dVtxCoord[ 3 * ind_vert + 3 * i + 2 ] = coords[ i + 2 * coord_size ];
//                 dVtxTag[ind_vert+i] = izone;
//             }
//         } else {
//             double *x_coords, *y_coords, *z_coords;
//             E_Int x_coord_size, y_coord_size, z_coord_size;
//             E_Int x_coord_nfld, y_coord_nfld, z_coord_nfld;
//             PyObject *py_x_coords = PyTuple_GetItem( zone.vertices_coordinates, 0 );
//             PyObject *py_y_coords = PyTuple_GetItem( zone.vertices_coordinates, 1 );
//             PyObject *py_z_coords = PyTuple_GetItem( zone.vertices_coordinates, 2 );
//             if ( not K_NUMPY::getFromNumpyArray( py_x_coords, x_coords, x_coord_size, x_coord_nfld ) ) {
//                 PyErr_SetString( PyExc_TypeError,
//                                  R"RAW(The X coordinates of vertices array must be a numpy array containing double values)RAW" );
//                 return NULL;
//             }
//             if ( not K_NUMPY::getFromNumpyArray( py_y_coords, y_coords, y_coord_size, y_coord_nfld ) ) {
//                 PyErr_SetString( PyExc_TypeError,
//                                  R"RAW(The Y coordinates of vertices array must be a numpy array containing double values)RAW" );
//                 return NULL;
//             }
//             if ( not K_NUMPY::getFromNumpyArray( py_z_coords, z_coords, z_coord_size, z_coord_nfld ) ) {
//                 PyErr_SetString( PyExc_TypeError,
//                                  R"RAW(The Z coordinates of vertices array must be a numpy array containing double values)RAW" );
//                 return NULL;
//             }
//             if ( ( x_coord_nfld != 1 ) or ( y_coord_nfld != 1 ) or ( z_coord_nfld != 1 ) ) {
//                 PyErr_SetString( PyExc_TypeError,
//                                  R"RAW(The coordinates of vertices arrays must be one dimensional numpy array)RAW" );
//                 return NULL;
//             }
//             if ( ( x_coord_size != y_coord_size ) or ( x_coord_size != z_coord_size ) ) {
//                 PyErr_SetString( PyExc_TypeError,
//                                  R"RAW(The coordinates of vertices arrays must have same number of values)RAW" );
//                 return NULL;
//             }
//             nb_loc_stored_verts += x_coord_size;
//             for ( E_Int i = 0; i < x_coord_size; ++i ) {
//                 dVtxCoord[ 3 * ind_vert + 3 * i + 0 ] = x_coords[ i ];
//                 dVtxCoord[ 3 * ind_vert + 3 * i + 1 ] = y_coords[ i ];
//                 dVtxCoord[ 3 * ind_vert + 3 * i + 2 ] = z_coords[ i ];
//                 dVtxTag[ind_vert+i] = izone;
//             }
//         }
//         ind_vert += zone.nb_vertices;
//         for ( auto &elts_type : zone.per_element_connectivities ) {
//             unordered_set<PDM_g_num_t> set_of_vertices;
//             set_of_vertices.insert( elts_type.elt2vertsglob.begin(), elts_type.elt2vertsglob.end() );
//             //===========================================================================
//             // int                        type_element, nb_verts_per_elt;
//             // K_MEMORY::vector_view<PDM_g_num_t> elt2vertsglob;
//             //===========================================================================
//             nb_loc_elements += elts_type.elt2vertsglob.size() / elts_type.nb_verts_per_elt;
//             nb_loc_vertices += set_of_vertices.size();
//             if ( elt2global_vertices.find( elts_type.type_element ) == elt2global_vertices.end() )
//             {
//                 elt2global_vertices[ elts_type.type_element ] =
//                     std::vector<PDM_g_num_t>( elts_type.elt2vertsglob.begin(), elts_type.elt2vertsglob.end() );
//                 elttag_global[ elts_type.type_element ] = std::vector<int>( nb_loc_elements, elts_type.tag );
//             }
//             else
//             {
//                 elt2global_vertices[ elts_type.type_element ].insert(
//                     elt2global_vertices[ elts_type.type_element ].end(), elts_type.elt2vertsglob.begin(),
//                     elts_type.elt2vertsglob.end() );
//                 elttag_global[ elts_type.type_element ].insert(elttag_global[ elts_type.type_element ].end(), 
//                                                                nb_loc_elements,
//                                                                elts_type.tag);
//             }
//         }
//         izone += 1;
//     }  // for ( auto& zone : domain )
// #if defined( SPLITTER_TRACE )
//     std::cout << "nb_tot_vertices : " << nb_tot_vertices << ", nb_tot_elements : " << nb_tot_elements << std::endl;
//     std::flush( std::cerr );
//     std::flush( std::cout );
// #endif
//     // Construction du elttag global :
//     std::vector<int> dCellTag;
//     dCellTag.reserve(nb_tot_elements);
//     for (const auto& tag : elttag_global )
//     {
//         dCellTag.insert(dCellTag.end(), tag.second.begin(), tag.second.end() );
//     }    
//     // Construction de la structure DMesh de Paradigma pour le découpage :
//     xcore::Ext_Communicator &ext_comm = comm.get_implementation();
//     PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm( &ext_comm );
//     int hdl_dmesh = PDM_DMesh_nodal_create( pdm_comm, nb_tot_vertices, nb_tot_elements );
// #   if defined( SPLITTER_TRACE )
//     std::cout << "Creation DMesh nodal fini" << std::endl;
//     std::flush(std::cout);
// #   endif
//     for ( auto &elt_type : elt2global_vertices ) {
//         int type_elt = elt_type.first;
//         int nb_vert_per_elt = nb_verts_for_element( type_elt );
//         int id_section = PDM_DMesh_nodal_section_add( hdl_dmesh, (PDM_Mesh_nodal_elt_t)type_elt );
//         PDM_DMesh_nodal_section_std_set( hdl_dmesh, id_section, elt_type.second.size() / nb_vert_per_elt,
//                                          elt_type.second.data() );
//     }
//     PDM_DMesh_nodal_coord_set( hdl_dmesh, nb_loc_stored_verts, dVtxCoord );
//     comm.barrier();
// #   if defined( SPLITTER_TRACE )
//     std::cout << "Coord DMesh nodal setted" << std::endl;
//     std::flush(std::cout);
// #   endif
//     PDM_DMesh_nodal_cell_face_compute( hdl_dmesh );
// #   if defined( SPLITTER_TRACE )
//     std::cout << "Cell face computed for DMesh nodal" << std::endl;
//     std::flush(std::cout);
// #   endif

//     int hdl_ppart;
//     const int nb_properties_for_cells = 3;
//     const int nCellPerCache = 1024, is_asynchrone = 1, is_vectorisation = 1;
//     int cell_properties[ nb_properties_for_cells ] = {nCellPerCache, is_asynchrone, is_vectorisation};

//     int nFaceGroupe = 0;  // Nombre de types de conditition limite dont les raccords...
//     int *dFaceGroupIdx = NULL, *dFaceGroup = NULL;

//     int have_dCellPart = 0;
//     // nb local d'éléments
//     std::vector<int> dCellPart( nb_tot_elements );

//     int nb_face_tot = PDM_DMesh_nodal_total_n_face_get( hdl_dmesh );
//     PDM_g_num_t *dFaceCell;
//     int nb_faces = PDM_DMesh_nodal_face_cell_get( hdl_dmesh, &dFaceCell );
//     int nb_vertices = PDM_DMesh_nodal_n_vtx_get( hdl_dmesh );
//     const double *d_coords = PDM_DMesh_nodal_vtx_get( hdl_dmesh );
//     int *dFaceVtxIdx;
//     PDM_g_num_t *dFaceVtx = NULL;
//     int nb_face_loc = PDM_DMesh_nodal_face_vtx_get( hdl_dmesh, &dFaceVtxIdx, &dFaceVtx );
//     PDM_part_create( &hdl_ppart, pdm_comm, PDM_PART_SPLIT_PTSCOTCH, "PDM_PART_RENUM_CELL_NONE",
//                      "PDM_PART_RENUM_FACE_NONE", nb_properties_for_cells, cell_properties, 0,
//                      NULL,  // Two last for face properties ( eq. to cell properties but for faces )
//                      nbparts, nb_loc_elements, nb_face_loc, nb_vertices, nFaceGroupe, NULL, NULL, dCellTag.data(), NULL,
//                      have_dCellPart, dCellPart.data(), dFaceCell, dFaceVtxIdx, dFaceVtx,
//                      NULL, d_coords, dVtxTag, dFaceGroupIdx, dFaceGroup );
// #   if defined( SPLITTER_TRACE )
//     std::cout << "DMesh nodal splitted !" << std::endl;
//     std::flush(std::cout);
// #   endif

//     PyObject *py_part_list = PyList_New( nbparts );
//     for ( int i = 0; i < nbparts; ++i )
//             PyList_SetItem( py_part_list, i, (PyObject *)build_splitted_zone( hdl_ppart, i ) );
//     delete[] dVtxCoord;
//     delete [] dVtxTag;
//     PDM_DMesh_nodal_free( hdl_dmesh );
//     return py_part_list;
// }  // Fin fonction splitElements
