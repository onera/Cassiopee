/*    
    Copyright 2013-2023 Onera.

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


# include <string>
# include <sstream> 
# include "intersector.h"
# include "Nuga/include/ngon_t.hxx"
# include "Nuga/include/Triangulator.h"
# include "Nuga/include/Splitter.h"
# include <memory>
//#include <iostream>

#include "Nuga/include/DynArray.h"
#include "Nuga/include/merge.h"
#include "Nuga/include/mesh_t.hxx"
#include "Nuga/include/BbTree.h"
#include "Nuga/include/medit.hxx"
#include "Nuga/include/sync_faces.hxx"
#include "Nuga/include/supermesh.hxx"

using namespace std;
using namespace NUGA;

using zmesh_t = NUGA::pg_smesh_t;
using ngon_type = ngon_t<K_FLD::IntArray>;

// update the point list according to a split expressed by oids
PyObject* K_INTERSECTOR::updatePointLists(PyObject* self, PyObject* args)
{
  PyObject *py_oids, *py_ptLists;

  if (!PyArg_ParseTuple(args, "OO", &py_oids, &py_ptLists)) return NULL;

  E_Int nb_bcs = PyList_Size(py_ptLists);

  E_Int sz{0}, r;
  E_Int* oids;
  E_Int res = K_NUMPY::getFromNumpyArray(py_oids, oids, sz, r, true/*shared*/);
  if (res != 1) return NULL;
  
  // WARNING : oids might have IDX_NONE (created entities, e.g. internal faces with adaptCells) and is 0-based 
  E_Int nb_pgs = 0;
  for (size_t i=0; i < sz; ++i)
  {
    if (oids[i] == IDX_NONE) continue;
    nb_pgs = std::max(nb_pgs, oids[i]+1);
  }

  ngon_unit split_graph;
  K_CONNECT::IdTool::reverse_indirection(nb_pgs,  oids, sz, split_graph);

  PyObject *l(PyList_New(0)), *tpl;
  std::vector<E_Int> new_ptl;

  for (E_Int i=0; i < nb_bcs; ++i)
  {
    new_ptl.clear();
    PyObject* o = PyList_GetItem(py_ptLists, i);

    E_Int ptl_sz;
    FldArrayI out;
    K_ARRAY::getFromList(o, out);
    ptl_sz = out.getSize();

    if (split_graph.size() != 0)
    {
      
    for (E_Int j=0; j < ptl_sz; ++j)
    {
      E_Int oid = out(j,1)-1;

      // need to check in case last face has disappeared
      // and size of split_graph is smaller than number of faces
      if ( oid < split_graph.size() )
      {
      
	E_Int nb_bits = split_graph.stride(oid);
	const E_Int* pbits = split_graph.get_facets_ptr(oid);
	// printf("pbits : %d \n", pbits );

	if (nb_bits == 1 && pbits[0] == E_IDX_NONE)  // gone
	  continue;
	else
	  for (E_Int u=0; u<nb_bits; ++u )
	  {
	    new_ptl.push_back(pbits[u]+1);
	  }
      }
    }
    }

    if (new_ptl.empty())
      tpl = Py_None;
    else
      tpl = K_NUMPY::buildNumpyArray(&new_ptl[0], new_ptl.size(), 1, 0);
    
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

   return l;
  
}

//=============================================================================
/* Triangulates exterior faces (any Polygon). */
//=============================================================================
PyObject* K_INTERSECTOR::triangulateExteriorFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Int int_or_ext(2); //0 : internals only, 1: external only, 2: both
  E_Int improve_qual(0);

  if (!PYPARSETUPLEI(args, "Oll", "Oii", &arr, &int_or_ext, &improve_qual)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  ngon_type ngi(cnt), ngo;
  NUGA::transfo_t qual_param;
  qual_param.improve_qual=improve_qual ? true : false;

  Splitter::triangulate_external_pgs<DELAUNAY::Triangulator>(crd, ngi, int_or_ext, qual_param, ngo);
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;
  
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* Triangulates BC. */
//=============================================================================
PyObject* K_INTERSECTOR::triangulateSpecifiedFaces(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_pgs;
  E_Int /*int_or_ext(2), */improve_qual(1); //0 : internals only, 1: external only, 2: both

  if (!PYPARSETUPLEI(args, "OOl", "OOi", &arr, &py_pgs, &improve_qual)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  E_Int res=0;
  E_Int* pgsList=NULL;
  E_Int size, nfld;
  if (py_pgs != Py_None)
    res = K_NUMPY::getFromNumpyArray(py_pgs, pgsList, size, nfld, true/*shared*/);

  if (res != 1) return NULL;
  
  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  ngon_type ngi(cnt), ngo;

  // enable history
  ngi.PGs._ancEs.resize((E_Int)2, (E_Int)ngi.PGs.size(), (E_Int)0);
  for (E_Int i=0; i < ngi.PGs.size(); ++i) ngi.PGs._ancEs(0,i)=i;

  NUGA::transfo_t qual_param;
  qual_param.improve_qual=improve_qual ? true : false;

  Splitter::triangulate_specified_pgs<DELAUNAY::Triangulator>(crd, ngi, pgsList, size, qual_param, ngo);
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);

  PyObject *l(PyList_New(0)), *tpl;
  
  // pushing out the mesh
  tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // pushing out history  
  E_Int sz = ngo.PGs.size();

  std::vector<E_Int> oids(sz);
  for (E_Int i = 0; i < sz; ++i) oids[i] = ngo.PGs._ancEs(0,i);

  tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);
  
  
  delete f; delete cn;
  return l;
}

PyObject* K_INTERSECTOR::triangulateNFaces(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_pgs;
  E_Int improve_qual(1), min_nvertices(4);

  if (!PYPARSETUPLEI(args, "OllO", "OiiO", &arr, &improve_qual, &min_nvertices, &py_pgs)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  E_Int* pgsList=NULL;
  E_Int size, nfld, res(1);
  if (py_pgs != Py_None)
    res = K_NUMPY::getFromNumpyArray(py_pgs, pgsList, size, nfld, true/*shared*/);

  if (res != 1) return NULL;
  
  ngon_type ngi(cnt), ngo;

  // enable history
  E_Int nb_pgs=(E_Int)ngi.PGs.size();
  ngi.PGs._ancEs.resize((E_Int)2, nb_pgs, (E_Int)0);
  for (E_Int i=0; i < nb_pgs; ++i) ngi.PGs._ancEs(0,i)=i;

  // discard pgs belonging to given BC
  std::vector<bool> addit(nb_pgs, true);
  if (pgsList != NULL)
    for (E_Int k = 0; k < size; ++k)
    {
     E_Int pg=pgsList[k];
      addit[pg] = false;
    }

  std::vector<E_Int> pgs;
  for (E_Int i=0; i < ngi.PGs.size(); ++i) if (ngi.PGs.stride(i) >= min_nvertices && addit[i]) pgs.push_back(i);

  NUGA::transfo_t qual_param;
  qual_param.improve_qual=improve_qual ? true : false;

  //std::cout << "nb of pgs to proceed over total : " << pgs.size() << "/" << ngi.PGs.size() << std::endl;

  Splitter::triangulate_specified_pgs<DELAUNAY::Triangulator>(crd, ngi, &pgs[0], pgs.size(), qual_param, ngo);
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);

  PyObject *l(PyList_New(0)), *tpl;

  //std::cout << "nb pgs before/after : " << ngi.PGs.size() << "/" << ngo.PGs.size() << std::endl;
  
  // pushing out the mesh
  tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // pushing out history  
  E_Int sz = ngo.PGs.size();

  std::vector<E_Int> oids(sz);
  for (E_Int i = 0; i < sz; ++i) oids[i] = ngo.PGs._ancEs(0,i);

  tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);
  
  
  delete f; delete cn;
  return l;
}

//=============================================================================
/* Convexify any concave polygon found in the input polyhedral mesh. */
//=============================================================================
PyObject* K_INTERSECTOR::convexifyFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float convexity_tol(1.e-8);

  if (!PYPARSETUPLEF(args, "Od", "Of", &arr, &convexity_tol)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  ngon_type ngi(cnt), ngo;
  Splitter::split_pgs(Splitter::convexify_pgs<DELAUNAY::Triangulator>, crd, ngi, convexity_tol, ngo);
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;
  
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX. */
//=============================================================================
PyObject* K_INTERSECTOR::superMesh(PyObject* self, PyObject* args)
{
  PyObject *arr1{nullptr}, *arr2{nullptr};
  E_Float rtol(1.e-6);

  if (!PYPARSETUPLEF(args, "OOd", "OOf", &arr1, &arr2, &rtol)) return NULL;

  K_FLD::FloatArray* f1(0);
  K_FLD::IntArray* cn1(0);
  char* varString1, *eltType1;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f1, cn1, varString1, eltType1);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd1 = *f1;
  K_FLD::IntArray & cnt1 = *cn1;

  K_FLD::FloatArray* f2(0);
  K_FLD::IntArray* cn2(0);
  char* varString2, *eltType2;
  // Check array # 2
  err = check_is_NGON(arr2, f2, cn2, varString2, eltType2);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd2 = *f2;
  K_FLD::IntArray & cnt2 = *cn2;

  /*std::cout << "crd1 : " << crd1.cols() << "/" << crd1.rows() << std::endl;
  std::cout << "cnt1 : " << cnt1.cols() << "/" << cnt1.rows() << std::endl;
  std::cout << "crd2 : " << crd2.cols() << "/" << crd2.rows() << std::endl;
  std::cout << "cnt2 : " << cnt2.cols() << "/" << cnt2.rows() << std::endl;*/


  // construction des structures de type "mesh" a pertir des ngon
  zmesh_t m0(crd1, cnt1);
  zmesh_t m1(crd2, cnt2);

  //std::cout << "m0/m1 cells : " << m0.ncells() << "/" << m1.ncells() << std::endl;

  // decoupage //////////////////////////////////////////

  //NUGA::chrono c;
  //c.start();

  zmesh_t xmesh;                 // maillage des morceaux polygonaux
  std::vector<E_Int> anc0, anc1; // anc0/anc1 : indice ancetre d'un polygone de xmesh dans m0/m1

  double ARTOL = ::sqrt(m0.Lref2())* rtol;
  /*std::cout << "Lref : " << ::sqrt(m0.Lref2()) << std::endl;
  std::cout << "RTOL : " << rtol << std::endl;
  std::cout << "ARTOL : " << ARTOL << std::endl;*/

  double ARTOL1 = ::sqrt(m1.Lref2())* rtol;
  //std::cout << "Lre : " << ::sqrt(m1.Lref2()) << std::endl;
  //std::cout << "ARTOL1 : " << ARTOL1 << std::endl;

  ARTOL = std::min(ARTOL, ARTOL1);

  NUGA::xmatch<zmesh_t>(m0, m1, ARTOL, anc0, anc1, xmesh);

  //std::cout << "xmesh cells : " << xmesh.ncells() << std::endl;

  PyObject *l(PyList_New(0));

  if (xmesh.ncells())
  {
    ngon_type ngo(xmesh.cnt, true);
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);

    PyObject* tpl = K_ARRAY::buildArray(xmesh.crd, varString1, cnto, -1, eltType1, false);

    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // pushing out arr1 history  
    tpl = K_NUMPY::buildNumpyArray(&anc0[0], anc0.size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // pushing out arr2 history  
    tpl = K_NUMPY::buildNumpyArray(&anc1[0], anc1.size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  delete f1; delete cn1;
  delete f2; delete cn2;
  return l;
}

//=============================================================================
/* XXX. */
//=============================================================================
PyObject* K_INTERSECTOR::replaceFaces(PyObject* self, PyObject* args)
{
  PyObject *arr{nullptr}, *arr_soup{nullptr}, *py_vfoid{nullptr};

  if (!PyArg_ParseTuple(args, "OOO", &arr, &arr_soup, &py_vfoid)) return NULL;

  K_FLD::FloatArray* f1(0);
  K_FLD::IntArray* cn1(0);
  char* varString1, *eltType1;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f1, cn1, varString1, eltType1);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd1 = *f1;
  K_FLD::IntArray & cnt1 = *cn1;

  K_FLD::FloatArray* f2(0);
  K_FLD::IntArray* cn2(0);
  char* varString2, *eltType2;
  // Check array # 2
  err = check_is_NGON(arr_soup, f2, cn2, varString2, eltType2);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd2 = *f2;
  K_FLD::IntArray & cnt2 = *cn2;

  E_Int sz{0}, r;
  E_Int* vfoid;
  E_Int res = K_NUMPY::getFromNumpyArray(py_vfoid, vfoid, sz, r, true/*shared*/);
  if (res != 1) return NULL;

  // construction des structures de type "mesh" a pertir des ngon
  NUGA::ph_mesh_t m0(crd1, cnt1);
  NUGA::pg_smesh_t xmesh(crd2, cnt2);

  //std::cout << "m0 faces : " << m0.cnt.PGs.size() << std::endl;
  // std::cout << "xmesh faces : " << xmesh.ncells() << std::endl;
  // std::cout << "vfoid sz : " << sz << std::endl;

  assert (sz == xmesh.ncells());

  // clean xmesh
  // K_FLD::ArrayAccessor<K_FLD::FloatArray> ca(xmesh.crd);
  // Vector_t<E_Int> nids;
  // E_Int nb_merges = ::merge(ca, 1.e-8, nids);
  // if (nb_merges)
  //   xmesh.cnt.change_indices(nids);
  // ngon_type::compact_to_used_nodes(xmesh.cnt, xmesh.crd);

  K_CONNECT::IdTool::init_inc(m0.cnt.PGs._type, m0.cnt.PGs.size()); // first, all faces assumed intact
  xmesh.cnt._type.resize(sz, IDX_NONE);
  for (size_t i=0; i < sz; ++i)
  {
    xmesh.cnt._type[i] = vfoid[i]; // assign history to faces to be appended referring to faces in m0
    assert (vfoid[i] < m0.cnt.PGs.size());
    //if (vfoid[i] < 0 || vfoid[i] >= m0.cnt.PGs.size()) std::cout << "error vfoid" << std::endl;
    m0.cnt.PGs._type[vfoid[i]] = IDX_NONE; // disable such face that will be replaced (to avoid to append both a face and its decomposition)
  }


  E_Int N0 = m0.crd.cols();
  m0.crd.pushBack(xmesh.crd);
  xmesh.cnt.shift(N0);
  m0.cnt.PGs.append(xmesh.cnt); //appending faces

  std::map<E_Int, std::vector<E_Int>> split_map;
  K_CONNECT::IdTool::reverse_indirection(&m0.cnt.PGs._type[0], m0.cnt.PGs._type.size(), split_map);

  double RTOL = 0.1;
  std::set<E_Int> modifiedPHs;
  replace_faces(m0, split_map, -RTOL, modifiedPHs);

  K_FLD::IntArray cnto;
  m0.cnt.export_to_array(cnto);

  PyObject* tpl = K_ARRAY::buildArray(m0.crd, varString1, cnto, -1, eltType1, false);

  PyObject *l(PyList_New(0));

  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // pushing out arr1 history  
  tpl = K_NUMPY::buildNumpyArray(&m0.cnt.PGs._type[0], m0.cnt.PGs._type.size(), 1, 0);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  //std::cout << "nb pgs : " << m0.cnt.PGs.size() << std::endl;
  //std::cout << "anc size : " << m0.cnt.PGs._type.size() << std::endl;

  // pushing out arr2 history  
  //pl = K_NUMPY::buildNumpyArray(&anc1[0], anc1.size(), 1, 0);
  //PyList_Append(l, tpl);
  //Py_DECREF(tpl);

  delete f1; delete cn1;
  delete f2; delete cn2;
  return l;
}

//============================================================================================================
/* Split (convexify, starify) some targeted polygons on targeted cells 
 * (typically bad polyhedra -concaves, non-centroid-star-shaped-)
 * to prepare the split of those bad cells.*/
//============================================================================================================
PyObject* K_INTERSECTOR::prepareCellsSplit(PyObject* self, PyObject* args)
{
  PyObject *arr(0);
  E_Float PH_conc_threshold(1./3.);
  E_Float PH_cvx_threshold(0.05);
  E_Float PG_cvx_threshold(1.e-8);
  E_Int PH_set(0); // 0 for concave cells or 1 for non-centroid-star_shaped cells
  E_Int split_policy (0); // 0 : convexify concave pgs on PH set. 1 : starify concave pgs on PH set. 2 : starify any pgs at concave-chains ends.

  if (!PYPARSETUPLE(args, 
                    "Ollddd", "Oiiddd", "Ollfff", "Oiifff",
                    &arr, &PH_set, &split_policy, &PH_conc_threshold, &PH_cvx_threshold, &PG_cvx_threshold)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  ngon_type ngi(cnt), ngo;
  
  transfo_t::ePHset phset = transfo_t::ePHset(PH_set);
  transfo_t::eSplitPolicy policy = transfo_t::eSplitPolicy(split_policy);
  
  err = Splitter::prepareCellsSplit<DELAUNAY::Triangulator>(crd, ngi, phset, policy, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold, ngo);

  PyObject* tpl = NULL;
  if (err)
  {
    if (err == 1)
      PyErr_SetString(PyExc_TypeError, "prepareCellsSplit : could not orient properly.");
    else if (err == 2)
      PyErr_SetString(PyExc_TypeError, "prepareCellsSplit : could not dectect correctly pathologies.");
    else if (err == 3)
      PyErr_SetString(PyExc_TypeError, "prepareCellsSplit : could not split polygons correctly.");
  }
  else
  {
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
  
    tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);
  }

  delete f; delete cn;
  return tpl;
}


///
PyObject* K_INTERSECTOR::syncMacthPeriodicFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Float center[3], axis[3], trans[3], artol(-0.01);

  if (!PYPARSETUPLEF(args, "O(ddd)(ddd)(ddd)d", "O(fff)(fff)(fff)f", &arr, 
                                                                     &center[0], &center[1], &center[2],
                                                                     &axis[0], &axis[1], &axis[2],
                                                                     &trans[0], &trans[1], &trans[2], &artol)) return nullptr;
  
  bool is_rot = false;
  double rot_angle = 0.;
  double * paxis{nullptr}, *ptrans{nullptr};

  if (axis[0] != 0. || axis[1] != 0. || axis[2] != 0.)
  {
    is_rot = true;
    rot_angle = ::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    paxis = &axis[0];
  }
  bool is_trans = false;
  if (trans[0] != 0. || trans[1] != 0. || trans[2] != 0.)
  {
    is_trans = true;
    ptrans = &trans[0];
  }

  if (is_trans && is_rot)
    return nullptr;
  if (!is_trans && !is_rot)
    return nullptr;

  K_FLD::FloatArray *f(0);
  K_FLD::IntArray *cn(0);
  char *varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return nullptr;


  std::unique_ptr<K_FLD::FloatArray> pf(f);   //for memory cleaning
  std::unique_ptr<K_FLD::IntArray> pcn(cn); //for memory cleaning

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  NUGA::ph_mesh_t m(crd, cnt);

  // for history
  m.cnt.PGs._ancEs.resize(2, m.cnt.PGs.size(), IDX_NONE);
  m.cnt.PHs._ancEs.resize(2, m.cnt.PHs.size(), IDX_NONE);
  K_CONNECT::IdTool::init_inc(m.cnt.PGs._ancEs, 0, m.cnt.PGs.size());
  K_CONNECT::IdTool::init_inc(m.cnt.PHs._ancEs, 0, m.cnt.PHs.size());
  
  bool carry_on = true;
  int itermax=10;
  int iter=0;
  do
  {
    std::map<E_Int, std::vector<E_Int>> glob_face_to_bits;
    detect_async_modified_faces(m, center, paxis, rot_angle, ptrans, artol, glob_face_to_bits);

    carry_on = !glob_face_to_bits.empty();
    if (!carry_on) break;

    duplicate_and_move_period_faces(m, center, paxis, rot_angle, ptrans, glob_face_to_bits);

    sync_faces(m, glob_face_to_bits, artol);

  } while (carry_on and iter++ < itermax);

  K_FLD::IntArray cnto;
  m.cnt.export_to_array(cnto);

  PyObject *l(PyList_New(0)), *tpl;
  
  // pushing out the mesh
  tpl = K_ARRAY::buildArray(m.crd, varString, cnto, -1, eltType, false);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // pushing out PG history  
  E_Int sz = m.cnt.PGs.size();

  std::vector<E_Int> oids(sz);
  for (E_Int i = 0; i < sz; ++i) oids[i] = m.cnt.PGs._ancEs(0,i);

  tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // pushing out PH history  
  sz = m.cnt.PHs.size();
  oids.clear(); oids.resize(sz);

  for (E_Int i = 0; i < sz; ++i) oids[i] = m.cnt.PHs._ancEs(0,i);

  tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
 
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  return l;
}


//=======================  Intersector/PolyMeshTools/split_faces.cpp ====================
