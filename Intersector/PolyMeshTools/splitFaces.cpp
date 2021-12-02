/*    
    Copyright 2013-2022 Onera.

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

using namespace std;
using namespace NUGA;


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
            
      E_Int nb_bits = split_graph.stride(oid);
      const E_Int* pbits = split_graph.get_facets_ptr(oid);

      if (nb_bits == 1 && pbits[0] == E_IDX_NONE)  // gone
        continue;
      else
        for (E_Int u=0; u<nb_bits; ++u )
          new_ptl.push_back(pbits[u]+1);
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
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
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
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
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

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
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
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ngi(cnt), ngo;
  Splitter::split_pgs(Splitter::convexify_pgs<DELAUNAY::Triangulator>, crd, ngi, convexity_tol, ngo);
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;
  
  
  delete f; delete cn;
  return tpl;
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
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
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


// match perio stuff ////
using ngon_type = ngon_t<K_FLD::IntArray>;
using crd_t = K_FLD::FloatArray;
using acrd_t = K_FLD::ArrayAccessor<crd_t>;

///
void detect_async_modified_faces(NUGA::ph_mesh_t& vmesh, const double* center, const double * axis, double angle, const double* translation,
  double ARTOL, std::map<int, std::vector<int>>& glob_face_to_bits)
{
  DELAUNAY::Triangulator dt;

  NUGA::pg_smesh_t m;
  std::vector<int> oids, ancestors;

  vmesh.get_boundary(m, oids, ancestors);

  std::map<int, std::vector<int>> loc_face_to_bits;

  glob_face_to_bits.clear();

  // 0. TOL computation
  m.get_nodal_metric2();

  if (ARTOL == 0.) ARTOL = -0.01;

  double TOL = ARTOL;
  if (ARTOL < 0.) // relative
  {
    E_Float Lmin, Lmax;
    E_Int imin, imax;
    ngon_type::edge_length_extrema(m.cnt, m.crd, Lmin, imin, Lmax, imax);
    TOL = -ARTOL * Lmin;
    assert(TOL > 0.); // no degen in m
  }

  E_Int npgs = m.cnt.size();

  // 1. compute centroids
  K_FLD::FloatArray Cs(3, npgs);
  for (E_Int i = 0; i < npgs; ++i)
  {
    const E_Int* nodes = m.cnt.get_facets_ptr(i);
    E_Int nb_nodes = m.cnt.stride(i);
    //K_MESH::Polygon::simplified_centroid(crd, nodes, nb_nodes, 1, Cs.col(i));
    //K_MESH::Polygon::centroid<3>(crd, nodes, nb_nodes, 1, Cs.col(i));
    K_MESH::Polygon::iso_barycenter(m.crd, nodes, nb_nodes, 1, Cs.col(i));
  }

  // 2. apply transfo
  K_FLD::FloatArray CsR(Cs);
  if (axis)
    NUGA::axial_rotate(CsR, center, axis, angle);
  else if (translation){
    for (size_t k = 0; k < CsR.cols(); ++k)
    {
      CsR(0, k) += translation[0];
      CsR(1, k) += translation[1];
      CsR(2, k) += translation[2];
    }
  }
  else
  {
    std::cout << "Error : no transfo defined" << std::endl;
    return;
  }

  ////////////////////////////////////////////////////////////////////////////

  acrd_t aCs(Cs);
  K_SEARCH::KdTree<> tree(aCs, EPSILON);

  std::vector<E_Int> nids(npgs, IDX_NONE);
  std::vector<double> d2s(npgs, IDX_NONE);

  ////////////////////////////////////////////////////////////////////////////

  // 3. append both
  K_FLD::FloatArray Cs_glog{ Cs };
  Cs_glog.pushBack(CsR);

#ifdef DEBUG_MATCH_PERIO
  {
    K_FLD::IntArray cntt(2, 1, 0);
    cntt(1, 0) = 1;
    medith::write("cloud", Cs_glog, cntt, "BAR");
  }
#endif

  // 4. Filtrate exact-matches

  K_FLD::ArrayAccessor<K_FLD::FloatArray> ca(Cs_glog);
  E_Int nmerges = ::merge(ca, TOL, nids);
  //std::cout << "nmerges : " << nmerges << std::endl;

  // check wrong auto_match
  for (size_t i = 0; i <npgs; ++i)
    assert(nids[i] == i);

  std::vector<E_Int> left, right, remain;
  for (size_t i = npgs; i <nids.size(); ++i)
  {
    if (nids[i] == i)
    {
      //remain.push_back(i-npgs);
      continue;
    }

    left.push_back(i - npgs);
    right.push_back(nids[i]);
    E_Int nid = nids[i];
    nids[i - npgs] = nid;
    nids[nid] = i - npgs;
  }

  nids.resize(npgs);

#ifdef DEBUG_MATCH_PERIO
  medith::write("left", m.crd, m.cnt, &left, 0);
  medith::write("right", m.crd, m.cnt, &right, 0);
  //medith::write("remain", crd, PGS, &remain, 0);
#endif

  // 5. Filtrate non-perio by boxes
  auto loc = m.get_localizer();
  std::vector<int> cands;

  K_FLD::FloatArray crdR(m.crd);
  if (axis)
    NUGA::axial_rotate(crdR, center, axis, angle);
  else if (translation)
    for (size_t k = 0; k < crdR.cols(); ++k)
    {
      crdR(0, k) += translation[0];
      crdR(1, k) += translation[1];
      crdR(2, k) += translation[2];
    }
  else
  {
    std::cout << "Error : no transfo defined" << std::endl;
    return;
  }

  acrd_t acrdR(crdR);
  acrd_t acrd(m.crd);

  for (size_t i = 0; i < npgs; ++i)
  {
    if (nids[i] != i) continue; // exact match already found

    auto face = m.element(i);

    K_SEARCH::BBox3D b;
    face.bbox(crdR, b); // compute box in rotated frame

    loc->get_candidates(b, cands);
    if (cands.empty()) continue;

    /*{
      K_FLD::FloatArray cc;
      K_FLD::IntArray ccnt;
      ngon_type ng;
      b.convert2NG(cc, ng);
      ng.export_to_array(ccnt);
      tp::write("D:\\slandier\\DATA\\tmp\\match_perio\\box.tp", cc, ccnt, "NGON");

      medith::write("face", crdR, m.cnt, i);
    }*/

    double n[3], nc[3];
    face.normal<K_FLD::FloatArray, 3>(crdR, n);

    double G[3];
    face.iso_barycenter<acrd_t, 3>(acrdR, G);

    for (size_t c = 0; c< cands.size(); ++c)
    {
      E_Int cand = cands[c];

      auto faceC = m.element(cand);
      faceC.normal<K_FLD::FloatArray, 3>(m.crd, nc);

      // true candidates must be colinear and opposedly oriented, so ps should near -1
      double ps = NUGA::dot<3>(n, nc);
      if (ps > -0.99) continue;

      // if G is in faceC, face is a piece of faceC
      bool is_in1;
      int err = faceC.fast_is_in_pred<DELAUNAY::Triangulator>(dt, m.crd, G, is_in1);

      //reciprocal test
      double GC[3];
      faceC.iso_barycenter<acrd_t, 3>(acrd, GC);
      bool  is_in2;
      err = face.fast_is_in_pred<DELAUNAY::Triangulator>(dt, crdR, GC, is_in2);

      if (!is_in1 && !is_in2) continue;

      if (is_in1 && !is_in2)
      {
        loc_face_to_bits[-(cand+1)].push_back(i); //neg storing to apply -Transfo
        continue;
      }

      if (is_in2 && !is_in1)
      {
        loc_face_to_bits[i].push_back(-(cand + 1)); //pos storing for i to apply +Transfo
        continue;
      }

      // if mutual inclusion : surface test
      double s1 = face.surface<3>(crdR);
      double s2 = faceC.surface<3>(m.crd);

      if (s1 < s2)
      {
        loc_face_to_bits[-(cand + 1)].push_back(i);
        continue;
      }
      else
      {
        loc_face_to_bits[i].push_back(-(cand + 1));
        continue;
      }
    }
  }

  // FILTER : ONLY ONE FACE PER PH PER ITER
  std::map<int, std::vector<int>> tmp;
  std::set<int> involved_phs;

  for (auto ii : loc_face_to_bits)
  {
    int id = ii.first;
    if (id < 0) id = -(id + 1);

    if (!involved_phs.insert(ancestors[id]).second) // already in
      continue;

    tmp.insert(std::make_pair(ii.first, ii.second));
  }

  loc_face_to_bits = tmp;

  // indirection to refer to volume mesh ids
  for (auto i : loc_face_to_bits)
  {
    int lface = i.first;
    int gface = 0;
    if (lface < 0)
    {
      lface = -(lface + 1);
      gface = oids[lface];
      gface = -(gface + 1);
    }
    else
      gface = oids[lface];

    auto gbits = i.second; //copy
    for (auto& k : gbits)
    {
      if (k < 0)
      {
        k = oids[-k - 1];
        k = -(k + 1);
      }
      else
        k = oids[k];
    }
    glob_face_to_bits[gface] = gbits;
  }

  //make it reciprocal (ONLY IF bits belong to the same cell)
  std::map<int, std::vector<int>> tmp1;
  for (auto i : glob_face_to_bits)
  {
    int master = i.first;
    auto& bits = i.second;
    tmp1[master] = bits;

    bool ancPH = ancestors[bits[0]];
    bool samePH = true;
    for (size_t k = 0; (k < bits.size() - 1) && samePH; ++k)
    {
      int bk = bits[k];
      int bkp1 = bits[k + 1];
      bk = (bk < 0) ? -(bk + 1) : bk;
      bkp1 = (bkp1 < 0) ? -(bkp1 + 1) : bkp1;

      samePH = (ancestors[bk] == ancestors[bkp1]);
    }
    if (!samePH) continue;

    for (auto& b : bits)tmp1[b].push_back(master);
  }
  glob_face_to_bits = tmp1;
}

///
void duplicate_and_move_period_faces
(NUGA::ph_mesh_t& m, const double* center, const double * axis, double angle, const double* translation,
 std::map<int, std::vector<int>>& face_to_bits)
{
  std::set<int> leftF; //apply +Transfo
  std::set<int> rightF;//apply -Transfo

  // 1. update ids (i.e. remove neg vals), separate ids in these two groups to apply appropriate transformation
  std::map<int, std::vector<int>> new_face_to_bits;
  for (auto i : face_to_bits)
  {
    int master = i.first;
    auto& bits = i.second;

    if (master < 0) {
      master = -master - 1;
      rightF.insert(master);
    }
    else
      leftF.insert(master);

    for (auto& b : bits)
    {
      if (b < 0) {
        b = -b - 1;
        rightF.insert(b);
      }
      else
        leftF.insert(b);
    }

    new_face_to_bits[master] = bits;                  // no more neg vals
  }

  face_to_bits = new_face_to_bits;//no more neg vals

  //2. extract these faces and apply Transfo

  std::vector<int> oids1;
  ngon_unit pgsL;
  std::vector<E_Int> lF(ALL(leftF));
  m.cnt.PGs.extract(lF, pgsL, oids1);

  //modify history such moved bits receive target location history (no sense for split master that will vanish anyway)
  //reset first
  pgsL._ancEs.clear();
  pgsL._ancEs.resize(2, pgsL.size(), IDX_NONE);
  K_CONNECT::IdTool::init_inc(pgsL._ancEs, 0, pgsL.size());
  int k = -1;
  for (auto i : lF)
  {
    ++k;
    auto tgtid = face_to_bits.find(i);
    assert(tgtid != face_to_bits.end());
    if (tgtid->second.size() != 1) continue; //split master
    pgsL._ancEs(0, k) = tgtid->second[0];
  }


  K_FLD::FloatArray crdL(m.crd);
  ngon_type::compact_to_used_nodes(pgsL, crdL);

  // reverse orient
  /*{
    NUGA::pg_smesh_t toto;
    toto.crd = crdL;
    toto.cnt = pgsL;
    toto.oriented = 0;
    toto.reverse_orient();
    pgsL = toto.cnt;
  }*/

  //
  if (axis)
    NUGA::axial_rotate(crdL, center, axis, angle);
  else if (translation)
    for (size_t k = 0; k < crdL.cols(); ++k)
    {
      crdL(0, k) += translation[0];
      crdL(1, k) += translation[1];
      crdL(2, k) += translation[2];
    }
  else
  {
    std::cout << "Error : no transfo defined" << std::endl;
    return;
  }


#ifdef DEBUG_MATCH_PERIO
  {
    ngon_type ngo(pgsL, false);
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    tp::write("D:\\slandier\\DATA\\tmp\\match_perio\\pgsL.tp", crdL, cnto, "NGON");
  }
#endif

  std::vector<int> oids2;
  ngon_unit pgsR;
  std::vector<E_Int> rF(ALL(rightF));
  m.cnt.PGs.extract(rF, pgsR, oids2);

  //modify history such moved bits receive target location history (no sense for split master that will vanish anyway)
  //reset first
  pgsR._ancEs.clear();
  pgsR._ancEs.resize(2, pgsR.size(), IDX_NONE);
  K_CONNECT::IdTool::init_inc(pgsR._ancEs, 0, pgsR.size());
  k = -1;
  for (auto i : rF)
  {
    ++k;
    auto tgtid = face_to_bits.find(i);
    assert(tgtid != face_to_bits.end());
    if (tgtid->second.size() != 1) continue; //split master
    pgsR._ancEs(0, k) = tgtid->second[0];
  }

  K_FLD::FloatArray crdR(m.crd);
  ngon_type::compact_to_used_nodes(pgsR, crdR);

  // reverse orient
  /*{
    NUGA::pg_smesh_t toto;
    toto.crd = crdR;
    toto.cnt = pgsR;
    toto.oriented = 0;
    toto.reverse_orient();
    pgsR = toto.cnt;
  }*/

  //
  if (axis)
    NUGA::axial_rotate(crdR, center, axis, -angle);
  else if (translation)
    for (size_t k = 0; k < crdR.cols(); ++k)
    {
      crdR(0, k) -= translation[0];
      crdR(1, k) -= translation[1];
      crdR(2, k) -= translation[2];
    }
  else
  {
    std::cout << "Error : no transfo defined" << std::endl;
    return;
  }

#ifdef DEBUG_MATCH_PERIO
  {
    ngon_type ngo(pgsR, false);
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    tp::write("D:\\slandier\\DATA\\tmp\\match_perio\\pgsR.tp", crdR, cnto, "NGON");
  }
#endif

  //3. append these faces (geometrically, with their own coordinates)

  std::vector<int> nids;
  K_CONNECT::IdTool::init_inc(nids, m.cnt.PGs.size());

  E_Int shift = m.cnt.PGs.size();

  // append L
  int pt_shift = m.crd.cols();
  pgsL.shift(pt_shift);
  m.crd.pushBack(crdL);
  m.cnt.PGs.append(pgsL);

  for (size_t i = 0; i < oids1.size(); ++i)
    nids[oids1[i]] = i + shift;

  shift = m.cnt.PGs.size();

  // append R
  pt_shift = m.crd.cols();
  pgsR.shift(pt_shift);
  m.crd.pushBack(crdR);
  m.cnt.PGs.append(pgsR);

  for (size_t i = 0; i < oids2.size(); ++i)
    nids[oids2[i]] = i + shift;

  //update face_to_bits
  for (auto& i : face_to_bits)
  {
    int master = i.first;
    auto& bits = i.second;
    for (auto& b : bits) b = nids[b];
  }
}

///
void sync_faces
(NUGA::ph_mesh_t& m, const std::map<int, std::vector<int>>& face_to_bits, double ARTOL)
{
  // 
  std::vector<int> molecPH;
  ngon_unit new_phs;

  // 1. replace faces and keep track of modified PHs
  std::set<E_Int> modifiedPHs;
  for (E_Int i = 0; i < m.cnt.PHs.size(); ++i)
  {
    const E_Int* pPGi = m.cnt.PHs.get_facets_ptr(i);
    E_Int nb_pgs = m.cnt.PHs.stride(i);

    molecPH.clear();

    for (E_Int p = 0; p < nb_pgs; ++p)
    {
      E_Int PGi = *(pPGi + p) - 1;

      auto itBits = face_to_bits.find(PGi);
      if (itBits == face_to_bits.end())
        molecPH.push_back(PGi + 1);
      else
      {
        modifiedPHs.insert(i);
        for (E_Int k = 0; k < itBits->second.size(); ++k)
          molecPH.push_back(itBits->second[k] + 1);
      }
    }

    new_phs.add(molecPH.size(), &molecPH[0]);
  }

  new_phs._type = m.cnt.PHs._type;  // hack to preserve flags (externality)
  new_phs._ancEs = m.cnt.PHs._ancEs;// hack

  m.cnt.PHs = new_phs;
  m.cnt.PHs.updateFacets();
  m.cnt.PHs.remove_duplicated(); //several occurence of the same face in each phs

                                 // 2. merge coincident nodes
  if (ARTOL == 0.) ARTOL = -0.01;

  double TOL = ARTOL;
  if (ARTOL < 0.) // relative
  {
    E_Float Lmin, Lmax;
    E_Int imin, imax;
    ngon_type::edge_length_extrema(m.cnt.PGs, m.crd, Lmin, imin, Lmax, imax);
    TOL = -ARTOL * Lmin;
    assert(TOL > 0.); // no degen in m
  }

  //std::cout << "TOL : " << TOL << std::endl;
  E_Int nb_merges = m.cnt.join_phs(m.crd, TOL);
  //std::cout << "nmerges : " << nb_merges << std::endl;

  //3. close_phs
  std::vector<E_Int> modPHs(ALL(modifiedPHs));
  ngon_type::close_phs(m.cnt, m.crd, &modPHs);

  //4. replace moved master by modified bits : i.e replace original bits by their modified version
  //4.a reverse face_to_bits
  std::map<int, std::vector<int>> face_to_bits_rev, tmp;
  for (auto i : face_to_bits)
  {
    auto bits = i.second;
    int f = i.first;
    for (auto k : bits)
      face_to_bits_rev[k].push_back(f);
  }
  //4.b keep only releavnt
  for (auto i : face_to_bits_rev)
  {
    if (i.second.size() == 1) continue;
    tmp[i.first] = i.second;
  }
  face_to_bits_rev = tmp;

  //4.c replace back
  if (!face_to_bits_rev.empty())
  {
    new_phs.clear();
    for (E_Int i = 0; i < m.cnt.PHs.size(); ++i)
    {
      const E_Int* pPGi = m.cnt.PHs.get_facets_ptr(i);
      E_Int nb_pgs = m.cnt.PHs.stride(i);

      molecPH.clear();

      for (E_Int p = 0; p < nb_pgs; ++p)
      {
        E_Int PGi = *(pPGi + p) - 1;

        auto itBits = face_to_bits_rev.find(PGi);
        if (itBits == face_to_bits_rev.end())
          molecPH.push_back(PGi + 1);
        else
        {
          for (E_Int k = 0; k < itBits->second.size(); ++k)
            molecPH.push_back(itBits->second[k] + 1);
        }
      }

      new_phs.add(molecPH.size(), &molecPH[0]);
    }

    new_phs._type = m.cnt.PHs._type;  // hack to preserve flags (externality)
    new_phs._ancEs = m.cnt.PHs._ancEs;// hack

    m.cnt.PHs = new_phs;
    m.cnt.PHs.updateFacets();

    ngon_type::close_phs(m.cnt, m.crd, &modPHs);
  }

  //5. clean
  std::vector<E_Int> pgnids, phnids;
  m.cnt.remove_unreferenced_pgs(pgnids, phnids);

  // assert closed
}

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
  do
  {
    std::map<int, std::vector<int>> glob_face_to_bits;
    detect_async_modified_faces(m, center, paxis, rot_angle, ptrans, artol, glob_face_to_bits);

    carry_on = !glob_face_to_bits.empty();
    if (!carry_on) break;

    duplicate_and_move_period_faces(m, center, paxis, rot_angle, ptrans, glob_face_to_bits);

    sync_faces(m, glob_face_to_bits, artol);

  } while (carry_on);

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
