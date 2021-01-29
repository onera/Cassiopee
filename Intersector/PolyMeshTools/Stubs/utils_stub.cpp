/*    
    Copyright 2013-2021 Onera.

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

# include "intersector.h"
# include "stub.h"

# include "Nuga/include/ngon_t.hxx"
# include "Nuga/include/Triangulator.h"


//=============================================================================
/* Creates 4 zones : 1) uncomputable polygons 2) uncomputable polyhedra 
   3) uncomputable polyhedra & neighbors 4) complementary of 3) */
//=============================================================================
PyObject* K_INTERSECTOR::extractUncomputables(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX*/
//=============================================================================
PyObject* K_INTERSECTOR::extractPathologicalCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}


//=============================================================================
/* Creates 2 zones : 1) outerlayer with firt neighborhoo 2) complementary */
//=============================================================================
PyObject* K_INTERSECTOR::extractOuterLayers(PyObject* self, PyObject* args)
{

  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extractNthCell(PyObject* self, PyObject* args)
{

  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extractBiggestCell(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::removeNthCell(PyObject* self, PyObject* args)
{

  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::detectIdenticalCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::detectOverConnectedFaces(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extractNthFace(PyObject* self, PyObject* args)
{

  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::checkCellsClosure(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  //K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  err = ngon_type::check_phs_closure(ngi);

  delete f; delete cn;

#ifdef E_DOUBLEINT
    return Py_BuildValue("l", long(err));
#else
    return Py_BuildValue("i", err);
#endif
}

PyObject* K_INTERSECTOR::checkCellsFlux(PyObject* self, PyObject* args)
{
  PyObject *arr, *PE;

  if (!PyArg_ParseTuple(args, "OO", &arr, &PE)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  // Check numpy (parentElement)
  FldArrayI* cFE;
  E_Int res = K_NUMPY::getFromNumpyArray(PE, cFE, true);

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (ngi.PGs.size() != cFE->getSize())
  {
    std::cout << "le ParentElment ne correpsond pas au nb de pgs" << std::endl;
    delete f; delete cn;
    return nullptr;
  }

  std::vector<E_Int> orient;
  E_Int imax=-1;
  E_Float fluxmax = -1;
  for (E_Int i=0; i < ngi.PHs.size(); ++i)
  {
    orient.clear();

    const E_Int* pF = ngi.PHs.get_facets_ptr(i);
    E_Int nbf = ngi.PHs.stride(i);
    orient.resize(nbf, 1);

    for (E_Int f = 0; f < nbf; ++f)
    {
      E_Int PGi = *(pF+f) - 1;
      //std::cout << "PGi bef wwong :" << PGi << std::endl;
      if ((*cFE)(PGi, 1) != i+1) orient[f] = -1;
      assert (((*cFE)(PGi, 1) == i+1) || ((*cFE)(PGi, 2) == i+1) );
    }

    //std::cout << "computing flux for PH : " << i << std::endl;
    K_MESH::Polyhedron<0> PH(ngi, i);
    E_Float flxVec[3];
    PH.flux(crd, &orient[0], flxVec);

    E_Float flux = ::sqrt(K_FUNC::sqrNorm<3>(flxVec));
    E_Float s = PH.surface(crd);

    flux /= s; // normalizing

    if (flux > fluxmax)
    {
      imax = i;
      fluxmax = flux;
    }
  }

  std::cout << "normalized max flux is : " << fluxmax << " reached at cell : " << imax << std::endl;

  delete f; delete cn;

  PyObject *l(PyList_New(0)), *tpl;

#ifdef E_DOUBLEINT
  tpl =  Py_BuildValue("l", long(imax));
   PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("i", imax);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(fluxmax));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", float(fluxmax);
  PyList_Append(l, tpl);
#endif
  
  return l;
}

PyObject* K_INTERSECTOR::checkCellsVolume(PyObject* self, PyObject* args)
{
  PyObject *arr, *PE;

  if (!PyArg_ParseTuple(args, "OO", &arr, &PE)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  // Check numpy (parentElement)
  FldArrayI* cFE;
  E_Int res = K_NUMPY::getFromNumpyArray(PE, cFE, true);

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (ngi.PGs.size() != cFE->getSize())
  {
    std::cout << "le ParentElment ne correpsond pas au nb de pgs" << std::endl;
    delete f; delete cn;
    return nullptr;
  }

  std::vector<E_Int> orient;
  E_Int imin=-1;
  E_Float vmin = NUGA::FLOAT_MAX;
  for (E_Int i=0; i < ngi.PHs.size(); ++i)
  {
    orient.clear();

    const E_Int* pF = ngi.PHs.get_facets_ptr(i);
    E_Int nbf = ngi.PHs.stride(i);
    orient.resize(nbf, 1);

    for (E_Int f = 0; f < nbf; ++f)
    {
      E_Int PGi = *(pF+f) - 1;
      //std::cout << "PGi bef wwong :" << PGi << std::endl;
      if ((*cFE)(PGi, 1) != i+1) orient[f] = -1;
      assert (((*cFE)(PGi, 1) == i+1) || ((*cFE)(PGi, 2) == i+1) );
    }

    //std::cout << "computing flux for PH : " << i << std::endl;
    K_MESH::Polyhedron<0> PH(ngi, i);
    double v;
    E_Int err = PH.volume<DELAUNAY::Triangulator>(crd, &orient[0], v);

    if (!err && v < vmin)
    {
      imin = i;
      vmin = v;
    }
    if (err)
    {
      std::cout << "error to triangulate cell " << i << "at face : " << err-1 << std::endl;
      //medith::write("badcell", crd, ngi, i);
      //medith::write("faultyPG", crd, ngi.PGs.get_facets_ptr(err-1), ngi.PGs.stride(err-1), 1);
    }
  }

  std::cout << "min vol is : " << vmin << " reached at cell : " << imin << std::endl;

  delete f; delete cn;

  PyObject *l(PyList_New(0)), *tpl;

#ifdef E_DOUBLEINT
  tpl =  Py_BuildValue("l", long(imin));
   PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("i", imin);
  PyList_Append(l, tpl);
#endif

#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", double(vmin));
  PyList_Append(l, tpl);
#else
  tpl =  Py_BuildValue("f", float(vmin);
  PyList_Append(l, tpl);
#endif
  
  return l;
}

PyObject* K_INTERSECTOR::volume(PyObject* self, PyObject* args)
{
  PyObject *arr, *axcelln;

  if (!PyArg_ParseTuple(args, "OO", &arr, &axcelln)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  // Check numpy (xcelln)
  bool use_xcelln = false;

  K_FLD::FloatArray* xcelln(nullptr);
  K_FLD::IntArray *cn1(0);
  char *varString1;
  E_Int ni, nj, nk;
  E_Int res = 0;
  if (axcelln != Py_None) res = K_ARRAY::getFromArray(axcelln, varString1, xcelln, ni, nj, nk, cn1, eltType);
  if (res == 1) use_xcelln = true;

  //std::cout << py_xcelln << std::endl;

  // E_Int res = 0;
  // if (py_xcelln != Py_None)
  // {
  //   std::cout << "get numpy " << std::endl;

  //   res = K_NUMPY::getFromNumpyArray(py_xcelln, xcelln, true);
  //   std::cout << xcelln << std::endl;
  //   if (res == 1) use_xcelln = true;
  // }

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  // std::cout << "use_xcelln ? " << use_xcelln << std::endl;
  // std::cout << "xcelln ? " << xcelln << std::endl;
  // std::cout << "res : " << res << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  if (use_xcelln && (xcelln != nullptr) && (ngi.PHs.size() != xcelln->getSize()))
  {
    std::cout << "le champ xcelln ne correpsond pas au nb de polyedres => pas pris en compte" << std::endl;
    std::cout << "nb phs : " << ngi.PHs.size() << std::endl;
    std::cout << "taille xcelln : " << xcelln->getSize() << std::endl;
    use_xcelln = false;
  }

  std::vector<E_Float> vols;
  ngon_type::volumes<DELAUNAY::Triangulator>(crd, ngi, vols, false/*not all cvx*/, true/*new algo*/);

  //std::cout << "use_xcelln ?" << use_xcelln << std::endl;
  E_Float V = 0.;
  if (use_xcelln)
  {
    for (size_t i = 0; i < vols.size(); ++i)
      V += vols[i] * (*xcelln)[i];
  }
  else
    for (size_t i = 0; i < vols.size(); ++i)
      V += vols[i];

  delete f; delete cn;

#ifdef E_DOUBLEREAL
  return Py_BuildValue("d", V);
#else
    return Py_BuildValue("f", float(V));
#endif
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::removeBaffles(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::checkForDegenCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::checkForBigCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::edgeLengthExtrema(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::computeAspectRatio(PyObject* self, PyObject* args)
{

  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extrudeBC(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extrudeSurf(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extrudeRevolSurf(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* reorient specified polygons. */
//=============================================================================
PyObject* K_INTERSECTOR::reorientSpecifiedFaces(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_pgs;
  E_Int dir(1); //1 : outward -1 : inward

  if (!PYPARSETUPLEI(args, "OOl", "OOi", &arr, &py_pgs, &dir)) return NULL;

  if (dir != -1 && dir != 1) dir = 1;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //std::cout << "before numpy" << std::endl;

  E_Int res=0;
  E_Int* pgsList=NULL;
  E_Int size, nfld;
  if (py_pgs != Py_None)
    res = K_NUMPY::getFromNumpyArray(py_pgs, pgsList, size, nfld, true/*shared*/, 0);

  //std::cout << "after numpy" << std::endl;

  std::cout << res << std::endl;

  if (res != 1) return NULL;
  
  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ngio(cnt);

  //std::cout << "after ngio construct" << std::endl;

  std::vector<E_Int> plist, oids;
  plist.insert(plist.end(), pgsList, pgsList+size);

  //std::cout << "after insert" << std::endl;
  //K_CONNECT::IdTool::shift(plist, -1);

  //std::cout << "min pg specified : " << *std::min_element(pgsList, pgsList+size) << std::endl;
  //std::cout << "max pg specified : " << *std::max_element(pgsList, pgsList+size) << std::endl;

  ngon_unit pgs;
  ngio.PGs.extract(plist, pgs, oids);

  std::vector<E_Int> orient;
  ngon_type::reorient_connex_PGs(pgs, (dir==-1), orient);

  // replace reverted polygons
  E_Int count(0);
  E_Int nb_pgs = pgs.size();
  for (E_Int i=0; i < nb_pgs; ++i)
  {
    if (orient[i] == 1) continue;

    ++count;
    
    E_Int PGi = oids[i];

    E_Int s = ngio.PGs.stride(PGi);
    E_Int* p = ngio.PGs.get_facets_ptr(PGi);
    std::reverse(p, p + s);
  }

  //std::cout << "nb of reoriented : "  << count  << " over " << size << " in pglist"<< std::endl;
    
  K_FLD::IntArray cnto;
  ngio.export_to_array(cnto);
  
  // pushing out the mesh
  PyObject *tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);   
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::reorient(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Int dir(1); //1 : outward -1 : inward

  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arr, &dir)) return NULL;

  if (dir != -1 && dir != 1) dir = 1;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ngio(cnt);

  ngio.flag_externals(INITIAL_SKIN);
  bool has_been_reversed;
  DELAUNAY::Triangulator t;
  err = ngon_type::reorient_skins(t, crd, ngio, has_been_reversed);

  if (dir == -1)
  {
    Vector_t<E_Int> oids;
    ngon_unit pg_ext;
    ngio.PGs.extract_of_type(INITIAL_SKIN, pg_ext, oids);

    E_Int nb_pgs = pg_ext.size();
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      E_Int PGi = oids[i];

      E_Int s = ngio.PGs.stride(PGi);
      E_Int* p = ngio.PGs.get_facets_ptr(PGi);
      std::reverse(p, p + s);
    }
  }
    
  K_FLD::IntArray cnto;
  ngio.export_to_array(cnto);
  
  // pushing out the mesh
  PyObject *tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);   
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::externalFaces(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::diffMesh(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::statsUncomputableFaces(PyObject* self, PyObject* args)
{

  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::statsSize(PyObject* self, PyObject* args)
{

  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::convert2Polyhedron(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::oneZonePerCell(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::oneZonePerFace(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::immerseNodes(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::closeCells(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

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
  ngon_type ngi(cnt);

  ngon_type::close_phs(ngi, crd);

  K_FLD::IntArray cnto;
  ngi.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* Converts a surfacic NGON from Cassiopee format to nuga format*/
//=============================================================================
PyObject* K_INTERSECTOR::convertNGON2DToNGON3D(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

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
  ngon_type ngi(cnt), outer, remaining;

  ngi.export_surfacic_view(cnt);

  ngon_unit pgs(&cnt[0]);

  PyObject* tpl = NULL;

  if (cnt.cols() != 0)
  {
  	ngon_type ngo(pgs, true);
  	K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  }

  delete f; delete cn;

  return tpl;
}

//=============================================================================
/* Converts a surfacic Basic to NGON nuga format*/
//=============================================================================
PyObject* K_INTERSECTOR::convertBasic2NGONFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_BASICF(arr, f, cn, varString, eltType);
  if (err)
  {
    std::cout << "convertBasic2NGONFaces : ERROR : " << err << std::endl;
    return nullptr;
  }

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  ngon_unit pgs;
  ngon_unit::convert_fixed_stride_to_ngon_unit(cnt, 1, pgs);

  PyObject* tpl = NULL;

  if (pgs.size() != 0)
  {
    ngon_type wNG(pgs, true);
    K_FLD::IntArray cnto;
    wNG.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  }

  delete f; delete cn;

  return tpl;
}

//=============================================================================
/* remove any cell contributing to a non-manifold boundary */
//=============================================================================
PyObject* K_INTERSECTOR::removeNonManifoldExternalCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* Computes centroids */
//=============================================================================
PyObject* K_INTERSECTOR::centroids(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  using ngon_type = ngon_t<K_FLD::IntArray>;
  ngon_type::eGEODIM geodim = ngon_type::get_ngon_geodim(cnt);

  //std::cout << "GEO dIM ? " << geodim << std::endl;

  if (geodim == ngon_type::eGEODIM::ERROR)
  {
    std::cout << "centroids : Input Error : mesh is corrupted." << std::endl;
    return nullptr;
  }
  if (geodim == ngon_type::eGEODIM::MIXED)
  {
    std::cout << "centroids : Input Error : mesh mixed elt types (lineic and/or surfacic and /or volumic." << std::endl;
    return nullptr;
  }
  if (geodim == ngon_type::eGEODIM::LINEIC)
  {
    std::cout << "centroids : Unsupported : lineic NGON are not handled." << std::endl;
    return nullptr;
  }

  // so either SURFACIC, SURFACIC_CASSIOPEE or VOLUMIC

  if (geodim == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
  {
    ngon_type ng(cnt);
    // convert to SURFACIC (NUGA)
    K_FLD::IntArray cnt1;
    ng.export_surfacic_view(cnt1);
    //std::cout << "exported" << std::endl;
    geodim = ngon_type::eGEODIM::SURFACIC;
    cnt=cnt1;
  }

  ngon_type ngi(cnt);

  K_FLD::FloatArray cents;

  if (geodim == ngon_type::eGEODIM::SURFACIC)
    ngon_type::centroids(ngi.PGs, crd, cents);
  else // (geodim == ngon_type::eGEODIM::VOLUMIC)
    ngon_type::centroids<DELAUNAY::Triangulator>(ngi, crd, cents);

  K_FLD::IntArray cnto;

  PyObject* tpl = K_ARRAY::buildArray(cents, varString, cnto, 0, "NODE", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* Computes volumes */
//=============================================================================
PyObject* K_INTERSECTOR::volumes(PyObject* self, PyObject* args)
{
  PyObject *arr;
  E_Int algo(0);
  E_Int all_pgs_cvx(0);

  if (!PYPARSETUPLEI(args, "Oll", "Oii", &arr, &algo, &all_pgs_cvx)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  std::vector<E_Float> vols;
  ngon_type::volumes<DELAUNAY::Triangulator>(crd, ngi, vols, (all_pgs_cvx == 1), (algo == 1));

  size_t sz = vols.size();
  FloatArray ar(1, sz);
  for (size_t i = 0; i < sz; ++i) 
    ar[i] = vols[i];

  PyObject* tpl = K_ARRAY::buildArray(ar, "volumes", *cn, -1, "NGON", true);

  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* retrieves any polygon that are overlapping */
//=============================================================================
PyObject* K_INTERSECTOR::getOverlappingFaces(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* retrieves any cells that are colliding */
//=============================================================================
PyObject* K_INTERSECTOR::getCollidingCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::getNthNeighborhood(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* retrieves any polygon that are connecting 2 aniso HEXA */
//=============================================================================
PyObject* K_INTERSECTOR::getAnisoInnerFaces(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;

}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::merge(PyObject* self, PyObject* args)
{

  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::oneph(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::concatenate(PyObject* self, PyObject* args)
{

  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::drawOrientation(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::getFaceIdsWithCentroids(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::getFaceIdsCollidingVertex(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::getCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=======================  Intersector/PolyMeshTools/utils.cpp ====================
