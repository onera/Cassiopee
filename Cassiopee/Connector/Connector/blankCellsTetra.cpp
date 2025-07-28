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

# include "connector.h"
#include "maskGen.h"
#include "Nuga/include/Collider.h"
#include "Nuga/include/ZoneExtractor.h"
#include <memory>
#include "Nuga/include/ngacc.hxx"
#include "Nuga/include/localizer.hxx"
#include "Nuga/include/collider.hxx"
#include "Nuga/include/Pentahedron.h"
#include "Nuga/include/Pyramid.h"
//#include <sstream>
//#include <fstream>
//#include <iostream>
//#include "chrono.h"

using namespace std;
using namespace K_FLD;

#define TETRA_HOOK_ID 9
#define TETRA_HOOK_DATA_NB 6
#define TRI_HOOK_ID 10
#define TRI_HOOK_DATA_NB 3

//============================================================================
/* Create a tet mask and returns a hook */
//============================================================================
PyObject* K_CONNECTOR::createTetraMask(PyObject* self, PyObject* args)
{
  PyObject* maskV, *maskS;
  E_Float tol;
  char *varStringV, *varStringS, *eltTypeV, *eltTypeS;
  E_Int ni, nj, nk;
  K_FLD::FldArrayF *fV(0), *fS(0);
  K_FLD::FldArrayI *cV(0), *cS(0);

  if (!PYPARSETUPLE_(args, OO_ R_, &maskV, &maskS, &tol)) return NULL;

  // Extraction des donnees
  // Volumic mask
  E_Int res = K_ARRAY::getFromArray3(maskV, varStringV, fV, ni, nj, nk, cV, eltTypeV);
  
  bool ok = (res == 2) && (strcmp(eltTypeV, "TETRA") == 0);
  if (!ok)
  {
    //std::cout << "res : " << res << std::endl;
    //std::cout << "eltTypeV : " << *eltTypeV << std::endl;
    
    PyErr_SetString(PyExc_TypeError,
		    "blankCellsTetra: input mask must be an unstructured Tetrahedral volume mesh.");
    return NULL;
  }
  
  //std::cout << "HOOK CREATION : mask size : " << cV->getNfld() << "/" << cV->getSize() << std::endl;
       
  E_Int posx = K_ARRAY::isCoordinateXPresent(varStringV);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varStringV);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varStringV);

  void* maskingV = new K_CONNECTOR::maskGen(*fV, posx+1, posy+1, posz+1, *cV, K_CONNECTOR::maskGen::TH4, tol);
  
  // Mask exterior faces
  res = K_ARRAY::getFromArray3(maskS, varStringS, fS, ni, nj, nk, cS, eltTypeS);
  ok = ( (res == 2) && (strcmp(eltTypeS, "TRI") == 0) );
  if (!ok)
  {
    //std::cout << "res : " << res << std::endl;
    //std::cout << "eltS : " << *eltS << std::endl;
    
    PyErr_SetString(PyExc_TypeError,
		    "blankCellsTetra: input mask skin must be triangular.");
    return NULL;
  }
       
  posx = K_ARRAY::isCoordinateXPresent(varStringS);
  posy = K_ARRAY::isCoordinateYPresent(varStringS);
  posz = K_ARRAY::isCoordinateZPresent(varStringS);

  void* maskingS = new K_CONNECTOR::maskGen(*fS, posx+1, posy+1, posz+1, *cS, K_CONNECTOR::maskGen::TRI, tol);
    
  PyObject* hook;
  E_Int* type = new E_Int; *type = TETRA_HOOK_ID;
  void** packet = new void* [TETRA_HOOK_DATA_NB+1];// +1 : adding id;
  packet[0] = type; // hook type
  packet[1] = maskingV; //to delete when releasing the hook
  packet[2] = fV; //to delete when releasing the hook
  packet[3] = cV; //to delete when releasing the hook
  packet[4] = maskingS; //to delete when releasing the hook
  packet[5] = fS; //to delete when releasing the hook
  packet[6] = cS; //to delete when releasing the hook
  
  //std::cout << "HOOK CREATE :: FMASK SZ : " <<  fskin->getSize() << std::endl;
  //std::cout << "HOOK CREATE :: cmask SZ : " <<  skin->getSize() << std::endl;
  
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif
  
  //RELEASESHAREDU(maskV, fV, cV);
  //RELEASESHAREDU(maskS, fS, cS);
  
  return hook;
}

//============================================================================
/* Create a tet mask and returns a hook */
//============================================================================
PyObject* K_CONNECTOR::createTriMask(PyObject* self, PyObject* args)
{
  PyObject* maskT;
  E_Float tol;
  char *varStringMsk, *eltTypeMsk;
  E_Int ni, nj, nk;
  K_FLD::FldArrayF *fmask(0);
  K_FLD::FldArrayI *cmask(0);
  
  if (!PYPARSETUPLE_(args, O_ R_, &maskT, &tol)) return NULL;
  

  /////////////////////////////////////////////////////////////////////////
  // Extraction des donnees
  
  E_Int res = K_ARRAY::getFromArray3(maskT, varStringMsk, fmask, ni, nj, nk, cmask, eltTypeMsk);
  bool err = (res != 2);
  err &=  (strcmp(eltTypeMsk, "TRI") != 0);
  if (err)
  {
    //std::cout << "res : " << res << std::endl;
    //std::cout << "eltTypeMsk : " << *eltTypeMsk << std::endl;
    
    PyErr_SetString(PyExc_TypeError,
		    "blankCellsTri: input mask must be an unstructured Triangulated surface mesh.");
    return NULL;
  }
  
  K_CONNECTOR::maskGen::eType ELT = K_CONNECTOR::maskGen::TRI;
       
  E_Int posx = K_ARRAY::isCoordinateXPresent(varStringMsk);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varStringMsk);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varStringMsk);

  void* maskingTool = new K_CONNECTOR::maskGen(*fmask, posx+1, posy+1, posz+1, *cmask, ELT, tol);
    
  PyObject* hook;
  E_Int* type = new E_Int; *type = TRI_HOOK_ID;
  void** packet = new void* [4];
  packet[0] = type; // hook type
  packet[1] = maskingTool;
  packet[2] = fmask; //to delete when releasing the hook
  packet[3] = cmask; //to delete when releasing the hook
  
  //std::cout << "HOOK CREATE :: FMASK SZ : " <<  fmask->getSize() << std::endl;
  //std::cout << "HOOK CREATE :: cmask SZ : " <<  cmask->getSize() << std::endl;
  
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif
  
  //RELEASESHAREDU(maskT, fmask, cmask);
  return hook;
}

//=============================================================================
/* 
   Free hook
 */
//=============================================================================
PyObject* K_CONNECTOR::deleteTetraMask(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PyArg_ParseTuple(args, "O", &hook)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0]; // type of hook

  assert (*type == TETRA_HOOK_ID); // same as creation for volume tree with a surface tree joined

  K_CONNECTOR::maskGen* ptMask = (K_CONNECTOR::maskGen*)packet[1];
  delete ptMask;
  K_FLD::FldArrayF* ptFldF = (K_FLD::FldArrayF*)packet[2];
  delete ptFldF;
  K_FLD::FldArrayI* ptFldI = (K_FLD::FldArrayI*)packet[3];
  delete ptFldI;
  ptMask = (K_CONNECTOR::maskGen*)packet[4];
  delete ptMask;
  ptFldF = (K_FLD::FldArrayF*)packet[5];
  delete ptFldF;
  ptFldI = (K_FLD::FldArrayI*)packet[6];
  delete ptFldI;
  delete type;
  delete [] packet;
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* 
   Free hook
 */
//=============================================================================
PyObject* K_CONNECTOR::deleteTriMask(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PyArg_ParseTuple(args, "O", &hook)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0]; // type of hook

  assert (*type == TRI_HOOK_ID); // K_CONNECTOR::maskT4Fld*

  K_CONNECTOR::maskGen* pt = (K_CONNECTOR::maskGen*)packet[1];
  delete pt;
  
  K_FLD::FldArrayF* ptfmask = (K_FLD::FldArrayF*)packet[2];
  delete ptfmask;
  K_FLD::FldArrayI* ptcmask = (K_FLD::FldArrayI*)packet[3];
  delete ptcmask;
  
  delete type; delete [] packet;
  Py_INCREF(Py_None);
  return Py_None;
}

//static int xcount=0;

struct MaskEntity
{
  K_CONNECTOR::maskGen* _main; //surfacic for TRI or volumic for tet
  K_CONNECTOR::maskGen* _auxiliary; //surfacic for tet mask
};

template <typename ELT_t>
E_Int do_the_blanking(E_Int blankingType, MaskEntity& maskE, const K_FLD::FldArrayF& fmesh, E_Int posx, E_Int posy, E_Int posz,
                      const K_FLD::FldArrayI* cmesh, E_Int CELLNVAL, E_Int overwrite, K_FLD::FldArrayI& cN);

//============================================================================
/* Blank cells defined in arrays by a Tetra mesh mask (as an input hook) */
//============================================================================
PyObject* K_CONNECTOR::blankCellsTetra(PyObject* self, PyObject* args)
{
  char* cellNName;
  PyObject* mesh;
  PyObject* celln;
  PyObject* maskHook;
  E_Int blankingType, CELLNVAL(0), overwrite(0);
  char *varString, *varStringC, *eltType(0), *eltTypeC;
  E_Int ni, nj, nk;
  K_FLD::FldArrayF *fmesh(0), *fC(0);
  K_FLD::FldArrayI *cmesh(0), *cC(0);
  
  if (!PYPARSETUPLE_(args, OOO_ III_ S_, 
                    &mesh, &celln, &maskHook, &blankingType, &CELLNVAL, &overwrite, &cellNName))
  {
    return NULL;
  }

  if (blankingType < 0 || blankingType > 1) 
  {
    printf("Warning: blankCellsTetra: blankingType is invalid. Set to default (1)\n");
    blankingType = 1;
  }
  
  if (maskHook == 0)
  {
    PyErr_SetString(PyExc_ValueError,
		    "blankCellsTetra: the hook for Tetra mask is missing.");
    return NULL;
  }

  if (CELLNVAL == 1) CELLNVAL=0; //1 is reserved as outside color so can't be inside color at the same time.
  
  /////////////////////////////////////////////////////////////////////////
  // Extraction des donnees
  
  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(maskHook);
#else
  packet = (void**) PyCapsule_GetPointer(maskHook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];

  if (*type != TETRA_HOOK_ID && *type != TRI_HOOK_ID)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "blankCellsTetra: this function requires a identify K_CONNECTOR::maskT4Fld hook.");
    return NULL;
  }

  E_Int res = K_ARRAY::getFromArray3(mesh, varString, fmesh, ni, nj, nk, cmesh, eltType);
  std::unique_ptr<K_FLD::FldArrayF> afmesh(fmesh); // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<K_FLD::FldArrayI> acmesh(cmesh); // to avoid to call explicit delete at several places in the code.
  
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                      "blankCellsTetra: input mesh must contain coordinates.");
    return NULL;
  }

  if (blankingType == 1 /*cell _intersect*/)
  {
    if ( (res != 2) || ((res == 2) && (strstr(eltType, "TETRA") == 0 && strstr(eltType, "PYRA") == 0 && strstr(eltType, "PENTA") == 0 && strstr(eltType, "HEXA") == 0 && strstr(eltType, "NGON") == 0) ) )
    {
     PyErr_SetString(PyExc_ValueError,
		   "blankCellsTetra: the mesh can currently only be of type TETRA, PYRA, PENTA, HEXA or NGON in cell_intersect mode.");
      return NULL;
    }
  }

  E_Int res2 = K_ARRAY::getFromArray3(celln, varStringC, fC, ni, nj, nk, cC, eltTypeC);
  std::unique_ptr<K_FLD::FldArrayF> afC(fC); // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<K_FLD::FldArrayI> acC(cC); // to avoid to call explicit delete at several places in the code.
   
  if (res2 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                      "blankCellsTetra: cellN must be specified upon entry.");
    return NULL;
  }
  bool struct_celln = (res2 == 1);

  // Blanking
  K_CONNECTOR::maskGen* mask = (K_CONNECTOR::maskGen*)(packet[1]);

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString) +1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString) +1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString) +1;
  
  if (mask == 0)
  {
    PyErr_SetString(PyExc_ValueError,
	     "blankCellsTetra: the input hook is invalid.");
    return NULL;
  }

  if (mask->nb_points() == 0 || mask->nb_elts() == 0)
  {
    //std::cout << "nb points/tets : " << mask->nb_points() << "/" << mask->nb_elts() << std::endl;
    PyErr_SetString(PyExc_ValueError,
	     "blankCellsTetra: the input mask is empty.");
    return NULL;
  }

  MaskEntity maske;
  maske._main=maske._auxiliary=mask;

  if (blankingType == 1 && *type == TETRA_HOOK_ID) //tetra mask
  {
    maske._auxiliary = (K_CONNECTOR::maskGen*)(packet[4]);
  }

  size_t sz = 0;
  if (blankingType == 0)
    sz = fmesh->getSize();
  else if (eltType && strstr(eltType, "NGON") == 0)
    sz = cmesh->getSize();
  else if (eltType && strstr(eltType, "NGON") != 0)
  {
    E_Int pos = (*cmesh)[1] + 2;
    sz = (*cmesh)[pos];
    //std::cout << "nb phs : " << sz << std::endl;
  }

  if (sz == 0)
  {
    PyErr_SetString(PyExc_ValueError,
       "blankCellsTetra: the input combination is not handled.");
    return NULL;
  }

  K_FLD::FldArrayI cN(sz);
  for (size_t i = 0; i < sz; ++i) cN[i] = E_Int((*fC)[i]);
  
  E_Int err = 0;
  if (eltType && strstr(eltType, "TETRA") != 0)
    err = do_the_blanking<K_MESH::Tetrahedron>(blankingType, maske, *fmesh, posx, posy, posz, cmesh, CELLNVAL, overwrite, cN);
  else if (eltType && strstr(eltType, "PYRA") != 0)
    err = do_the_blanking<K_MESH::Pyramid>(blankingType, maske, *fmesh, posx, posy, posz, cmesh, CELLNVAL, overwrite, cN);
  else if (eltType && strstr(eltType, "PENTA") != 0)
    err = do_the_blanking<K_MESH::Pentahedron>(blankingType, maske, *fmesh, posx, posy, posz, cmesh, CELLNVAL, overwrite, cN);
  else if (eltType && strstr(eltType, "HEXA") != 0)
    err = do_the_blanking<K_MESH::Hexahedron>(blankingType, maske, *fmesh, posx, posy, posz, cmesh, CELLNVAL, overwrite, cN);
  else if (eltType && strstr(eltType, "NGON") != 0)
    err = do_the_blanking<K_MESH::Polyhedron<UNKNOWN> >(blankingType, maske, *fmesh, posx, posy, posz, cmesh, CELLNVAL, overwrite, cN);
  else // STRUCTURED
  {
  	do_the_blanking<K_MESH::Hexahedron/*dummy*/>(blankingType, maske, *fmesh, posx, posy, posz, cmesh, CELLNVAL, overwrite, cN);
  }

  if (err)
  {
    PyErr_SetString(PyExc_TypeError, "blankCellsTetra: failed.");
    return NULL;
  }
  
  K_FLD::FldArrayF cellnout(sz);
  for (size_t i = 0; i < sz; ++i) cellnout[i] = E_Float(cN[i]);
  
  PyObject* tpl = NULL;
  if (struct_celln)
    tpl = K_ARRAY::buildArray(cellnout, cellNName, ni, nj, nk);
  else
    tpl = K_ARRAY::buildArray(cellnout, cellNName, *cC, -1, eltTypeC, false);
  
  //RELEASESHAREDB(res, mesh, fmesh, cmesh);
  //RELEASESHAREDB(res2, celln, fC, cC);
  return tpl;
}

///
inline E_Int __do_the_node_blanking
(K_CONNECTOR::maskGen* mask, const K_FLD::FldArrayF& fmesh, E_Int posx, E_Int posy, E_Int posz, E_Int CELLNVAL, E_Int overwrite, K_FLD::FldArrayI& cN)
{
  // node_in or cell_in
  
  //Preconditioning : mask
  K_CONNECT::FldZoneExtractor extractor;      
  
  //mask bbox
  const K_SEARCH::BBox3D* bb = mask->getTree()->getGlobalBox();
  std::vector<E_Int> subset;
  K_FLD::ArrayAccessor<K_FLD::FldArrayF > acrd(fmesh, posx, posy, posz);
  extractor.getInBox<3>(bb->minB, bb->maxB, E_EPSILON, acrd, subset);
  
  if (subset.empty()) return 0;
  size_t sz = subset.size();
  
  K_FLD::FldArrayF fzoom;
  fzoom.resize(sz, 3);
  E_Float P[3];
  for (size_t i=0; i < sz; ++i)
  {
    acrd.getEntry(subset[i], P);
    for (E_Int k=0; k < 3; ++k) fzoom(i, k+1) = P[k];
  }
  
  K_FLD::FldArrayI cellNzoom;
  E_Int err = mask->blank(fzoom, 1,2,3, cellNzoom, CELLNVAL, overwrite);
  if (err)
    return err;
  
  // SET CELLN ////////////////////////////////////////////////////////////////
  if (cN.getSize() < fmesh.getSize())
    cN.resize(1, fmesh.getSize(), VISIBLE);

  if (overwrite)
  {
    for (size_t i=0; i < sz; ++i)
    {
      cN[subset[i]] = cellNzoom[i];
    }
  }
  else //append only blanked values
  {
    for (size_t i=0; i < sz; ++i)
    {
      if (cellNzoom[i]==CELLNVAL) cN[subset[i]] = CELLNVAL;
    }
  }
  return 0;
}

/// 
template <typename ELT_t>
inline E_Int __do_the_cell_blanking
(MaskEntity& maskE, const K_FLD::FldArrayF& fmesh, E_Int posx, E_Int posy, E_Int posz,
 const K_FLD::FldArrayI& cmesh,  E_Int CELLNVAL, E_Int overwrite, K_FLD::FldArrayI& cN)
{  
  // Preconditioning
  // Mask bbox
  const K_SEARCH::BBox3D* bb = maskE._main->getTree()->getGlobalBox();
  K_FLD::ArrayAccessor<K_FLD::FldArrayF > acrd(fmesh, posx, posy, posz);
  
  K_CONNECT::FldZoneExtractor extractor;
  K_FLD::ArrayAccessor<K_FLD::FldArrayI > acnt(cmesh, -1);
  
  std::vector<E_Int> subset;
  extractor.getInBox<3, ELT_t>(bb->minB, bb->maxB, E_EPSILON, acrd, acnt, subset);
  
  // CB : essai OMP
  //std::vector<E_Int> subset;
  //extractor.getInBox_omp<3>(bb->minB, bb->maxB, E_EPSILON, acrd, acnt, subset);
  
  if (subset.empty()) return 0;
  
  //std::cout << "nb of extracted : " << subset.size() << std::endl;
    
  size_t sz1 = subset.size();
  E_Int ROWS(cmesh.getNfld());
  K_FLD::FldArrayI connectZoom(sz1, ROWS);
  
  for (E_Int j = 0; j < ROWS; ++j)
  {
    E_Int* pZ = connectZoom.begin(j+1);
    const E_Int* pM = cmesh.begin(j+1);
    for (size_t i = 0; i < sz1; ++i) pZ[i] = pM[subset[i]];
  }
    
  // FIRST PASS : CELL IN
  typedef K_FLD::FldArrayF Coordinate_t;
  typedef K_FLD::FldArrayI Connectivity_t;

  Coordinate_t ZoomCentroids;
  ZoomCentroids.reserve(3, connectZoom.getSize());
  E_Float G[3];
  K_FLD::ArrayAccessor<Coordinate_t> ac(fmesh, posx, posy, posz);
  K_FLD::ArrayAccessor<Connectivity_t> acvz(connectZoom, -1);
  K_FLD::IntArray e;
  e.reserve(1,8);
  E_Int*pt = e.begin();
  for (size_t i=0; i < sz1; ++i)
  {    
    acvz.getEntry(i, pt);
    ELT_t Hi(pt);
    Hi.template iso_barycenter< K_FLD::ArrayAccessor<Coordinate_t> >(ac, G);
    ZoomCentroids.pushBack(G, G+3);
  }

  Connectivity_t cellNzoom(ZoomCentroids.getSize(), 1);
  cellNzoom.setAllValuesAt(1);
  E_Int err = maskE._main->blank(ZoomCentroids, 1,2,3, cellNzoom, CELLNVAL, overwrite);
  if (err) return err;
        
  // reduce the work set to those inside the zoom but not masked yet
  E_Int sz2=0, k=0;
  for (E_Int i = 0; i < cellNzoom.getSize(); ++i)
    if (cellNzoom[i] != CELLNVAL) ++sz2;

  Connectivity_t connectZoom2(sz2, ROWS);
  std::vector<E_Int> oldIds2;
  oldIds2.reserve(sz2);     
  for (E_Int i = 0; i < cellNzoom.getSize(); ++i)
  {
    if (cellNzoom[i]!=CELLNVAL)
    {
      for (E_Int j = 0; j < ROWS; ++j)
        connectZoom2(k,j+1) = connectZoom(i,j+1);
      oldIds2.push_back(i);
      ++k;
    }
  }
  
  //std::cout << " remaining unmasked in zoom :" << k << std::endl;

  // SECOND PASS : CELL INTERSECT
  std::vector<E_Int> is_colliding;
  if (!oldIds2.empty())
  {    
    using loc_t = NUGA::localizer<K_SEARCH::BbTree3D>;
    loc_t loc(*maskE._auxiliary->getTree(), maskE._auxiliary->tolerance());

    std::vector<E_Int> ids1;
    const K_FLD::ArrayAccessor<K_FLD::FldArrayF>& acrd1 = maskE._auxiliary->coords();
    const Connectivity_t& cnt1 = maskE._auxiliary->connect().array();
    
    const Coordinate_t& crd2 = fmesh;
    const Connectivity_t& cnt2 = connectZoom2;
                      
    NUGA::COLLIDE::compute<Coordinate_t, Connectivity_t, Connectivity_t, loc_t, K_MESH::Triangle, ELT_t, 3>(acrd1, cnt1, loc,
                                                        crd2, posx, posy, posz, cnt2, 
                                                        ids1, is_colliding, maskE._auxiliary->tolerance());
  }
  
  // SET CELLN ////////////////////////////////////////////////////////////////
  if (cN.getSize() < cmesh.getSize())
    cN.resize(1, cmesh.getSize(), VISIBLE);
  //std::set<int> uids;
  if (overwrite)
  {
    for (E_Int i=0; i < cellNzoom.getSize(); ++i){
      cN[subset[i]] = cellNzoom[i];
    }
  }
  else //append only blanked values
  {
    for (E_Int i=0; i < cellNzoom.getSize(); ++i){
      if (cellNzoom[i]==CELLNVAL)cN[subset[i]] = CELLNVAL;
    }
  }
  
  //cell_intersect
  for (size_t i=0; i < is_colliding.size(); ++i){
    if (is_colliding[i]){
     cN[subset[oldIds2[i]]] = CELLNVAL;
     //uids.insert(subset[oldIds2[i]]);
    }
  }
  
  //int minid = *std::min_element(ids.begin(), ids.end());
  //int maxid = *std::max_element(ids.begin(), ids.end());
  //std::cout << "ids stat : size/min/max : " << ids.size() << " " << minid << " " << maxid << std::endl;
  //std::cout << "cN size : " << cN.getSize() << std::endl;
  
  /*std::ostringstream o;
  o << "/home/slandier/tmp/mask_tree/mask_" << xcount++ << ".txt";
  std::ofstream of(o.str().c_str());
  for (std::set<int>::const_iterator it=uids.begin(); it != uids.end(); ++it)
    of <<*it << std::endl;
  of.close();*/
  
  return 0;
  
}

/// 
template <>
inline E_Int __do_the_cell_blanking<K_MESH::Polyhedron<UNKNOWN> >
(MaskEntity& maskE, const K_FLD::FldArrayF& fmesh, E_Int posx, E_Int posy, E_Int posz,
 const K_FLD::FldArrayI& cmesh,  E_Int CELLNVAL, E_Int overwrite, K_FLD::FldArrayI& cN)
{  
  using ELT_t = K_MESH::Polyhedron<UNKNOWN>;
  // Preconditioning
  //std::cout << "cell blanking PH" << std::endl;
  // Mask bbox
  const K_SEARCH::BBox3D* bb = maskE._main->getTree()->getGlobalBox();
  std::vector<E_Int> subset;
  K_FLD::ArrayAccessor<K_FLD::FldArrayF > acrd(fmesh, posx, posy, posz);
  
  using ngon_type = ngon_t<K_FLD::FldArrayI>;
    
  K_CONNECT::ZoneExtractor<K_FLD::FldArrayF, ngon_type> extractor;
  ngon_type ng(cmesh);
  K_FLD::ArrayAccessor< ngon_type > acnt(ng, -1);
  extractor.getInBox<3, ELT_t>(bb->minB, bb->maxB, E_EPSILON, acrd, acnt, subset);
    
  if (subset.empty()) return 0;
  
  //std::cout << "nb of extracted : " << subset.size() << std::endl;
  size_t sz1 = subset.size();
  ngon_t<K_FLD::FldArrayI> ngZoom;
  ngZoom.PGs = ng.PGs;
  //
  for (size_t i = 0; i < sz1; ++i)
    ngZoom.addPH(ng, subset[i]);
  ngZoom.PHs.updateFacets();

  // FIRST PASS : CELL IN
  typedef K_FLD::FldArrayF Coordinate_t;
  typedef K_FLD::FldArrayI Connectivity_t;

  Coordinate_t ZoomCentroids;
  ZoomCentroids.resize(ngZoom.PHs.size(), 3);
  E_Float G[3];
  
  K_FLD::ArrayAccessor<ngon_type> acvz(ngZoom, -1);
  
  ELT_t e;
  for (size_t i=0; i < sz1; ++i)
  {    
    acvz.getEntry(i, e);
    e.iso_barycenter<K_FLD::ArrayAccessor<Coordinate_t> >(acrd, G, 1);

    for (E_Int k=0; k < 3; ++k)ZoomCentroids(i, k+1) = G[k]; //hack because dynamic interface of FldArray is broken with _rake/_stride logic.
  }
 
  Connectivity_t cellNzoom(ZoomCentroids.getSize(), 1);
  cellNzoom.setAllValuesAt(1);
  E_Int err = maskE._main->blank(ZoomCentroids, 1,2,3, cellNzoom, CELLNVAL, overwrite);
  if (err) return err;
        
  // reduce the work set to those inside the zoom but not masked yet
  E_Int sz2=0;
  for (E_Int i = 0; i < cellNzoom.getSize(); ++i)
    if (cellNzoom[i] !=CELLNVAL) ++sz2;

  ngon_type ngZoom2;
  ngZoom2.PGs = ngZoom.PGs;
  std::vector<E_Int> oldIds2;
  oldIds2.reserve(sz2);     
  for (E_Int i = 0; i < cellNzoom.getSize(); ++i)
  {
    if (cellNzoom[i]!=CELLNVAL)
    {
      ngZoom2.addPH(ngZoom, i);
      oldIds2.push_back(i);
    }
  }
  
  ngZoom2.PHs.updateFacets();
  
  //std::cout << " remaining unmasked in zoom :" << k << std::endl;

  // SECOND PASS : CELL INTERSECT
  std::vector<E_Int> is_colliding;
  if (!oldIds2.empty())
  {
    using loc_t = NUGA::localizer<K_SEARCH::BbTree3D>;
    loc_t loc(*maskE._auxiliary->getTree(), maskE._auxiliary->tolerance());

    std::vector<E_Int> ids1;
    const K_FLD::ArrayAccessor<K_FLD::FldArrayF>& acrd1 = maskE._auxiliary->coords();
    const Connectivity_t& cnt1 = maskE._auxiliary->connect().array();
    
    const Coordinate_t& crd2 = fmesh;
    const ngon_type& ng2 = ngZoom2;
                      
    NUGA::COLLIDE::compute<Coordinate_t, Connectivity_t, ngon_type, loc_t, K_MESH::Triangle, ELT_t, 3>(acrd1, cnt1, loc, crd2, posx, posy, posz, ng2, ids1, is_colliding, maskE._auxiliary->tolerance());
  }
  
  // SET CELLN ////////////////////////////////////////////////////////////////
  if (cN.getSize() < cmesh.getSize())
    cN.resize(1, cmesh.getSize(), VISIBLE);
  //std::set<int> uids;
  if (overwrite)
  {
    for (E_Int i=0; i < cellNzoom.getSize(); ++i){
      cN[subset[i]] = cellNzoom[i];
    }
  }
  else //append only blanked values
  {
    for (E_Int i=0; i < cellNzoom.getSize(); ++i){
      if (cellNzoom[i]==CELLNVAL)cN[subset[i]] = CELLNVAL;
    }
  }
  
  //cell_intersect
  for (size_t i=0; i < is_colliding.size(); ++i){
    if (is_colliding[i]){
     cN[subset[oldIds2[i]]] = CELLNVAL;
     //uids.insert(subset[oldIds2[i]]);
    }
  }
  
  //int minid = *std::min_element(ids.begin(), ids.end());
  //int maxid = *std::max_element(ids.begin(), ids.end());
  //std::cout << "ids stat : size/min/max : " << ids.size() << " " << minid << " " << maxid << std::endl;
  //std::cout << "cN size : " << cN.getSize() << std::endl;
  
  /*std::ostringstream o;
  o << "/home/slandier/tmp/mask_tree/mask_" << xcount++ << ".txt";
  std::ofstream of(o.str().c_str());
  for (std::set<int>::const_iterator it=uids.begin(); it != uids.end(); ++it)
    of <<*it << std::endl;
  of.close();*/
  
  return 0;
  
}

///
template <typename ELT_t>
inline E_Int do_the_blanking
(E_Int blankingType, MaskEntity& maskE, const K_FLD::FldArrayF& fmesh, E_Int posx, E_Int posy, E_Int posz,
 const K_FLD::FldArrayI* cmesh, E_Int CELLNVAL, E_Int overwrite, K_FLD::FldArrayI& cN)
{
  // node_in or cell_in
  if (blankingType != 1)   
    return __do_the_node_blanking(maskE._main, fmesh, posx, posy, posz, CELLNVAL, overwrite, cN);
  
  return __do_the_cell_blanking<ELT_t>(maskE, fmesh, posx, posy, posz, *cmesh, CELLNVAL, overwrite, cN);
}
