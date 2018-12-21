/*    
    Copyright 2013-2019 Onera.

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

// #define FLAG_STEP

# include "intersector.h"
#include "Nuga/Boolean/NGON_BooleanOperator.h"
#include <memory>
//#include <sstream>
//#include <fstream>
//#include <iostream>
#ifdef FLAG_STEP
#include "chrono.h"
E_Int chrono::verbose=1;
#endif

using namespace std;
using namespace K_FLD;
using namespace NUGA;

//============================================================================
/* Blank cells defined in arrays by a Tetra mesh mask (as an input hook) */
//============================================================================
PyObject* K_INTERSECTOR::XcellN(PyObject* self, PyObject* args)
{
  //
  typedef FloatArray crd_t;
  typedef IntArray cnt_t;
  typedef NGON_BooleanOperator<crd_t, cnt_t> boolean_t;
  
  PyObject* bgm;
  PyObject* celln;
  PyObject* maskingMesh;
  PyObject* wall_pgl; // list of body walls polygons to set correctly the cellN (to 0) inside bodies.
  PyObject* ghost_pgl; // list of boundary polygons to extrude to avoid unecessary X computations.
  char *varString1, *varString2, *varString3, *eltType1, *eltType2, *eltType3;
  E_Int ni, nj, nk;
  crd_t *crd(0), *crdMask(0), *fCelln(0);
  cnt_t *cnt(0), *cntMask(0), *cCelln(0);
  
  if (!PyArg_ParseTuple(args, "OOOOO", &bgm, &celln, &maskingMesh, &wall_pgl, &ghost_pgl)) return NULL;
    
  /////////////////////////////////////////////////////////////////////////
  // Extraction des donnees
  
  E_Int res = K_ARRAY::getFromArray(bgm, varString1, crd, ni, nj, nk, cnt, eltType1);

  if (res != 2 || strcmp(eltType1, "NGON") != 0 )
  {
    PyErr_SetString(PyExc_ValueError,
		   "XcellN: the input mesh must be NGON.");
      return NULL;
  }
  
  std::unique_ptr<crd_t> afmesh(crd); // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<cnt_t> acmesh(cnt); // to avoid to call explicit delete at several places in the code.
  
  //std::cout << "bgm : " << crd->cols() << "/" << cnt->cols() << std::endl;
  //std::cout << "res : " << res << std::endl;
  
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                      "XcellN: input mesh must contain coordinates.");
    return NULL;
  }

  res = K_ARRAY::getFromArray(celln, varString2, fCelln, ni, nj, nk, cCelln, eltType2);
  
  std::unique_ptr<crd_t> af2(fCelln); // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<cnt_t> ac2(cCelln); // to avoid to call explicit delete at several places in the code.
   
  if (res == -1 || strcmp(eltType2, "NGON*") != 0 )
  {
    PyErr_SetString(PyExc_TypeError,
                      "XcellN: cellN must be specified at centers.");
    return NULL;
  }

  E_Int posc = K_ARRAY::isCellNatureField2Present(varString2);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                      "XcellN: celln variable not found for one structured array.");
    return NULL;
  }
  res = K_ARRAY::getFromArray(maskingMesh, varString3, crdMask, ni, nj, nk, cntMask, eltType3);

  std::unique_ptr<crd_t> af3(crdMask); // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<cnt_t> ac3(cntMask); // to avoid to call explicit delete at several places in the code.

  if (res != 2 || strcmp(eltType3, "NGON") != 0 )
  {
    PyErr_SetString(PyExc_ValueError,
		   "XcellN: the input mesh must be NGON.");
      return NULL;
  }

  boolean_t::eAggregation agglo_pol = boolean_t::CONVEX;
  assert (agglo_pol == boolean_t::CONVEX); // volume is now based on cvx_triangulation that require convex PGs
  boolean_t oper (*crd, *cnt, *crdMask, *cntMask, 0., agglo_pol);

  // For checking PG list consistency
  E_Int maxPGid=0;
  {
    ngon_t<K_FLD::IntArray> ng(*cntMask);
    maxPGid = ng.PGs.size();
    //std::cout << "mask PGs : " << maxPGid << std::endl;
  }

  // Passing the specified wall pgs to the boolean to ignore cells that fall inside bodies
  {

    FldArrayI* wall_ids=NULL;
    if (wall_pgl != Py_None)
      res = K_NUMPY::getFromNumpyArray(wall_pgl, wall_ids, true);

    std::unique_ptr<FldArrayI> pL(wall_ids); // to avoid to call explicit delete at several places in the code.
  
    //std::cout << "result for NUMPY (WALL) is : " << res << std::endl;
    if ((res == 1) && (wall_ids != NULL)  && (wall_ids->getSize() != 0))
    {
      E_Int nb_special_pgs = wall_ids->getSize();
      E_Int minid(INT_MAX), maxid(-1);
      std::vector<E_Int> pgsList(nb_special_pgs);
      for (E_Int i = 0; i < nb_special_pgs; ++i) 
      {
        pgsList[i]=(*wall_ids)[i]-1;
        //std::cout << pgsList[i] << std::endl;
        minid = std::min(minid, pgsList[i]);
        maxid = std::max(maxid, pgsList[i]);
      }

      if (nb_special_pgs)
      {
        if (minid >  maxPGid || maxid > maxPGid)
          printf("Warning: XcellN: the BCWall polygon list is invalid. Skipped...\n");
        else
          oper.passPGs(0, pgsList);
      }
    }
  }
  // Passing the specified ghost pgs to the boolean to avoid some X computations (connectNoMatch situation)
  {

    FldArrayI* ghost_ids=NULL;
    if (ghost_pgl != Py_None)
      res = K_NUMPY::getFromNumpyArray(ghost_pgl, ghost_ids, true);

    std::unique_ptr<FldArrayI> pL(ghost_ids); // to avoid to call explicit delete at several places in the code.
  
    //std::cout << "result for NUMPY (GHOST) is : " << res << std::endl;
    if ((res == 1) && (ghost_ids != NULL)  && (ghost_ids->getSize() != 0))
    {
      E_Int nb_special_pgs = ghost_ids->getSize();
      E_Int minid(INT_MAX), maxid(-1);
      std::vector<E_Int> pgsList(nb_special_pgs);
      for (E_Int i = 0; i < nb_special_pgs; ++i) 
      {
        pgsList[i]=(*ghost_ids)[i]-1;
        //std::cout << pgsList[i] << std::endl;
        minid = std::min(minid, pgsList[i]);
        maxid = std::max(maxid, pgsList[i]);
      }

      if (nb_special_pgs)
      {
        if (minid >  maxPGid || maxid > maxPGid)
          printf("Warning: XcellN: the BCWall polygon list is invalid. Skipped...\n");
        else
          oper.passPGs(1, pgsList);
      }
    }
  }

  
  E_Int sz = fCelln->getSize();
  std::vector<E_Float> cN(sz);
  for (E_Int i = 0; i < sz; ++i) cN[i]=(*fCelln)[i];
  
  E_Int err = oper.XcellN(cN);

  if (err)
  {
    PyErr_SetString(PyExc_TypeError, "XcellN: failed.");
    return NULL;
  }

  //assert (cN.size() == sz);
  
  FloatArray cellnout(1, sz);
  for (E_Int i = 0; i < sz; ++i) cellnout[i] = E_Float(cN[i]);
 
  return K_ARRAY::buildArray(cellnout, "cellN", *cCelln, -1, "NGON", true); 
}
