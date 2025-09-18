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
#include <unordered_set>
#include <algorithm> // for std::copy

#include "kcore.h"
#include "converter.h"
using namespace K_FLD;

//=============================================================================
/* Adapt a face point list numpy to a vertex point list numpy for unstructured 
   arrays */
//=============================================================================
PyObject* K_CONVERTER::adaptBCFacePL2VertexPL(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* arrayFPL;
  if (!PYPARSETUPLE_(args, OO_, &array, &arrayFPL)) return NULL;

  // Check numpy (face point list)
  FldArrayI* facePL;
  E_Int res = K_NUMPY::getFromPointList(arrayFPL, facePL);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "adaptBCFacePL2VertexPL: facePL numpy is invalid.");
    return NULL;
  }
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  
  res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  
  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "adaptBCFacePL2VertexPL: for unstructured arrays only.");
    RELEASESHAREDN(arrayFPL, facePL);
    RELEASESHAREDS(array, f);
    return NULL;
  }
  else if (res == 2)
  { 
    PyObject* tpl = NULL;
    if (K_STRING::cmp(eltType, "NGON") == 0 || K_STRING::cmp(eltType, "NGON*") == 0)
    {
      tpl = adaptBCFacePL2VertexPL_NGON(cn, facePL);
    }
    else
    {
      //tpl = adaptBCFacePL2VertexPL_ME(cn, facePL);
      PyErr_SetString(PyExc_TypeError,
                      "adaptBCFacePL2VertexPL: not implemented yet for ME arrays.");
    }

    RELEASESHAREDN(arrayFPL, facePL);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "adaptBCFacePL2VertexPL: unrecognised type of array.");
    return NULL;
  }
}

PyObject* K_CONVERTER::adaptBCFacePL2VertexPL_NGON(FldArrayI* cn, FldArrayI* fpl)
{
  E_Int nv, npts, fidx;
  E_Int* ngon = cn->getNGon();
  E_Int* indPG = cn->getIndPG();
  E_Int nfaces = fpl->getSize();
  E_Int* fplp = fpl->begin();
  
  std::unordered_set<E_Int> vertexSet;
  
  // Loop over all faces of the face point list and fill vertex set
  for (E_Int i = 0; i < nfaces; i++)
  {
    fidx = fplp[i]-1;
    E_Int* face = cn->getFace(fidx, nv, ngon, indPG);
    for (E_Int p = 0; p < nv; p++) vertexSet.insert(face[p]);
  }
  npts = vertexSet.size();
  
  PyObject* tpl = K_NUMPY::buildNumpyArray(npts, 1, 1);
  E_Int* vertexPL = K_NUMPY::getNumpyPtrI(tpl);
  std::copy(vertexSet.begin(), vertexSet.end(), vertexPL);
  return tpl;
}

PyObject* K_CONVERTER::adaptBCFacePL2VertexPL_ME(FldArrayI* cn, FldArrayI* fpl)
{
  return NULL;
}
