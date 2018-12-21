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

#include "post.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Enforce the indicator field to 0 for the finest level */
//=============================================================================
PyObject* K_POST::enforceIndicatorForFinestLevel(PyObject* self, 
                                                 PyObject* args)
{
  PyObject *indicator, *octree;
  if (!PyArg_ParseTuple(args, "OO", &indicator, &octree)) return NULL;
  // Verif octree HEXA/QUAD
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(octree, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "enforceIndicatorForFinestLevel: array must be unstructured.");
    RELEASESHAREDB(res, octree, f, cn); return NULL;   
  }
  //E_Int dim = 3;
  //if (strcmp(eltType,"HEXA") == 0) dim = 3;
  //else if (strcmp(eltType,"QUAD") == 0) dim = 2;
  //else
  //{
  //  PyErr_SetString(PyExc_TypeError, 
  //                  "enforceIndicatorForFinestLevel: array must be HEXA or QUAD.");
  //  RELEASESHAREDU(octree, f, cn); return NULL;
  //}
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforceIndicatorForFinestLevel: coordinates not found in array.");
    RELEASESHAREDU(octree, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  // Verif indicator
  E_Int nii, nji, nki;
  FldArrayF* fi; FldArrayI* cni;
  char* varStringi; char* eltTypei;
  E_Int resi = K_ARRAY::getFromArray(indicator, varStringi, fi, 
                                     nii, nji, nki, cni, eltTypei, true);
  if (resi != 1 && resi != 2) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforceIndicatorForFinestLevel: indic array must be structured.");
    RELEASESHAREDU(octree, f, cn); return NULL;
  }
  E_Int posi = K_ARRAY::isNamePresent("indicator", varStringi);
  if (posi == -1) 
  {
    RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: enforceIndicatorForFinestLevel: no refinement indicator given. Nothing done."); 
    return indicator;
  }
  posi++;
  if (fi->getSize() != cn->getSize()) 
  {
    RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: enforceIndicatorForFinestLevel: refinement indicator size must be equal to the number of elements. Nothing done."); 
    return indicator;
  }
  /*-----------FIN DES VERIFS ------------------*/
  E_Float* xt = f->begin(posx);
  E_Float* indict = fi->begin(posi);
  E_Int nelts = cn->getSize();
  E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
  E_Float dhmin = K_CONST::E_MAX_FLOAT;
  FldArrayF dht(nelts);
  E_Float* dhtp = dht.begin();
  E_Int ind1, ind2;
  for (E_Int et = 0; et < nelts; et++)
  {
    ind1 = cn1[et]-1; ind2 = cn2[et]-1;
    dhtp[et] = K_FUNC::E_abs(xt[ind2]-xt[ind1]);   
    dhmin = K_FUNC::E_min(dhtp[et], dhmin);
  }
  E_Float eps = 1.e-10;
  for (E_Int et = 0; et < nelts; et++)
  {
    if (K_FUNC::fEqualZero(dhtp[et]-dhmin, eps) == true) 
      indict[et] = -2000.;
  }
  /*-----------CONSTRUCTION ARRAY DE SORTIE ------------------*/
  PyObject* tpl;
  if (resi == 1) 
    tpl = K_ARRAY::buildArray(*fi, varStringi, nii, nji, nki);
  else 
    tpl = K_ARRAY::buildArray(*fi, varStringi, *cni, -1, eltTypei, 
                              false);
  RELEASESHAREDB(resi, indicator, fi, cni);
  return tpl;
}
