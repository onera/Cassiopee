/*    
    Copyright 2013-2024 Onera.

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

// Elliptic generators

# include "generator.h"

using namespace std;
using namespace K_FLD;

extern "C"
{
  void k6coeffttm_(E_Float*x, E_Float*y, 
                   E_Int& n, E_Int& m, 
                   E_Float* a11, E_Float* a12, E_Float* a22, 
                   E_Float* b11, E_Float* b12, E_Float* b22,
                   E_Float* c11, E_Float* c12, E_Float* c22);
  void k6sor_(E_Float* x, E_Float* y, E_Int& n, E_Int& m,
              E_Float* a11, E_Float* a12, E_Float*a22,
              E_Float* b11, E_Float* b12, E_Float* b22,
              E_Float* c11, E_Float* c12, E_Float* c22,
              E_Float* rhsx, E_Float* rhsy, E_Float& err);
}

// ============================================================================
/* TTM */
// ============================================================================
PyObject* K_GENERATOR::TTMMesh(PyObject* self, PyObject* args)
{
  E_Int nit;
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &nit)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int posx, posy, posz;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, 
                                    im, jm, km, cn, eltType);
  
  // Check
  if (res != 1)
  {
    delete f;
    if (res == 2) delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "TTM: array must structured.");
    return NULL;
  }
  if (km != 1)
  {
    printf("Warning: TTM: applied only on k=1 zones.\n"); 
    delete f;
    return array;
  }
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "TTM: can't find coordinates in array.");
    return NULL;    
  }
  posx++; posy++; posz++;

  E_Int nj = jm;
  E_Int ni = im;
  
  FldArrayF* coord = new FldArrayF(*f);

  E_Int it = 0;
  FldArrayF a11(ni*nj); a11.setAllValuesAtNull();
  FldArrayF a12(ni*nj); a12.setAllValuesAtNull();
  FldArrayF a22(ni*nj); a22.setAllValuesAtNull();
  FldArrayF b11(ni*nj); b11.setAllValuesAtNull();
  FldArrayF b12(ni*nj); b12.setAllValuesAtNull();
  FldArrayF b22(ni*nj); b22.setAllValuesAtNull();
  FldArrayF c11(ni*nj); c11.setAllValuesAtNull();
  FldArrayF c12(ni*nj); c12.setAllValuesAtNull();
  FldArrayF c22(ni*nj); c22.setAllValuesAtNull();
  FldArrayF rhsx(ni*nj); rhsx.setAllValuesAtNull();
  FldArrayF rhsy(ni*nj); rhsy.setAllValuesAtNull();
  E_Float err = 1000.;

  while (it < nit && err > 1.e-13)
  {
    k6coeffttm_(coord->begin(posx), coord->begin(posy), 
                ni, nj, a11.begin(), a12.begin(), a22.begin(), 
                b11.begin(), b12.begin(), b22.begin(), 
                c11.begin(), c12.begin(), c22.begin());
    k6sor_(coord->begin(posx), coord->begin(posy),
           ni, nj, a11.begin(), a12.begin(), a22.begin(), 
           b11.begin(), b12.begin(), b22.begin(), 
           c11.begin(), c12.begin(), c22.begin(), 
           rhsx.begin(), rhsy.begin(), err);
    it++;
  }
 
  // build array
  delete f;
  PyObject* tpl = K_ARRAY::buildArray(*coord, varString, ni, nj, 1);
  delete coord;
  return tpl;
}
