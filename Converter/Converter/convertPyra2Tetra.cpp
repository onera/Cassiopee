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

// convertit un maillage pyramidal (pyra) en maillage tetraedrique
// voir article Dompierre et al 1999
// How to subdivide Pyramids, Prisms, and Hexahedra into Tetrahedra

# include "converter.h"
# include "kcore.h"
# include <string.h>
# include <stdio.h>

using namespace K_FLD;
using namespace std;
using namespace K_FUNC;

//=============================================================================
/* Conversion du maillage pyramidal en maillage tetraedrique. Chaque pyramide
   est decomposee en 2 tetraedres. Pour chaque pyramide, on determine les 2 
   diagonales de la face quadrilatere. La division en deux tetraedres se fait 
   suivant la plus petite de ces 2 diagonales. */
//=============================================================================
PyObject* K_CONVERTER::convertPyra2Tetra(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, 
                                    eltType, true);

  // Test non structure ?
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertPyra2Tetra: invalid array.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "convertPyra2Tetra: input array must be unstructured.");
    return NULL;
  }

  // Test pyra type ?
  if (strcmp(eltType, "PYRA") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "convertPyra2Tetra: unstructured array must be pentahedrical.");
    return NULL;
  }
  
  // tableau des coordonnees 
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertPyra2Tetra: coord must be present in array.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  // Chaque pyra se decompose en 2 tetraedres
  E_Int neltsp = cn->getSize();
  E_Int nelts = 2*neltsp;
  FldArrayI& cn0 = *cn;
  E_Int eltt = 4;//TETRA 
  PyObject* tpl = K_ARRAY::buildArray(f->getNfld(), varString, f->getSize(), nelts, eltt, NULL);
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(f->getSize(), f->getNfld(), fnp, true); fn = *f;
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  FldArrayI ct(nelts, 4, cnnp, true);

  E_Int* ct1 = ct.begin(1);
  E_Int* ct2 = ct.begin(2);
  E_Int* ct3 = ct.begin(3);
  E_Int* ct4 = ct.begin(4);
  E_Int* cn01 = cn0.begin(1);
  E_Int* cn02 = cn0.begin(2);
  E_Int* cn03 = cn0.begin(3);
  E_Int* cn04 = cn0.begin(4);
  E_Int* cn05 = cn0.begin(5);

#pragma omp parallel default(shared)
  {
  E_Int cnt;
  E_Int squareDiag1, squareDiag2;
  E_Int i1, i2, i3, i4;

#pragma omp for
  for (E_Int elt = 0; elt < neltsp; elt++)
  {
    /* determination des indices de l'element */
    i1 = cn01[elt]-1;
    i2 = cn02[elt]-1;
    i3 = cn03[elt]-1;
    i4 = cn04[elt]-1;
    //i5 = cn05[elt]-1;

    /* determination de la plus petite diagonale de la face quad i1i2i3i4 */
    // on retient les 2 sommets de la diag min
    squareDiag1 = ((x[i3]-x[i1])*(x[i3]-x[i1])+(y[i3]-y[i1])*(y[i3]-y[i1])+(z[i3]-z[i1])*(z[i3]-z[i1]));
    squareDiag2 = ((x[i4]-x[i2])*(x[i4]-x[i2])+(y[i4]-y[i2])*(y[i4]-y[i2])+(z[i4]-z[i2])*(z[i4]-z[i2]));

    /* construction des elements tetra */
    if (squareDiag1 <= squareDiag2)
    {
      // build tetras: I1I2I3I5,I1I3I4I5
      // t1: I1I2I3I5
      cnt = 2*elt;
      ct1[cnt] = cn01[elt];
      ct2[cnt] = cn02[elt];
      ct3[cnt] = cn03[elt];
      ct4[cnt] = cn05[elt];

      // t2: I1I3I4I5
      cnt = 2*elt+1;
      ct1[cnt] = cn01[elt];
      ct2[cnt] = cn03[elt];
      ct3[cnt] = cn04[elt];
      ct4[cnt] = cn05[elt];
    }
    else
    {
      // build tetras: I2I3I4I5, I2I4I1I5
      // t1: I1I2I3I5
      cnt = 2*elt;
      ct1[cnt] = cn02[elt];
      ct2[cnt] = cn03[elt];
      ct3[cnt] = cn04[elt];
      ct4[cnt] = cn05[elt];

      // t2: I1I5I3I6
      cnt = 2*elt+1;
      ct1[cnt] = cn02[elt];
      ct2[cnt] = cn04[elt];
      ct3[cnt] = cn01[elt];
      ct4[cnt] = cn05[elt];
      cnt++;
    }
  }
  }

  RELEASESHAREDU(array, f, cn);
  return tpl;
}
