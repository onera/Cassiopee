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
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;
  
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f,
                                     ni, nj, nk, cn, eltType);

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

  // Build new connectivity and fields
  FldArrayI& cm = *(cn->getConnect(0));
  E_Int neltsp = cm.getSize();
  // Chaque pyra se decompose en 2 tetraedres
  E_Int nelts = 2*neltsp;
  E_Int npts = f->getSize(), api = f->getApi(), nfld = f->getNfld();

  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts,
                                       "TETRA", false, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);
  FldArrayI& cm2 = *(cn2->getConnect(0));

#pragma omp parallel default(shared)
  {
    E_Int cnt;
    E_Int squareDiag1, squareDiag2;
    E_Int i1, i2, i3, i4;

#pragma omp for
    for (E_Int elt = 0; elt < neltsp; elt++)
    {
      /* determination des indices de l'element */
      i1 = cm(elt,1)-1;
      i2 = cm(elt,2)-1;
      i3 = cm(elt,3)-1;
      i4 = cm(elt,4)-1;
      //i5 = cm(elt,5)-1;

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
        cm2(cnt,1) = cm(elt,1);
        cm2(cnt,2) = cm(elt,2);
        cm2(cnt,3) = cm(elt,3);
        cm2(cnt,4) = cm(elt,5);

        // t2: I1I3I4I5
        cnt = 2*elt+1;
        cm2(cnt,1) = cm(elt,1);
        cm2(cnt,2) = cm(elt,3);
        cm2(cnt,3) = cm(elt,4);
        cm2(cnt,4) = cm(elt,5);
      }
      else
      {
        // build tetras: I2I3I4I5, I2I4I1I5
        // t1: I1I2I3I5
        cnt = 2*elt;
        cm2(cnt,1) = cm(elt,2);
        cm2(cnt,2) = cm(elt,3);
        cm2(cnt,3) = cm(elt,4);
        cm2(cnt,4) = cm(elt,5);

        // t2: I1I5I3I6
        cnt = 2*elt+1;
        cm2(cnt,1) = cm(elt,2);
        cm2(cnt,2) = cm(elt,4);
        cm2(cnt,3) = cm(elt,1);
        cm2(cnt,4) = cm(elt,5);
      }
    }

    // Copy fields to f2
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* f2p = f2->begin(n);
#pragma omp for
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }

  RELEASESHAREDU(array, f, cn);
  RELEASESHAREDU(tpl, f2, cn2);
  return tpl;
}