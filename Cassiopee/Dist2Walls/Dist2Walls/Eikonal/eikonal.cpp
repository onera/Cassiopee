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

# include "dist2walls.h"
# include "eikonalSolver.h"
# include "eikonalFMMSolver.h"
# include "eikonalFIMSolver.h"
# include <iostream>
# include <ctime>

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Eikonal solver (array) 
   IN: cartesian grid */
// ============================================================================
PyObject* K_DIST2WALLS::eikonal(PyObject* self, PyObject* args)
{
  PyObject* array;
  //struct timespec beg, end;
  //clock_gettime(CLOCK_REALTIME, &beg);
  int algo = 0;
  if (!PyArg_ParseTuple(args, "O|i", &array, &algo)) return NULL;
  
  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl, 
                               cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "eikonal: invalid array.");
    return NULL;
  }
  if (res == 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "eikonal: only for structured grids.");
    RELEASESHAREDB(res, array, f, cn);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "eikonal: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
    
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();

  // Construit l'array resultat et l'initialise par copie
  PyObject* tpl;
  tpl = K_ARRAY::buildArray(nfld, varString, nil, njl, nkl);
  
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(npts, nfld, fnp, true);
  fn.setAllValuesAt(*f);

  // Get the pointer on phi
  E_Int pos = K_ARRAY::isNamePresent("Phi", varString);
  if (pos == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "eikonal: cannot find Phi in array.");
    return NULL;
  }
  E_Float* phi = fn.begin(pos+1);
  E_Float max_float = 0.;
  for ( int i = 0; i < nil*njl*nkl; ++i )
    max_float = std::max(max_float,phi[i]);

  // Get the pointer on v (speed)
  E_Int posv = K_ARRAY::isNamePresent("speed", varString);
  if (posv == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "eikonal: cannot find front speed in array.");
    return NULL;
  }

  E_Float* v = fn.begin(posv+1);
  E_Float dh = x[1]-x[0];
  //E_Int nbSubIter = 5;

  //clock_gettime(CLOCK_REALTIME, &end);
  //double seconds = (double)((end.tv_sec+end.tv_nsec*1.E-9) - (beg.tv_sec+beg.tv_nsec*1.E-9));
  //std::cout << "Temps passé en C avant appel Eikonal solver : " << seconds << "secondes" << std::endl;  
  //E_Int nt = __NUMTHREADS__;  
  //nt = 0; // pas de multithread pour l'instant
  //clock_gettime(CLOCK_REALTIME, &beg);
  if (algo == 0 ) // Algorithme d'origine FMM
  {
    Eikonal::FMM::solveOnIsotropGrid( nil, njl, nkl, x[0], y[0], z[0], dh, phi, v);
  }
  if (algo == 1 ) // Variant de l'algorithme d'origine : FIM
  {
    Eikonal::FIM::solveOnIsotropGrid( nil, njl, nkl, x[0], y[0], z[0], dh, phi, v, max_float);
  }
  if (algo == 2 )// Si on a choisit l'algorithme FIM
  //if (nt == 0) 
  {
    solveEikonalOnIsotropGrid(nil, njl, nkl,
                              x[0], y[0], z[0],
                              dh, v, phi);
  }
  // else 
  // {
  //   blockFIM( nil, njl, nkl, x[0], y[0], z[0], dh, niBlk, njBlk, nkBlk, nbSubIter, v, phi);
  // }
  //clock_gettime(CLOCK_REALTIME, &end);
  //seconds = (double)((end.tv_sec+end.tv_nsec*1.E-9) - (beg.tv_sec+beg.tv_nsec*1.E-9));
  //std::cout << "Temps passé pour Eikonal solver : " << seconds << "secondes" << std::endl;  
  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}
