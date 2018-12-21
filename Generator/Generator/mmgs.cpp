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

// Remaillage surfacique avec mmgs

#include "generator.h"
#include "MMGS/mmgs.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC; 

// ============================================================================
/* MMGS
   IN: maillage TRI
   IN: eventuellement metric ou solution
   IN: 
   IN: 
   OUT: maillage TRI remaille. */
// ============================================================================
PyObject* K_GENERATOR::mmgs(PyObject* self, PyObject* args)
{
  E_Float ridgeAngle=30.; // angle detection
  E_Int anisotropy = 0;
  E_Float hmin = 0.;
  E_Float hmax = 0.;
  E_Float hausd = 0.;
  E_Float hgrad = -1.;
  PyObject* array;
  if (!PYPARSETUPLE(args, 
                    "Odddddl", "Odddddi", 
                    "Offfffl", "Offfffi", 
                    &array, &ridgeAngle, &hmin, &hmax, &hausd, &hgrad, &anisotropy))
  {
    return NULL;
  }

  // Check data

  
  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
                                     cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "mmgs: invalid array.");
    return NULL;
  }

  /* update Info */
  MMG5_Info info;
  info.imprim = -99;
  info.ddebug = 0;
  info.mem = -1;
  
  /* ridge angle */
  info.dhd = ridgeAngle;
  info.dhd = max(0.0, min(180.0,info.dhd));
  info.dhd = cos(info.dhd*3.14159265359/180.0);
  
  info.hmax = hmax;
  info.hmin = hmin;
  info.hausd = hausd;
  if (hgrad < 0.0) hgrad = -1;
  else hgrad = log(hgrad);
  info.hgrad = hgrad;

  /* load mesh */
  
  
  /* load met */

  /* Analysis */
  //_MMG5_analys(&mesh);
  
  /* Main call */
  //_MMG5_mmgs1(MMG5_pMesh mesh, MMG5_pSol met);

  /* Export */
  //unscaleMesh(&mesh,&met)

  RELEASESHAREDB(res, array, f, cn);
  return NULL;
}
