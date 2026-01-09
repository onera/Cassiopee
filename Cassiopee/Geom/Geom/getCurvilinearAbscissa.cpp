/*
    Copyright 2013-2026 ONERA.

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

# include "geom.h"
using namespace K_FUNC;
using namespace K_FLD;
using namespace std;
using namespace K_CONST;

// ============================================================================
/* Return the curvilinear abscissa of a 1D array defining a mesh */
// ============================================================================
PyObject* K_GEOM::getCurvilinearAbscissa(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  E_Int api = f->getApi();
  E_Int posx, posy, posz;
  E_Float length = 0.;
  E_Float l, dx, dy, dz;
  PyObject* tpl = NULL;
  FldArrayF* ab;

  if (res == 1)
  {
    if (jm != 1 || km != 1)
      printf("Warning: getCurvilinearAbscissa: only line j=1, k=1 is taken into account.\n");

    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "getCurvilinearAbscissa: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;


    tpl = K_ARRAY::buildArray3(1, "s", im, jm, km, api);
    K_ARRAY::getFromArray3(tpl, ab);
    E_Float* abp = ab->begin();

    abp[0] = 0.;
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);

    for (E_Int i = 1; i < im; i++)
    {
      dx = xt[i] - xt[i-1];
      dy = yt[i] - yt[i-1];
      dz = zt[i] - zt[i-1];

      l = sqrt(dx*dx + dy*dy + dz*dz);
      length += l;
      abp[i] = abp[i-1] + l;
    }

    E_Float inv = 1./length;
    for (E_Int i = 1; i < im; i++) abp[i] *= inv;

    RELEASESHAREDS(tpl, ab);
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2)
  {
    if (strcmp(eltType, "BAR") != 0)
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "getCurvilinearAbscissa: only for BAR-array.");
      return NULL;
    }
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "getCurvilinearAbscissa: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int npts = f->getSize();
    tpl = K_ARRAY::buildArray3(1, "s", npts, *cn, eltType, false, api, true);
    K_ARRAY::getFromArray3(tpl, ab);
    E_Float* abp = ab->begin();

    abp[0] = 0.;
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    
    // pt de depart : 0
    E_Int nelts = cn->getSize();
    E_Int ind2;
    FldArrayIS dejaVu(nelts); dejaVu.setAllValuesAtNull();
    short* dejaVup = dejaVu.begin();
    for (E_Int ind1 = 0; ind1 < npts-1; ind1++)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        if (dejaVup[et] == 0 && (*cn)(et,1)-1 == ind1)
        {
          ind2 = (*cn)(et,2)-1;
          dx = xt[ind2] - xt[ind1];
          dy = yt[ind2] - yt[ind1];
          dz = zt[ind2] - zt[ind1];

          l = sqrt(dx*dx + dy*dy + dz*dz);
          length += l;
          abp[ind2] = abp[ind1] + l;
          dejaVup[et] = 1;
          break;
        }
      }
    }

    E_Float inv = 1./length;
    for (E_Int i = 1; i < npts; i++) abp[i] *= inv;

    RELEASESHAREDS(tpl, ab);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "getCurvilinearAbscissa: invalid array.");
    return NULL;
  }
}
