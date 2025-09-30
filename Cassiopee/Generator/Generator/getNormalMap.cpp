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

// getNormalMap

# include "generator.h"

using namespace K_CONST;
using namespace K_FLD;

// ============================================================================
/* Return normals map of a surface array */
// ============================================================================
PyObject* K_GENERATOR::getNormalMapOfMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int posx, posy, posz;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res == 1 || res == 2)
  {
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      PyErr_SetString(PyExc_ValueError,
                      "getNormalMap: can't find coordinates in array.");
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }
    posx++; posy++; posz++;

    E_Int api = f->getApi();
    E_Int npts = f->getSize();

    if (res == 1) // cas structure
    {
      E_Int im1 = im-1;
      E_Int jm1 = jm-1;
      E_Int km1 = km-1;

      if ((im == 1 && jm == 1) ||
          (im == 1 && km == 1) ||
          (jm == 1 && km == 1))
      {
        PyErr_SetString(PyExc_TypeError,
                        "getNormalMap: a surface array is required.");
        RELEASESHAREDS(array, f); return NULL;
      }

      if (im == 1)
      {
        im = jm; jm = km; im1 = 1;
      }
      else if (jm == 1)
      {
        jm = km; jm1 = 1;
      }
      else if (km == 1)
      {
        km1 = 1;
      }
      else
      {
        PyErr_SetString(PyExc_TypeError,
                        "getNormalMap: a surface array is required.");
        RELEASESHAREDS(array, f); return NULL;
      }

      PyObject* tpl = K_ARRAY::buildArray3(3, "sx,sy,sz", im1, jm1, km1, api);
      FldArrayF* nsurf;
        K_ARRAY::getFromArray3(tpl, nsurf);
      K_METRIC::compNormStructSurf(
        im, jm, f->begin(posx), f->begin(posy), f->begin(posz),
        nsurf->begin(1), nsurf->begin(2), nsurf->begin(3)
      );
      RELEASESHAREDS(tpl, nsurf);
      RELEASESHAREDS(array, f);
      return tpl;
    }
    else // if (res == 2) // cas non structure
    {
      if (strcmp(eltType, "NGON") == 0)
      {
        // on verifie que le NGON est surfacique
        E_Int dim = cn->getDim();
        if (dim != 2)
        {
          PyErr_SetString(PyExc_TypeError,
                          "getNormalMap: NGON array must be a surface.");
          RELEASESHAREDU(array, f, cn); return NULL;
        }

        // Build array contenant la surface
        PyObject* tpl = K_ARRAY::buildArray3(3, "sx,sy,sz", npts,
                                             *cn, eltType, true, api, true);
        FldArrayF* nsurf;
        K_ARRAY::getFromArray3(tpl, nsurf);
        E_Int ierr = K_METRIC::compSurfNGon(
          f->begin(posx), f->begin(posy), f->begin(posz), *cn,
          nsurf->begin(1), nsurf->begin(2), nsurf->begin(3)
        );

        // sortie si une erreur a ete trouvee
        if (ierr == 1)
        {
          PyErr_SetString(PyExc_TypeError,
                          "getNormalMap: only valid for surface NGons.");
          RELEASESHAREDS(tpl, nsurf);
          RELEASESHAREDU(array, f, cn);
          return NULL;
        }
        RELEASESHAREDS(tpl, nsurf);
        RELEASESHAREDU(array, f, cn);
        return tpl;
        
      }
      else  // cas Elements basiques
      {
        E_Int nelts = cn->getSize(); // nb d'elements
        // elts doivent etre tri ou quad : 2D
        if (strcmp(eltType, "TRI") != 0 && strcmp(eltType, "QUAD") != 0)
        {
          PyErr_SetString(PyExc_TypeError,
                          "getNormalMap: elements must be TRI or QUAD.");
          RELEASESHAREDU(array, f, cn); return NULL;
        }

        PyObject* tpl = K_ARRAY::buildArray3(3, "sx,sy,sz", npts,
                                             *cn, eltType, true, api, true);
        FldArrayF* nsurf;
        K_ARRAY::getFromArray3(tpl, nsurf);

        FldArrayF surf(nelts, 1);

        K_METRIC::compSurfUnstruct(
          *cn, eltType,
          f->begin(posx), f->begin(posy), f->begin(posz),
          nsurf->begin(1), nsurf->begin(2), nsurf->begin(3), surf.begin());
        RELEASESHAREDS(tpl, nsurf);
        RELEASESHAREDU(array, f, cn);
        return tpl;
      }
    }
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "getNormalMap: unknown type of array.");
    return NULL;
  }
}
