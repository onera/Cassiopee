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

# include "transform.h"
using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;

//=============================================================================
/* Déclarations et Fonctions */
//=============================================================================

struct Point
{
  E_Float x;
  E_Float y;
  E_Float z;
};

inline Point vecSub(Point v1, Point v2)
{
  return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

inline Point vecScale(Point v, E_Float scalaire)
{
  return {v.x * scalaire, v.y * scalaire, v.z * scalaire};
}

inline Point vecAdd(Point v1, Point v2)
{
  return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

inline E_Float vecDot(Point v1, Point v2)
{
  return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

inline Point vecRot90(Point v) /* {x}***/
{
  return {-v.y, v.x, 0.0};
}

inline Point vecCross (Point v1, Point v2)
{
  return {v1.y * v2.z - v1.z * v2.y, v1.z * v2.x- v1.x * v2.z, v1.x * v2.y - v1.y * v2.x} ;
}

inline bool isSamePoint(Point p1, Point p2) 
{
  return (std::abs(p1.x - p2.x) < 1e-10 && std::abs(p1.y - p2.y) < 1e-10);
}

// ============================================================================
/* Conservative smoothing */
// ============================================================================
PyObject* K_TRANSFORM::consSmooth(PyObject* self, PyObject* args)
{
  PyObject* array; // E_Int sweeps;

  E_Int sweeps = 1;
  
  if (!PYPARSETUPLE_(args, O_ I_, &array, &sweeps))
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =  K_ARRAY::getFromArray3(array, varString,
                                      f, im, jm, km, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: array is invalid.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // Build output
  PyObject* tpl = K_ARRAY::buildArray3(*f, varString, im, jm, km, f->getApi());
  FldArrayF* fo;
  // K_ARRAY::getFromArray3(array, fo);
  K_ARRAY::getFromArray3(tpl, fo);

  E_Float* a = fo->begin(posx);
  E_Float* b = fo->begin(posy);
  E_Float* c = fo->begin(posz);

  E_Int npts = im;
  E_Int nUnique = npts-1;

  Point P0 = {a[0], b[0], c[0]};
  Point Pn = {a[nUnique], b[nUnique], c[nUnique]};

  E_Int ouvert;
  E_Int start;
  E_Int end;

  if (isSamePoint(P0, Pn))
  {
    ouvert = 0;
    // printf ("La surface est fermée");
    start = 0;
    end = nUnique;
  }
  else 
  {
    ouvert = 1;
    //printf ("La surface est ouverte ");
    start = 1;
    end = npts-4;
  }
  
  for (E_Int k = 0; k < sweeps; k++)
  {
    for (E_Int i = start; i < end; i++)
    {
      /* On repart à 0 si on dépasse la fin du tableau */
      E_Int idx0 = i % (nUnique);
      E_Int idx1 = (i + 1) % (nUnique);
      E_Int idx2 = (i + 2) % (nUnique);
      E_Int idx3 = (i + 3) % (nUnique);

      /* On récupère les coordonnées des points i, i+1, i+2, i+3 */
      Point Pi   = {a[idx0], b[idx0], c[idx0]};
      Point Pip1 = {a[idx1], b[idx1], c[idx1]};
      Point Pip2 = {a[idx2], b[idx2], c[idx2]};
      Point Pip3 = {a[idx3], b[idx3], c[idx3]};

      /* On calcule les différences de points qui nous interessent (i+3 - i), (i+2 - i) et (i+1 - i) */
      Point dv3 = vecSub(Pip3, Pi); /* xi+3 - xi */
      Point dv2 = vecSub(Pip2, Pi);
      Point dv1 = vecSub(Pip1, Pi);

      /* On définit uNormal = unit normal to baseline (i+3;i) */
      E_Float normeDv3 = vecDot(dv3, dv3); /*  ||{xi+3-xi}**||² = ||xi+3-xi||² */
      if (normeDv3 < 1e-12) continue;
      E_Float divNorme = 1.0 / normeDv3 ;/*  1 / ||xi+3-xi||² */

      Point uNormal = vecScale(vecRot90 (dv3), divNorme); /* {xi+3-xi}** / ||xi+3-xi||²*/

      /* Calcul de l'aire signée */
      E_Float aire = 0.5 * vecDot( vecCross (dv3,dv2), {0.0, 0.0, 1.0}) + \
                    0.5 * vecDot( vecCross (dv2,dv1), {0.0, 0.0, 1.0}) ;

      E_Float h = 1.5 * aire;

      Point newPip1 = vecAdd(vecAdd(vecScale(Pi, 2.0/3.0), vecScale(Pip3, 1.0/3.0)), vecScale(uNormal, h));
      Point newPip2 = vecAdd(vecAdd(vecScale(Pi, 1.0/3.0), vecScale(Pip3, 2.0/3.0)), vecScale(uNormal, h));

      a[idx1] = newPip1.x; b[idx1] = newPip1.y; c[idx1] = newPip1.z;
      a[idx2] = newPip2.x; b[idx2] = newPip2.y; c[idx2] = newPip2.z;

      if (ouvert == 0)
      {
        if (idx1 == 0)
        {
          a[nUnique] = newPip1.x; b[nUnique] = newPip1.y; c[nUnique] = newPip1.z;
        }
        if (idx2 == 0 )
        {
          a[nUnique] = newPip2.x; b[nUnique] = newPip2.y; c[nUnique] = newPip2.z;
        }
      }



    }
  
  }

  RELEASESHAREDB(res, array, f, cn); 
  RELEASESHAREDS(tpl, fo); 
  return tpl;

}

