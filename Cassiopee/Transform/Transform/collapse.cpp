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

# include "transform.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
// Transforme un tri en bar: recolle les sommets de distance
// minimale ou de hauteur minimale
//=============================================================================
PyObject* K_TRANSFORM::collapse(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString,
                                     f, im, jm, km, cn, eltType);
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "collapse: array must be unstructured.");
    if (res == 1) delete f;
    return NULL;
  }
  if (K_STRING::cmp(eltType, "TRI") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "collapse: array must be TRI.");
    delete f; delete cn; return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "collapse: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);
  E_Int nelts = cn->getSize();
  E_Int api = f->getApi();

  // Construction de la BAR ou TRI
  FldArrayI* cn2 = new FldArrayI();
  E_Float eps = 1.e-12;
  char eltType2[256];
  if (K_STRING::cmp(eltType, "TRI") == 0) //TRI->BAR
  {
    collapseMinVertexInTriangle(*cn, xt, yt, zt);
    cn2->malloc(nelts, 2); strcpy(eltType2, "BAR");
    // elimination des elements degeneres dans la connectivite TRI
    for (E_Int et = 0; et < nelts; et++)
    {
      E_Int ind1 = (*cn)(et,1) - 1;
      E_Int ind2 = (*cn)(et,2) - 1;
      E_Int ind3 = (*cn)(et,3) - 1;
      E_Float dx = xt[ind1] - xt[ind2];
      E_Float dy = yt[ind1] - yt[ind2];
      E_Float dz = zt[ind1] - zt[ind2];
      E_Float dist = dx*dx + dy*dy + dz*dz;
      if (dist <= eps*eps) // doublon
      {
        (*cn2)(et,1) = ind1+1; (*cn2)(et,2) = ind3+1;
        goto nexttri;
      }
      dx = xt[ind1] - xt[ind3];
      dy = yt[ind1] - yt[ind3];
      dz = zt[ind1] - zt[ind3];
      dist = dx*dx + dy*dy + dz*dz;
      if (dist <= eps*eps) // doublon
      {
        (*cn2)(et,1) = ind1+1; (*cn2)(et,2) = ind2+1;
        goto nexttri;
      }
      dx = xt[ind3] - xt[ind2];
      dy = yt[ind3] - yt[ind2];
      dz = zt[ind3] - zt[ind2];
      dist = dx*dx + dy*dy + dz*dz;
      if ( dist <= eps*eps ) // doublon
      {
        (*cn2)(et,1) = ind3+1; (*cn2)(et,2) = ind1+1;
        goto nexttri;
      }
      PyErr_SetString(PyExc_TypeError,
                      "collapse: triangle not degenerated.");
      delete f; delete cn; delete cn2;
      nexttri:;
    }
  }

  K_CONNECT::cleanConnectivity(posx, posy, posz, eps, eltType2, *f, *cn2);
  PyObject* tpl = K_ARRAY::buildArray3(*f, varString, *cn2, eltType2, api);
  RELEASESHAREDU(array, f, cn);
  delete cn2;
  return tpl;
}

//=============================================================================
/* Collapse un vertex du triangle sur un autre
   On calcul les longeurs de chaque edge et la hauteur de chaque point.
   Le point de hauteur minimum est supprime si sa hauteur est inferieure
   a la longueur des edges. Sinon, l'edge le plus petit est contracte. */
//=============================================================================
void
K_TRANSFORM::collapseMinVertexInTriangle(FldArrayI& cn,
                                         E_Float* xt, E_Float* yt, E_Float* zt)
{
  E_Int nov1, nov2, nov3;
  E_Float dx, dy, dz;
  E_Int vert1, vert2;
  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2); E_Int* cn3 = cn.begin(3);
  E_Float p1[3], p2[3], p3[3];
  E_Float xp, yp, zp;
  E_Bool in;
  E_Float h1, h2, h3, d12, d13, d23, h, d;
  E_Int ind1, ind2;

  for ( E_Int et = 0; et < cn.getSize(); et++)
  {
    nov1 = cn1[et]-1; nov2 = cn2[et]-1; nov3 = cn3[et]-1;

    // Calcul des hauteurs pour chaque point
    p1[0] = xt[nov1]; p1[1] = yt[nov1]; p1[2] = zt[nov1];
    p2[0] = xt[nov2]; p2[1] = yt[nov2]; p2[2] = zt[nov2];
    p3[0] = xt[nov3]; p3[1] = yt[nov3]; p3[2] = zt[nov3];

    K_COMPGEOM::distanceToBar(p2, p3, p1, 0, xp, yp, zp, in, h1);
    K_COMPGEOM::distanceToBar(p1, p3, p2, 0, xp, yp, zp, in, h2);
    K_COMPGEOM::distanceToBar(p1, p3, p3, 0, xp, yp, zp, in, h3);

    // Distances des edges
    dx = xt[nov1]-xt[nov2]; dy = yt[nov1]-yt[nov2]; dz = zt[nov1]-zt[nov2];
    d12 = dx*dx + dy*dy + dz*dz;

    dx = xt[nov1]-xt[nov3]; dy = yt[nov1]-yt[nov3]; dz = zt[nov1]-zt[nov3];
    d13 = dx*dx + dy*dy + dz*dz;

    dx = xt[nov2]-xt[nov3]; dy = yt[nov2]-yt[nov3]; dz = zt[nov2]-zt[nov3];
    d23 = dx*dx + dy*dy + dz*dz;

    h = K_FUNC::E_min(h1, h2, h3);
    d = K_FUNC::E_min(d12, d13, d23);
    if (h < d)
    {
      if (h == h1)
      {
        if (d12 < d13)
        {vert1 = 1; vert2 = 2;}
        else
        {vert1 = 1; vert2 = 3;}
      }
      else if (h == h2)
      {
        if (d12 < d23)
        {vert1 = 2; vert2 = 1;}
        else
        {vert1 = 2; vert2 = 3;}
      }
      else // h == h3
      {
        if (d13 < d23)
        {vert1 = 3; vert2 = 1;}
        else
        {vert1 = 3; vert2 = 2;}
      }
    }
    else
    {
      if (d == d12) {vert1 = 1; vert2 = 2;}
      else if (d == d13) {vert1 = 1; vert2 = 3;}
      else {vert1 = 2; vert2 = 3;}
    }

    ind1 = cn(et, vert1)-1; ind2 = cn(et, vert2)-1;
    for (E_Int vert = 1; vert <= 3; vert++)
    {
      for (E_Int et2 = 0; et2 < cn.getSize(); et2++)
        if ( cn(et2,vert) == ind2+1 ) cn(et2,vert) = ind1+1;
    }
  }
}
