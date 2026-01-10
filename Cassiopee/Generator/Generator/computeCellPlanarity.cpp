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

# include "generator.h"

using namespace K_FLD;
using namespace K_CONST;

//=============================================================================
/* 
   Calcul une distance mesurant la non-platitude des mailles.
   IN: un maillage surfacique
   OUT: cette mesure.
*/
//=============================================================================
PyObject*
K_GENERATOR::computeCellPlanarity( PyObject* self, PyObject* args )
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check arrays
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);

  // Check data
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeCellPlanarity: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    if (ni != 1 && nj != 1 && nk != 1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "computeCellPlanarity: array must be a surface array.");
      return NULL;
    }
  }
  if (res == 2)
  {
    if (strcmp(eltType, "TRI") != 0 && strcmp(eltType, "QUAD") != 0)
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "computeCellPlanarity: array must be a surface array.");
      return NULL;
    }
  }
  
  E_Int api = f->getApi();
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "computeCellPlanarity: coordinates must be present in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  if (res == 1)
  {
    E_Int ni1, nj1, nk1;
    FldArrayF* dist = new FldArrayF();
    cellPlanarityStructured(ni, nj, nk,
                            *f,
                            posx, posy, posz,
                            *dist, ni1, nj1, nk1);
    PyObject* tpl = K_ARRAY::buildArray3(*dist, "dist", ni1, nj1, nk1, api);
    delete dist; 
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else
  {
    FldArrayF* dist = new FldArrayF();
    FldArrayI* connect = new FldArrayI();
    (*connect) = (*cn);

    if (strcmp(eltType, "TRI") == 0)
    {
      dist->malloc(cn->getSize());
      dist->setAllValuesAtNull(); // les triangles sont forcement planaires
    }
    else
    {
      // QUADS
      cellPlanarityUnstructured(*f, *cn,
                                posx, posy, posz,
                                *dist);
    }
    PyObject* tpl = K_ARRAY::buildArray3(*dist, "dist", *connect, eltType, api);
    delete dist; delete connect;
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
}

//=============================================================================
void K_GENERATOR::cellPlanarityStructured(
  E_Int ni, E_Int nj, E_Int nk,
  FldArrayF& f,
  E_Int posx, E_Int posy, E_Int posz,
  FldArrayF& dist, E_Int& ni1, E_Int& nj1, E_Int& nk1)
{
  E_Float sigma0, sigma1;

  ni1 = K_FUNC::E_max(1, ni-1);
  nj1 = K_FUNC::E_max(1, nj-1);
  nk1 = K_FUNC::E_max(1, nk-1);

  dist.malloc(ni1*nj1*nk1);
  
  E_Int dir = 1;
  if (ni == 1) dir = 1;
  else if (nj == 1) dir = 2;
  else if (nk == 1) dir = 3;

  E_Int inc1 = 0;
  E_Int inc2 = 0;

  switch(dir)
  {
    case 1:
      inc1 = ni;
      inc2 = ni*nj;
      break;

    case 2:
      inc1 = 1;
      inc2 = ni*nj;
      break;

    case 3:
      inc1 = 1;
      inc2 = ni;
      break;
  }

  E_Int ind, ind1, ind2, ind3, ind4;
  E_Float xp, yp, zp, distLoc;
  E_Bool in;
  E_Float pQ[3];
  E_Float p1[3];
  E_Float p2[3];
  E_Float p3[3];
  E_Float* xt = f.begin(posx);
  E_Float* yt = f.begin(posy);
  E_Float* zt = f.begin(posz);

  for (E_Int k = 0; k < nk1; k++)
    for (E_Int j = 0; j < nj1; j++)
      for (E_Int i = 0; i < ni1; i++)
      {
        ind = i + j*ni1 + k*ni1*nj1;
        ind1 = i + j*ni + k*ni*nj;
        ind2 = i + inc1 + j*ni + k*ni*nj;
        ind3 = i + inc2 + j*ni + k*ni*nj;
        ind4 = i + inc1 + inc2 + j*ni + k*ni*nj;
        pQ[0] = ONE_FOURTH*(xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4]);
        pQ[1] = ONE_FOURTH*(yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4]); 
        pQ[2] = ONE_FOURTH*(zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4]); 

        dist[ind] = 0.;
        distLoc = 0.;
        p1[0] = xt[ind1];
        p1[1] = yt[ind1];
        p1[2] = zt[ind1];
        p2[0] = xt[ind2];
        p2[1] = yt[ind2];
        p2[2] = zt[ind2];
        p3[0] = xt[ind3];
        p3[1] = yt[ind3];
        p3[2] = zt[ind3];

        //pQ[0] = f(ind4, posx);
        //pQ[1] = f(ind4, posy);
        //pQ[2] = f(ind4, posz);

        K_COMPGEOM::distanceToTriangle(
          p1, p2, p3,
          pQ, 1,
          distLoc, in, 
          xp, yp, zp, sigma0, sigma1 );
        //if (in == true)
        dist[ind] = K_FUNC::E_max(distLoc, dist[ind]);
          //printf("  Q: %f %f %f %12.8g\n", pQ[0], pQ[1], pQ[2], distLoc);
          //printf("  1: %f %f %f %f \n", p1[0], p1[1], p1[2]);
          //printf("  2: %f %f %f %f \n", p2[0], p2[1], p2[2]);
          //printf("  3: %f %f %f %f \n", p3[0], p3[1], p3[2]);

        p1[0] = xt[ind4];
        p1[1] = yt[ind4];
        p1[2] = zt[ind4];

        //pQ[0] = f(ind1, posx);
        //pQ[1] = f(ind1, posy);
        //pQ[2] = f(ind1, posz);
        K_COMPGEOM::distanceToTriangle(
          p1, p2, p3,
          pQ, 1,
          distLoc, in, 
          xp, yp, zp, sigma0, sigma1);
        //if (in == true)
        dist[ind] = K_FUNC::E_max(distLoc, dist[ind]);
          //printf("  Q: %f %f %f %12.8g\n", pQ[0], pQ[1], pQ[2], distLoc);
          //printf("  1: %f %f %f %f \n", p1[0], p1[1], p1[2]);
          //printf("  2: %f %f %f %f \n", p2[0], p2[1], p2[2]);
          //printf("  3: %f %f %f %f \n", p3[0], p3[1], p3[2]);

        p1[0] = xt[ind1];
        p1[1] = yt[ind1];
        p1[2] = zt[ind1];
        p2[0] = xt[ind4];
        p2[1] = yt[ind4];
        p2[2] = zt[ind4];
        p3[0] = xt[ind3];
        p3[1] = yt[ind3];
        p3[2] = zt[ind3];
        //pQ[0] = f(ind2, posx);
        //pQ[1] = f(ind2, posy);
        //pQ[2] = f(ind2, posz);

        K_COMPGEOM::distanceToTriangle(
          p1, p2, p3,
          pQ, 1,
          distLoc, in, 
          xp, yp, zp, sigma0, sigma1 );
        //if (in == true)
        dist[ind] = K_FUNC::E_max(distLoc, dist[ind]);
          //printf("  Q: %f %f %f %12.8g\n", pQ[0], pQ[1], pQ[2], distLoc);
          //printf("  1: %f %f %f %f \n", p1[0], p1[1], p1[2]);
          //printf("  2: %f %f %f %f \n", p2[0], p2[1], p2[2]);
          //printf("  3: %f %f %f %f \n", p3[0], p3[1], p3[2]);

        p3[0] = xt[ind2];
        p3[1] = yt[ind2];
        p3[2] = zt[ind2];
        //pQ[0] = f(ind3, posx);
        //pQ[1] = f(ind3, posy);
        //pQ[2] = f(ind3, posz);

        K_COMPGEOM::distanceToTriangle(
          p1, p2, p3,
          pQ, 1,
          distLoc, in, 
          xp, yp, zp, sigma0, sigma1 );
        //if (in == true)
        dist[ind] = K_FUNC::E_max(distLoc, dist[ind]);
        //printf("  Q: %f %f %f %12.8g \n", pQ[0], pQ[1], pQ[2], dist[ind]);
        //printf("  1: %f %f %f %f \n", p1[0], p1[1], p1[2]);
        //printf("  2: %f %f %f %f \n", p2[0], p2[1], p2[2]);
        //printf("  3: %f %f %f %f \n", p3[0], p3[1], p3[2]);
      }
  dist.sqrt();
}

//=============================================================================
void K_GENERATOR::cellPlanarityUnstructured(
  FldArrayF& f, FldArrayI& cn,
  E_Int posx, E_Int posy, E_Int posz,
  FldArrayF& dist)
{
  E_Float sigma0, sigma1;

  // les elements sont supposes etre des QUAD
  dist.malloc(cn.getSize());
  E_Int ind, ind1, ind2, ind3, ind4;
  E_Float xp, yp, zp, distLoc;
  E_Bool in;
  E_Float pQ[3];
  E_Float p1[3];
  E_Float p2[3];
  E_Float p3[3];

  E_Int ne = cn.getSize();
  E_Float* xt = f.begin(posx);
  E_Float* yt = f.begin(posy);
  E_Float* zt = f.begin(posz);

  for (E_Int i = 0; i < ne; i++)
  {
    ind = i;
    ind1 = cn(i,1)-1;
    ind2 = cn(i,2)-1;
    ind3 = cn(i,3)-1;
    ind4 = cn(i,4)-1;

    pQ[0] = ONE_FOURTH*(xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4]);
    pQ[1] = ONE_FOURTH*(yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4]);
    pQ[2] = ONE_FOURTH*(zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4]);

    dist[ind] = 0;
    distLoc = 0.;
    p1[0] = xt[ind1];
    p1[1] = yt[ind1];
    p1[2] = zt[ind1];
    p2[0] = xt[ind2];
    p2[1] = yt[ind2];
    p2[2] = zt[ind2];
    p3[0] = xt[ind3];
    p3[1] = yt[ind3];
    p3[2] = zt[ind3];

    K_COMPGEOM::distanceToTriangle(
      p1, p2, p3,
      pQ, 1,
      distLoc, in, 
      xp, yp, zp, sigma0, sigma1 );
    //if (in == true)
    dist[ind] = K_FUNC::E_max(distLoc, dist[ind]);
      //printf("- %f %f %f %12.8g - \n", pQ[0], pQ[1], pQ[2], distLoc);
      //printf("  1- %f %f %f %f \n", p1[0], p1[1], p1[2]);
      //printf("  2- %f %f %f %f \n", p2[0], p2[1], p2[2]);
      //printf("  3- %f %f %f %f \n", p3[0], p3[1], p3[2]);

    p1[0] = xt[ind4];
    p1[1] = yt[ind4];
    p1[2] = zt[ind4];
    
    K_COMPGEOM::distanceToTriangle(
      p1, p2, p3,
      pQ, 1,
      distLoc, in, 
      xp, yp, zp, sigma0, sigma1 );
    //if (in == true)
    dist[ind] = K_FUNC::E_max(distLoc, dist[ind]);
      //printf("- %f %f %f %12.8g - \n", pQ[0], pQ[1], pQ[2], distLoc);
      //printf("  1- %f %f %f %f \n", p1[0], p1[1], p1[2]);
      //printf("  2- %f %f %f %f \n", p2[0], p2[1], p2[2]);
      //printf("  3- %f %f %f %f \n", p3[0], p3[1], p3[2]);

    p1[0] = xt[ind1];
    p1[1] = yt[ind1];
    p1[2] = zt[ind1];
    p2[0] = xt[ind4];
    p2[1] = yt[ind4];
    p2[2] = zt[ind4];
    p3[0] = xt[ind3];
    p3[1] = yt[ind3];
    p3[2] = zt[ind3];
        
    K_COMPGEOM::distanceToTriangle(
      p1, p2, p3,
      pQ, 1,
      distLoc, in, 
      xp, yp, zp, sigma0, sigma1 );
    //if (in == true)
      dist[ind] = K_FUNC::E_max(distLoc, dist[ind]);
      //printf("- %f %f %f %12.8g - \n", pQ[0], pQ[1], pQ[2], distLoc);
      //printf("  1- %f %f %f %f \n", p1[0], p1[1], p1[2]);
      //printf("  2- %f %f %f %f \n", p2[0], p2[1], p2[2]);
      //printf("  3- %f %f %f %f \n", p3[0], p3[1], p3[2]);

    p3[0] = xt[ind2];
    p3[1] = yt[ind2];
    p3[2] = zt[ind2];
    K_COMPGEOM::distanceToTriangle(
      p1, p2, p3,
      pQ, 1,
      distLoc, in, 
      xp, yp, zp, sigma0, sigma1 );
    //if (in == true)
      dist[ind] = K_FUNC::E_max(distLoc, dist[ind]);
      //printf("- %f %f %f %12.8g - \n", pQ[0], pQ[1], pQ[2], distLoc);
      //printf("  1- %f %f %f %f \n", p1[0], p1[1], p1[2]);
      //printf("  2- %f %f %f %f \n", p2[0], p2[1], p2[2]);
      //printf("  3- %f %f %f %f \n", p3[0], p3[1], p3[2]);
  }
  dist.sqrt();
}
