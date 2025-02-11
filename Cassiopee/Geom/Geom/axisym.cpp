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

#include "geom.h"

using namespace K_FLD;

extern "C"
{
  void k6axisym_(const E_Int& npts,
                 const E_Float* x, const E_Float* y, const E_Float* z,
                 const E_Float& xc, const E_Float& yc, const E_Float& zc,
                 const E_Float& nx, const E_Float& ny, const E_Float& nz,
                 const E_Float& teta, const E_Float& rmod,
                 E_Float* xo, E_Float* yo, E_Float* zo);
}

//=========================================================================
/* Generation de maillage axisymetrique
   maillage structure (ni,nj) -> (ni,nj,nk)
   maillage non structure TRI -> maillage non structure PENTA
   maillage non structure QUAD -> maillage non structure HEXA */
//=========================================================================
PyObject* K_GEOM::axisym(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float xc, yc, zc, nx, ny, nz, teta;
  E_Int nteta;
  PyObject* arrayR;

  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ R_ I_ O_,
                    &array, &xc, &yc, &zc, &nx, &ny, &nz, &teta, &nteta, &arrayR))
  {
      return NULL;
  }
  // check axis vector
  if (K_FUNC::fEqualZero(nx) == true && 
      K_FUNC::fEqualZero(ny) == true &&
      K_FUNC::fEqualZero(nz) == true)
  {
    PyErr_SetString(PyExc_ValueError,
                    "axisym: vector has null norm.");
    return NULL; 
  }
  // check teta
  if (K_FUNC::fEqualZero(teta) == true)
  {
    PyErr_SetString(PyExc_ValueError,
                    "axisym: teta must not be null.");
    return NULL;
  }
  if (nteta < 2)
  {
    PyErr_SetString(PyExc_ValueError,
                    "axisym: number of points in the azimuthal direction must be greater than 1.");
    return NULL;  
  }
  // check arrayR : si existant remplace teta et nteta
  // Cet array servira a definir un coefficient pour r en fonction de teta
  // Il n'est pas encore utilise. Il faudrait calculer r = f(teta) en fonction de
  // arrayR et de l'axe de symetrie.
  E_Int useR = 0;
  FldArrayF* fR; FldArrayI* cnR;
  E_Float* xr=NULL, *yr=NULL, *zr=NULL;
  E_Float xrC=0., yrC=0., zrC=0.;
  E_Float om0x=0., om0y=0., om0z=0., om1x=0., om1y=0., om1z=0.;
  if (arrayR != Py_None)
  {
    E_Int niR, njR, nkR;
    char* varStringR; char* eltTypeR;    
    useR = K_ARRAY::getFromArray(arrayR, varStringR, fR, 
                                 niR, njR, nkR, cnR, eltTypeR, true);

    if (useR != 1) { useR = 0; goto next; }

    // barycentre de la courbe
    E_Int posxr = K_ARRAY::isCoordinateXPresent(varStringR);
    E_Int posyr = K_ARRAY::isCoordinateYPresent(varStringR);
    E_Int poszr = K_ARRAY::isCoordinateZPresent(varStringR);
    if (posxr == -1 || posyr == -1 || poszr == -1)
    {
      useR = 0; goto next;
    }
    else
    {
      xr = fR->begin(posxr+1);
      yr = fR->begin(posyr+1);
      zr = fR->begin(poszr+1);
      // On ne prend pas le dernier (structure bouclant)
      for (E_Int i = 0; i < niR*njR*nkR-1; i++) 
      {
        xrC += xr[i]; yrC += yr[i]; zrC += zr[i];
      }
      xrC = xrC / (niR*njR*nkR);
      yrC = yrC / (niR*njR*nkR);
      zrC = zrC / (niR*njR*nkR);
      E_Float vx, vy, vz, v;
      E_Int nteta0;
      nteta = niR*njR*nkR;
      om0x = xr[0]-xrC; om0y = yr[0]-yrC; om0z = zr[0]-zrC;
      nteta0 = nteta/8;
      if (nteta0 == 0) nteta0 = 1;
      om1x = xr[nteta0]-xrC;
      om1y = yr[nteta0]-yrC;
      om1z = zr[nteta0]-zrC;
      vx = om0y*om1z-om0z*om1y;
      vy = om0z*om1x-om0x*om1z;
      vz = om0x*om1y-om0y*om1x;
      om1x = -(om0y*vz-om0z*vy);
      om1y = -(om0z*vx-om0x*vz);
      om1z = -(om0x*vy-om0y*vx);
      v = sqrt(om0x*om0x+om0y*om0y+om0z*om0z);
      om0x = om0x / v; om0y = om0y / v; om0z = om0z / v;
      v = sqrt(om1x*om1x+om1y*om1y+om1z*om1z);
      om1x = om1x / v; om1y = om1y / v; om1z = om1z / v;
    }
  }

  next: ;
  // Check array
  E_Int ni0, nj0, nk0;
  FldArrayF* f0; FldArrayI* cn0;
  char* varString0; char* eltType0;
  E_Int res = K_ARRAY::getFromArray(array, varString0, f0, 
                                    ni0, nj0, nk0, cn0, eltType0, true); 
  switch (res)
  {
    case 1:
      // check 2D i,j
      if (nk0 != 1)
      {
        RELEASESHAREDS(array, f0);
        PyErr_SetString(PyExc_TypeError,
                        "axisym: structured array must be an i-array or a i,j-array.");
        return NULL; 
      }
      break;

    case 2:
      if (strcmp(eltType0, "TRI") != 0 
          && strcmp(eltType0, "QUAD") != 0
          && strcmp(eltType0, "BAR") != 0)
      {
        RELEASESHAREDU(array, f0, cn0);
        PyErr_SetString(PyExc_TypeError,
                        "axisym: unstructured array must be a BAR, a TRI or QUAD-array.");
        return NULL;
      }
      break;
      
    default:
      PyErr_SetString(PyExc_TypeError,
                      "axisym: invalid array.");
      return NULL;       
  }

  // coordonnees
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString0);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString0);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString0);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f0, cn0);
    PyErr_SetString(PyExc_TypeError,
                    "axisym: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  
  // Transformation de teta en radians
  E_Float pi = 4*atan(1.);
  teta = teta*pi/180.;
  E_Float dteta = teta/(nteta-1);

  E_Int npts0 = f0->getSize();
  E_Int nfld0 = f0->getNfld();
  E_Float* xt0 = f0->begin(posx);
  E_Float* yt0 = f0->begin(posy);
  E_Float* zt0 = f0->begin(posz);
  
  E_Int npts = npts0*nteta;
  FldArrayF* f = new FldArrayF(npts, nfld0);
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  FldArrayF tmp(npts0, nfld0);
  E_Float* xo = tmp.begin(1);
  E_Float* yo = tmp.begin(2);
  E_Float* zo = tmp.begin(3);

  // Init
  for (E_Int eq = 1; eq <= nfld0; eq++)
  {
    E_Float* fp = f->begin(eq);
    E_Float* fp0 = f0->begin(eq);
    for (E_Int ind = 0; ind < npts0; ind++)
    { 
      fp[ind] = fp0[ind];
    }
  }

  // Calcul des coordonnees
  /*
  for (E_Int ind = 0; ind < npts0; ind++)
  { xt[ind] = xt0[ind]; yt[ind] = yt0[ind]; zt[ind] = zt0[ind]; }
  */
  /*
  if (userR == true)
  {
    omx = xr[0]-xrC;
    omy = yr[0]-yrC;
    omz = zr[0]-zrC;
    rm = sqrt(omx*omx+omy*omy+omz*omz);
    for (E_Int ind = 0; ind < npts0; ind++)
    { 
      xt[ind] = xt0[ind]; yt[ind] = yt0[ind]; zt[ind] = zt0[ind]; 
    }
  }
  */

  E_Float tetai, rm=1.;
  //printf("xc,yc,zc %f %f %f\n", xrC, yrC, zrC);
  //printf("axis0 %f %f %f\n", om0x,om0y,om0z);
  //printf("axis1 %f %f %f\n", om1x,om1y,om1z);
  E_Float omx, omy, omz, X, Y;
  for (E_Int k = 0; k < nteta; k++)
  { 
    if (useR == true)
    {
      omx = xr[k]-xrC; omy = yr[k]-yrC; omz = zr[k]-zrC;
      X = omx*om0x+omy*om0y+omz*om0z;
      Y = omx*om1x+omy*om1y+omz*om1z;
      tetai = atan2(Y, X);
      if (tetai < 0) tetai = 2*pi+tetai;
      rm = sqrt(omx*omx+omy*omy+omz*omz); // r(teta) factor
      //printf("%d -  tetai=%f, r=%f\n", k, tetai, rm);
      //tetai = k*dteta;
    }
    else tetai = k*dteta;
    k6axisym_(npts0, xt0, yt0, zt0, xc, yc, zc, nx, ny, nz, tetai, rm,
              xo, yo, zo);
    
    for (E_Int eq = 1; eq <= nfld0; eq++)
    {
      E_Float* fp = f->begin(eq);
      E_Float* fp0 = f0->begin(eq);

      for (E_Int ind = 0; ind < npts0; ind++)
      { 
        E_Int indn = ind + k*npts0;
        fp[indn] = fp0[ind];
      }
    }
    for (E_Int ind = 0; ind < npts0; ind++)
    { 
      E_Int indn = ind + k*npts0;
      xt[indn] = xo[ind]; yt[indn] = yo[ind]; zt[indn] = zo[ind]; 
    }
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f0);
    if (useR != 0) RELEASESHAREDB(useR, arrayR, fR, cnR);
    PyObject* tpl = K_ARRAY::buildArray(*f, varString0, ni0, nj0, nteta);
    delete f;
    return tpl;
  }
  else
  {
    char eltType[20]; E_Int nvert = 0;
    if (strcmp(eltType0, "TRI") == 0) 
    { strcpy(eltType, "PENTA"); nvert = 6; }
    else if (strcmp(eltType0, "QUAD") == 0)
    { strcpy(eltType, "HEXA"); nvert = 8; }
    else if (strcmp(eltType0, "BAR") == 0)
    { strcpy(eltType, "QUAD"); nvert = 4; }

    E_Int nelts0 = cn0->getSize(); E_Int nvert0 = cn0->getNfld();
    E_Int nelts = nelts0*(nteta-1); 
    
    FldArrayI* cn = new FldArrayI(nelts, nvert);
    cn->setAllValuesAt(-1);

    /* calcul de la connectivite */
    for (E_Int v = 1; v <= nvert0; v++)
    {
      E_Int* cnv = cn->begin(v); 
      E_Int* cnv0 = cn0->begin(v);
      
      for (E_Int k = 0; k < nteta; k++)
        for (E_Int et = 0; et < nelts0; et++)
        {
          E_Int eti = et + k*nelts0;
          cnv[eti] = cnv0[et] + k *npts0;    
        }
    }

    if (nvert == 4)
    {
      E_Int* cnv = cn->begin(3); 
      E_Int* cnvp = cn->begin(2);
      for (E_Int et = 0; et < nelts; et++) cnv[et] = cnvp[et]+npts0;
      cnv = cn->begin(4); 
      cnvp = cn->begin(1);
      for (E_Int et = 0; et < nelts; et++) cnv[et] = cnvp[et]+npts0;
    }
    else
    {
      for (E_Int v = nvert0+1; v <= nvert; v++)
      {
        E_Int* cnv = cn->begin(v); 
        E_Int* cnvp = cn->begin(v-nvert0);
        for (E_Int et = 0; et < nelts; et++) cnv[et] = cnvp[et]+npts0;
      }
    }

    RELEASESHAREDU(array, f0, cn0);
    if (useR != 0) RELEASESHAREDB(useR, arrayR, fR, cnR);
    PyObject* tpl = K_ARRAY::buildArray(*f, varString0, *cn, -1, eltType);
    delete f; delete cn;
    return tpl;
  }
}
