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

// Analytical grids creation

#include "generator.h"
#include <math.h>
using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

// ============================================================================
/* Create a regular cylindrical mesh of nixnjxnk points */
// ============================================================================
PyObject* K_GENERATOR::cylinderMesh(PyObject* self, PyObject* args)
{
  E_Int ni, nj, nk;
  E_Float xo, yo, zo;
  E_Float tetas, tetae, R1, R2, H;
  if (!PYPARSETUPLE_(args, TRRR_ RRRR_ R_ TIII_, 
                    &xo, &yo, &zo, &R1, &R2, &tetas, &tetae, 
                    &H, &ni, &nj, &nk)) return NULL;

  E_Float pi = 4*atan(1.);
  E_Float t1 = tetas*pi/180.;
  E_Float t2 = tetae*pi/180.;
  E_Float hk;

  // check ni, nj, nk
  if (ni < 1 || nj < 1 || nk < 1)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cylinder: ni, nj, nk must be >=1.");
    return NULL;
  }

  // Check R1, R2  
  if (R1 < 0. || R2 < 0.)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cylinder: R1 and R2 must be >=0.");
    return NULL;
  }
  
  // Check ni, nj, nk
  if (nk == 1) hk = 0.;
  else hk = H/(nk-1);

  E_Float hR, delta;
  if (nj == 1) hR = 0.;
  else hR = (R2-R1)/(nj-1);
  if (ni == 1) delta = 0.;
  else delta =  1./(ni-1.);
  E_Float dt2mt1 = delta* (t2-t1);

  // Create a portion of cylinder
  E_Int ninj = ni*nj;
  FldArrayF* coord = new FldArrayF(ninj*nk, 3);

  E_Float alpha, calpha, salpha, coef;
  E_Int ind;
  E_Float* xt = coord->begin(1);
  E_Float* yt = coord->begin(2);
  E_Float* zt = coord->begin(3);

  for (E_Int i = 0; i < ni; i++)
  {
    alpha = t1+i*dt2mt1;
    calpha = cos(alpha);
    salpha = sin(alpha);

    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
      {
        ind = i+j*ni+k*ninj;
        coef = j*hR+R1;
        xt[ind] = xo + coef*calpha;
        yt[ind] = yo + coef*salpha;
        zt[ind] = zo + k*hk;
      }  
  }
  // Build array
  PyObject* tpl = K_ARRAY::buildArray(*coord, "x,y,z", 
                                      ni, nj, nk);
  delete coord;
  return tpl;
}
// ============================================================================
/* Create a cylindrical mesh of nixnjxnk points from distributions,
   all distrib in [0,1]. */
// ============================================================================
PyObject* K_GENERATOR::cylinderMesh2(PyObject* self, PyObject* args)
{
  E_Float xo, yo, zo;
  E_Float tetas, tetae, R1, R2, H;
  PyObject* arrayR; PyObject* arrayT; PyObject* arrayZ;

#ifdef E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "(ddd)dddddOOO", 
                        &xo, &yo, &zo, &R1, &R2, &tetas, &tetae, 
                        &H, 
                        &arrayR, &arrayT, &arrayZ))   
#else
    if (!PyArg_ParseTuple(args, "(fff)fffffOOO",
                          &xo, &yo, &zo, &R1, &R2, &tetas, &tetae, 
                          &H, 
                          &arrayR, &arrayT, &arrayZ))   
#endif
  {
    PyErr_SetString(PyExc_TypeError, 
                    "cylinder2: wrong arguments.");
    return NULL;
  }

  // Check R1, R2
  if (R1 < 0. || R2 < 0.)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cylinder2: R1 and R2 must be >=0.");
    return NULL;
  }

  // Check array
  E_Int niR, njR, nkR, niT, njT, nkT, niZ, njZ, nkZ;
  FldArrayF* fR; FldArrayF* fT; FldArrayF* fZ;
  FldArrayI* cnR; FldArrayI* cnT; FldArrayI* cnZ;
  char* varStringR; char* varStringT; char* varStringZ;
  char* eltTypeR; char* eltTypeT; char* eltTypeZ;
  
  E_Int resR = K_ARRAY::getFromArray(arrayR, varStringR, fR, niR, njR, nkR, 
                                     cnR, eltTypeR); 
  E_Int resT = K_ARRAY::getFromArray(arrayT, varStringT, fT, niT, njT, nkT, 
                                     cnT, eltTypeT);
  E_Int resZ = K_ARRAY::getFromArray(arrayZ, varStringZ, fZ, niZ, njZ, nkZ, 
                                     cnZ, eltTypeZ);
  
  if (resR == 1 && resT == 1 && resZ == 1)
  {
    E_Float pi = 4*atan(1.);
    E_Float t1 = tetas*pi/180.;
    E_Float t2 = tetae*pi/180.;
    E_Int ni = niT;
    E_Int nj = niR;
    E_Int nk = niZ;
    E_Int posxr = K_ARRAY::isCoordinateXPresent(varStringR);
    E_Int posxt = K_ARRAY::isCoordinateXPresent(varStringT);
    E_Int posxz = K_ARRAY::isCoordinateXPresent(varStringZ);

    if (posxr == -1 || posxt == -1 || posxz == -1)
    {
      delete fR; delete fT; delete fZ;
      PyErr_SetString(PyExc_TypeError,
                      "cylinder2: can't find coordinates in array.");
      return NULL;
    }
    posxr++; posxt++; posxz++;

    if (nj == 1 || ni == 1)
    {
      delete fR; delete fT; delete fZ;
      PyErr_SetString(PyExc_TypeError, 
                      "cylinder2: creation of cylinder with ni = 1 or nj = 1 is impossible.");
      return NULL;
    }
    E_Float hR = (R2-R1);
    FldArrayF& dr = *fR;
    FldArrayF& dt = *fT;
    FldArrayF& dz = *fZ;

    // Create a portion of cylinder
    FldArrayF* coord = new FldArrayF(ni*nj*nk, 3);
    
    E_Float alpha, cosalpha, sinalpha;
    E_Int ind;
    E_Int ninj = ni*nj;
    E_Float t2mt1 = t2-t1;
    E_Float* xt = coord->begin(1);
    E_Float* yt = coord->begin(2);
    E_Float* zt = coord->begin(3);
    E_Float* drx = dr.begin(posxr);
    E_Float coef;

    for (E_Int i = 0; i < ni; i++)
    {
      alpha = t1 + t2mt1 * dt(i,posxt);
      cosalpha = cos(alpha);
      sinalpha = sin(alpha);
      
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj; j++)
        {        
          ind = i + j*ni + k*ninj;
          coef = hR * drx[j] + R1;
          xt[ind] = xo + coef * cosalpha;
          yt[ind] = yo + coef * sinalpha;
          zt[ind] = zo + H * dz(k,posxz);
        }
    }
    delete fR; delete fT; delete fZ;

    // Build array 
    PyObject* tpl = K_ARRAY::buildArray(*coord, "x,y,z", ni, nj, nk);
    delete coord;
    return tpl;
  }
  else if (resR == 2 || resT == 2 || resZ == 2)
  {
    delete fR; delete fT; delete fZ;
    if (resR == 2) delete cnR;
    if (resT == 2) delete cnT;
    if (resZ == 2) delete cnZ;
    PyErr_SetString(PyExc_TypeError,
                    "cylinder2: cannot be used with unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "cylinder2: one array is invalid.");
    return NULL;
  }
}

// ============================================================================
/* Create a cylindrical mesh of nixnjxnk points from distributions,
   all distrib in [0,1]. */
// ============================================================================
PyObject* K_GENERATOR::cylinderMesh3( PyObject* self,
                                      PyObject* args )
{
  E_Float tetas, tetae;
  PyObject* arrayXZ; PyObject* arrayT;
  if (!PYPARSETUPLE_(args, O_ RR_ O_,
                     &arrayXZ, &tetas, &tetae, &arrayT))
  {
    PyErr_SetString(PyExc_TypeError, 
                    "cylinder3: wrong arguments.");
    return NULL;
  }

  // Check array
  E_Int im1, jm1, km1, im2, jm2, km2;
  FldArrayF* fxz; FldArrayI* cnxz;
  char* varStringxz; char* eltTypexz;
  FldArrayF* ft; FldArrayI* cnt;
  char* varStringt; char* eltTypet;

  E_Int resxz = 
    K_ARRAY::getFromArray(arrayXZ, varStringxz, fxz, im1, jm1, km1, 
                          cnxz, eltTypexz); 
  E_Int rest = 
    K_ARRAY::getFromArray(arrayT, varStringt, ft, im2, jm2, km2, 
                          cnt, eltTypet);

  if (resxz == 1 && rest == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varStringxz);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varStringxz);
    E_Int post = K_ARRAY::isCoordinateXPresent(varStringt); 
    if (posx == -1 || posz == -1 || post == -1)
    {
      delete fxz; delete ft;
      PyErr_SetString(PyExc_TypeError,
                      "cylinder3: can't find coordinates in array.");
      return NULL;
    }
    posx++; posz++; post++;

    E_Int ni, nj, nk;
    E_Int ind, indp;
    E_Float alpha, cosalpha, sinalpha;
    E_Float pi = 4*atan(1.);
    E_Float t1 = tetas*pi/180.;
    E_Float t2 = tetae*pi/180.;
    
    FldArrayF& dt = *ft;
    FldArrayF& dxz = *fxz;

    ni = im2;
    
    if (im1 == 1)
    {
      nj = jm1;
      nk = km1;
    }
    else if ( jm1 == 1)
    {
      nj = im1;
      nk = km1;
    }
    else
    {
      nj = im1;
      nk = jm1;
    }

    // Create a portion of cylinder
    FldArrayF* coord = new FldArrayF(ni*nj*nk, 3);
    E_Int ninj = ni*nj;
    E_Float t2mt1 = t2-t1;
    E_Float* xt = coord->begin(1);
    E_Float* yt = coord->begin(2);
    E_Float* zt = coord->begin(3);
    E_Float* dxp = dxz.begin(posx);
    E_Float* dzp = dxz.begin(posz);
    E_Float* dtp = dt.begin(post);
 
    for (E_Int i = 0; i < ni; i++)
    {
      alpha = t1 + t2mt1 * dtp[i];
      cosalpha = cos(alpha);
      sinalpha = sin(alpha);
      
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj; j++)
        {        
          ind = i + j*ni + k*ninj;
          indp = j + k*nj;
          xt[ind] = dxp[indp] * cosalpha;
          yt[ind] = dxp[indp] * sinalpha;
          zt[ind] = dzp[indp];
        }
    }
    delete fxz; delete ft;

    // Build array
    PyObject* tpl =
      K_ARRAY::buildArray(*coord, "x,y,z", ni, nj, nk);
    delete coord;
    return tpl;
  }
  else if (resxz == 2 || rest == 2)
  {
    delete fxz; delete ft;
    if (resxz == 2) delete cnxz; 
    if (rest == 2) delete cnt;
    PyErr_SetString(PyExc_TypeError,
                    "cylinder3: can not be used on unstructured arrays.");
    return NULL;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                 "cylinder3: unknown type of array.");
    return NULL;
  }
}
