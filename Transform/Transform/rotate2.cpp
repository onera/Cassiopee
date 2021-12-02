/*    
    Copyright 2013-2022 Onera.

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

// Rotate functions in place only on coordinates
# include "transform.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Rotate an array describing a mesh - array2 in place */
// ============================================================================
PyObject* K_TRANSFORM::_rotateA1(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float xc, yc, zc;
  E_Float nx, ny, nz, teta;
  if (!PYPARSETUPLEF(args,
                    "O(ddd)(ddd)d", "O(fff)(fff)f",
                    &array, &xc, &yc, &zc, &nx, &ny, &nz, &teta))
  {
      return NULL;
  }

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;

  // Preliminary check
  if (K_FUNC::fEqualZero(nx) == true && 
      K_FUNC::fEqualZero(ny) == true &&
      K_FUNC::fEqualZero(nz) == true)
  {
    PyErr_SetString(PyExc_ValueError,
                    "rotate: vector has null norm.");
    return NULL; 
  }

  res = K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
                              cn, eltType); 
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "rotate: invalid array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "rotate: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // Transformation en radians
  E_Float pi = 4.*atan(1.);
  teta = teta*pi/180.;
    
  E_Int npts = f->getSize();
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  // rotate
  E_Float norm, unx, uny, unz;
  E_Float sinteta, sinteta5;
  E_Float e0,e1,e2,e3,a1;

  norm = nx*nx+ny*ny+nz*nz;
  norm = K_FUNC::E_max(norm, 1.e-20);  
  norm = 1./sqrt(norm);
  unx = nx*norm;
  uny = ny*norm;
  unz = nz*norm;

  sinteta = sin(teta);
  sinteta5 = sin(teta*0.5);

  // quaternions
  e0 = cos(teta*0.5);
  e1 = -unx*sinteta5;
  e2 = -uny*sinteta5;
  e3 = -unz*sinteta5;
  a1 = e0*e0-e1*e1-e2*e2-e3*e3;
  
#pragma omp parallel default(shared)
  {
    E_Float rx,ry,rz,a2,px,py,pz;
#pragma omp for 
    for (E_Int ind = 0; ind < npts; ind++)
    {
      rx = xt[ind]-xc;
      ry = yt[ind]-yc;
      rz = zt[ind]-zc;
      a2 = e1*rx+e2*ry+e3*rz;
      px = a1*rx+2*e1*a2-(ry*unz-rz*uny)*sinteta;
      py = a1*ry+2*e2*a2-(rz*unx-rx*unz)*sinteta;
      pz = a1*rz+2*e3*a2-(rx*uny-ry*unx)*sinteta;
      xt[ind] = xc+px;
      yt[ind] = yc+py;
      zt[ind] = zc+pz;
    }
  }

  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}

// ============================================================================
/* Rotate an array describing a mesh, input is two frames */
// ============================================================================
PyObject* K_TRANSFORM::_rotateA2(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float xc, yc, zc;
  E_Float e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z;
  E_Float f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z;

  if (!PYPARSETUPLEF(args,
                     "O(ddd)((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))", 
                     "O(fff)((fff)(fff)(fff))((fff)(fff)(fff))",
                     &array,
                     &xc, &yc, &zc, 
                     &e1x, &e1y, &e1z,
                     &e2x, &e2y, &e2z,
                     &e3x, &e3y, &e3z,
                     &f1x, &f1y, &f1z,
                     &f2x, &f2y, &f2z,
                     &f3x, &f3y, &f3z))
  {
      return NULL;
  }

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;

  // Make unity
  E_Float e1 = sqrt(e1x*e1x+e1y*e1y+e1z*e1z);
  E_Float e2 = sqrt(e2x*e2x+e2y*e2y+e2z*e2z);
  E_Float e3 = sqrt(e3x*e3x+e3y*e3y+e3z*e3z);
  E_Float f1 = sqrt(f1x*f1x+f1y*f1y+f1z*f1z);
  E_Float f2 = sqrt(f2x*f2x+f2y*f2y+f2z*f2z);
  E_Float f3 = sqrt(f3x*f3x+f3y*f3y+f3z*f3z);  

  if (K_FUNC::fEqualZero(e1) == true || 
      K_FUNC::fEqualZero(e2) == true ||
      K_FUNC::fEqualZero(e3) == true ||
      K_FUNC::fEqualZero(f1) == true || 
      K_FUNC::fEqualZero(f2) == true ||
      K_FUNC::fEqualZero(f3) == true)
  {
    PyErr_SetString(PyExc_ValueError,
                    "rotate: frame vectors must have non null norm.");
    return NULL;
  }

  e1 = 1./e1; e2 = 1./e2; e3 = 1./e3;
  f1 = 1./f1; f2 = 1./f2; f3 = 1./f3;
  e1x = e1x * e1; e1y = e1y * e1; e1z = e1z * e1;
  e2x = e2x * e2; e2y = e2y * e2; e2z = e2z * e2;
  e3x = e3x * e3; e3y = e3y * e3; e3z = e3z * e3;
  f1x = f1x * f1; f1y = f1y * f1; f1z = f1z * f1;
  f2x = f2x * f2; f2y = f2y * f2; f2z = f2z * f2;
  f3x = f3x * f3; f3y = f3y * f3; f3z = f3z * f3;

  res = K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
                               cn, eltType); 
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "rotate: invalid array.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "rotate: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  
  E_Int npts = f->getSize();
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);
  
  // rotate
  E_Float m11 = f1x*e1x+f2x*e2x+f3x*e3x;
  E_Float m21 = f1y*e1x+f2y*e2x+f3y*e3x;
  E_Float m31 = f1z*e1x+f2z*e2x+f3z*e3x;
  E_Float m12 = f1x*e1y+f2x*e2y+f3x*e3y;
  E_Float m22 = f1y*e1y+f2y*e2y+f3y*e3y;
  E_Float m32 = f1z*e1y+f2z*e2y+f3z*e3y;
  E_Float m13 = f1x*e1z+f2x*e2z+f3x*e3z;
  E_Float m23 = f1y*e1z+f2y*e2z+f3y*e3z;
  E_Float m33 = f1z*e1z+f2z*e2z+f3z*e3z;

#pragma omp parallel default(shared)
  {
    E_Float x,y,z;
#pragma omp for nowait
    for (E_Int i = 0; i < npts; i++)
    {
      x = xt[i]-xc;
      y = yt[i]-yc;
      z = zt[i]-zc;
      xt[i] = xc + m11*x + m12*y + m13*z;
      yt[i] = yc + m21*x + m22*y + m23*z;
      zt[i] = zc + m31*x + m32*y + m33*z;      
    }
  }
  
  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}

// ============================================================================
/* Rotate an array describing a mesh */
// ============================================================================
PyObject* K_TRANSFORM::_rotateA3(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float xc, yc, zc;
  E_Float alpha, beta, gamma;
  if (!PYPARSETUPLEF(args,
                    "O(ddd)(ddd)", "O(fff)(fff)",
                    &array, &xc, &yc, &zc, &alpha, &beta, &gamma))
  {
      return NULL;
  }

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
                               cn, eltType); 
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "rotate: invalid array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "rotate: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  
  // Transformation en radians
  E_Float pi = 4*atan(1.);
  alpha = alpha*pi/180.;
  beta = beta*pi/180.;
  gamma = gamma*pi/180.;

  E_Int npts = f->getSize();

  // rotate
  E_Float calpha = cos(alpha);
  E_Float salpha = sin(alpha);
  E_Float cbeta = cos(beta);
  E_Float sbeta = sin(beta);
  E_Float cgamma = cos(gamma);
  E_Float sgamma = sin(gamma);
    
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  
#pragma omp parallel default(shared)
  {
    E_Float dx, dy, dz;
    E_Float x1, y1, z1, x2, y2, z2;
#pragma omp for
    for (E_Int i = 0; i < npts; i++)
    {
      // Rotation autour de Oz (Ox->Ox1, Oy->Oy1, Oz->Oz)
      dx = x[i]-xc; dy = y[i]-yc; dz = z[i]-zc;
      x1 = xc + cgamma*dx - sgamma*dy;
      y1 = yc + sgamma*dx + cgamma*dy;
      z1 = zc + dz;
      
      // Rotation autour de Oy1 (Ox1->Ox2, Oy1->Oy1, Oz1->Oz2)
      dx = x1-xc; dy = y1-yc; dz = z1-zc;
      x2 = xc + cbeta*dx - sbeta*dz;
      y2 = yc + dy;
      z2 = zc + sbeta*dx + cbeta*dz;
      
      // Rotation autour de Oz2 (Ox2->Ox3, Oy2->Oy3, Oz2->Oz2)
      dx = x2-xc; dy = y2-yc; dz = z2-zc;
      x[i] = xc + dx;
      y[i] = yc + calpha*dy - salpha*dz;
      z[i] = zc + salpha*dy + calpha*dz;
    } 
  }
  
  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}
