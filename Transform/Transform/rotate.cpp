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

// Rotate functions
# include "transform.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

extern "C"
{
  void k6rotate_(const E_Int& npts,
                 const E_Float* x, const E_Float* y, const E_Float* z,
                 const E_Float& xc, const E_Float& yc, const E_Float& zc,
                 const E_Float& nx, const E_Float& ny, const E_Float& nz,
                 const E_Float& teta,
                 E_Float* xo, E_Float* yo, E_Float* zo);
}

//==============================================================================
E_Int K_TRANSFORM::extractVectorComponents(char* varString, 
                                           PyObject* listOfFieldVectors, 
                                           vector<E_Int>& posvx, vector<E_Int>& posvy, vector<E_Int>& posvz)
{
  if (PyList_Check(listOfFieldVectors) == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "rotate: last argument must be a list.");
    return -1;
  }
  E_Int nvectors = PyList_Size(listOfFieldVectors);
  if (nvectors == 0) return 0;

  for (E_Int novar = 0; novar < nvectors; novar++)
  {
    PyObject* tpl0 = PyList_GetItem(listOfFieldVectors, novar);
    if (PyList_Check(tpl0) == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "rotate: vector fields must be a list.");
      return -1;
    }
    E_Int sizeVect = PyList_Size(tpl0);
    if (sizeVect != 3)
    {
      PyErr_SetString(PyExc_TypeError,
                      "rotate: vector fields must be defined by 3 components.");
      return -1;       
    }
    // Check if each component is a string
    vector<char*> vars;
    for (E_Int noc = 0; noc < 3; noc++)
    {
      PyObject* tpl1 = PyList_GetItem(tpl0, noc);
      char* vect;
      if (PyString_Check(tpl1))
      {
        vect = PyString_AsString(tpl1); 
      }
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(tpl1))
      {
        vect = PyBytes_AsString(PyUnicode_AsUTF8String(tpl1)); 
      }
#endif 
      else 
      {
        PyErr_SetString(PyExc_TypeError,
                        "rotate: vector component name must be a string.");
        return -1;
      }
      vars.push_back(vect);
    }
    E_Int posu = K_ARRAY::isNamePresent(vars[0], varString);
    E_Int posv = K_ARRAY::isNamePresent(vars[1], varString);
    E_Int posw = K_ARRAY::isNamePresent(vars[2], varString);
    if (posu == -1 || posv == -1 || posw == -1) 
    {
      // printf("Warning: rotate: vector field (%s,%s,%s) not found in array.\n",vars[0],vars[1],vars[2]);
      ;
    }
    else 
    {
      posvx.push_back(posu+1); posvy.push_back(posv+1); posvz.push_back(posw+1); 
    }
  }
  return 0;
}
// ============================================================================
/* Rotate an array describing a mesh */
// ============================================================================
PyObject* K_TRANSFORM::rotateA1(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* listOfFieldVectors;
  E_Float xc, yc, zc;
  E_Float nx, ny, nz, teta;
  if (!PYPARSETUPLEF(args,
                    "O(ddd)(ddd)dO", "O(fff)(fff)fO",
                    &array, &xc, &yc, &zc, &nx, &ny, &nz, &teta, &listOfFieldVectors))
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

  res = K_ARRAY::getFromArray(array, varString, f, nil, njl, nkl, 
                              cn, eltType, true); 
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "rotate: invalid array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  posx++; posy++; posz++;
  vector<E_Int> posvx; vector<E_Int> posvy; vector<E_Int> posvz;
  E_Int ok = extractVectorComponents(varString, listOfFieldVectors, posvx, posvy, posvz);
  if (ok == -1) 
  {
    RELEASESHAREDB(res, array, f, cn); 
    return NULL;
  }

  // Transformation en radians
  E_Float pi = 4*atan(1.);
  teta = teta*pi/180.;
    
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();

  // Construit l'array resultat et l'initialise par copie
  PyObject* tpl;
  if (res == 1) //structured
  {
    tpl = K_ARRAY::buildArray(nfld, varString, 
                              nil, njl, nkl);
  } 
  else //unstructured 
  {
    E_Int csize = cn->getSize()*cn->getNfld(); 
    tpl = K_ARRAY::buildArray(nfld, varString,
                              npts, cn->getSize(),
                              -1, eltType, false, csize);
  }
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(npts, nfld, fnp, true);
  fn.setAllValuesAt(*f);

  if (res == 2)
  {
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
  }

  // rotate
  E_Int nvect = posvx.size();
  if (posx>0 && posy>0 && posz>0) //coord + champs
  {
    k6rotate_(npts, f->begin(posx), f->begin(posy), f->begin(posz),
              xc, yc, zc, nx, ny, nz, teta, 
              fn.begin(posx), fn.begin(posy), fn.begin(posz));          
  }
  for (E_Int nov = 0; nov < nvect; nov++)
  {
    posx = posvx[nov]; posy = posvy[nov]; posz = posvz[nov];
    k6rotate_(npts, f->begin(posx), f->begin(posy), f->begin(posz),
              0., 0., 0., nx, ny, nz, teta, 
              fn.begin(posx), fn.begin(posy), fn.begin(posz));          
  }

  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}

// ============================================================================
/* Rotate an array describing a mesh, input is two frames */
// ============================================================================
PyObject* K_TRANSFORM::rotateA2(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float xc, yc, zc;
  E_Float e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z;
  E_Float f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z;
  PyObject* listOfFieldVectors;

  if (!PYPARSETUPLEF(args,
                    "O(ddd)((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))O", "O(fff)((fff)(fff)(fff))((fff)(fff)(fff))O",
                    &array,
                    &xc, &yc, &zc, 
                    &e1x, &e1y, &e1z,
                    &e2x, &e2y, &e2z,
                    &e3x, &e3y, &e3z,
                    &f1x, &f1y, &f1z,
                    &f2x, &f2y, &f2z,
                    &f3x, &f3y, &f3z, 
                    &listOfFieldVectors))
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

  res = K_ARRAY::getFromArray(array, varString, f, nil, njl, nkl, 
                              cn, eltType, true); 
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "rotate: invalid array.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  posx++; posy++; posz++;
  vector<E_Int> posvx; vector<E_Int> posvy; vector<E_Int> posvz;
  E_Int ok = extractVectorComponents(varString, listOfFieldVectors, posvx, posvy, posvz);
  if (ok == -1) 
  {
    RELEASESHAREDB(res, array, f, cn); 
    return NULL;
  }
    
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();

  // Construit l'array resultat et l'initialise par copie
  PyObject* tpl;
  if (res == 1) //structured
  {
    tpl = K_ARRAY::buildArray(nfld, varString, 
                              nil, njl, nkl);
  } 
  else //unstructured 
  {
    E_Int csize = cn->getSize()*cn->getNfld(); 
    tpl = K_ARRAY::buildArray(nfld, varString,
                              npts, cn->getSize(),
                              -1, eltType, false, csize);
  }
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(npts, nfld, fnp, true);
  fn.setAllValuesAt(*f);
  if (res == 2)
  {
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
  }

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

  if (posx>0 && posy>0 && posz>0) //coord + champs
  {
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    E_Float* xp = fn.begin(posx);
    E_Float* yp = fn.begin(posy);
    E_Float* zp = fn.begin(posz);

#pragma omp parallel shared (npts, xt, yt, zt, xp, yp, zp, xc, yc, zc, m11, m12, m13, m21, m22, m23, m31, m32, m33) if (npts > 100)
    {
      E_Float x, y, z;
#pragma omp for nowait
      for (E_Int i = 0; i < npts; i++)
      {
        x = (xt[i]-xc);
        y = (yt[i]-yc);
        z = (zt[i]-zc);
        xp[i] = xc + m11*x + m12*y + m13*z;
        yp[i] = yc + m21*x + m22*y + m23*z;
        zp[i] = zc + m31*x + m32*y + m33*z;      
      }
    }
  }
  E_Int nvect = posvx.size();
  for (E_Int nov = 0; nov < nvect; nov++)
  {
    E_Int posx = posvx[nov];
    E_Int posy = posvy[nov];
    E_Int posz = posvz[nov];
    E_Float* fxt = f->begin(posx);
    E_Float* fyt = f->begin(posy);
    E_Float* fzt = f->begin(posz);
    E_Float* fxp = fn.begin(posx);
    E_Float* fyp = fn.begin(posy);
    E_Float* fzp = fn.begin(posz);
#pragma omp parallel shared (npts, fxt, fyt, fzt, fxp, fyp, fzp, m11, m12, m13, m21, m22, m23, m31, m32, m33) if (npts > 100)
    {
#pragma omp for nowait
      for (E_Int i = 0; i < npts; i++)
      {
        fxp[i] = m11*fxt[i] + m12*fyt[i] + m13*fzt[i];
        fyp[i] = m21*fxt[i] + m22*fyt[i] + m23*fzt[i];
        fzp[i] = m31*fxt[i] + m32*fyt[i] + m33*fzt[i];      
      }
    }
  }
  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}

// ============================================================================
/* Rotate an array describing a mesh */
// ============================================================================
PyObject* K_TRANSFORM::rotateA3(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float xc, yc, zc;
  E_Float alpha, beta, gamma;
  PyObject* listOfFieldVectors;

  if (!PYPARSETUPLEF(args,
                    "O(ddd)(ddd)O", "O(fff)(fff)O",
                    &array, &xc, &yc, &zc, &alpha, &beta, &gamma, &listOfFieldVectors))
  {
      return NULL;
  }

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray(array, varString, f, nil, njl, nkl, 
                              cn, eltType, true); 
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "rotate: invalid array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  posx++; posy++; posz++;
  vector<E_Int> posvx; vector<E_Int> posvy; vector<E_Int> posvz;
  E_Int ok = extractVectorComponents(varString, listOfFieldVectors, posvx, posvy, posvz);
  if (ok == -1) 
  {
    RELEASESHAREDB(res, array, f, cn); 
    return NULL;
  }
  // Transformation en radians
  E_Float pi = 4*atan(1.);
  alpha = alpha*pi/180.;
  beta = beta*pi/180.;
  gamma = gamma*pi/180.;

  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();

  // Construit l'array resultat et l'initialise par copie
  PyObject* tpl;
  if (res == 1) //structured
  {
    tpl = K_ARRAY::buildArray(nfld, varString, 
                              nil, njl, nkl);
  } 
  else //unstructured 
  {
    E_Int csize = cn->getSize()*cn->getNfld(); 
    tpl = K_ARRAY::buildArray(nfld, varString,
                              npts, cn->getSize(),
                              -1, eltType, false, csize);
  }
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(npts, nfld, fnp, true);
  fn.setAllValuesAt(*f);

  if (res == 2)
  {
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
  }

  // rotate
  E_Float calpha = cos(alpha);
  E_Float salpha = sin(alpha);
  E_Float cbeta = cos(beta);
  E_Float sbeta = sin(beta);
  E_Float cgamma = cos(gamma);
  E_Float sgamma = sin(gamma);
  if (posx>0 && posy>0 && posz>0) //coord + champs
  {
    E_Float* xp = fn.begin(posx);
    E_Float* yp = fn.begin(posy);
    E_Float* zp = fn.begin(posz);
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
        xp[i] = xc + dx;
        yp[i] = yc + calpha*dy - salpha*dz;
        zp[i] = zc + salpha*dy + calpha*dz;
      }
    } 
  }
  
  E_Int nvect = posvx.size();
  for (E_Int nov = 0; nov < nvect; nov++)
  {
    E_Int posx = posvx[nov];
    E_Int posy = posvy[nov];
    E_Int posz = posvz[nov];
    E_Float* fxt = f->begin(posx);
    E_Float* fyt = f->begin(posy);
    E_Float* fzt = f->begin(posz);
    E_Float* fxp = fn.begin(posx);
    E_Float* fyp = fn.begin(posy);
    E_Float* fzp = fn.begin(posz);
//#pragma omp parallel default(shared)
    {
      E_Float dx, dy, dz;
      E_Float x1, y1, z1, x2, y2, z2;
//#pragma omp for
      for (E_Int i = 0; i < npts; i++)
      {
        // Rotation autour de Oz (Ox->Ox1, Oy->Oy1, Oz->Oz)
        dx = fxt[i]; dy = fyt[i]; dz = fzt[i];
        x1 = cgamma*dx - sgamma*dy;
        y1 = sgamma*dx + cgamma*dy;
        z1 = dz;
      
        // Rotation autour de Oy1 (Ox1->Ox2, Oy1->Oy1, Oz1->Oz2)
        dx = x1; dy = y1; dz = z1;
        x2 = cbeta*dx - sbeta*dz;
        y2 = dy;
        z2 = sbeta*dx + cbeta*dz;
        
        // Rotation autour de Oz2 (Ox2->Ox3, Oy2->Oy3, Oz2->Oz2)
        dx = x2; dy = y2; dz = z2;
        fxp[i] = dx;
        fyp[i] = calpha*dy - salpha*dz;
        fzp[i] = salpha*dy + calpha*dz;
      }
    }
  }

  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}
