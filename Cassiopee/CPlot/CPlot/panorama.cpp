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
#include "cplot.h"
#include "Data.h"
#include <math.h>

#define M_PI 3.1415926535897932384626433832795

// interpolation when using images
void interp(E_Int ind,
            E_Float* final1, E_Float* final2, E_Float* final3, E_Float* final4,
            E_Float* im1, E_Float* im2, E_Float* im3, E_Float* im4, 
            E_Float px, E_Float py, E_Int ni1, E_Int nj1)
{
  E_Int i1 = E_Int(px*ni1);
  E_Int i1p1 = std::min(i1+1, ni1);
  E_Int j1 = E_Int(py*nj1);
  E_Int j1p1 = std::min(j1+1, nj1);
  E_Int ind1 = i1+j1*(ni1+1);
  E_Int ind2 = i1p1+j1*(ni1+1);
  E_Int ind3 = i1+j1p1*(ni1+1);
  E_Int ind4 = i1p1+j1p1*(ni1+1);
  E_Float ax = px*ni1-i1;
  E_Float bx = 1.-ax;
  E_Float ay = py*nj1-j1;
  E_Float by = 1.-ay;
  //printf("%g %g // %g %d\n", ax, ay, py, nj1);
  final1[ind] = bx*by*im1[ind1]+ax*by*im1[ind2]+bx*ay*im1[ind3]+ax*ay*im1[ind4];
  final2[ind] = bx*by*im2[ind1]+ax*by*im2[ind2]+bx*ay*im2[ind3]+ax*ay*im2[ind4];
  final3[ind] = bx*by*im3[ind1]+ax*by*im3[ind2]+bx*ay*im3[ind3]+ax*ay*im3[ind4];
  final4[ind] = bx*by*im4[ind1]+ax*by*im4[ind2]+bx*ay*im4[ind3]+ax*ay*im4[ind4];
}

// interpolation when using slit buffers
void interp2(E_Int ind,
             char* final, char* im, 
             E_Float px, E_Float py, E_Int ni1, E_Int nj1)
{
  E_Int i1 = E_Int(px*ni1);
  E_Int i1p1 = std::min(i1+1, ni1);
  E_Int j1 = E_Int(py*nj1);
  E_Int j1p1 = std::min(j1+1, nj1);

  E_Int ind1 = i1+j1*(ni1+1);
  E_Int ind2 = i1p1+j1*(ni1+1);
  E_Int ind3 = i1+j1p1*(ni1+1);
  E_Int ind4 = i1p1+j1p1*(ni1+1);
  
  E_Float ax = px*ni1-i1;
  E_Float bx = 1.-ax;
  E_Float ay = py*nj1-j1;
  E_Float by = 1.-ay;

  ax = 0.5; bx = 0.5;
  ay = 0.5; by = 0.5;
  
  //printf("%g %g // %g %d\n", ax, ay, py, nj1);
  //printf("%g %g\n", ax, bx);

  E_Float a1, a2, a3, a4;
  uint8_t ret;

  //char c = (uint8_t)im[3*ind1];
  //ret = (int)c;

  a1 = bx*by*(uint8_t)im[3*ind1];
  a2 = ax*by*(uint8_t)im[3*ind2];
  a3 = bx*ay*(uint8_t)im[3*ind3];
  a4 = ax*ay*(uint8_t)im[3*ind4];
  ret = std::round(a1 + a2 + a3 + a4);
  //printf("ret=%d // %d %d\n", ret, im[3*ind1], im[3*ind2]);
  final[3*ind] = (char)ret;
  a1 = bx*by*(uint8_t)im[3*ind1+1];
  a2 = ax*by*(uint8_t)im[3*ind2+1];
  a3 = bx*ay*(uint8_t)im[3*ind3+1];
  a4 = ax*ay*(uint8_t)im[3*ind4+1];
  ret = std::round(a1 + a2 + a3 + a4);
  final[3*ind+1] = (char)(ret);
  a1 = bx*by*(uint8_t)im[3*ind1+2];
  a2 = ax*by*(uint8_t)im[3*ind2+2];
  a3 = bx*ay*(uint8_t)im[3*ind3+2];
  a4 = ax*ay*(uint8_t)im[3*ind4+2];
  ret = std::round(a1 + a2 + a3 + a4);
  final[3*ind+2] = (char)(ret);

  //final[3*ind] = E_Int(im[3*ind4]);
  //final[3*ind+1] = E_Int(im[3*ind4+1]);
  //final[3*ind+2] = E_Int(im[3*ind4+2]);
  
  //final[3*ind]   = bx*by*im[3*ind1]  +ax*by*im[3*ind2]  +bx*ay*im[3*ind3]  +ax*ay*im[3*ind4];
  //final[3*ind+1] = bx*by*im[3*ind1+1]+ax*by*im[3*ind2+1]+bx*ay*im[3*ind3+1]+ax*ay*im[3*ind4+1];
  //final[3*ind+2] = bx*by*im[3*ind1+2]+ax*by*im[3*ind2+2]+bx*ay*im[3*ind3+2]+ax*ay*im[3*ind4+2];
}

// perform the stitching (identical to panorama.frag but on the cpu)
// it doesnt have the texture size limit
// si type360=0 -> 360 deg
// si type360=1 -> 180 deg
// shift: eye shift
// fov2: fov enlargement
PyObject* K_CPLOT::panorama(PyObject* self, PyObject* args)
{
  // Get the 4 arrays of cube images (left, right, bottom, top, back, front)
  PyObject* leftArray; PyObject* rightArray;
  PyObject* bottomArray; PyObject* topArray;
  PyObject* backArray; PyObject* frontArray;
  PyObject* finalArray;
  E_Int type360;
  if (!PYPARSETUPLE_(args, OOOO_ OOO_ I_, &leftArray, &rightArray,
    &bottomArray, &topArray, &backArray, &frontArray, &finalArray, 
    &type360))
  {
    return NULL;
  }
  char* varString;
  E_Int ni, nj, nk, res;
  FldArrayF* left; FldArrayI* cn; char* eltType;
  res = K_ARRAY::getFromArray3(leftArray, varString, left, 
                               ni, nj, nk, cn, eltType);
  if (res != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "panorama: requires a structured array (left).");
    return NULL;
  }

  FldArrayF* right;
  res = K_ARRAY::getFromArray3(rightArray, varString, right, 
                               ni, nj, nk, cn, eltType);
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "panorama: requires a structured array (right).");
    return NULL;
  }

  FldArrayF* bottom;
  res = K_ARRAY::getFromArray3(bottomArray, varString, bottom, 
                               ni, nj, nk, cn, eltType);
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "panorama: requires a structured array (bottom).");
    return NULL;
  }

  FldArrayF* top;
  res = K_ARRAY::getFromArray3(topArray, varString, top, 
                               ni, nj, nk, cn, eltType);
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "panorama: requires a structured array (top).");
    return NULL;
  }

  FldArrayF* back;
  res = K_ARRAY::getFromArray3(backArray, varString, back, 
                               ni, nj, nk, cn, eltType);
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "panorama: requires a structured array (back).");
    return NULL;
  }

  FldArrayF* front;
  res = K_ARRAY::getFromArray3(frontArray, varString, front, 
                               ni, nj, nk, cn, eltType);
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "panorama: requires a structured array (front).");
    return NULL;
  }

  FldArrayF* final;
  E_Int nil, njl, nkl;
  res = K_ARRAY::getFromArray3(finalArray, varString, final, 
                               nil, njl, nkl, cn, eltType);
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "panorama: requires a structured array (final).");
    return NULL;
  }

  // code
  E_Int nijl = nil*njl;
  E_Int nil1 = nil-1;
  E_Int njl1 = njl-1;
  E_Int ni1 = ni-1;
  E_Int nj1 = nj-1;

  E_Float* left1 = left->begin(4); // r,g,b,a
  E_Float* left2 = left->begin(5);
  E_Float* left3 = left->begin(6);
  E_Float* left4 = left->begin(7);
  
  E_Float* right1 = right->begin(4); // r,g,b,a
  E_Float* right2 = right->begin(5);
  E_Float* right3 = right->begin(6);
  E_Float* right4 = right->begin(7);

  E_Float* bottom1 = bottom->begin(4); // r,g,b,a
  E_Float* bottom2 = bottom->begin(5);
  E_Float* bottom3 = bottom->begin(6);
  E_Float* bottom4 = bottom->begin(7);

  E_Float* top1 = top->begin(4); // r,g,b,a
  E_Float* top2 = top->begin(5);
  E_Float* top3 = top->begin(6);
  E_Float* top4 = top->begin(7);
  
  E_Float* back1 = back->begin(4); // r,g,b,a
  E_Float* back2 = back->begin(5);
  E_Float* back3 = back->begin(6);
  E_Float* back4 = back->begin(7);
  
  E_Float* front1 = front->begin(4); // r,g,b,a
  E_Float* front2 = front->begin(5);
  E_Float* front3 = front->begin(6);
  E_Float* front4 = front->begin(7);
  
  E_Float* final1 = final->begin(4); // r,g,b,a
  E_Float* final2 = final->begin(5);
  E_Float* final3 = final->begin(6);
  E_Float* final4 = final->begin(7);
  E_Float tinf, tsup;

  if (type360 == 0) { tinf = -M_PI; tsup = 2*M_PI; } // 360
  else  { tinf = -M_PI/2.; tsup = M_PI; } // 180

  // fov of each image
  //E_Float fov = 90.;
  //printf("fov=%g\n", fov);
  
  #pragma omp parallel
  {
    E_Int ii, jj;
    E_Float tx, ty, px, py;
    E_Float theta, phi, x, y, z;
    E_Float scale;
    # pragma omp for
    for (E_Int ind = 0; ind < nijl; ind++)
    {
      jj = E_Int(ind/nil);
      ii = ind - jj*nil;
      tx = (1.*ii) / nil1;
      ty = (1.*jj) / njl1;

      theta = tinf + tx * tsup; // between -pi and pi
      phi = -M_PI/2. + ty * M_PI; // between pi/2 and -pi/2

      x = cos(phi) * sin(theta);
      y = sin(phi);
      z = cos(phi) * cos(theta);
      //printf("%d %d -tx=%g %g %g %g %g\n", ii,jj,tx,ty,x,y,z);
    
      if (std::abs(x) >= std::abs(y) && std::abs(x) >= std::abs(z)) 
      {
        if (x < 0.0)
        {
          scale = -1.0 / x;
          px = ( z*scale + 1.0) / 2.0;
          py = ( y*scale + 1.0) / 2.0;
          interp(ind, final1, final2, final3, final4,
                 left1, left2, left3, left4,
                 px, py, ni1, nj1);
        }
        else
        {
          scale = 1.0 / x;
          px = (-z*scale + 1.0) / 2.0;
          py = ( y*scale + 1.0) / 2.0;
          interp(ind, final1, final2, final3, final4,
                 right1, right2, right3, right4, 
                 px, py, ni1, nj1);
        }
      }
      else if (std::abs(y) >= std::abs(z))
      {
        if (y > 0.0)
        {
          scale = -1.0 / y;
          px = (-x*scale + 1.0) / 2.0;
          py = ( z*scale + 1.0) / 2.0;
          interp(ind, final1, final2, final3, final4,
                 bottom1, bottom2, bottom3, bottom4, 
                 px, py, ni1, nj1);
        }
        else 
        {
          scale = 1.0 / y;
          px = (-x*scale + 1.0) / 2.0;
          py = (-z*scale + 1.0) / 2.0;
          interp(ind, final1, final2, final3, final4,
                 top1, top2, top3, top4,
                 px, py, ni1, nj1);
        }
      }
      else 
      {
        if (z < 0.0) 
        {
          scale = -1.0 / z;
          px = (-x*scale + 1.0) / 2.0;
          py = ( y*scale + 1.0) / 2.0;
          interp(ind, final1, final2, final3, final4,
                 back1, back2, back3, back4, 
                 px, py, ni1, nj1);
        }
        else
        {
          scale = 1.0 / z;
          px = ( x*scale + 1.0) / 2.0;
          py = ( y*scale + 1.0) / 2.0;
          interp(ind, final1, final2, final3, final4,
                 front1, front2, front3, front4, 
                 px, py, ni1, nj1);
        }
      }
    }
  }

  RELEASESHAREDS(frontArray, front);
  RELEASESHAREDS(backArray, back);
  RELEASESHAREDS(leftArray, left);
  RELEASESHAREDS(rightArray, right);
  RELEASESHAREDS(topArray, top);
  RELEASESHAREDS(bottomArray, bottom);
  RELEASESHAREDS(finalArray, final);
  
  return Py_None;
}

//===========================================================================

PyObject* K_CPLOT::panoramaODS(PyObject* self, PyObject* args)
{
  // Get the n arrays of cube images
  PyObject* front;
  PyObject* top;
  PyObject* bottom;
  PyObject* finalArray;
  E_Int type360;
  if (!PYPARSETUPLE_(args, OOOO_ I_, &front, &top, &bottom, &finalArray, &type360))
  {
    return NULL;
  }

  E_Int nangles = PyList_Size(front);
  char* varString;
  E_Int ni, nj, nk, res;
  FldArrayI* cn; char* eltType;
  std::vector<FldArrayF*> frontImages(nangles);
  for (E_Int i = 0; i < nangles; i++)
  {
    res = K_ARRAY::getFromArray3(PyList_GetItem(front, i), varString, frontImages[i], 
                                 ni, nj, nk, cn, eltType);
    if (res != 1) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "panorama: requires a structured array.");
      return NULL;
    }
  }

  std::vector<FldArrayF*> topImages(nangles);
  for (E_Int i = 0; i < nangles; i++)
  {
    res = K_ARRAY::getFromArray3(PyList_GetItem(top, i), varString, topImages[i], 
                                 ni, nj, nk, cn, eltType);
    if (res != 1) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "panorama: requires a structured array.");
      return NULL;
    }
  }

  std::vector<FldArrayF*> botImages(nangles);
  for (E_Int i = 0; i < nangles; i++)
  {
    res = K_ARRAY::getFromArray3(PyList_GetItem(bottom, i), varString, botImages[i], 
                                 ni, nj, nk, cn, eltType);
    if (res != 1) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "panorama: requires a structured array.");
      return NULL;
    }
  }

  FldArrayF* final;
  E_Int nil, njl, nkl;
  res = K_ARRAY::getFromArray3(finalArray, varString, final, 
                               nil, njl, nkl, cn, eltType);
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "panorama: requires a structured array (final).");
    return NULL;
  }

  assert(njl == nj);

  // slit image pointers
  std::vector<E_Float*> imf1(nangles);
  std::vector<E_Float*> imf2(nangles);
  std::vector<E_Float*> imf3(nangles);
  std::vector<E_Float*> imf4(nangles);
  for (E_Int i = 0; i < nangles; i++)
  {
    imf1[i] = frontImages[i]->begin(4);
    imf2[i] = frontImages[i]->begin(5);
    imf3[i] = frontImages[i]->begin(6);
    imf4[i] = frontImages[i]->begin(7);
  }
  std::vector<E_Float*> imt1(nangles);
  std::vector<E_Float*> imt2(nangles);
  std::vector<E_Float*> imt3(nangles);
  std::vector<E_Float*> imt4(nangles);
  for (E_Int i = 0; i < nangles; i++)
  {
    imt1[i] = topImages[i]->begin(4);
    imt2[i] = topImages[i]->begin(5);
    imt3[i] = topImages[i]->begin(6);
    imt4[i] = topImages[i]->begin(7);
  }
  std::vector<E_Float*> imb1(nangles);
  std::vector<E_Float*> imb2(nangles);
  std::vector<E_Float*> imb3(nangles);
  std::vector<E_Float*> imb4(nangles);
  for (E_Int i = 0; i < nangles; i++)
  {
    imb1[i] = botImages[i]->begin(4);
    imb2[i] = botImages[i]->begin(5);
    imb3[i] = botImages[i]->begin(6);
    imb4[i] = botImages[i]->begin(7);
  }

  // final image pointers
  E_Float* final1 = final->begin(4); // r,g,b,a
  E_Float* final2 = final->begin(5);
  E_Float* final3 = final->begin(6);
  E_Float* final4 = final->begin(7);

  //E_Float tinf, tsup;
  //if (type360 == 0) { tinf = -M_PI; tsup = 2*M_PI; } // 360
  //else  { tinf = -M_PI/2.; tsup = M_PI; } // 180

  E_Int ni1 = ni-1; // cube image
  E_Int nj1 = nj-1;
  
  E_Int ind;
  
  // direct sliting
  /*
  E_Int mid = ni/2;
  for (E_Int i = 0; i < nangles; i++)
  {
    for (E_Int j = 0; j < njl; j++)
    {
      ind = i + j*nil;
      final1[ind] = imb1[i][mid+j*ni];
      final2[ind] = imb2[i][mid+j*ni];
      final3[ind] = imb3[i][mid+j*ni];
      final4[ind] = imb4[i][mid+j*ni];
    }
  }
  */
  //return Py_None;

  // transformation
  E_Float y, z, scale, py, phi, ty;

  for (E_Int i = 0; i < nangles; i++)
  {
    // 1D interpolation
    for (E_Int j = 0; j < njl; j++)
    {
      ty = (1.*j)/nj1;
      phi = - M_PI/2. + ty * M_PI; // between pi/2 and -pi/2
  
      ind = i + j*nil;
      
      y = sin(phi); // entre -1 et 1
      z = cos(phi); // entre 0 et 0

      if (phi > -M_PI/4. && phi < M_PI/4.)
      {
        scale = 1.0 / z;
        py = ( y*scale + 1.0) / 2.0;
        //printf("%g %g\n", py, phi);
        interp(ind, final1, final2, final3, final4,
              imf1[i], imf2[i], imf3[i], imf4[i],
              0.5, py, ni1, nj1);
      }
      else if (phi >= M_PI/4.) // bottom
      {
        scale = -1.0 / y;
        py = ( z*scale + 1.0) / 2.0;
        //printf("%g %g\n", py, phi);
        interp(ind, final1, final2, final3, final4,
              imb1[i], imb2[i], imb3[i], imb4[i],
              0.5, py, ni1, nj1);
      }
      else if (phi <= -M_PI/4.) // top
      {
        scale = -1.0 / y;
        py = ( z*scale + 1.0) / 2.0;
        //printf("%g %g\n", py, phi);
        interp(ind, final1, final2, final3, final4,
              imt1[i], imt2[i], imt3[i], imt4[i],
              0.5, py, ni1, nj1);
      }
    }
  }

  for (E_Int i = 0; i < nangles; i++)
  {
    RELEASESHAREDS(PyList_GetItem(front, i), frontImages[i]);
    RELEASESHAREDS(PyList_GetItem(top, i), topImages[i]);
    RELEASESHAREDS(PyList_GetItem(bottom, i), botImages[i]);
  } 
  
  RELEASESHAREDS(finalArray, final);
  
  return Py_None;
}

// acumulate one slit in image at place i
void accumulateSlit(E_Int ni, E_Int nj, char* imf, char* imt, char* imb, 
                    E_Int i, E_Int nil, E_Int njl, char* imOut)
{
  E_Int ni1 = ni-1; // slit image
  E_Int nj1 = nj-1;
  
  {
    E_Int ind;
    E_Float scale, ty, py, phi, y, z;

    // 1D interpolation
    for (E_Int j = 0; j < njl; j++)
    {
      ty = (1.*j)/nj1;
      phi = M_PI/2. - ty * M_PI; // between pi/2 and -pi/2
      ind = i + (njl-j-1)*nil;

      y = sin(phi); // entre -1 et 1
      z = cos(phi); // entre 0 et 0

      if (phi > -M_PI/4. && phi < M_PI/4.)
      {
        scale = 1.0 / z;
        py = ( y*scale + 1.0) / 2.0;
        interp2(ind, imOut, imf,
                0.5, py, ni1, nj1);
      }
      else if (phi >= M_PI/4.) // top
      {
        scale = -1.0 / y;
        py = ( z*scale + 1.0) / 2.0;
        interp2(ind, imOut, imt,
                0.5, py, ni1, nj1);
      }
      else if (phi <= -M_PI/4.) // bottom
      {
        scale = -1.0 / y;
        py = ( z*scale + 1.0) / 2.0;
        interp2(ind, imOut, imb,
                0.5, py, ni1, nj1);
      }
    }
  } 
}