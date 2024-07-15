/*    
    Copyright 2013-2024 Onera.

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
  //printf("ind=%d %g\n", ind, final1[ind]);
  final2[ind] = bx*by*im2[ind1]+ax*by*im2[ind2]+bx*ay*im2[ind3]+ax*ay*im2[ind4];
  final3[ind] = bx*by*im3[ind1]+ax*by*im3[ind2]+bx*ay*im3[ind3]+ax*ay*im3[ind4];
  final4[ind] = bx*by*im4[ind1]+ax*by*im4[ind2]+bx*ay*im4[ind3]+ax*ay*im4[ind4];
}

// perform the stitching (identical to panorama.frag but on the cpu)
// it doesnt have the texture size limit
// si type360=0 -> 360 deg
// si type360=1 -> 180 deg
PyObject* K_CPLOT::panorama(PyObject* self, PyObject* args)
{
  // Get the 4 arrays of cube images (left, right, bottom, top, back, front)
  PyObject* leftArray; PyObject* rightArray;
  PyObject* bottomArray; PyObject* topArray;
  PyObject* backArray; PyObject* frontArray;
  PyObject* finalArray;
  E_Int type360;
  if (!PYPARSETUPLE_(args, OOOO_ OOO_ I_, &leftArray, &rightArray,
    &bottomArray, &topArray, &backArray, &frontArray, &finalArray, &type360))
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
#define M_PI 3.1415926535897932384626433832795

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
      phi = -M_PI/2. + ty * M_PI; // between -pi/2 and pi/2

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
          //printf("1px=%g %g\n", px, py);
          interp(ind, final1, final2, final3, final4,
                 left1, left2, left3, left4,
                 px, py, ni1, nj1);
        }
        else
        {
          scale = 1.0 / x;
          px = (-z*scale + 1.0) / 2.0;
          py = ( y*scale + 1.0) / 2.0;
          //printf("2px=%g %g\n", px, py);
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
          px = (-x*scale + 1.0) / 2.0; // a shifter a droite
          py = ( z*scale + 1.0) / 2.0;
          //printf("3px=%g %g\n", px, py);
          interp(ind, final1, final2, final3, final4,
                 bottom1, bottom2, bottom3, bottom4, 
                 px, py, ni1, nj1);
        }
        else 
        {
          scale = 1.0 / y;
          px = (-x*scale + 1.0) / 2.0;
          py = (-z*scale + 1.0) / 2.0;
          //printf("4px=%g %g\n", px, py);
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
          //printf("5px=%g %g\n", px, py);
          interp(ind, final1, final2, final3, final4,
                 back1, back2, back3, back4, 
                 px, py, ni1, nj1);
        }
        else
        {
          scale = 1.0 / z;
          px = ( x*scale + 1.0) / 2.0;      
          py = ( y*scale + 1.0) / 2.0;
          //printf("6px=%g %g\n", px, py);
          interp(ind, final1, final2, final3, final4,
                 front1, front2, front3, front4, 
                 px, py, ni1, nj1);
        }
      }
    }
  }
  return Py_None;
}

