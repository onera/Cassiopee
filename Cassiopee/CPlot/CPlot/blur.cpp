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
#include <math.h>

#define M_PI 3.1415926535897932384626433832795

//=======================================================================
// Create 1D half gaussian kernel coefficients
// IN: sigma: in pixels
// IN: n: half size of kernel
// OUT: c: half kernel coefficients
//=======================================================================
void createGaussFilter(E_Float sigma, E_Int n, E_Float* c)
{
  E_Float sigma2 = (sigma*sigma);
  for (E_Int i = 0; i <= n; i++)
  {
    c[i] = exp( -(i*i) / (2.*sigma2)) / (sigma*sqrt(2.*M_PI));
  }

  /*
  printf("gaussian coefficients: ");
  for (E_Int i = 0; i <= n; i++) printf("%g ", c[i]);
  printf("\n");
  E_Float sum = 0.;
  sum = c[0];
  for (E_Int i = 1; i <= n; i++)
  {
    sum += c[i];
    sum += c[i];
  }
  printf("sum=%g\n", sum);
  */
}

//=======================================================================
// Apply gaussian blur to in (single scalar)
// IN: in: color array
// IN: ni,nj: image size
// IN: c: kernel coef
// IN: n: kernel half size
// OUT: out: color array (already allocated)
// Apply 1D filter in both directions
//=======================================================================
void gaussianBlur2(E_Float* in, E_Int ni, E_Int nj, E_Float* c, E_Int n, E_Float* out)
{
  E_Int ind;

  // filter
  for (E_Int j = 0; j < nj; j++)
  for (E_Int i = 0; i < ni; i++)
  {
    ind = i+j*ni;
    out[ind] = in[ind]*c[0];
  }

  // filter en i
  for (E_Int j = 0; j < nj; j++)
  for (E_Int i = n; i < ni-n; i++)
  {
    ind = i+j*ni;
    for (E_Int k = 1; k <= n; k++)
    {
      out[ind] += in[ind-k]*c[k]; // right
      out[ind] += in[ind+k]*c[k]; // left
    }
  }

  // filter en j
  for (E_Int j = n; j < nj-n; j++)
  for (E_Int i = 0; i < ni; i++)
  {
    ind = i+j*ni;
    for (E_Int k = 1; k <= n; k++)
    {
      out[ind] += in[ind-k*ni]*c[k];
      out[ind] += in[ind+k*ni]*c[k];
    }
  }
}

//=======================================================================
// blur array, python function
//=======================================================================
PyObject* K_CPLOT::blur(PyObject* self, PyObject* args)
{
  // Get the 4 arrays of cube images (left, right, bottom, top, back, front)
  PyObject* array; E_Float blurSigma;
  if (!PYPARSETUPLE_(args, O_ R_, &array, &blurSigma))
  {
    return NULL;
  }
  char* varString;
  E_Int ni, nj, nk, res;
  FldArrayF* im; FldArrayI* cn; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString, im, 
                               ni, nj, nk, cn, eltType);
  if (res != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "blur: requires a structured array.");
    return NULL;
  }

  E_Float* im1 = im->begin(4); // r,g,b,a
  E_Float* im2 = im->begin(5);
  E_Float* im3 = im->begin(6);
  //E_Float* im4 = im->begin(7);


  // gaussian blur needed for headsets
  if (blurSigma > 0.)
  {
    printf("bluring...\n");
    E_Int n = 3; // half kernel size
    E_Float* c = new E_Float [n+1]; // kernel coeff
    createGaussFilter(blurSigma, n, c);
    FldArrayF out(ni*nj, 3);
    E_Float* out1 = out.begin(1);
    E_Float* out2 = out.begin(2);
    E_Float* out3 = out.begin(3);

    gaussianBlur2(im1, ni, nj, c, n, out1);
    gaussianBlur2(im2, ni, nj, c, n, out2);
    gaussianBlur2(im3, ni, nj, c, n, out3);
    
    #pragma omp parallel
    {
      #pragma omp for
      for (E_Int ind = 0; ind < ni*nj; ind++) 
      {  
        im1[ind] = out1[ind];
        im2[ind] = out2[ind];
        im3[ind] = out3[ind];
      }
    }
  
    delete [] c;
  }

  RELEASESHAREDS(array, im);
  return Py_None;
}
