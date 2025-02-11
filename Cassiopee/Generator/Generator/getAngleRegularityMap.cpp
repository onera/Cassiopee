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

// getAngleRegularityMap

# include "generator.h"

using namespace K_CONST;
using namespace K_FLD;
using namespace K_FUNC;

inline E_Float computeAngle(E_Float x1,E_Float y1,E_Float z1,E_Float x2,E_Float y2,E_Float z2,E_Float x3,E_Float y3,E_Float z3) {
  E_Float a2, b2, c2, a, b, d;
  E_Float pi = 4*atan(1.);
  E_Float degconst = 180.0 / pi;
  a2 = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
  b2 = (x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2);
  c2 = (x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+(z3-z1)*(z3-z1);
  if ((a2 == 0) || (b2==0) || (c2==0))
  {
    return 0.;
  }
  else
  {
    a = sqrt(a2);
    b = sqrt(b2);
    d = E_max(E_min((a2+b2-c2)/(2.*a*b),1.),-1.);       
    return E_abs(acos(d)*degconst - 180.);
  }
}

// ============================================================================
/* Return angle regularity map */
// ============================================================================
PyObject* K_GENERATOR::getAngleRegularityMap(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int posx, posy, posz;
  E_Int res;
  res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, 
                              eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getAngleOrthogonalityMap: unknown type of array.");
    return NULL;
  }

  PyObject* tpl = NULL;
  
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "getAngleOrthogonalityMap: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // pointeurs associes aux coordonnees
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  
  if (res == 1) // cas structure
  {
    // Dimension du tableau
    E_Int dim = 3;
    E_Int im1 = im-1;
    E_Int jm1 = jm-1;
    E_Int km1 = km-1;
    if (im == 1) im1 = 1;
    if (jm == 1) jm1 = 1;
    if (km == 1) km1 = 1;
    // Direction constante (en 2D)
    // 0: 3D
    // 1: 2D i constant
    // 2: 2D j constant
    // 3: 2D k constant
    E_Int dir = 0;
    if (im == 1)
    {
      dim = 2; dir=1;
      if ((jm ==1)||(km == 1)) dim = 1;
    }
    else if (jm == 1)
    {
      dim = 2; dir = 2;
      if ((im ==1)||(km == 1)) dim = 1;
    } 
    else if (km == 1)
    {
      dim = 2; dir = 3;
      if ((im ==1)||(jm == 1)) dim = 1;
    }
    
    // Construction du tableau numpy stockant les angles 
    // definissant l'orthogonalite
    tpl = K_ARRAY::buildArray(1, "regularityAngle", im1, jm1, km1);
    // pointeur sur le tableau d'angle
    E_Float* alphamax = K_ARRAY::getFieldPtr(tpl);
    //E_Int ncells = im1*jm1*km1;
    // FldArrayF alphamax(ncells, 1, alpha, true);
    
    // calcul de l'orthogonalite globale
    #pragma omp parallel
    {
      E_Int ithread = __CURRENT_THREAD__;
      if (dim == 1) // dimension = 1D
      {
        E_Int ni = E_max(im, E_max(jm, km));
        E_Int ni1 = ni - 1;

        E_Int ind1, ind2, ind3;

        // Initialisation
        #pragma omp for schedule(static)
        for (E_Int i = 0; i < ni1; i++)
        {
          alphamax[i] = -42.;
        }

        // Corners
        if (ithread == 0)
        {
          // imin --------------------------------------------------
          E_Int i = 0;
          // calcul de l'angle | i->(i+1)
          ind1 = i;
          ind2 = i+1;
          ind3 = i+2;
          alphamax[i] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[i]);

          // imax --------------------------------------------------
          i = (ni1-1);
          // calcul de l'angle | (i-1)->i
          ind1 = i-1;
          ind2 = i;
          ind3 = i+1;
          alphamax[i] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[i]);
        }

        // Global loop
        #pragma omp for schedule(static)
        for (E_Int i = 1; i < ni1-1; i++)
        {
          // calcul de l'angle | i->(i+1)
          ind1 = i;
          ind2 = i+1;
          ind3 = i+2;
          alphamax[i] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[i]);
          
          // calcul de l'angle | (i-1)->i
          ind1 = i-1;
          ind2 = i;
          ind3 = i+1;
          alphamax[i] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[i]);
        }
        
      }
      else if (dim == 2) // dimension = 2D
      {
        // (ni,nj) : dimensions of 2D-grid 
        E_Int ni=1, nj=1;
        switch (dir)
        {
          case 3: // k constant
            ni = im; nj = jm;
            break;
          case 2: // j constant
            ni = im; nj = km;
            break;
          case 1: // i constant
            ni = jm; nj = km;
            break;
        }

        E_Int ni1 = ni - 1; E_Int nj1 = nj - 1;
        if (ni == 1) ni1 = 1;
        if (nj == 1) nj1 = 1;

        E_Int ind, indn, ind1, ind2, ind3;

        // Initialisation
        for (E_Int j = 0; j < nj1; j++)
        {
          #pragma omp for schedule(static)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind = j*ni1+i;
            alphamax[ind] = -42.;
          }
        }

        // Corners
        if (ithread == 0)
        {
          // imin, jmin --------------------------------------------------
          ind  = 0*ni1+0;
          indn = 0*ni+0;
          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax, jmin --------------------------------------------------
          ind  = 0*ni1+(ni1-1);
          indn = 0*ni+(ni1-1);
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imin, jmax --------------------------------------------------
          ind  = (nj1-1)*ni1+0;
          indn = (nj1-1)*ni+0;
          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax, jmax --------------------------------------------------
          ind  = (nj1-1)*ni1+(ni1-1);
          indn = (nj1-1)*ni+(ni1-1);
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn; 
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
        }

        // Edges
        #pragma omp for schedule(static)
        for (E_Int j = 1; j < nj1-1; j++) // imin | imax
        {
          // imin --------------------------------------------------
          ind  = j*ni1+0;
          indn = j*ni+0;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax --------------------------------------------------
          ind  = j*ni1+(ni1-1);
          indn = j*ni+(ni1-1);

          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
        }

        #pragma omp for schedule(static)
        for (E_Int i = 1; i < ni1-1; i++) // jmin | jmax
        {
          // jmin --------------------------------------------------
          ind  = 0*ni1+i;
          indn = 0*ni+i;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // jmax --------------------------------------------------
          ind  = (nj1-1)*ni1+i;
          indn = (nj1-1)*ni+i;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
        }

        // Global loop
        for (E_Int j = 1; j < nj1-1; j++)
        {
          #pragma omp for schedule(static)
          for (E_Int i = 1; i < ni1-1; i++)
          {
            ind  = j*ni1+i;
            indn = j*ni+i;

            // calcul de l'angle | i->(i+1)
            ind1 = indn;
            ind2 = indn+1;
            ind3 = indn+2;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
            
            // calcul de l'angle | (i-1)->i
            ind1 = indn-1;
            ind2 = indn;
            ind3 = indn+1;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | j->(j+1)
            ind1 = indn;
            ind2 = indn+ni;
            ind3 = indn+2*ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (j-1)->j
            ind1 = indn-ni;
            ind2 = indn;
            ind3 = indn+ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          }
        }
      }
      else  // dimension = 3D
      {
        // (ni,nj,nk) : dimensions of 3D-grid 
        E_Int ni = im;
        E_Int nj = jm;
        E_Int nk = km;
        E_Int ni1 = ni - 1; E_Int nj1 = nj - 1; E_Int nk1 = nk - 1;
        if (ni == 1) ni1 = 1;
        if (nj == 1) nj1 = 1;
        if (nk == 1) nk1 = 1;

        E_Int ind, indn, ind1, ind2, ind3;

        // Initialisation
        for (E_Int k = 0; k < nk1; k++)
        {
          for (E_Int j = 0; j < nj1; j++)
          {
            #pragma omp for schedule(static)
            for (E_Int i = 0; i < ni1; i++)
            {
              ind = k*ni1*nj1+j*ni1+i;
              alphamax[ind] = -42.;
            }
          }
        }

        if (ithread == 0)
        {
          // Corners
          // imin, jmin, kmin --------------------------------------------------
          ind  = 0*ni1*nj1+0*ni1+0;
          indn = 0*ni*nj+0*ni+0;
          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imin, jmin, kmax --------------------------------------------------
          ind  = (nk1-1)*ni1*nj1+0*ni1+0;
          indn = (nk1-1)*ni*nj+0*ni+0;
          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imin, jmax, kmin --------------------------------------------------
          ind  = 0*ni1*nj1+(nj1-1)*ni1+0;
          indn = 0*ni*nj+(nj1-1)*ni+0;
          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imin, jmax, kmax --------------------------------------------------
          ind  = (nk1-1)*ni1*nj1+(nj1-1)*ni1+0;
          indn = (nk1-1)*ni*nj+(nj1-1)*ni+0;
          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax, jmin, kmin --------------------------------------------------
          ind  = 0*ni1*nj1+0*ni1+(ni1-1);
          indn = 0*ni*nj+0*ni+(ni1-1);      
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax, jmin, kmax --------------------------------------------------
          ind  = (nk1-1)*ni1*nj1+0*ni1+(ni1-1);
          indn = (nk1-1)*ni*nj+0*ni+(ni1-1);      
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax, jmax, kmin --------------------------------------------------
          ind  = 0*ni1*nj1+(nj1-1)*ni1+(ni1-1);
          indn = 0*ni*nj+(nj1-1)*ni+(ni1-1);
          
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax, jmax, kmax
          ind  = (nk1-1)*ni1*nj1+(nj1-1)*ni1+(ni1-1);
          indn = (nk1-1)*ni*nj+(nj1-1)*ni+(ni1-1);      
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
        }

        // Edges
        #pragma omp for schedule(static)
        for (E_Int k = 1; k < nk1-1; k++) // imin | imax and jmin | jmax
        {
          // imin and jmin --------------------------------------------------
          ind  = k*ni1*nj1+0*ni1+0;
          indn = k*ni*nj+0*ni+0;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax and jmin --------------------------------------------------
          ind  = k*ni1*nj1+0*ni1+(ni1-1);
          indn = k*ni*nj+0*ni+(ni1-1);
          
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imin and jmax --------------------------------------------------
          ind  = k*ni1*nj1+(nj1-1)*ni1+0;
          indn = k*ni*nj+(nj1-1)*ni+0;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax and jmax --------------------------------------------------
          ind  = k*ni1*nj1+(nj1-1)*ni1+(ni1-1);
          indn = k*ni*nj+(nj1-1)*ni+(ni1-1);
          
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
        }
        
        #pragma omp for schedule(static)
        for (E_Int j = 1; j < nj1-1; j++) // imin | imax and kmin | kmax
        {
          // imin and kmin --------------------------------------------------
          ind  = 0*ni1*nj1+j*ni1+0;
          indn = 0*ni*nj+j*ni+0;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax and kmin --------------------------------------------------
          ind  = 0*ni1*nj1+j*ni1+(ni1-1);
          indn = 0*ni*nj+j*ni+(ni1-1);

          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imin and kmax --------------------------------------------------
          ind  = (nk1-1)*ni1*nj1+j*ni1+0;
          indn = (nk1-1)*ni*nj+j*ni+0;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // imax and kmax --------------------------------------------------
          ind  = (nk1-1)*ni1*nj1+j*ni1+(ni1-1);
          indn = (nk1-1)*ni*nj+j*ni+(ni1-1);

          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
        }

        #pragma omp for schedule(static)
        for (E_Int i = 1; i < ni1-1; i++) // jmin | jmax and kmin | kmax
        {
          // jmin and kmin --------------------------------------------------
          ind  = 0*ni1*nj1+0*ni1+i;
          indn = 0*ni*nj+0*ni+i;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // jmax and kmin --------------------------------------------------
          ind  = 0*ni1*nj1+(nj1-1)*ni1+i;
          indn = 0*ni*nj+(nj1-1)*ni+i;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = indn+ni*nj;
          ind3 = indn+2*ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // jmin and kmax --------------------------------------------------
          ind  = (nk1-1)*ni1*nj1+0*ni1+i;
          indn = (nk1-1)*ni*nj+0*ni+i;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = indn+ni;
          ind3 = indn+2*ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // jmax and kmax --------------------------------------------------
          ind  = (nk1-1)*ni1*nj1+(nj1-1)*ni1+i;
          indn = (nk1-1)*ni*nj+(nj1-1)*ni+i;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = indn+1;
          ind3 = indn+2;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          
          // calcul de l'angle | (i-1)->i
          ind1 = indn-1;
          ind2 = indn;
          ind3 = indn+1;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = indn-ni;
          ind2 = indn;
          ind3 = indn+ni;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = indn-ni*nj;
          ind2 = indn;
          ind3 = indn+ni*nj;
          alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
        }

        // Faces
        for (E_Int k = 1; k < nk1-1; k++) // imin | imax
        {
          #pragma omp for schedule(static)
          for (E_Int j = 1; j < nj1-1; j++)
          {
            // imin --------------------------------------------------
            ind  = k*ni1*nj1+j*ni1+0;
            indn = k*ni*nj+j*ni+0;

            // calcul de l'angle | i->(i+1)
            ind1 = indn;
            ind2 = indn+1;
            ind3 = indn+2;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | j->(j+1)
            ind1 = indn;
            ind2 = indn+ni;
            ind3 = indn+2*ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (j-1)->j
            ind1 = indn-ni;
            ind2 = indn;
            ind3 = indn+ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | k->(k+1)
            ind1 = indn;
            ind2 = indn+ni*nj;
            ind3 = indn+2*ni*nj;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (k-1)->k
            ind1 = indn-ni*nj;
            ind2 = indn;
            ind3 = indn+ni*nj;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // imax --------------------------------------------------
            ind  = k*ni1*nj1+j*ni1+(ni1-1);
            indn = k*ni*nj+j*ni+(ni1-1);

            // calcul de l'angle | (i-1)->i
            ind1 = indn-1;
            ind2 = indn;
            ind3 = indn+1;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | j->(j+1)
            ind1 = indn;
            ind2 = indn+ni;
            ind3 = indn+2*ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (j-1)->j
            ind1 = indn-ni;
            ind2 = indn;
            ind3 = indn+ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | k->(k+1)
            ind1 = indn;
            ind2 = indn+ni*nj;
            ind3 = indn+2*ni*nj;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (k-1)->k
            ind1 = indn-ni*nj;
            ind2 = indn;
            ind3 = indn+ni*nj;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          }
        }

        for (E_Int k = 1; k < nk1-1; k++) // jmin | jmax
        {
          #pragma omp for schedule(static)
          for (E_Int i = 1; i < ni1-1; i++)
          {
            // jmin --------------------------------------------------
            ind  = k*ni1*nj1+0*ni1+i;
            indn = k*ni*nj+0*ni+i;

            // calcul de l'angle | i->(i+1)
            ind1 = indn;
            ind2 = indn+1;
            ind3 = indn+2;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (i-1)->i
            ind1 = indn-1;
            ind2 = indn;
            ind3 = indn+1;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | j->(j+1)
            ind1 = indn;
            ind2 = indn+ni;
            ind3 = indn+2*ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | k->(k+1)
            ind1 = indn;
            ind2 = indn+ni*nj;
            ind3 = indn+2*ni*nj;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (k-1)->k
            ind1 = indn-ni*nj;
            ind2 = indn;
            ind3 = indn+ni*nj;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // jmax --------------------------------------------------
            ind  = k*ni1*nj1+(nj1-1)*ni1+i;
            indn = k*ni*nj+(nj1-1)*ni+i;

            // calcul de l'angle | i->(i+1)
            ind1 = indn;
            ind2 = indn+1;
            ind3 = indn+2;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (i-1)->i
            ind1 = indn-1;
            ind2 = indn;
            ind3 = indn+1;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (j-1)->j
            ind1 = indn-ni;
            ind2 = indn;
            ind3 = indn+ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | k->(k+1)
            ind1 = indn;
            ind2 = indn+ni*nj;
            ind3 = indn+2*ni*nj;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (k-1)->k
            ind1 = indn-ni*nj;
            ind2 = indn;
            ind3 = indn+ni*nj;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          }
        }

        for (E_Int j = 1; j < nj1-1; j++) // kmin | kmax
        {
          #pragma omp for schedule(static)
          for (E_Int i = 1; i < ni1-1; i++)
          {
            // kmin --------------------------------------------------
            ind  = 0*ni1*nj1+j*ni1+i;
            indn = 0*ni*nj+j*ni+i;

            // calcul de l'angle | i->(i+1)
            ind1 = indn;
            ind2 = indn+1;
            ind3 = indn+2;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (i-1)->i
            ind1 = indn-1;
            ind2 = indn;
            ind3 = indn+1;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | j->(j+1)
            ind1 = indn;
            ind2 = indn+ni;
            ind3 = indn+2*ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (j-1)->j
            ind1 = indn-ni;
            ind2 = indn;
            ind3 = indn+ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | k->(k+1)
            ind1 = indn;
            ind2 = indn+ni*nj;
            ind3 = indn+2*ni*nj;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // kmax --------------------------------------------------
            ind  = (nk1-1)*ni1*nj1+j*ni1+i;
            indn = (nk1-1)*ni*nj+j*ni+i;

            // calcul de l'angle | i->(i+1)
            ind1 = indn;
            ind2 = indn+1;
            ind3 = indn+2;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (i-1)->i
            ind1 = indn-1;
            ind2 = indn;
            ind3 = indn+1;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | j->(j+1)
            ind1 = indn;
            ind2 = indn+ni;
            ind3 = indn+2*ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (j-1)->j
            ind1 = indn-ni;
            ind2 = indn;
            ind3 = indn+ni;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

            // calcul de l'angle | (k-1)->k
            ind1 = indn-ni*nj;
            ind2 = indn;
            ind3 = indn+ni*nj;
            alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
          }
        }

        // Global loop
        for (E_Int k = 1; k < nk1-1; k++)
        {
          for (E_Int j = 1; j < nj1-1; j++)
          {
            #pragma omp for schedule(static)
            for (E_Int i = 1; i < ni1-1; i++)
            {
              ind  = k*ni1*nj1+j*ni1+i;
              indn = k*ni*nj+j*ni+i;

              // calcul de l'angle | i->(i+1)
              ind1 = indn;
              ind2 = indn+1;
              ind3 = indn+2;
              alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
              
              // calcul de l'angle | (i-1)->i
              ind1 = indn-1;
              ind2 = indn;
              ind3 = indn+1;
              alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

              // calcul de l'angle | j->(j+1)
              ind1 = indn;
              ind2 = indn+ni;
              ind3 = indn+2*ni;
              alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

              // calcul de l'angle | (j-1)->j
              ind1 = indn-ni;
              ind2 = indn;
              ind3 = indn+ni;
              alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

              // calcul de l'angle | k->(k+1)
              ind1 = indn;
              ind2 = indn+ni*nj;
              ind3 = indn+2*ni*nj;
              alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);

              // calcul de l'angle | (k-1)->k
              ind1 = indn-ni*nj;
              ind2 = indn;
              ind3 = indn+ni*nj;
              alphamax[ind] = E_max(computeAngle(x[ind1],y[ind1],z[ind1],x[ind2],y[ind2],z[ind2],x[ind3],y[ind3],z[ind3]), alphamax[ind]);
            }
          }
        }
      }
    }
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else // if (res == 2)
  {
    if (strcmp(eltType, "NGON") != 0) // Elements basiques
    { 
      E_Int nelts = cn->getSize(); // nb d'elements
      tpl = K_ARRAY::buildArray(1, "regularityAngle", f->getSize(), nelts, -1, eltType, true);
      E_Float* fieldp = K_ARRAY::getFieldPtr(tpl);
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      FldArrayI cnn(nelts, cn->getNfld(), cnnp, true); cnn = *cn;
      for (E_Int i = 0; i < nelts; i++) fieldp[i] = 0.;
    }
    else
    {
      E_Int* cnp = cn->begin(); // pointeur sur la connectivite NGon
      E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
      E_Int nelts = cnp[sizeFN+2];  // nombre total d elements
      E_Int npts = f->getSize();
      tpl = K_ARRAY::buildArray(1, "regularityAngle",
                                npts, nelts,
                                -1, eltType, true, cn->getSize());
      E_Float* fieldp = K_ARRAY::getFieldPtr(tpl);
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      FldArrayI cnn(cn->getSize(), 1, cnnp, true); cnn = *cn;
      for (E_Int i = 0; i < nelts; i++) fieldp[i] = 0.;
    } 
    
    return tpl;
    //RELEASESHAREDB(res, array, f, cn);
    //PyErr_SetString(PyExc_ValueError,
    //                "getAngleOrthogonalityMap: not for unstructured arrays.");
    //return NULL;
  }
}
