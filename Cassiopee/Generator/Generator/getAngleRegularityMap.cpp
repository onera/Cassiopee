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

// getAngleRegularityMap

# include "generator.h"

using namespace K_CONST;
using namespace K_FLD;
using namespace K_FUNC;

inline E_Float computeAngle2(
  const E_Float* x, const E_Float* y, const E_Float* z,
  E_Int ind1, E_Int ind2, E_Int ind3
)
{
  E_Float x1 = x[ind1], x2 = x[ind2], x3 = x[ind3];
  E_Float y1 = y[ind1], y2 = y[ind2], y3 = y[ind3];
  E_Float z1 = z[ind1], z2 = z[ind2], z3 = z[ind3];

  E_Float a2 = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
  E_Float b2 = (x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+(z3-z1)*(z3-z1);
  E_Float c2 = (x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2);

  if (a2 < E_GEOM_CUTOFF || b2 < E_GEOM_CUTOFF) // security check
  { 
    return 0.;
  }

  E_Float a = sqrt(a2);
  E_Float b = sqrt(b2);
  E_Float cosalpha = E_max(E_min((a2+b2-c2)/(2.*a*b),1.),-1.); // law of cosines
  E_Float degconst = 180.0 / K_CONST::E_PI;

  return E_abs(acos(cosalpha)*degconst - 180.);
}

// ============================================================================
/* Return angle regularity map */
// ============================================================================
PyObject* K_GENERATOR::getAngleRegularityMap(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int posx, posy, posz;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, 
                               eltType);

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

  E_Int api = f->getApi();
  E_Int nvertex = f->getSize();
  
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
    E_Int ncells = im1*jm1*km1;

    E_Int ni = 1, nj = 1, nk = 1;
    if (im == 1)
    {
      if (jm == 1) {dim = 1; ni = km;}
      else if (km == 1) {dim = 1; ni = jm;}
      else {dim = 2; ni = jm; nj = km;}
    }
    else if (jm == 1)
    {
      if (im == 1) {dim = 1; ni = km;}
      else if (km == 1) {dim = 1; ni = im;}
      else {dim = 2; ni = im; nj = km;}
    } 
    else if (km == 1)
    {
      if (im == 1) {dim = 1; ni = jm;}
      else if (jm == 1) {dim = 1; ni = im;}
      else {dim = 2; ni = im; nj = jm;}
    }
    else
    {
      ni = im; nj = jm; nk = km;
    }
    E_Int ni1 = E_max(1, E_Int(ni)-1);
    E_Int nj1 = E_max(1, E_Int(nj)-1);
    E_Int nk1 = E_max(1, E_Int(nk)-1);
    E_Int ni1nj1 = ni1*nj1;
    E_Int ninj = ni*nj;
    
    // Construction du tableau numpy stockant les angles 
    // definissant l'orthogonalite
    tpl = K_ARRAY::buildArray3(1, "regularityAngle", im1, jm1, km1, api);
    // pointeur sur le tableau d'angle
    FldArrayF* f2;
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* alphamax = f2->begin(1);

    if (ncells == 1) // mono cell mesh -> early exit
    {
      alphamax[0] = 0.;
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDS(array, f);
      return tpl;
    }
    
    // calcul de l'orthogonalite globale
    #pragma omp parallel
    {
      E_Int ithread = __CURRENT_THREAD__;
      E_Int ind, indn, ind1, ind2, ind3;
      E_Int iprev, jprev, kprev;
      E_Int inext, jnext, knext;
      E_Int inext2, jnext2, knext2;

      if (dim == 1) // dimension = 1D
      {
        #pragma omp for
        for (E_Int i = 0; i < ni1; i++)
        {
          iprev = K_FUNC::E_max(i-1, 0);
          inext = i+1;
          inext2 = K_FUNC::E_min(i+2, ni-1);
          
          alphamax[i] = 0.;
          alphamax[i] = E_max(computeAngle2(x, y, z, i, iprev, inext), alphamax[i]);
          alphamax[i] = E_max(computeAngle2(x, y, z, inext, i, inext2), alphamax[i]);
        }        
      }
      else if (dim == 2) // dimension = 2D
      {
        #pragma omp for collapse(2)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          iprev = K_FUNC::E_max(i-1, 0);
          inext = i+1;
          inext2 = K_FUNC::E_min(i+2, ni-1);

          jprev = K_FUNC::E_max(j-1, 0);
          jnext = j+1;
          jnext2 = K_FUNC::E_min(j+2, nj-1);

          ind  = j*ni1+i;
          indn = j*ni+i;
          alphamax[ind] = 0.;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = j*ni+inext;
          ind3 = j*ni+inext2;
          alphamax[ind] = E_max(computeAngle2(x, y, z, ind2, ind1, ind3), alphamax[ind]);
          
          // calcul de l'angle | (i-1)->i
          ind1 = j*ni+iprev;
          ind2 = indn;
          ind3 = j*ni+inext;
          alphamax[ind] = E_max(computeAngle2(x, y, z, ind2, ind1, ind3), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = jnext*ni+i;
          ind3 = jnext2*ni+i;
          alphamax[ind] = E_max(computeAngle2(x, y, z, ind2, ind1, ind3), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = jprev*ni+i;
          ind2 = indn;
          ind3 = jnext*ni+i;
          alphamax[ind] = E_max(computeAngle2(x, y, z, ind2, ind1, ind3), alphamax[ind]);
        }
      }
      else  // dimension = 3D
      {
        #pragma omp for collapse(3)
        for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          iprev = K_FUNC::E_max(i-1, 0);
          inext = i+1;
          inext2 = K_FUNC::E_min(i+2, ni-1);

          jprev = K_FUNC::E_max(j-1, 0);
          jnext = j+1;
          jnext2 = K_FUNC::E_min(j+2, nj-1);

          kprev = K_FUNC::E_max(k-1, 0);
          knext = k+1;
          knext2 = K_FUNC::E_min(k+2, nk-1);

          ind  = k*ni1nj1+j*ni1+i;
          indn = k*ninj+j*ni+i;
          alphamax[ind] = 0.;

          // calcul de l'angle | i->(i+1)
          ind1 = indn;
          ind2 = k*ninj+j*ni+inext;
          ind3 = k*ninj+j*ni+inext2;
          alphamax[ind] = E_max(computeAngle2(x, y, z, ind2, ind1, ind3), alphamax[ind]);
          
          // calcul de l'angle | (i-1)->i
          ind1 = k*ninj+j*ni+iprev;
          ind2 = indn;
          ind3 = k*ninj+j*ni+inext;
          alphamax[ind] = E_max(computeAngle2(x, y, z, ind2, ind1, ind3), alphamax[ind]);

          // calcul de l'angle | j->(j+1)
          ind1 = indn;
          ind2 = k*ninj+jnext*ni+i;
          ind3 = k*ninj+jnext2*ni+i;
          alphamax[ind] = E_max(computeAngle2(x, y, z, ind2, ind1, ind3), alphamax[ind]);

          // calcul de l'angle | (j-1)->j
          ind1 = k*ninj+jprev*ni+i;
          ind2 = indn;
          ind3 = k*ninj+jnext*ni+i;
          alphamax[ind] = E_max(computeAngle2(x, y, z, ind2, ind1, ind3), alphamax[ind]);

          // calcul de l'angle | k->(k+1)
          ind1 = indn;
          ind2 = knext*ninj+j*ni+i;
          ind3 = knext2*ninj+j*ni+i;
          alphamax[ind] = E_max(computeAngle2(x, y, z, ind2, ind1, ind3), alphamax[ind]);

          // calcul de l'angle | (k-1)->k
          ind1 = kprev*ninj+j*ni+i;
          ind2 = indn;
          ind3 = knext*ninj+j*ni+i;
          alphamax[ind] = E_max(computeAngle2(x, y, z, ind2, ind1, ind3), alphamax[ind]);
        }
      }
    }
    RELEASESHAREDS(tpl, f2);
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else // if (res == 2)
  {
    PyObject* tpl = K_ARRAY::buildArray3(1, "regularityAngle", nvertex, *cn, eltType, 1, api, true);
    FldArrayF* f2;
    K_ARRAY::getFromArray3(tpl, f2);
    f2->setAllValuesAtNull();
    RELEASESHAREDS(tpl, f2);
    return tpl;
  }
}
