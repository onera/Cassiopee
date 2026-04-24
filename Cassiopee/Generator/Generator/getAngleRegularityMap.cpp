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

using namespace K_FLD;
using namespace K_FUNC;

// ============================================================================
/* Return angle regularity map */
// ============================================================================
PyObject* K_GENERATOR::getAngleRegularityMap(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int normalized;

  if (!PYPARSETUPLE_(args, O_ I_, &array, &normalized)) return NULL;
  
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

  E_Int api = f->getApi();
  E_Int npts = f->getSize();
  PyObject* tpl = NULL;
  FldArrayF* f2;
  
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
      //E_Int ithread = __CURRENT_THREAD__;
      E_Int ind, indn, ind1, ind2, ind3;
      E_Int iprev, jprev, kprev;
      E_Int inext, jnext, knext;
      E_Int inext2, jnext2, knext2;

      if (dim == 1) // dimension = 1D
      {
        #pragma omp for
        for (E_Int i = 0; i < ni1; i++)
        {
          iprev = E_max(i-1, 0);
          inext = i+1;
          inext2 = E_min(i+2, ni-1);
          
          alphamax[i] = 0.;
          alphamax[i] = E_max(K_COMPGEOM::computeSkewness(x, y, z, iprev, i, inext, 180., normalized), alphamax[i]);
          alphamax[i] = E_max(K_COMPGEOM::computeSkewness(x, y, z, i, inext, inext2, 180., normalized), alphamax[i]);
        }        
      }
      else if (dim == 2) // dimension = 2D
      {
        #pragma omp for collapse(2)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          iprev = E_max(i-1, 0);
          inext = i+1;
          inext2 = E_min(i+2, ni-1);

          jprev = E_max(j-1, 0);
          jnext = j+1;
          jnext2 = E_min(j+2, nj-1);

          ind  = j*ni1+i;
          indn = j*ni+i;
          alphamax[ind] = 0.;

          // angle i->(i+1)
          ind1 = indn;
          ind2 = j*ni+inext;
          ind3 = j*ni+inext2;
          alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 180., normalized), alphamax[ind]);
          
          // angle (i-1)->i
          ind1 = j*ni+iprev;
          ind2 = indn;
          ind3 = j*ni+inext;
          alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 180., normalized), alphamax[ind]);

          // angle j->(j+1)
          ind1 = indn;
          ind2 = jnext*ni+i;
          ind3 = jnext2*ni+i;
          alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 180., normalized), alphamax[ind]);

          // angle (j-1)->j
          ind1 = jprev*ni+i;
          ind2 = indn;
          ind3 = jnext*ni+i;
          alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 180., normalized), alphamax[ind]);
        }
      }
      else  // dimension = 3D
      {
        #pragma omp for collapse(3)
        for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          iprev = E_max(i-1, 0);
          inext = i+1;
          inext2 = E_min(i+2, ni-1);

          jprev = E_max(j-1, 0);
          jnext = j+1;
          jnext2 = E_min(j+2, nj-1);

          kprev = E_max(k-1, 0);
          knext = k+1;
          knext2 = E_min(k+2, nk-1);

          ind  = k*ni1nj1+j*ni1+i;
          indn = k*ninj+j*ni+i;
          alphamax[ind] = 0.;

          // angle i->(i+1)
          ind1 = indn;
          ind2 = k*ninj+j*ni+inext;
          ind3 = k*ninj+j*ni+inext2;
          alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 180., normalized), alphamax[ind]);
          
          // angle (i-1)->i
          ind1 = k*ninj+j*ni+iprev;
          ind2 = indn;
          ind3 = k*ninj+j*ni+inext;
          alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 180., normalized), alphamax[ind]);

          // angle j->(j+1)
          ind1 = indn;
          ind2 = k*ninj+jnext*ni+i;
          ind3 = k*ninj+jnext2*ni+i;
          alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 180., normalized), alphamax[ind]);

          // angle (j-1)->j
          ind1 = k*ninj+jprev*ni+i;
          ind2 = indn;
          ind3 = k*ninj+jnext*ni+i;
          alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 180., normalized), alphamax[ind]);

          // angle k->(k+1)
          ind1 = indn;
          ind2 = knext*ninj+j*ni+i;
          ind3 = knext2*ninj+j*ni+i;
          alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 180., normalized), alphamax[ind]);

          // angle (k-1)->k
          ind1 = kprev*ninj+j*ni+i;
          ind2 = indn;
          ind3 = knext*ninj+j*ni+i;
          alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 180., normalized), alphamax[ind]);
        }
      }
    }
    RELEASESHAREDS(tpl, f2);
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else // if (res == 2)
  {
    PyObject* tpl = K_ARRAY::buildArray3(
      1, "regularityAngle", npts, *cn, eltType, true, api, true
    );
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* alphamax = f2->begin(1); // pointeur sur le tableau d'angle

    if (strcmp(eltType, "NGON") == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getAngleRegumarityMap: not implemented for NGON arrays.");
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }

    E_Int nc = cn->getNConnect();

    // Number of facets per element
    std::vector<E_Int> nfpe;
    E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, true);
    if (ierr != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getAngleRegularityMap: Error computing nfpe.");
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }

    // Total number of elements/facets
    E_Int ntotFacets = 0;
    E_Int ntotElts = 0;
    std::vector<E_Int> nepc(nc+1), nfpc(nc+1);
    nepc[0] = 0; nfpc[0] = 0;
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      E_Int nelts = cm.getSize();
      nepc[ic+1] = nepc[ic] + nelts;
      nfpc[ic+1] = nfpc[ic] + nfpe[ic]*nelts;  // number of facets per connectivity
      ntotElts += nelts;
      ntotFacets += nfpe[ic]*nelts;
    }

    // Calcul de la connectivite face -> elements
    FldArrayI cFE(ntotFacets,2);
    ierr = K_CONNECT::connectEV2FE(eltType, *cn, cFE, true);
    if (ierr == 1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getRegularityMap: Error computing cFE.");
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }

    // Compute center of facets
    FldArrayF fxint(ntotFacets), fyint(ntotFacets), fzint(ntotFacets);
    K_METRIC::compUnstructCenterInt(*cn, eltType, x, y, z,
      fxint.begin(), fyint.begin(), fzint.begin(), true
    );
    E_Float* xint = fxint.begin(1);
    E_Float* yint = fyint.begin(1);
    E_Float* zint = fzint.begin(1);

    // Compute center of elements
    FldArrayF fxb(ntotElts), fyb(ntotElts), fzb(ntotElts);
    K_METRIC::compUnstructCellCenter(*cn, x, y, z,
      fxb.begin(), fyb.begin(), fzb.begin()
    );
    E_Float* xb = fxb.begin(1);
    E_Float* yb = fyb.begin(1);
    E_Float* zb = fzb.begin(1);

    // calcul de la regularite
    #pragma omp parallel
    {
      E_Int nelts, ind, pos;
      E_Int ind2;
      E_Int elOffset, fctOffset;
      E_Float x1, y1, z1;
      E_Float x2, y2, z2;
      E_Float x3, y3, z3;
      // loop over all connectivities
      for (E_Int ic = 0; ic < nc; ic++)
      {
        FldArrayI& cm = *(cn->getConnect(ic));
        nelts = cm.getSize();
        elOffset = nepc[ic];
        fctOffset = nfpc[ic];
        
        // loop over all elements of connectivity cm
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind = i + elOffset; // true element index
          x1 = xb[ind]; y1 = yb[ind]; z1 = zb[ind]; // cell center
          alphamax[ind] = 0.; // initialization

          // loop over all faces of element i
          for (E_Int f = 0; f < nfpe[ic]; f++)
          {
            pos = f + i*nfpe[ic] + fctOffset;
            //ind1 = cFE(pos, 1) - 1;
            ind2 = cFE(pos, 2) - 1;
            if (ind2 < 0) continue; // facet has only one neighbor element
            x2 = xint[pos]; y2 = yint[pos]; z2 = zint[pos]; // facet center
            x3 = xb[ind2]; y3 = yb[ind2]; z3 = zb[ind2]; // neighbor element center

            alphamax[ind] = E_max(K_COMPGEOM::computeSkewness(x1,y1,z1,x2,y2,z2,x3,y3,z3,180.,normalized), alphamax[ind]);
          }
        }
      }
    }
    RELEASESHAREDS(tpl, f2);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
}
