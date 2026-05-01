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

// getNonOrthogonalityMap

# include "generator.h"
# include <math.h>

using namespace K_FLD;
using namespace K_FUNC;

// ============================================================================
/* Return angle regularity map */
// ============================================================================
PyObject* K_GENERATOR::getNonOrthogonalityMap(PyObject* self, PyObject* args)
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
                    "getNonOrthogonalityMap: unknown type of array.");
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
                    "getNonOrthogonalityMap: can't find coordinates in array.");
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
    tpl = K_ARRAY::buildArray3(1, "nonOrthogonality", im1, jm1, km1, api);
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
          alphamax[i] = 0.;
        }        
      }
      else if (dim == 2) // dimension = 2D
      {
        #pragma omp for collapse(2)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          alphamax[i] = 0.;
        }
      }
      else  // dimension = 3D
      {
        #pragma omp for collapse(3)
        for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          alphamax[i] = 0.;
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
      1, "nonOrthogonality", npts, *cn, eltType, true, api, true
    );
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* alphamax = f2->begin(1); // pointeur sur le tableau d'angle

    if (strcmp(eltType, "NGON") == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getNonOrthogonalityMap: not implemented for NGON arrays.");
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
                      "getNonOrthogonalityMap: Error computing nfpe.");
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

    // Get dimensionality
    E_Int dim = K_CONNECT::getDimME(eltType);

    if (dim == 1) // dim 1 -> early exit
    {
      #pragma omp parallel
      {
        #pragma omp for
        for (E_Int i = 0; i < ntotElts; i++)
        {
          alphamax[i] = 0.;
        }
      }
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
      return tpl;
    } 

    // Calcul de la connectivite face -> elements
    FldArrayI cFE(ntotFacets,2);
    ierr = K_CONNECT::connectEV2FE(eltType, *cn, cFE, true);
    if (ierr == 1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getNonOrthogonalityMap: Error computing cFE.");
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }

    // Compute center of elements
    FldArrayF fxb(ntotElts), fyb(ntotElts), fzb(ntotElts);
    K_METRIC::compUnstructCellCenter(*cn, x, y, z,
      fxb.begin(), fyb.begin(), fzb.begin()
    );
    E_Float* xb = fxb.begin(1);
    E_Float* yb = fyb.begin(1);
    E_Float* zb = fzb.begin(1);

    if (dim == 3)
    {
      // Compute facet normals
      FldArrayF fsnx(ntotFacets), fsny(ntotFacets), fsnz(ntotFacets), surf(ntotFacets);
      K_METRIC::compSurfUnstruct(*cn, eltType, x, y, z,
        fsnx.begin(), fsny.begin(), fsnz.begin(), surf.begin()
      );
      E_Float* snx = fsnx.begin(1);
      E_Float* sny = fsny.begin(1);
      E_Float* snz = fsnz.begin(1);

      E_Float degconst = 180.0 / K_CONST::E_PI;

      // calcul de la regularite
      #pragma omp parallel
      {
        E_Int nelts, ind, pos;
        E_Int ind2;
        E_Int elOffset, fctOffset;
        E_Float x1, y1, z1;
        E_Float x2, y2, z2;
        E_Float ux, uy, uz;
        E_Float vx, vy, vz;
        E_Float normu, normv;
        E_Float cosalpha, alpha, skewness;
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
              ind2 = cFE(pos, 2) - 1;
              if (ind2 < 0) continue; // facet has only one neighbor element
              x2 = xb[ind2]; y2 = yb[ind2]; z2 = zb[ind2]; // neighbor element center
              
              ux = snx[pos]; uy = sny[pos]; uz = snz[pos]; // face normal
              vx = x2-x1; vy = y2-y1; vz = z2-z1; // centroid to centroid vector

              normu = sqrt(ux*ux + uy*uy + uz*uz);
              normv = sqrt(vx*vx + vy*vy + vz*vz);

              cosalpha = K_FUNC::E_max(K_FUNC::E_min((ux*vx + uy*vy + uz*vz)/(normu*normv),1.),-1.);
              alpha = acos(cosalpha);
              skewness = K_FUNC::E_abs(alpha*degconst);

              if (false)
              {
                printf("center 1: %f/%f/%f\n", x1,y1,z1);
                printf("center 2: %f/%f/%f\n", x2,y2,z2);
                printf("u : %f/%f/%f\n", ux,uy,uz);
                printf("v : %f/%f/%f\n", vx,vy,vz);
                printf("alpha : %f (degree: %f)\n", alpha, alpha*degconst);
              }

              alphamax[ind] = E_max(skewness, alphamax[ind]);
            }
          }
        }
      }
    }
    else // dim == 2
    {
      // Compute QUAD/TRI element normals
      FldArrayF fsnx(ntotElts), fsny(ntotElts), fsnz(ntotElts), surf(ntotElts);
      K_METRIC::compSurfUnstruct(*cn, eltType, x, y, z,
        fsnx.begin(), fsny.begin(), fsnz.begin(), surf.begin()
      );
      E_Float* snx = fsnx.begin(1);
      E_Float* sny = fsny.begin(1);
      E_Float* snz = fsnz.begin(1);

      E_Float degconst = 180.0 / K_CONST::E_PI;

      // get eltTypes
      std::vector<char*> eltTypes;
      K_ARRAY::extractVars(eltType, eltTypes);

      // calcul de la regularite
      #pragma omp parallel
      {
        E_Int nelts, ind, pos;
        E_Int ind1, ind2;
        E_Int iver1, iver2;
        E_Int elOffset, fctOffset;
        E_Float x1, y1, z1;
        E_Float x2, y2, z2;
        E_Float xver1, yver1, zver1;
        E_Float xver2, yver2, zver2;
        E_Float ux, uy, uz;
        E_Float vx, vy, vz;
        E_Float nx, ny, nz;
        E_Float ex, ey, ez;
        E_Float normu, normv;
        E_Float cosalpha, alpha, skewness;
        std::vector<std::vector<E_Int> > facets;
        // loop over all connectivities
        for (E_Int ic = 0; ic < nc; ic++)
        {
          FldArrayI& cm = *(cn->getConnect(ic));
          nelts = cm.getSize();
          elOffset = nepc[ic];
          fctOffset = nfpc[ic];
          K_CONNECT::getEVFacets(facets, eltTypes[ic], false, true);
          
          // loop over all elements of connectivity cm
          #pragma omp for
          for (E_Int i = 0; i < nelts; i++)
          {
            ind = i + elOffset; // true element index
            x1 = xb[ind]; y1 = yb[ind]; z1 = zb[ind]; // cell center
            nx = snx[ind]; ny = sny[ind]; nz = snz[ind]; // face normal
            alphamax[ind] = 0.; // initialization

            // loop over all faces of element i
            for (E_Int f = 0; f < nfpe[ic]; f++)
            {
              pos = f + i*nfpe[ic] + fctOffset;
              ind2 = cFE(pos, 2) - 1;
              if (ind2 < 0) continue; // facet has only one neighbor element
              x2 = xb[ind2]; y2 = yb[ind2]; z2 = zb[ind2]; // neighbor element center

              iver1 = cm(i,facets[f][0])-1;
              iver2 = cm(i,facets[f][1])-1;
              xver1 = x[iver1]; yver1 = y[iver1]; zver1 = z[iver1]; // first edge vertex
              xver2 = x[iver2]; yver2 = y[iver2]; zver2 = z[iver2]; // second edge vertex

              ex = xver2-xver1; ey = yver2-yver1; ez = zver2-zver1; // edge vector
              ux = ey*nz-ez*ny; uy = ez*nx-ex*nz; uz = ex*ny-ey*nx; // edge normal vector = e x n
              vx = x2-x1; vy = y2-y1; vz = z2-z1; // centroid to centroid vector

              normu = sqrt(ux*ux + uy*uy + uz*uz);
              normv = sqrt(vx*vx + vy*vy + vz*vz);

              cosalpha = K_FUNC::E_max(K_FUNC::E_min((ux*vx + uy*vy + uz*vz)/(normu*normv),1.),-1.);
              alpha = acos(cosalpha);
              skewness = K_FUNC::E_abs(alpha*degconst);
              
              if (false)
              {
                printf("center 1: %f/%f/%f\n", x1,y1,z1);
                printf("center 2: %f/%f/%f\n", x2,y2,z2);
                printf("n : %f/%f/%f\n", nx,ny,nz);
                printf("e : %f/%f/%f\n", ex,ey,ez);
                printf("u : %f/%f/%f\n", ux,uy,uz);
                printf("v : %f/%f/%f\n", vx,vy,vz);
                printf("alpha : %f (degree: %f)\n", alpha, alpha*degconst);
              }

              alphamax[ind] = E_max(skewness, alphamax[ind]);
            }
          }
        }
      }
      for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
    }
    
    RELEASESHAREDS(tpl, f2);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
}
