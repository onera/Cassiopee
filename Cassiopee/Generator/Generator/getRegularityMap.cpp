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
# include "generator.h"

// ============================================================================
// Inline functions
// ============================================================================
inline E_Float ratioMax1(E_Float v, E_Float v1)
{
  return K_FUNC::E_abs(v1 - v)/(K_FUNC::E_max(v, K_CONST::E_GEOM_CUTOFF));
}

inline E_Float ratioMax2(E_Float v, E_Float v1, E_Float v2)
{
  return K_FUNC::E_max(ratioMax1(v, v1), ratioMax1(v, v2));
}

inline E_Float ratioMax3(E_Float v, E_Float v1, E_Float v2, E_Float v3)
{
  return K_FUNC::E_max(ratioMax2(v, v1, v2), ratioMax1(v, v3));
}

inline E_Float ratioMax4(E_Float v, E_Float v1, E_Float v2, E_Float v3, E_Float v4)
{
  return K_FUNC::E_max(ratioMax2(v, v1, v2), ratioMax2(v, v3, v4));
}

inline E_Float ratioMax5(E_Float v, E_Float v1, E_Float v2, E_Float v3, E_Float v4, E_Float v5)
{
  return K_FUNC::E_max(ratioMax3(v, v1, v2, v3), ratioMax2(v, v4, v5));
}

inline E_Float ratioMax6(E_Float v, E_Float v1, E_Float v2, E_Float v3, E_Float v4, E_Float v5, E_Float v6)
{
  return K_FUNC::E_max(ratioMax3(v, v1, v2, v3), ratioMax3(v, v4, v5, v6));
}

// ============================================================================
/* Return regularity map */
// ============================================================================
PyObject* K_GENERATOR::getRegularityMap(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int posx, posy, posz;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getRegularityMap: unknown type of array.");
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
                    "getRegularityMap: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // pointeurs associes aux coordonnees
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  
  if (res == 1) // cas structure
  {
    // Calcul et affectation de :
    //    - la dimension du tableau (dim)
    //    - les directions de maillage I,J,K definies (dir)
    //    - les noms des champs de sortie
    E_Int dim = 3;
    E_Int dimC = 3;
    E_Int im1 = im-1;
    E_Int jm1 = jm-1;
    E_Int km1 = km-1;
    E_Int dirI = 2;
    E_Int dirJ = 3;
    E_Int dirK = 4;
    if (im == 1) dirI = 1;
    if (jm == 1) dirJ = 1;
    if (km == 1) dirK = 1;
    E_Int dir = dirI*dirJ*dirK;
    E_Int ni, nj, nk;
    ni = im; nj = jm; nk = km;
    switch (dir)
    {
      case 2: // dim 1 - dir I
        ni = im; dim = 1; dimC = 1;
        break;
      case 3: // dim 1 - dir J
        ni = jm; dim = 1; dimC = 1;
        break;
      case 4: // dim 1 - dir K
        ni = km; dim = 1; dimC = 1;
        break;
      case 6: // dim 2 - dir IJ
        ni = im;
        nj = jm;
        dim = 2;
        dimC = 2;
        if (im == 2) { dimC = 1; ni = jm; nj = 1; }
        if (jm == 2) { dimC = 1; ni = im; nj = 1; }
        break;
      case 8: // dim 2 - dir IK
        ni = im;
        nj = km;
        dim = 2;
        dimC = 2;
        if (im == 2) { dimC = 1; ni = km; nj = 1; }
        if (km == 2) { dimC = 1; ni = im; nj = 1; }
        break;
      case 12: // dim 2 - dir JK
        ni = jm;
        nj = km;
        dim = 2;
        dimC = 2;
        if (im == 2) { dimC = 1; ni = km; nj = 1; }
        if (km == 2) { dimC = 1; ni = jm; nj = 1; }
        break;
      default:
        if (im == 2) { dimC = 2; ni = jm; nj = km; }
        if (jm == 2) { dimC = 2; ni = im; nj = km; }
        if (km == 2) { dimC = 2; ni = im; nj = jm; }
    }

    if (im == 1) im1 = 1;
    if (jm == 1) jm1 = 1;
    if (km == 1) km1 = 1;
    E_Int ncells = im1*jm1*km1;
    E_Int ninti = im*jm1*km1;
    E_Int nintj = im1*jm*km1;
    E_Int nintk = im1*jm1*km;
    E_Int nint =  ninti + nintj + nintk;
    E_Int ni1 = ni-1; E_Int nj1 = nj-1; E_Int nk1 = nk-1;
    if (ni == 1) ni1 = 1;
    if (nj == 1) nj1 = 1;
    if (nk == 1) nk1 = 1;

    // Construction du tableau numpy stockant les champs 
    // definissant la regularite
    tpl = K_ARRAY::buildArray3(1, "regularity", im1, jm1, km1, api);
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* reg = f2->begin(1);
      
    // calcul du volume
    FldArrayF vol(ncells);
    if (dim == 1)
      K_METRIC::compSurfStruct1D(im, jm, km, xp, yp, zp, vol.begin());
    else if (dim == 2)
      K_METRIC::compSurfStruct2D(im, jm, km, xp, yp, zp, vol.begin());
    else
    {
      FldArrayF surf(nint,3);
      FldArrayF snorm(nint);
      FldArrayF centerInt(nint, 3);
      K_METRIC::compMetricStruct(
        im, jm, km, ninti, nintj, nintk,
        xp, yp, zp,
        vol.begin(), surf.begin(1), surf.begin(2), surf.begin(3), 
        snorm.begin(), 
        centerInt.begin(1), centerInt.begin(2), centerInt.begin(3));
    }

    if (ncells == 1) // mono cell mesh -> early exit
    {
      reg[0] = 0.;
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDS(array, f);
      return tpl;
    }

    // calcul de la regularite
    #pragma omp parallel
    {
      E_Int ithread = __CURRENT_THREAD__;
      // variables locales pour les indices
      E_Int ind, ind1, ind2, ind3, ind4, ind5, ind6;
      if (dimC == 1) // dimension 1D
      {
        // Aux frontieres, traitement degenere.
        if (ithread == 0)
        {
          reg[0] = ratioMax1(vol[0], vol[1]);
          reg[ni1-1] = ratioMax1(vol[ni1-1], vol[ni1-2]);
        }

        // Boucle sur les indices
        #pragma omp for
        for (E_Int i = 1; i < ni1-1; i++)
        {
          reg[i] = ratioMax2(vol[i], vol[i-1], vol[i+1]);
        }
      }
      else if (dimC == 2) // dimension = 2D
      {
        E_Float etVol;
        // Aux coins, traitement degenere.
        if (ithread == 0)
        {
          // imin, jmin
          reg[0] = ratioMax2(vol[0], vol[1], vol[ni1]);
          // imax, jmin
          ind  = ni1-1;
          ind1 = ind-1;
          ind2 = ind+ni1;
          reg[ind] = ratioMax2(vol[ind], vol[ind1], vol[ind2]);
          // imin, jmax
          ind  = (nj1-1)*ni1;
          ind1 = ind+1;
          ind2 = ind-ni1;
          reg[ind] = ratioMax2(vol[ind], vol[ind1], vol[ind2]);
          // imax, jmax
          ind  = nj1*ni1-1;
          ind1 = ind-1;
          ind2 = ind-ni1;
          reg[ind] = ratioMax2(vol[ind], vol[ind1], vol[ind2]);
        }
        
        // Aux aretes, traitement degenere.
        #pragma omp for nowait
        for (E_Int i = 1; i < ni1-1; i++)
        {
          // jmin
          ind  = i;
          ind1 = ind-1;
          ind2 = ind+1;
          ind3 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol, vol[ind1], vol[ind2], vol[ind3]);
          // jmax
          ind  = (nj1-1)*ni1 + i;
          ind1 = ind-1;
          ind2 = ind+1;
          ind3 = ind - ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol, vol[ind1], vol[ind2], vol[ind3]);
        }

        #pragma omp for nowait
        for (E_Int j = 1; j < nj1-1; j++)
        {
          // imin
          ind  = j*ni1;
          ind1 = ind + 1;
          ind2 = ind - ni1;
          ind3 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol, vol[ind1], vol[ind2], vol[ind3]);
          // imax
          ind  = j*ni1 + ni1-1;
          ind1 = ind - 1;
          ind2 = ind - ni1;
          ind3 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol, vol[ind1], vol[ind2], vol[ind3]);
        }

        // Boucle generale sur les indices des cellules interieures
        #pragma omp for collapse(2)
        for (E_Int j = 1; j < nj1-1; j++)
        for (E_Int i = 1; i < ni1-1; i++)
        {
          ind  = j*ni1 + i;
          ind1 = ind - 1;
          ind2 = ind + 1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
        }
      }
      else if (dimC == 3)  // dimension = 3D
      {	      
        E_Int ni1nj1 = ni1*nj1;
        E_Float etVol;
        // Aux coins, traitement degenere.
        if (ithread == 0)
        {
          // imin, jmin, kmin
          ind  = 0;
          ind1 = ind + 1;
          ind2 = ind + ni1;
          ind3 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imax, jmin, kmin
          ind  = ni1 - 1;
          ind1 = ind - 1;
          ind2 = ind + ni1;
          ind3 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imin, jmax, kmin
          ind  = (nj1-1)*ni1;
          ind1 = ind + 1;
          ind2 = ind - ni1;
          ind3 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imax, jmax, kmin
          ind  = ni1nj1 - 1;
          ind1 = ind - 1;
          ind2 = ind - ni1;
          ind3 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imin, jmin, kmax
          ind  = ni1nj1*(nk1-1);
          ind1 = ind + 1;
          ind2 = ind + ni1;
          ind3 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imax, jmin, kmax
          ind  = ni1nj1*(nk1-1)+ni1-1;
          ind1 = ind - 1;
          ind2 = ind + ni1;
          ind3 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imin, jmax, kmax
          ind  = ni1nj1*(nk1-1)+(nj1-1)*ni1;
          ind1 = ind + 1;
          ind2 = ind - ni1;
          ind3 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imax, jmax, kmax
          ind  = ni1nj1*(nk1-1)+ni1nj1 - 1;
          ind1 = ind - 1;
          ind2 = ind - ni1;
          ind3 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax3(etVol,vol[ind1],vol[ind2],vol[ind3]);
        }
    
        // Aux aretes, traitement degenere.
        #pragma omp for nowait
        for (E_Int i=1; i<ni1-1; i++)
        {
          // jmin, kmin
          ind  = i;
          ind1 = ind + 1;
          ind2 = ind - 1;
          ind3 = ind + ni1nj1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // jmax, kmin
          ind  = (nj1-1)*ni1+i;
          ind1 = ind + 1;
          ind2 = ind - 1;
          ind3 = ind + ni1nj1;
          ind4 = ind - ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // jmin, kmax
          ind  = ni1nj1*(nk1-1)+i;
          ind1 = ind + 1;
          ind2 = ind - 1;
          ind3 = ind - ni1nj1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // jmax, kmax
          ind  = ni1nj1*(nk1-1)+(nj1-1)*ni1+i;
          ind1 = ind + 1;
          ind2 = ind - 1;
          ind3 = ind - ni1nj1;
          ind4 = ind - ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
        }

        #pragma omp for nowait
        for (E_Int j=1; j< nj1-1; j++)
        {
          // imin, kmin
          ind  = j*ni1;
          ind1 = ind + 1;
          ind2 = ind + ni1nj1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imax, kmin
          ind  = j*ni1 + ni1-1;
          ind1 = ind - 1;
          ind2 = ind + ni1nj1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imin, kmax
          ind  = ni1nj1*(nk1-1)+j*ni1;
          ind1 = ind + 1;
          ind2 = ind - ni1nj1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imax, kmax
          ind  = ni1nj1*(nk1-1)+j*ni1 + ni1-1;
          ind1 = ind - 1;
          ind2 = ind - ni1nj1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
        }

        #pragma omp for nowait
        for (E_Int k = 1; k < nk1-1; k++)
        {
          // imin, jmin
          ind  = k*ni1nj1;
          ind1 = ind + 1;
          ind2 = ind + ni1;
          ind3 = ind - ni1nj1;
          ind4 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imax, jmin
          ind = k*ni1nj1 + ni1-1;
          ind1 = ind - 1;
          ind2 = ind + ni1;
          ind3 = ind - ni1nj1;
          ind4 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imin, jmax
          ind = k*ni1nj1+ni1*(nj1-1);
          ind1 = ind + 1;
          ind2 = ind - ni1;
          ind3 = ind - ni1nj1;
          ind4 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imax, jmax
          ind = k*ni1nj1+ni1nj1 -1;
          ind1 = ind - 1;
          ind2 = ind - ni1;
          ind3 = ind - ni1nj1;
          ind4 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
        }

        // Aux faces, traitement degenere.
        #pragma omp for nowait collapse(2)
        for (E_Int k=1; k < nk1-1; k++)
        for (E_Int j=1; j < nj1-1; j++)
        {
          // face imin
          ind  = k*ni1nj1 + j*ni1;
          ind1 = ind + 1;
          ind2 = ind + ni1;
          ind3 = ind - ni1;
          ind4 = ind + ni1nj1;
          ind5 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);

          // face imax
          ind  = k*ni1nj1 + j*ni1 + ni1 - 1;
          ind1 = ind - 1;
          ind2 = ind + ni1;
          ind3 = ind - ni1;
          ind4 = ind + ni1nj1;
          ind5 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);
        }

        #pragma omp for nowait collapse(2)
        for (E_Int k = 1; k < nk1-1; k++)
        for (E_Int i = 1; i < ni1-1; i++)
        {
          // face jmin
          ind  = k*ni1nj1 + i;
          ind1 = ind - 1;
          ind2 = ind + 1;
          ind3 = ind + ni1;
          ind4 = ind - ni1nj1;
          ind5 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);

          // face jmax
          ind  = k*ni1nj1 + (nj1-1)*ni1 + i;
          ind1 = ind - 1;
          ind2 = ind + 1;
          ind3 = ind - ni1;
          ind4 = ind - ni1nj1;
          ind5 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);
        }

        #pragma omp for nowait collapse(2)
        for (E_Int j = 1; j < nj1-1; j++)
        for (E_Int i = 1; i < ni1-1; i++)
        {
          // face kmin
          ind  =  j*ni1 + i;
          ind1 = ind - 1;
          ind2 = ind + 1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          ind5 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);

          // face kmax
          ind  = (nk1-1)*ni1nj1 + j*ni1 + i;
          ind1 = ind - 1;
          ind2 = ind + 1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          ind5 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);
        }

        // Boucle generale sur les indices des cellules interieures
        #pragma omp for collapse(3)
        for (E_Int k = 1; k < nk1-1; k++)
        for (E_Int j = 1; j < nj1-1; j++)
        for (E_Int i = 1; i < ni1-1; i++)
        {
          ind  = k*ni1nj1 + j*ni1 + i;
          ind1 = ind - 1;
          ind2 = ind + 1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          ind5 = ind - ni1nj1;
          ind6 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = ratioMax6(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5],vol[ind6]);
        }
      }
    }

    RELEASESHAREDS(tpl, f2);
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else // cas non structure
  {
    // Construction du tableau numpy stockant le ratio max de volumes entre
    // elements voisins, definissant la regularite
    PyObject* tpl = K_ARRAY::buildArray3(
      1, "regularity", npts, *cn, eltType, 1, api, true
    );
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* reg = f2->begin(1);

    if (strcmp(eltType, "NGON") == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getRegularityMap: not implemented for NGON arrays.");
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }

    // Pre-compute total number of elements and facets
    E_Int ierr;
    E_Int ntotFacets = 0, ntotElts = 0;
    E_Int nc = cn->getNConnect();
    std::vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);

    // Number of facets per element
    std::vector<E_Int> nfpe;
    ierr = K_CONNECT::getNFPE(nfpe, eltType, false);
    if (ierr != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getRegularityMap: Error computing nfpe.");
      for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }

    // Total number of elements/facets
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
      E_Int nelts = cm.getSize();
      ntotFacets += nfpe[ic]*nelts;
      ntotElts += nelts;
    }

    // Calcul de la connectivite element -> elements voisins
    std::vector<std::vector<E_Int> > cEEN(ntotElts);
    ierr = K_CONNECT::connectEV2EENbrs(eltType, npts, *cn, cEEN);
    if (ierr == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getRegularityMap: Error computing cEEN.");
      for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }

    // Get dimensionality
    E_Int dim = K_CONNECT::getDimME(eltType);
    
    K_FLD::FldArrayF snx(ntotFacets), sny(ntotFacets), snz(ntotFacets);
    K_FLD::FldArrayF surf(ntotFacets);
    K_FLD::FldArrayF vol(ntotElts);

    if (dim == 3)
    {
      K_METRIC::compMetricUnstruct(
        *cn, eltType, xp, yp, zp,
        snx.begin(), sny.begin(), snz.begin(), surf.begin(), vol.begin()
      );
    }
    else // 1D/2D
    {
      K_METRIC::compSurfUnstruct(
        *cn, eltType, xp, yp, zp,
        snx.begin(), sny.begin(), snz.begin(), vol.begin()
      );
    }

    // Calcul du ratio maximum 
    // entre les volumes des elements voisins et celui de l'element courant
    #pragma omp parallel
    {
      E_Int nneis, inei;
      E_Float maxr, voli;

      #pragma omp for
      for (E_Int i = 0; i < ntotElts; i++)
      {
        maxr = 0.;
        voli = vol[i];
        // Loop over neighbouring elements
        const std::vector<E_Int>& cEENi = cEEN[i];
        nneis = cEENi.size(); // number of neighboring cells
        for (E_Int n = 0; n < nneis; n++)
        {
          inei = cEENi[n]; // neighbor index
          maxr = K_FUNC::E_max(maxr, ratioMax1(voli, vol[inei]));
        }
        reg[i] = maxr;
      }
    }     

    for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
    RELEASESHAREDS(tpl, f2);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
}
