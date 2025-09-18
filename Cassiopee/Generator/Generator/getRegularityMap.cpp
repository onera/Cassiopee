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

// getRegularityMap

# include "generator.h"

using namespace std;
using namespace K_CONST;
using namespace K_FLD;
using namespace K_FUNC;

#define RATIOMAX2(v,v1,v2)              E_max(E_abs(v1-v)/E_max(v,E_GEOM_CUTOFF),E_abs(v2-v)/E_max(v,E_GEOM_CUTOFF))
#define RATIOMAX3(v,v1,v2,v3)           E_max(E_max(E_abs(v1-v)/E_max(v,E_GEOM_CUTOFF),E_abs(v2-v)/E_max(v,E_GEOM_CUTOFF)),E_abs(v3-v)/E_max(v,E_GEOM_CUTOFF))
#define RATIOMAX4(v,v1,v2,v3,v4)        E_max(E_max(E_abs(v1-v)/E_max(v,E_GEOM_CUTOFF),E_abs(v2-v)/E_max(v,E_GEOM_CUTOFF)),E_max(E_abs(v3-v)/E_max(v,E_GEOM_CUTOFF),E_abs(v4-v)/E_max(v,E_GEOM_CUTOFF)))
#define RATIOMAX5(v,v1,v2,v3,v4,v5)     E_max(E_max(E_max(E_abs(v1-v)/E_max(v,E_GEOM_CUTOFF),E_abs(v2-v)/E_max(v,E_GEOM_CUTOFF)),E_max(E_abs(v3-v)/E_max(v,E_GEOM_CUTOFF),E_abs(v4-v)/E_max(v,E_GEOM_CUTOFF))),E_abs(v5-v)/E_max(v,E_GEOM_CUTOFF))
#define RATIOMAX6(v,v1,v2,v3,v4,v5,v6)  E_max(E_max(E_max(E_abs(v1-v)/E_max(v,E_GEOM_CUTOFF),E_abs(v2-v)/E_max(v,E_GEOM_CUTOFF)),E_abs(v3-v)/E_max(v,E_GEOM_CUTOFF)),E_max(E_max(E_abs(v4-v)/E_max(v,E_GEOM_CUTOFF),E_abs(v5-v)/E_max(v,E_GEOM_CUTOFF)),E_abs(v6-v)/E_max(v,E_GEOM_CUTOFF)))

extern "C"
{
  void k6compunstrmetric_(E_Int& npts, E_Int& nelts, E_Int& nedges, 
                          E_Int& nnodes, E_Int* cn, 
                          E_Float* coordx, E_Float* coordy, E_Float* coordz, 
                          E_Float* xint, E_Float* yint, E_Float* zint, 
                          E_Float* snx, E_Float* sny, E_Float* snz, 
                          E_Float* surf, E_Float* vol);
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

  PyObject* tpl = NULL;
  
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
    E_Int dirI=2;
    E_Int dirJ=3;
    E_Int dirK=4;
    if (im == 1) dirI=1;
    if (jm == 1) dirJ=1;
    if (km == 1) dirK=1;
    E_Int dir = dirI*dirJ*dirK;
    E_Int ni,nj,nk;
    ni = im; nj = jm; nk = km;
    switch (dir)
    {
      case  2: // dim 1 - dir I
          ni = im; dim = 1; dimC = 1;
          break;
      case  3: // dim 1 - dir J
          ni = jm; dim = 1; dimC = 1;
          break;
      case  4: // dim 1 - dir K
          ni = km; dim = 1; dimC = 1;
          break;
      case  6: // dim 2 - dir IJ
          ni = im;
          nj = jm;
          dim = 2;
          dimC = 2;
          if (im == 2) { dimC = 1; ni = jm; nj = 1; }
          if (jm == 2) { dimC = 1; ni = im; nj = 1; }
          break;
      case  8: // dim 2 - dir IK
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
    tpl = K_ARRAY::buildArray(1, "regularity", im1, jm1, km1);
    // pointeur sur le tableau
    E_Float* reg = K_ARRAY::getFieldPtr(tpl);
      
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

    // calcul de la regularite	
    #pragma omp parallel
    {
      E_Int ithread = __CURRENT_THREAD__;
      // variables locales pour les indices
      E_Int ind,ind1,ind2,ind3,ind4,ind5,ind6;
      if (dimC == 1) // dimension 1D
      {
        // Aux frontieres, traitement degenere.
        if (ithread == 0)
        {
          reg[0] = E_abs(vol[1]-vol[0])/E_max(vol[0],E_GEOM_CUTOFF);
          reg[ni1-1] = E_abs(vol[ni1-2]-vol[ni1-1])/E_max(vol[ni1-1],E_GEOM_CUTOFF);
        }

        // Boucle sur les indices
        #pragma omp for schedule(static)
        for (E_Int i = 1; i < ni1-1; i++)
        {
          reg[i] = RATIOMAX2(vol[i],vol[i-1],vol[i+1]);
        }
      }
      else if (dimC == 2) // dimension = 2D
      {
        E_Float etVol;
        // Aux coins, traitement degenere.
        if (ithread == 0)
        {
          // imin, jmin
          reg[0] = RATIOMAX2(vol[0],vol[1],vol[ni1]);
          // imax, jmin
          ind  = ni1-1;
          ind1 = ind-1;
          ind2 = ind+ni1;
          reg[ind] = RATIOMAX2(vol[ind],vol[ind1],vol[ind2]);
          // imin, jmax
          ind  = (nj1-1)*ni1;
          ind1 = ind+1;
          ind2 = ind-ni1;
          reg[ind] = RATIOMAX2(vol[ind],vol[ind1],vol[ind2]);
          // imax, jmax
          ind  = nj1*ni1-1;
          ind1 = ind-1;
          ind2 = ind-ni1;
          reg[ind] = RATIOMAX2(vol[ind],vol[ind1],vol[ind2]);
        }
        
        // Aux aretes, traitement degenere.
        #pragma omp for schedule(static)
        for (E_Int i = 1; i < ni1-1; i++)
        {
          // jmin
          ind  = i;
          ind1 = ind-1;
          ind2 = ind+1;
          ind3 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // jmax
          ind  = (nj1-1)*ni1 + i;
          ind1 = ind-1;
          ind2 = ind+1;
          ind3 = ind - ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
        }

        #pragma omp for schedule(static)
        for (E_Int j = 1; j < nj1-1; j++)
        {
          // imin
          ind  = j*ni1;
          ind1 = ind + 1;
          ind2 = ind - ni1;
          ind3 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imax
          ind  = j*ni1 + ni1-1;
          ind1 = ind - 1;
          ind2 = ind - ni1;
          ind3 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
        }

        // Boucle generale sur les indices des cellules interieures
        for (E_Int j = 1; j < nj1-1; j++)
        {
          #pragma omp for schedule(static)
          for (E_Int i = 1; i < ni1-1; i++)
          {
            ind  = j*ni1 + i;
            ind1 = ind - 1;
            ind2 = ind + 1;
            ind3 = ind - ni1;
            ind4 = ind + ni1;
            etVol = vol[ind];
            reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          }
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
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imax, jmin, kmin
          ind  = ni1 - 1;
          ind1 = ind - 1;
          ind2 = ind + ni1;
          ind3 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imin, jmax, kmin
          ind  = (nj1-1)*ni1;
          ind1 = ind + 1;
          ind2 = ind - ni1;
          ind3 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imax, jmax, kmin
          ind  = ni1nj1 - 1;
          ind1 = ind - 1;
          ind2 = ind - ni1;
          ind3 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imin, jmin, kmax
          ind  = ni1nj1*(nk1-1);
          ind1 = ind + 1;
          ind2 = ind + ni1;
          ind3 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imax, jmin, kmax
          ind  = ni1nj1*(nk1-1)+ni1-1;
          ind1 = ind - 1;
          ind2 = ind + ni1;
          ind3 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imin, jmax, kmax
          ind  = ni1nj1*(nk1-1)+(nj1-1)*ni1;
          ind1 = ind + 1;
          ind2 = ind - ni1;
          ind3 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
          // imax, jmax, kmax
          ind  = ni1nj1*(nk1-1)+ni1nj1 - 1;
          ind1 = ind - 1;
          ind2 = ind - ni1;
          ind3 = ind - ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX3(etVol,vol[ind1],vol[ind2],vol[ind3]);
        }
    
        // Aux aretes, traitement degenere.
        #pragma omp for schedule(static)
        for (E_Int i=1; i<ni1-1; i++)
        {
          // jmin, kmin
          ind  = i;
          ind1 = ind + 1;
          ind2 = ind - 1;
          ind3 = ind + ni1nj1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // jmax, kmin
          ind  = (nj1-1)*ni1+i;
          ind1 = ind + 1;
          ind2 = ind - 1;
          ind3 = ind + ni1nj1;
          ind4 = ind - ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // jmin, kmax
          ind  = ni1nj1*(nk1-1)+i;
          ind1 = ind + 1;
          ind2 = ind - 1;
          ind3 = ind - ni1nj1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // jmax, kmax
          ind  = ni1nj1*(nk1-1)+(nj1-1)*ni1+i;
          ind1 = ind + 1;
          ind2 = ind - 1;
          ind3 = ind - ni1nj1;
          ind4 = ind - ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
        }

        #pragma omp for schedule(static)
        for (E_Int j=1; j< nj1-1; j++)
        {
          // imin, kmin
          ind  = j*ni1;
          ind1 = ind + 1;
          ind2 = ind + ni1nj1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imax, kmin
          ind  = j*ni1 + ni1-1;
          ind1 = ind - 1;
          ind2 = ind + ni1nj1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imin, kmax
          ind  = ni1nj1*(nk1-1)+j*ni1;
          ind1 = ind + 1;
          ind2 = ind - ni1nj1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imax, kmax
          ind  = ni1nj1*(nk1-1)+j*ni1 + ni1-1;
          ind1 = ind - 1;
          ind2 = ind - ni1nj1;
          ind3 = ind - ni1;
          ind4 = ind + ni1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
        }

        #pragma omp for schedule(static)
        for (E_Int k = 1; k < nk1-1; k++)
        {
          // imin, jmin
          ind  = k*ni1nj1;
          ind1 = ind + 1;
          ind2 = ind + ni1;
          ind3 = ind - ni1nj1;
          ind4 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imax, jmin
          ind = k*ni1nj1 + ni1-1;
          ind1 = ind - 1;
          ind2 = ind + ni1;
          ind3 = ind - ni1nj1;
          ind4 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imin, jmax
          ind = k*ni1nj1+ni1*(nj1-1);
          ind1 = ind + 1;
          ind2 = ind - ni1;
          ind3 = ind - ni1nj1;
          ind4 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
          // imax, jmax
          ind = k*ni1nj1+ni1nj1 -1;
          ind1 = ind - 1;
          ind2 = ind - ni1;
          ind3 = ind - ni1nj1;
          ind4 = ind + ni1nj1;
          etVol = vol[ind];
          reg[ind] = RATIOMAX4(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4]);
        }

        // Aux faces, traitement degenere.
        for (E_Int k=1; k < nk1-1; k++)
        {
          #pragma omp for schedule(static)
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
            reg[ind] = RATIOMAX5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);

            // face imax
            ind  = k*ni1nj1 + j*ni1 + ni1 - 1;
            ind1 = ind - 1;
            ind2 = ind + ni1;
            ind3 = ind - ni1;
            ind4 = ind + ni1nj1;
            ind5 = ind - ni1nj1;
            etVol = vol[ind];
            reg[ind] = RATIOMAX5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);
          }
        }

        for (E_Int k=1;k<nk1-1;k++)
        {
          #pragma omp for schedule(static)
          for (E_Int i=1;i<ni1-1;i++)
          {
            // face jmin
            ind  = k*ni1nj1 + i;
            ind1 = ind - 1;
            ind2 = ind + 1;
            ind3 = ind + ni1;
            ind4 = ind - ni1nj1;
            ind5 = ind + ni1nj1;
            etVol = vol[ind];
            reg[ind] = RATIOMAX5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);

            // face jmax
            ind  = k*ni1nj1 + (nj1-1)*ni1 + i;
            ind1 = ind - 1;
            ind2 = ind + 1;
            ind3 = ind - ni1;
            ind4 = ind - ni1nj1;
            ind5 = ind + ni1nj1;
            etVol = vol[ind];
            reg[ind] = RATIOMAX5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);
          }
        }

        for (E_Int j=1;j<nj1-1;j++)
        {
          #pragma omp for schedule(static)
          for (E_Int i=1;i<ni1-1;i++)
          {
            // face kmin
            ind  =  j*ni1 + i;
            ind1 = ind - 1;
            ind2 = ind + 1;
            ind3 = ind - ni1;
            ind4 = ind + ni1;
            ind5 = ind + ni1nj1;
            etVol = vol[ind];
            reg[ind] = RATIOMAX5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);

            // face kmax
            ind  = (nk1-1)*ni1nj1 + j*ni1 + i;
            ind1 = ind - 1;
            ind2 = ind + 1;
            ind3 = ind - ni1;
            ind4 = ind + ni1;
            ind5 = ind - ni1nj1;
            etVol = vol[ind];
            reg[ind] = RATIOMAX5(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5]);
          }
        }

        // Boucle generale sur les indices des cellules interieures
        for (E_Int k=1;k<nk1-1;k++)
        {
          for (E_Int j=1;j<nj1-1;j++)
          {
            #pragma omp for schedule(static)
            for (E_Int i=1;i<ni1-1;i++)
            {
              ind  = k*ni1nj1 + j*ni1 + i;
              ind1 = ind - 1;
              ind2 = ind + 1;
              ind3 = ind - ni1;
              ind4 = ind + ni1;
              ind5 = ind - ni1nj1;
              ind6 = ind + ni1nj1;
              etVol = vol[ind];
              reg[ind] = RATIOMAX6(etVol,vol[ind1],vol[ind2],vol[ind3],vol[ind4],vol[ind5],vol[ind6]);
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
    // cas non structure
    E_Int nelts = cn->getSize();
    E_Int nnodes = cn->getNfld(); // nb de noeuds ds 1 element
    E_Int npts = f->getSize();
    E_Float etVol; E_Int eti; E_Int indi;

    // Construction du tableau numpy stockant le ratio max de volumes entre elements voisins, 
    // definissant la regularite
    PyObject* tpl = K_ARRAY::buildArray(1, "regularity", nelts, nelts, -1, eltType, true);
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), nelts*nnodes);
    E_Float* regp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF reg(nelts,1, regp, true);

    // Calcul de la connectivite vertex->elements
    vector< vector<E_Int> > cVE(npts); 
    K_CONNECT::connectEV2VE(*cn, cVE);
  
    // Rapport MAX de volumes entre un element et ses voisins. Resultat au "centre"
    if (strcmp(eltType, "TRI") == 0)
    {
      E_Float maxratio;
      // Calcul du volume des elements
      FldArrayF snx(nelts);
      FldArrayF sny(nelts);
      FldArrayF snz(nelts);
      FldArrayF vol(nelts);
      FldArrayF volDummy(nelts);
      E_Int nedges = 1;
      //tableau local au fortran
      FldArrayF xint(nelts,nedges);
      FldArrayF yint(nelts,nedges);
      FldArrayF zint(nelts,nedges);
      //
      k6compunstrmetric_(npts, nelts, nedges, nnodes, cn->begin(), 
                         xp, yp, zp, xint.begin(), yint.begin(), zint.begin(),
                         snx.begin(), sny.begin(), 
                         snz.begin(), vol.begin(), volDummy.begin());
      
      // Calcul du ratio maximum 
      // entre les volumes des elements voisins et celui de l'element courant
      for (E_Int et = 0; et < nelts; et++)
      {
        etVol = vol[et];
        maxratio = 0;
        for (E_Int i = 0; i < nedges; i++)
        {
          E_Int* cni = cn->begin(i+1);
          indi = cni[et]-1;
          vector<E_Int>& cVEi = cVE[indi]; E_Int sizei = cVEi.size();
          for (E_Int noeti = 0; noeti < sizei; noeti++)
          {
            eti = cVEi[noeti];
            maxratio = E_max(maxratio,E_abs(vol[eti]-etVol)/E_max(etVol,E_GEOM_CUTOFF));
          }
        }
        reg[et] = maxratio;
      }
    }
    else if (strcmp(eltType, "QUAD") == 0)
    {
      E_Float maxratio;
      // Calcul du volume des elements
      FldArrayF snx(nelts);
      FldArrayF sny(nelts);
      FldArrayF snz(nelts);
      FldArrayF vol(nelts);
      FldArrayF volDummy(nelts);
      E_Int nedges = 1;
      // tableau local au fortran
      FldArrayF xint(nelts,nedges);
      FldArrayF yint(nelts,nedges);
      FldArrayF zint(nelts,nedges);
      //
      k6compunstrmetric_(npts, nelts, nedges, nnodes, cn->begin(), 
                         xp, yp, zp,
                         xint.begin(), yint.begin(), zint.begin(),
                         snx.begin(), sny.begin(), 
                         snz.begin(), vol.begin(), volDummy.begin());

      // Calcul du ratio maximum 
      // entre les volumes des elements voisins et celui de l'element courant
      for (E_Int et = 0; et < nelts; et++)
      {
        etVol = vol[et];
        maxratio = 0;
        for (E_Int i = 0; i < nedges; i++)
        {
          E_Int* cni = cn->begin(i+1);
          indi = cni[et]-1;
          vector<E_Int>& cVEi = cVE[indi]; E_Int sizei = cVEi.size();
          for (E_Int noeti = 0; noeti < sizei; noeti++)
          {
            eti = cVEi[noeti];
            maxratio = E_max(maxratio,E_abs(vol[eti]-etVol)/E_max(etVol,E_GEOM_CUTOFF));
          }
        }
        reg[et] = maxratio;
      }
   }
    else if (strcmp(eltType, "TETRA") == 0)
    {
      E_Int nedges = 4;
      E_Float maxratio;
      // Calcul du volume des elements
      FldArrayF snx(nelts, nedges);
      FldArrayF sny(nelts, nedges);
      FldArrayF snz(nelts, nedges);
      FldArrayF surf(nelts, nedges);
      FldArrayF vol(nelts);
      //tableau local au fortran
      FldArrayF xint(nelts,nedges);
      FldArrayF yint(nelts,nedges);
      FldArrayF zint(nelts,nedges);
      //
      k6compunstrmetric_(npts, nelts, nedges, nnodes, cn->begin(), 
                         xp, yp, zp,
                         xint.begin(), yint.begin(), zint.begin(),
                         snx.begin(), sny.begin(), 
                         snz.begin(), surf.begin(), vol.begin());
      
      // Calcul du ratio maximum 
      // entre les volumes des elements voisins et celui de l'element courant
      for (E_Int et = 0; et < nelts; et++)
      {
        etVol = vol[et];
        maxratio = 0;
        for (E_Int i = 0; i < nedges; i++)
        {
          E_Int* cni = cn->begin(i+1);
          indi = cni[et]-1;
          vector<E_Int>& cVEi = cVE[indi]; E_Int sizei = cVEi.size();
          for (E_Int noeti = 0; noeti < sizei; noeti++)
          {
            eti = cVEi[noeti];
            maxratio = E_max(maxratio,E_abs(vol[eti]-etVol)/E_max(etVol,E_GEOM_CUTOFF));
          }
        }
        reg[et] = maxratio;
      }
    }
    else if (strcmp(eltType, "PYRA") == 0)
    {
      // volume computation not implemented for PYRA
      PyErr_SetString(PyExc_TypeError,
                      "getRegularityMap: not yet implemented for PYRA.");
      return NULL;
    }
    else if (strcmp(eltType, "PENTA") == 0)
    {
      E_Int nedges = 5;
      E_Float maxratio;

      // Calcul du volume des elements
      FldArrayF snx(nelts, nedges);
      FldArrayF sny(nelts, nedges);
      FldArrayF snz(nelts, nedges);
      FldArrayF surf(nelts, nedges);
      FldArrayF vol(nelts);
      //tableau local au fortran
      FldArrayF xint(nelts,nedges);
      FldArrayF yint(nelts,nedges);
      FldArrayF zint(nelts,nedges);
      //
      k6compunstrmetric_(npts, nelts, nedges, nnodes, cn->begin(), 
                         xp, yp, zp,
                         xint.begin(), yint.begin(), zint.begin(),
                         snx.begin(), sny.begin(), 
                         snz.begin(), surf.begin(), vol.begin());

      // Calcul du ratio maximum 
      // entre les volumes des elements voisins et celui de l'element courant
      for (E_Int et = 0; et < nelts; et++)
      {
        etVol = vol[et];
        maxratio = 0;
        for (E_Int i = 0; i < nedges; i++)
        {
          E_Int* cni = cn->begin(i+1);
          indi = cni[et]-1;
          vector<E_Int>& cVEi = cVE[indi]; E_Int sizei = cVEi.size();
          for (E_Int noeti = 0; noeti < sizei; noeti++)
          {
            eti = cVEi[noeti];
            maxratio = E_max(maxratio,E_abs(vol[eti]-etVol)/E_max(etVol,E_GEOM_CUTOFF));
          }
        }
        reg[et] = maxratio;
      }
   }
    else if (strcmp(eltType, "HEXA") == 0)
    {
      E_Int nedges = 6;
      E_Float maxratio;
      // Calcul du volume des elements
      FldArrayF snx(nelts, nedges);
      FldArrayF sny(nelts, nedges);
      FldArrayF snz(nelts, nedges);
      FldArrayF surf(nelts, nedges);
      FldArrayF vol(nelts);
      //tableau local au fortran
      FldArrayF xint(nelts,nedges);
      FldArrayF yint(nelts,nedges);
      FldArrayF zint(nelts,nedges);
      //
      k6compunstrmetric_(npts, nelts, nedges, nnodes, cn->begin(), 
                         xp, yp, zp,
                         xint.begin(), yint.begin(), zint.begin(),
                         snx.begin(), sny.begin(), 
                         snz.begin(), surf.begin(), vol.begin());
      
      // Calcul du ratio maximum 
      // entre les volumes des elements voisins et celui de l'element courant
      for (E_Int et = 0; et < nelts; et++)
      {
        etVol = vol[et];
        maxratio = 0;
        for (E_Int i = 0; i < nedges; i++)
        {
          E_Int* cni = cn->begin(i+1);
          indi = cni[et]-1;
          vector<E_Int>& cVEi = cVE[indi]; E_Int sizei = cVEi.size();
          for (E_Int noeti = 0; noeti < sizei; noeti++)
          {
            eti = cVEi[noeti];
            maxratio = E_max(maxratio,E_abs(vol[eti]-etVol)/E_max(etVol,E_GEOM_CUTOFF));
          }
        }
        reg[et] = maxratio;
      }
    }
    else if (strcmp(eltType, "BAR") == 0)
    {
      // Calcul du volume des elements
      E_Int* cn1 = cn->begin(1);
      E_Int* cn2 = cn->begin(2);
      FldArrayF vol(nelts);
      E_Int ind1, ind2;
      E_Float dx, dy, dz;
      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cn1[i]-1; 
        ind2 = cn2[i]-1;
        dx = xp[ind2]-xp[ind1];
        dy = yp[ind2]-yp[ind1];
        dz = zp[ind2]-zp[ind1];
        vol[i] = sqrt(dx*dx + dy*dy + dz*dz);
      } 
    }
    else
    {
        PyErr_SetString(PyExc_TypeError,
                        "getRegularityMap: unknown type of element.");
        RELEASESHAREDU(array, f, cn);
        return NULL;
    }
    
    RELEASESHAREDU(array, f, cn); 
    return tpl;
  }
}