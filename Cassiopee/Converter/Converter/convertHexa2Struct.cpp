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
# include "converter.h"
# include <stack>
using namespace std;
using namespace K_FLD;

#define FINDNEXT { next = -1;                   \
  for (E_Int i = 0; i < nv; i++)                \
  { if (no[4*i] == 0) { next = i; break; } } }

//=============================================================================
// Numerote suivant i en partant de ind (i,j,k)
// IN: P1, P2: direction en ind
// IN: ind: indice du point de depart
// IN: i,j,k: indices corr. a ind
// IN/OUT: no: numerotation des pts
// OUT: nii: numerotation inverse nii[i,j,k] = ind
// IN: si: ?
// IN: cVN: connectivite vertex/vertex voisins
// IN: x,y,z: coords
// IN: maxAngle: angle max de deviation tolere
// Retourne le nbre de pts indices en i (ni).
//=============================================================================
E_Int numeroteI(E_Float* P1, E_Float* P2, E_Int ind, 
                E_Int i, E_Int j, E_Int k, E_Int* no, E_Int* nii, E_Int si,
                vector< vector<E_Int> >& cVN, 
                E_Float* x, E_Float* y, E_Float* z, E_Float maxAngle)
{
  E_Int indj;
  E_Float Pt[3]; E_Float* t;
  E_Float* P3 = Pt;
  // Cherche parmi les voisins celui qui a la deviance minimum par rapport
  // a P1, P2
  E_Float bestDev;
  E_Int indBest = 0;
  E_Float dx1, dy1, dz1, dx2, dy2, dz2, alpha, dev;
  E_Float dirVect[3];

  while (0 != 1)
  {
    // indice i,j,k
    no[4*ind] = 1; no[4*ind+1] = i; no[4*ind+2] = j; no[4*ind+3] = k;
    nii[i+j*si+k*2*si] = ind;
    vector<E_Int>& voisins = cVN[ind];
    E_Int n = voisins.size();
    bestDev = 1.e6;
    indBest = 0;
    //printf("n=" SF_D_ "\n", n);
    //if (n != 4) { no[4*ind] = 2; goto end; }

    for (E_Int j = 0; j < n; j++)
    {
      indj = voisins[j]-1;
      P3[0] = x[indj]; P3[1] = y[indj]; P3[2] = z[indj];
      dx1 = P2[0]-P1[0]; dy1 = P2[1]-P1[1]; dz1 = P2[2]-P1[2];
      dx2 = P3[0]-P1[0]; dy2 = P3[1]-P1[1]; dz2 = P3[2]-P1[2];
      dirVect[0] = dy1*dz2-dz1*dy2;
      dirVect[1] = dz1*dx2-dx1*dz2;
      dirVect[2] = dx1*dy2-dy1*dx2;
      alpha = K_COMPGEOM::getAlphaAngleBetweenBars(P1, P2, P1, P3, dirVect);
      dev = fabs(alpha-180.);
      if (dev < bestDev) { bestDev = dev; indBest = indj; }
    }
    printf("P1: %f %f %f\n", P1[0], P1[1], P1[2]);
    printf("P2: %f %f %f\n", P2[0], P2[1], P2[2]);
    P3[0] = x[indBest]; P3[1] = y[indBest]; P3[2] = z[indBest];
    printf("P3: %f %f %f\n\n", P3[0], P3[1], P3[2]);

    // verifie que le point pressenti n'est pas deja indice
    if (no[4*indBest] == 0 && bestDev < maxAngle)
    {
      i++;
      ind = indBest; t = P2; P2 = P1; P1 = P3; P3 = t;
    }
    else { no[4*indBest] = 2; goto end; }
  }

  end:
  return i;
}

//=============================================================================
PyObject* K_CONVERTER::convertHexa2Struct(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res > 2) return NULL;
  if (res == 1) { RELEASESHAREDS(array, f); return array; }
  
  if (K_STRING::cmp(eltType, "QUAD") != 0)
  {
    PyErr_SetString(PyExc_TypeError, "convertHexa2Struct: array must be of QUAD type.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if (posx == 0 || posy == 0 || posz == 0)
  {
    PyErr_SetString(PyExc_TypeError, "convertHexa2Struct: array must contain coordinates.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  // Construit les lignes i a partir de ce point de depart
  E_Int start = 0;
  
  E_Int nv = f->getSize();
  E_Int nfld = f->getNfld();
  vector< vector<E_Int> > cVN(nv);
  K_CONNECT::connectEV2VNbrs(*cn, cVN, 0);
  
  // Numerotation pour chaque noeud
  // no[ind] = 0 (non numerote), 1 (numerote), 2 (missed)
  // no[ind+1] = i
  // no[ind+2] = j
  // no[ind+3] = k
  FldArrayI no(4*nv); no.setAllValuesAtNull();
  FldArrayI ii(3*nv); ii.setAllValuesAt(-1);
  E_Float P1[3], P2[3];
  E_Int next; E_Float* fp;
  E_Int ind1, ind2;

  // Numerotation du premier noeud et d'un de ses voisins
  E_Int ind = start;

  ni = 1; nj = 1; nk = 1;

  // 0,0,0
  no[4*ind] = 1; no[4*ind+1] = 0; no[4*ind+2] = 0; no[4*ind+3] = 0;
  ii[0] = ind;
  vector<E_Int>& voisins = cVN[ind];

  E_Int istart = 0; E_Int jstart = 0; E_Int kstart = 0; // indice debut ligne en i
  //printf("Ligne j=" SF_D_ "\n", jstart);
  next = voisins[0]-1;
  P1[0] = x[next]; P1[1] = y[next]; P1[2] = z[next];
  P2[0] = x[ind]; P2[1] = y[ind]; P2[2] = z[ind];

  ni = numeroteI(P1, P2, next,
                 1, jstart, kstart, no.begin(), ii.begin(), nv,
                 cVN, x, y, z, 45.);

  // next j - recupere les P1,P2 de la ligne precedente
  next = voisins[1]-1;
  ind1 = ii[(istart+1)+jstart*nv+kstart*2*nv];
  P1[0] = x[ind1]; P1[1] = y[ind1]; P1[2] = z[ind1];
  ind2 = ii[istart+jstart*nv+kstart*2*nv];
  P2[0] = x[ind2]; P2[1] = y[ind2]; P2[2] = z[ind2];
  jstart++;
  //printf("Ligne j=" SF_D_ "\n", jstart);
  ni = numeroteI(P1, P2, next,
                 istart, jstart, kstart, no.begin(), ii.begin(), nv,
                 cVN, x, y, z, 45.);
  
  //printf("ni=" SF_D_ ", nj=" SF_D_ ", nk=" SF_D_ "\n", ni, nj, nk);  
  
  /* Reconstruit le maillage structure */
  E_Int api = f->getApi();
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, ni, nj, nk, api);
  FldArrayF* fn;
  K_ARRAY::getFromArray3(tpl, fn);

  for (E_Int n = 1; n <= nfld; n++) 
  {
    E_Float* fnp = fn->begin(n);
    fp = f->begin(n);
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          printf(SF_D3_ " - " SF_D_ "\n", i,j,k, ii[i+nv*j+2*nv*k]);
          fnp[i+j*ni+k*ni*nj] = fp[ii[i+nv*j+2*nv*k]];
        }
  }
  RELEASESHAREDS(tpl, fn);

  return tpl;
}

