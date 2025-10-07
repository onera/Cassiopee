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

// Operations on distributions

# include "generator.h"
# include <vector>

using namespace std;
using namespace K_CONST;
using namespace K_FUNC;
using namespace K_FLD;
# define MSGX printf("Warning: enforceX (x0=%f, eh=%f): ", x0, eh);
# define MSGPLUSX printf("Warning: enforcePlusX (eh=%f): ", eh);
# define MSGMOINSX printf("Warning: enforceMoinsX (eh=%f): ", eh);
# define MSGY printf("Warning: enforceY (x0=%f, eh=%f): ", y0, eh);
# define MSGPLUSY printf("Warning: enforcePlusY (eh=%f): ", eh);
# define MSGMOINSY printf("Warning: enforceMoinsY (eh=%f): ", eh);

extern "C"
{
  void k6stretch_(const E_Float& t1, const E_Float& t2,
                  E_Float* sn,
                  const E_Int& nbp, const E_Float& dsm,
                  const E_Float& dsp,
                  const E_Int& ityp, const E_Int& inewt);
}

// ============================================================================
/* Enforce a x line in a distribution */
// ============================================================================
// Renforce une distribution [0;ni-1] en 0 entre istart et iend autour de x0
// Remplace 2*supp points (ou moins si probleme) : de istart (non modifie) a
// iend (non modifie)
// Par 2*supp+2*add points (ou + ou - si probleme)
// ============================================================================
PyObject* K_GENERATOR::enforceXMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int supp, suppl, suppr; // supressed points, left, right
  E_Int add;                // added points
  E_Float eh, ehr;          // enforce length
  E_Float x0;
  E_Int verbose;

  if (!PYPARSETUPLE_(args, O_ RR_ TII_ I_,
                     &array, &x0, &eh, &supp, &add, &verbose)) return NULL;

  bool isVerbose = (verbose == 1) ? true : false;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);

  FldArrayF snr; FldArrayF snl; FldArrayF snl2;
  E_Int jc, il, i, k, istart, iend;
  E_Float hl;
  E_Bool pb;

  if (res == 1)
  {
    FldArrayF* an = new FldArrayF();
    FldArrayF& coord2 = *an;

    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "enforceX: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    pb = false;

    FldArrayF& coord = *f;
    // Determination de l'indice de modification
    jc = 0;
    hl = x0;

    E_Int api = f->getApi();
    E_Float* xt = coord.begin(posx);

    // Determination de la cellule de distrib contenant hl
    if (hl <= xt[0] || hl >= xt[ni-1])
    {
      char msg[256];
      sprintf(msg,
              "enforceX (x0=" SF_F_ ", eh=" SF_F_ "): cannot find hl in distrib, stopped.",
              x0, eh);
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError, msg);
      return NULL;
    }
    if (hl <= xt[1])
    {
      char msg[256];
      sprintf(msg,
              "enforceX (x0=" SF_F_ ", eh=" SF_F_ "): line to enforce is in the first cell of distribution. \n=> use enforcePlusX to enforce this mesh.", x0, eh);
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError, msg);
      return NULL;
    }
    if (hl >=xt[ni-2])
    {
      char msg[256];
      sprintf(msg,
              "enforceX (x0=" SF_F_ ", eh=" SF_F_ "): line to enforce is in the last cell of distribution. \n=> use enforceMoinsX to enforce this mesh.", x0, eh);
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError, msg);
      return NULL;
    }
    if (hl <= xt[1] + eh/2.)
    {
      char msg[256];
      sprintf(msg,
              "enforceX (x0=" SF_F_ ", eh=" SF_F_ "): overlapping around node 1.\n=> decrease eh or increase hl or enforce initial mesh.", x0, eh);
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError, msg);
      return NULL;
    }
    if (hl >= xt[ni-2] - eh/2.)
    {
      char msg[256];
      sprintf(msg,
              "enforceX (x0=" SF_F_ ", eh=" SF_F_ "): overlapping around node ni-2.\n=> decrease eh or decrease hl.\n", x0, eh);
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      msg);
      return NULL;
    }

    E_Float h ;
    E_Int ninj = ni*nj;
    for (i = 0; i < ni; i++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = i + jc*ni + k*ninj;
        h = xt[ind];
        if (h > hl) goto exit;
      }
    exit:

    il = i-1;

    suppl = supp;
    suppr = supp;

    // Verifications des parametres supp et add et modification si necessaire
    if (supp > il)
    {
      suppl = il;
      if (isVerbose)
      {
        MSGX;
        printf("cannot suppress " SF_D_ " points on left side.\n", supp);
        printf("supp set to " SF_D_ " on left side.\n", suppl);
      }
    }
    if (supp > ni-1-il-1)
    {
      suppr = ni-1-il-1;
      if (isVerbose)
      {
        MSGX;
        printf("cannot suppress " SF_D_ " points on right side.\n", supp);
        printf("supp set to " SF_D_ " on right side.\n", suppr);
      }
    }
    while ((suppl+add < 2) || (suppr+add < 2))
    {
      add++;
      pb = true;
    }
    if (pb == true && isVerbose)
    {
      MSGX;
      printf("supp+add < 2.\n");
      printf("   => add set to: " SF_D_ ".\n", add);
    }
    pb = false;
    if (add == 0) // add = 0 => risque de probleme sur la croissance ou decroissance des tailles de maille
    {
      add++;
      if (isVerbose)
      {
        MSGX;
        printf("add = 0 => problems will occur on step-size.\n");
        printf("   => add set to: " SF_D_ ".\n", add);
      }
    }

    istart = il - suppl;
    iend = il + suppr + 1;

    // Verification : comparaison entre eh et deltal, eh et deltar
    E_Float pt1 = xt[istart+jc*ni];
    E_Float pt1a = xt[istart+jc*ni+1];
    E_Float pt2a = xt[iend+jc*ni-1];
    E_Float pt2 = xt[iend+jc*ni];
    E_Float pt3 = hl-eh/2.;
    E_Float pt4 = hl+eh/2.;

    E_Float deltal = pt1a - pt1;
    E_Float deltar = pt2 - pt2a;

    if (add > 0 && isVerbose)
    {
      if (eh > deltal)
      {
        MSGX;
        printf("add = " SF_D_ " > 0  while eh = " SF_F_ " ",add,eh);
        printf("is greater than the first step on left side = " SF_F_ ".\n", deltal);
      }
      if (eh > deltar)
      {
        MSGX;
        printf("add = " SF_D_ " > 0  while eh = " SF_F_ " ",add,eh);
        printf("is greater than the first step on right side = " SF_F_ ".\n", deltar);
      }
    }
    if (add < 0 && isVerbose)
    {
      if (eh < deltal)
      {
        MSGX;
        printf("add = " SF_D_ " < 0  while eh = " SF_F_ " ",add,eh);
        printf("is lower than the first step on left side = " SF_F_ ".\n", deltal);
      }
      if (eh < deltar)
      {
        MSGX;
        printf("add = " SF_D_ " < 0  while eh = " SF_F_ " ",add,eh);
        printf("is lower than the last step on right side = " SF_F_ ".\n", deltar);
      }
    }

    // Distribution a gauche
    E_Int npl = suppl + add + 2 ;
    snl.malloc(npl);
    snl2.malloc(npl);

    k6stretch_(pt1, pt4, snl2.begin(), npl, eh, deltal, 2, 1);
    for (E_Int pp = 0 ; pp < npl; pp++)
    {
      snl[pp] = -snl2[npl-1-pp]+pt4+pt1;
    }
    // Verification de la distribution : decroissance des tailles de maille a gauche si add > 0
    pb = false;
    if (add > 0)
    {
      for (i = 0 ; i < npl-2 ; i++)
      {
        if (snl[i+2]-snl[i+1] >= snl[i+1]-snl[i])
        {
          pb = true;
        }
      }
      if (pb == true && isVerbose)
      {
        MSGX;
        printf("non-decreasing step-size in left-distribution.\n");
        printf("   => Please change add = " SF_D_ " or eh  = " SF_F_ ".\n", add, eh);
      }
    }
    // Verification de la distribution : croissance des tailles de maille a gauche si add < 0
    pb = false;
    if (add < 0)
    {
      for (i = 0 ; i < npl-2 ; i++)
      {
        if (snl[i+2]-snl[i+1] <= snl[i+1]-snl[i])
        {
          pb = true;
        }
      }
      if (pb == true && isVerbose)
      {
        MSGX;
        printf("non-increasing step-size in left-distribution.\n");
        printf("   => Please change add = " SF_D_ " or eh = " SF_F_ ".\n", add, eh);
      }
    }
    pb = false;

    // Distribution a droite
    E_Int npr = suppr + add + 2;
    snr.malloc(npr) ;
    k6stretch_(pt3, pt2, snr.begin(), npr, eh, deltar, 2, 1);
    // Verification de la distribution : croissance des tailles de maille a droite si add > 0
    pb = false;
    if (add > 0)
    {
      for (i = 0 ; i < npr-2 ; i++)
      {
        if (snr[i+2]-snr[i+1] <= snr[i+1]-snr[i])
        {
          pb = true;
        }
      }
      if (pb == true && isVerbose)
      {
        MSGX;
        printf("non-increasing step-size in right-distribution.\n");
        printf("   => Please change add = " SF_D_ " or eh = " SF_F_ ".\n", add, eh);
      }
    }
    // Verification de la distribution : decroissance des tailles de maille a droite si add < 0
    pb = false;
    if (add < 0)
    {
      for (i = 0 ; i < npr-2 ; i++)
      {
        if (snr[i+2]-snr[i+1] >= snr[i+1]-snr[i])
        {
          pb = true;
        }
      }
      if (pb == true && isVerbose)
      {
        MSGX;
        printf("non-decreasing step-size in right-distribution.\n");
        printf("   => Please change add = " SF_D_ " or eh  = " SF_F_ ".\n", add, eh);
      }
    }
    pb = false;

    // Distribution finale
    E_Int np = ni+2*add;
    coord2.malloc(np*nj*nk, 3);
    coord2.setAllValuesAtNull();
    E_Float* yt = coord.begin(posy);
    E_Float* zt = coord.begin(posz);
    E_Float* xt2 = coord2.begin(posx);
    E_Float* yt2 = coord2.begin(posy);
    E_Float* zt2 = coord2.begin(posz);

    for (i = 0 ; i < istart ; i++)
      for (k = 0 ; k < nk ; k++)
      {
        E_Int ind0 = i+jc*ni+k*ninj;
        E_Int ind = i + jc*np + k*np*nj;
        xt2[ind] = xt[ind0];
        yt2[ind] = yt[ind0];
        zt2[ind] = zt[ind0];
      }

    for (i = istart ; i < istart+npl-1 ; i++)
      for (k = 0 ; k < nk ; k++)
      {
        E_Int ind0 = jc*ni+k*ninj;
        E_Int ind = i + jc*np + k*np*nj;
        xt2[ind] = snl[i-istart];
        yt2[ind] = yt[ind0];
        zt2[ind] = zt[ind0];
      }

    for (i = istart+npl-1 ; i < istart+npl-1+npr-1 ; i++)
      for (k = 0 ; k < nk ; k++)
      {
        E_Int ind0 = jc*ni+k*ninj;
        E_Int ind = i + jc*np + k*np*nj;
        xt2[ind] = snr[i-istart-npl+2];
        yt2[ind] = yt[ind0];
        zt2[ind] = zt[ind0];
      }

    for (i = istart+npl-1+npr-1 ; i < ni+2*add ; i++)
      for (k = 0 ; k < nk ; k++)
      {
        E_Int ind0 = i-2*add+jc*ni+k*ninj;
        E_Int ind = i + jc*np + k*np*nj;
        xt2[ind] = xt[ind0];
        yt2[ind] = yt[ind0];
        zt2[ind] = zt[ind0];
      }

    // Verification de la nouvelle distribution globale : croissance des abscisses
    pb = false;
    for (i = 0 ; i < ni+2*add-1 ; i++)
      for (k = 0 ; k < nk ; k++)
      {
        E_Int ind = i + jc*np + k*np*nj;
        E_Int ind2 = i+1 + jc*np + k*np*nj;
        if ( xt2[ind] >= xt2[ind2] )
        {
          if (isVerbose)
          {
            MSGX;
            printf("X coordinates are not increasing:\n");
            printf("   X(" SF_D_ ") = " SF_F_ " and X(" SF_D_ ") = " SF_F_ ".\n",
                   i, coord2(ind,1), i+1, coord2(ind2,1));
          }
          pb = true;
        }
      }
    if (pb == true && isVerbose)
    {
      printf("      => Please change add = " SF_D_ " or supp = " SF_D_ " or eh  = " SF_F_ ".\n",
             add, supp, eh);
    }
    pb = false;

    // Verification de la nouvelle distribution globale : eh
    ehr = xt2[istart+npl-1] - xt2[istart+npl-2] - eh;
    if ((ehr > eh/1000. || ehr < -eh/1000.) && isVerbose)
    {
      MSGX;
      printf("eh different from the specified value.\n");
      printf("in new distribution, eh = " SF_F_ ".\n", ehr+eh);
    }

    jc++;
    // For debug
    // In reality, we must parametrize the curve
    while (jc < nj)
    {
      for (i = 0; i < np; i++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = i + jc*np + k*np*nj;
          E_Int ind0 = i + k*np*nj;
          E_Int ind1 = 0 + jc*ni + k*ni*nj;
          xt2[ind] = xt2[ind0];
          yt2[ind] = yt[ind1];
          zt2[ind] = zt[ind1];
        }
      jc++;
    }

    RELEASESHAREDS(array, f);  
    PyObject* tpl = K_ARRAY::buildArray3(*an, varString, np, nj, nk, api);
    delete an;
    return tpl;
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);  
    PyErr_SetString(PyExc_TypeError,
                    "enforceX: not for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforceX: invalid array.");
    return NULL;
  }
}

// ============================================================================
/* Enforce a x line in a distribution */
// ============================================================================
// Renforce une distribution [0;ni-1] en 0 entre istart = 0 et
// iend = istart + supp + 1
// Remplace supp points (ou ni-2 si supp > ni-2) : de 1 (= istart+1) a
// supp (= iend-1)
// Par supp+add points
// ============================================================================
PyObject* K_GENERATOR::enforcePlusXMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int supp;      // supressed points
  E_Int add;       // added points
  E_Float eh, ehr;    // enforce length
  E_Int verbose;

  if (!PYPARSETUPLE_(args, O_ R_ TII_ I_,
                    &array, &eh, &supp, &add, &verbose)) return NULL;

  bool isVerbose = (verbose == 1) ? true : false;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  FldArrayF coord; FldArrayF sn1;
  E_Int jc, np, np1;
  E_Float delta1;

  E_Int api = f->getApi();

  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "enforcePlusX: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int i, k;
    E_Int istart, iend;
    E_Bool pb = false;

    // Verification des parametres supp et add et modification si necessaire
    i = 0;
    istart = 0;
    if (supp <= 0)
    {
      supp = 1;
      if (isVerbose)
      {
        MSGPLUSX;
        printf("supp <= 0.\n");
        printf("   =>supp set to: " SF_D_ ".\n", supp);
      }
    }
    if (supp > ni-2)
    {
      supp = ni-2;
      if (isVerbose)
      {
        MSGPLUSX;
        printf("cannot suppress " SF_D_ " points.\n", supp);
        printf("   =>supp set to: " SF_D_ ".\n", supp);
      }
    }
    iend = istart + supp + 1 ;
    if (supp+add < 2) // supp+add < 2 => probleme sur-determine
    {
      while (supp+add < 2) {add++;}
      if (isVerbose)
      {
        MSGPLUSX;
        printf("supp+add=" SF_D_ " < 2.\n", supp+add);
        printf("   => add set to: " SF_D_ ".\n", add);
      }
    }
    if (add == 0) // add=0 => risque de probleme sur la croissance ou decroissance des tailles de maille
    {
      add++;
      if (isVerbose)
      {
        MSGPLUSX;
        printf("add=0 => problems will occur on step-size.\n");
        printf("   => add set to: " SF_D_ ".\n", add);
      }
    }

    // Distribution finale
    np = ni+add;
    coord.malloc(np*nj*nk, 3);
    coord.setAllValuesAtNull();

    // Nouvelle distribution sur la zone modifiee : istart <= i <= iend (istart et iend non modifies)
    np1 = supp+add+2;
    sn1.malloc(np1);

    // Distributions pour j = 0
    jc = 0;

    // Verification : comparaison entre eh et delta1
    delta1 = (*f)(iend+jc*ni, posx) - (*f)(iend-1+jc*ni, posx);

    if (isVerbose)
    {
      if (add > 0)
      {
        if (eh > delta1)
        {
          MSGPLUSX;
          printf("   add=" SF_D_ " >= 0 while eh=" SF_F_ "", add, eh);
          printf(" is greater than the last step=" SF_F_ ".\n", delta1);
        }
      }
      if (add < 0)
      {
        if (eh < delta1)
        {
          MSGPLUSX;
          printf("   add=" SF_D_ " <= 0 while eh=" SF_F_ "", add, eh);
          printf(" is lower than the last step=" SF_F_ ".\n",  delta1);
        }
      }
    }

    // Taille de la premiere maille=eh ; taille de la derniere maille inchangee = x(iend) - x(iend-1)
    E_Float pt1 = (*f)(istart+jc*ni, posx);
    E_Float pt2 = (*f)(iend+jc*ni, posx);
    k6stretch_(pt1, pt2, sn1.begin(), np1, eh, delta1, 2, 1);

    // Verification de la distribution : croissance des points
    for (i = 0 ; i < np1-1 ; i++)
    {
      if (sn1[i+1] <= sn1[i])
      {
        pb = true;
        if (isVerbose)
        {
          MSGPLUSX;
          printf("X coordinates are not increasing:\n");
          printf("   X(" SF_D_ ") = " SF_F_ "", i, sn1[i]);
          printf("   X(" SF_D_ ") = " SF_F_ ".\n", i+1, sn1[i+1]);
        }
      }
    }

    if (pb == true && isVerbose)
    {
      printf("      => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
    }

    // Verification de la distribution : croissance des tailles de maille si add > 0
    pb = false;
    if (add > 0)
    {
      for (i = 0 ; i < np1-2 ; i++)
      {
        if (sn1[i+2]-sn1[i+1] <= sn1[i+1]-sn1[i])
        {
          pb = true;
        }
      }
      if (pb == true && isVerbose)
      {
        MSGPLUSX;
        printf("non-increasing step-size in distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }

    // Verification de la distribution : decroissance des tailles de maille si add < 0
    pb = false;
    if (add < 0)
    {
      for (i = 0 ; i < np1-2 ; i++)
      {
        if (sn1[i+2]-sn1[i+1] >= sn1[i+1]-sn1[i])
        {
          pb = true;
        }
      }
      if (pb == true && isVerbose)
      {
        MSGPLUSX;
        printf("non-decreasing step-size in distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }

    // Creation de la nouvelle distribution
    for (i = 0 ; i < np1 ; i++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = i + jc*np + k*np*nj;
        coord(ind,1) = sn1[i];
        coord(ind,2) = (*f)(jc*ni+k*ni*nj, posy);
        coord(ind,3) = (*f)(jc*ni+k*ni*nj, posz);
      }

    for (i = np1 ; i < np ; i++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = i + jc*np + k*np*nj;
        coord(ind,1) = (*f)(i-add+jc*ni+k*ni*nj, posx);
        coord(ind,2) = (*f)(i-add+jc*ni+k*ni*nj, posy);
        coord(ind,3) = (*f)(i-add+jc*ni+k*ni*nj, posz);
      }
    // Verification de la nouvelle distribution globale : croissance des abscisses
    pb = false ;
    for (i = 0 ; i < np-1 ; i++)
      for (k = 0 ; k < nk ; k++)
      {
        E_Int ind = i + jc*np + k*np*nj;
        E_Int ind2 = i+1 + jc*np + k*np*nj;
        if (coord(ind,1) >= coord(ind2,1))
        {
          pb = true;
          if (isVerbose)
          {
            MSGPLUSX;
            printf("X coordinates are not increasing:\n");
            printf("   X(" SF_D_ ") = " SF_F_ "", i, coord(ind,1));
            printf("   X(" SF_D_ ") = " SF_F_ ".\n", i+1, coord(ind2,1));
          }
        }
      }
    if (pb == true and isVerbose)
    {
      printf("      => Please change add=" SF_D_ " or supp=" SF_D_ " or eh=" SF_F_ ".\n", add, supp, eh);
    }
    pb = false;

    // Verification de la nouvelle distribution globale : eh
    ehr = coord(1,1)-coord(0,1)-eh ;
    if (isVerbose && (ehr > eh/1000. || ehr < -eh/1000.))
    {
      MSGPLUSX;
      printf("eh different from the specified value.\n");
      printf("   in new distribution, eh=" SF_F_ ".\n", coord(1,1)-coord(0,1));
    }

    jc++;

    while (jc < nj)
    {
      for (i = 0; i < np; i++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = i + jc*np + k*np*nj;
          E_Int ind0 = i + k*np*nj;
          E_Int ind1 = 0 + jc*ni + k*ni*nj;
          coord(ind,1) = coord(ind0,1);
          coord(ind,2) = (*f)(ind1, posy);
          coord(ind,3) = (*f)(ind1,posz);
        }
      jc++;
    }
    RELEASESHAREDS(array, f);  
    // Build array
    PyObject* tpl = K_ARRAY::buildArray3(coord, varString, np, nj, nk, api);
    return tpl;
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "enforcePlusX: not used for unstructured arrays.");
    return NULL;

  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforcePlusX: unknown type of array.");
    return NULL;
  }
}

// ============================================================================
/* Enforce a x line in a distribution */
// ============================================================================
// Renforce une distribution [0;ni-1] en 0 entre istart = iend - supp - 1
// et iend = ni
// Remplace supp points (ou ni-2 si supp > ni-2) : de 1 (= istart+1) a
// supp (= iend-1) par supp+add points
// ============================================================================
PyObject* K_GENERATOR::enforceMoinsXMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int supp;      // supressed points
  E_Int add;       // added points
  E_Float eh, ehr;    // enforce height
  E_Int verbose;
  if (!PYPARSETUPLE_(args, O_ R_ TII_ I_,
                    &array, &eh, &supp, &add, &verbose)) return NULL;

  bool isVerbose = (verbose == 1) ? true : false;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  FldArrayF coord;
  FldArrayF sn1, sn2;

  E_Int api = f->getApi();

  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "enforceMoinsX: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int i, k;
    E_Int istart, iend;
    E_Bool pb = false;

    // Verification des parametres supp et add et modification si necessaire
    iend = ni-1 ;
    if (supp <= 0)
    {
      supp = 1;
      if (isVerbose)
      {
        MSGMOINSX;
        printf("supp <= 0.\n");
        printf("   =>supp set to: supp=" SF_D_ ".\n", supp);
      }
    }
    if (supp > ni-2)
    {
      supp = ni-2;
      if (isVerbose)
      {
        MSGMOINSX;
        printf("cannot suppress " SF_D_ " points.\n", supp);
        printf("   =>supp set to: supp=" SF_D_ ".\n", supp);
      }
    }
    istart = iend - supp - 1;
    if (supp+add < 2) // supp+add < 2 => probleme sur-determine
    {
      while (supp+add < 2) {add++;}
      if (isVerbose)
      {
        MSGMOINSX;
        printf("supp+add=" SF_D_ " < 2.\n", supp+add);
        printf("   => add set to: add=" SF_D_ ".\n", add);
      }
    }
    if (add == 0) // add = 0 => risque de probleme sur la croissance ou decroissance des tailles de maille
    {
      add++;
      if (isVerbose)
      {
        MSGMOINSX;
        printf("add=0 => problems will occur on step-size.\n");
        printf("   => add set to: add=" SF_D_ ".\n", add);
      }
    }

    // Distribution finale
    E_Int np = ni+add ;
    coord.malloc(np*nj*nk, 3);
    coord.setAllValuesAtNull();

    // Nouvelle distribution sur la zone modifiee : istart <= i <= iend (istart et iend non modifies)
    E_Int np1 = supp+add+2;
    sn1.malloc(np1);
    sn2.malloc(np1);

    // Distributions pour j = 0
    E_Int jc = 0;

    // Verification : comparaison entre eh et delta1
    E_Float delta1 =
      (*f)(istart+1+jc*ni,posx) - (*f)(istart+jc*ni+jc*ni,posx);

    if (isVerbose)
    {
      if (add > 0)
      {
        if (eh > delta1)
        {
          MSGMOINSX;
          printf("   add=" SF_D_ " > 0 while eh=" SF_F_ " is greater than the last step=" SF_F_ ".\n", add, eh, delta1);
        }
      }
      if (add < 0)
      {
        if (eh < delta1)
        {
          MSGMOINSX;
          printf("   add=" SF_D_ " < 0  while eh=" SF_F_ " is lower than the last step=" SF_F_ ".\n", add, eh, delta1);
        }
      }
    }

    // Taille de la premiere maille = delta1 ; taille de la derniere maille inchangee = x(iend) - x(iend-1)
    E_Float pt1 = (*f)(istart+jc*ni,posx);
    E_Float pt2 = (*f)(iend+jc*ni,posx);
    k6stretch_(pt1, pt2, sn2.begin(), np1, eh, delta1, 2, 1);
    for (i = 0; i < np1; i++)
      sn1[i] = -sn2[np1-1-i]+pt2+pt1;

    // Verification de la distribution : decroissance des tailles de maille si add > 0
    pb = false;
    if (add > 0)
    {
      for (i = 0 ; i < np1-2 ; i++)
      {
        if (sn1[i+2]-sn1[i+1] >= sn1[i+1]-sn1[i])
        {
          pb = true; break;
        }
      }
      if (isVerbose && pb == true)
      {
        MSGMOINSX;
        printf("Non-decreasing step-size in distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }

    // Verification de la distribution : decroissance des tailles de maille si add < 0
    pb = false;
    if (add < 0)
    {
      for (i = 0 ; i < np1-2 ; i++)
      {
        if (sn1[i+2]-sn1[i+1] <= sn1[i+1]-sn1[i])
        {
          pb = true; break;
        }
      }
      if (isVerbose && pb == true)
      {
        MSGMOINSX;
        printf("non-increasing step-size in distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }

    // Creation de la nouvelle distribution
    for (i = 0 ; i < istart ; i++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = i + jc*np + k*np*nj;
        coord(ind,1) = (*f)(i+jc*ni+k*ni*nj, posx);
        coord(ind,2) = (*f)(i+jc*ni+k*ni*nj, posy);
        coord(ind,3) = (*f)(i+jc*ni+k*ni*nj, posz);
      }

    for (i = istart ; i < np ; i++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = i + jc*np + k*np*nj;
        coord(ind,1) = sn1[i-istart];
        coord(ind,2) = (*f)(jc*ni+k*ni*nj, posy);
        coord(ind,3) = (*f)(jc*ni+k*ni*nj, posz);
      }

    // Verification de la nouvelle distribution globale : croissance des abscisses
    pb = false;
    for (i = 0 ; i < np-1 ; i++)
      for (k = 0 ; k < nk ; k++)
      {
        E_Int ind = i + jc*np + k*np*nj;
        E_Int ind2 = i+1 + jc*np + k*np*nj;
        if (coord(ind,1) >= coord(ind2,1))
        {
          if (isVerbose)
          {
            MSGMOINSX;
            printf("X coordinates are not increasing:\n");
            printf("   X(" SF_D_ ") = " SF_F_ "", i, coord(ind,1));
            printf("   X(" SF_D_ ") = " SF_F_ ".\n", i+1, coord(ind2,1));
          }
          pb = true;
        }
      }
    if (isVerbose && pb == true)
    {
      printf("      => Please change add=" SF_D_ ".\n", add);
      printf("      => or supp=" SF_D_ ".\n", supp);
      printf("      => or eh=" SF_F_ ".\n", eh);
    }
    pb = false;

    // Verification de la nouvelle distribution globale : eh
    ehr = coord(np-1,1)-coord(np-2,1)-eh ;
    if (isVerbose && (ehr > eh/1000. || ehr < -eh/1000.))
    {
      MSGMOINSX;
      printf("eh different from the specified value.\n");
      printf("   in new distribution, eh=" SF_F_ ".\n", coord(1,1)-coord(0,1));
    }

    jc++;
    // For debug
    // In reality, we must parametrize the curve
    while (jc < nj)
    {
      for (i = 0; i < np; i++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = i + jc*np + k*np*nj;
          E_Int ind0 = i + k*np*nj;
          E_Int ind1 = 0 + jc*ni + k*ni*nj;
          coord(ind,1) = coord(ind0,1);
          coord(ind,2) = (*f)(ind1, posy);
          coord(ind,3) = (*f)(ind1, posz);
        }
      jc++;
    }

    RELEASESHAREDS(array, f);  
    PyObject* tpl = K_ARRAY::buildArray3(coord, varString, np, nj, nk, api);
    return tpl;
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "enforceMoinsX: not used for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforceMoinsX: invalid type of array.");
    return NULL;
  }
}

// ============================================================================
/* Enforce a j line in a distribution */
// ============================================================================
PyObject* K_GENERATOR::enforceYMesh(PyObject* self, PyObject* args)
{
  E_Int supp, suppl, suppr; // supressed points
  E_Int add;                // added points
  E_Float eh, ehr;          // enforce height
  E_Float y0;
  PyObject* array;
  E_Int verbose;
  if (!PYPARSETUPLE_(args, O_ RR_ TII_ I_,
                    &array, &y0, &eh, &supp, &add, &verbose)) return NULL;

  bool isVerbose = (verbose == 1) ? true : false;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  FldArrayF coord;
  FldArrayF snr; FldArrayF snl; FldArrayF snl2;

  E_Int api = f->getApi();

  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "enforceY: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Float hl;
    E_Int jl, j, k;
    E_Int jstart, jend;
    E_Bool pb = false;

    // Determination de l'indice de modification
    E_Int ic = 0;
    hl = y0;

    // Determination de la cellule de distrib contenant hl
    if ( (hl <= (*f)(0, posy)) || (hl >= (*f)((nj-1)*ni, posy) ))
    {
      char msg[256];
      sprintf(msg,
              "enforceY (y0=" SF_F_ ", eh=" SF_F_ "): cannot find hl in distrib, stopped.",
              y0, eh);
      PyErr_SetString(PyExc_TypeError, msg);
      return NULL;
    }
    if ( hl <= (*f)(ni,posy) )
    {
      char msg[256];
      sprintf(msg,
              "enforceY (y0=" SF_F_ ", eh=" SF_F_ "): line to enforce is in the first cell of distribution. \n=> use enforcePlusY to enforce this mesh.", y0, eh);
      PyErr_SetString(PyExc_TypeError, msg);
      return NULL;
    }
    if ( hl >= (*f)((nj-2)*ni, posy) )
    {
      char msg[256];
      sprintf(msg,
              "enforceY (y0=" SF_F_ ", eh=" SF_F_ "): line to enforce is in the last cell of distribution. \n=> use enforceMoinsY to enforce this mesh.", y0, eh);
      PyErr_SetString(PyExc_TypeError, msg);
      return NULL;
    }

    if ( hl <= (*f)(ni,posy) + eh/2. )
    {
      char msg[256];
      sprintf(msg,
              "enforceY (y0=" SF_F_ ", eh=" SF_F_ "): overlapping around node 1.\n=> decrease eh or increase hl or enforce initial mesh.", y0, eh);
      PyErr_SetString(PyExc_TypeError, msg);
      return NULL;
    }
    if ( hl >= (*f)((nj-2)*ni, posy)-eh/2. )
    {
      char msg[256];
      sprintf(msg,
              "enforceY (y0=" SF_F_ ", eh=" SF_F_ "): overlapping around node ni-2.\n=> decrease eh or increase hl or enforce initial mesh.", y0, eh);
      PyErr_SetString(PyExc_TypeError, msg);
      return NULL;
    }

    E_Float h ;
    for (j = 0; j < nj; j++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = ic + j*ni + k*ni*nj;
        h = (*f)(ind, posy);
        if (h > hl) goto exit;
      }
    exit:

    jl = j-1;
    suppl = supp ;
    suppr = supp ;

    // Verifications des parametres supp et add et modification si necessaire
    if (supp > jl)
    {
      suppl = jl;
      if (isVerbose)
      {
        MSGY;
        printf("cannot suppress " SF_D_ " points on left side.\n", supp);
        printf("supp set to " SF_D_ " on left side.\n", suppl);
      }
    }
    if (supp > nj-1-jl-1)
    {
      suppr = nj-1-jl-1;
      if (isVerbose)
      {
        MSGY;
        printf("cannot suppress " SF_D_ " points on right side.\n", supp);
        printf("supp set to " SF_D_ " on right side.\n", suppr);
      }
    }
    while ( ( suppl+add < 2 ) || ( suppr+add < 2 ) )
    {
      add++;
      pb = true;
    }
    if (pb == true && isVerbose)
    {
      MSGY;
      printf("supp+add < 2.\n");
      printf("   => add set to: " SF_D_ ".\n", add);
    }
    pb = false;
    if (add == 0) // add = 0 => risque de probleme sur la croissance ou decroissance des tailles de maille
    {
      add++;
      if (isVerbose)
      {
        MSGY;
        printf("add=0 => problems will occur on step-size.\n");
        printf("   => add set to: " SF_D_ ".\n", add);
      }
    }

    jstart = jl - suppl;
    jend = jl + suppr + 1;
    ic = 0;

    // Verification : comparaison entre eh et deltal, eh et deltar
    E_Float pt1 = (*f)(jstart*ni, posy);
    E_Float pt1a = (*f)((jstart+1)*ni, posy);
    E_Float pt2a = (*f)((jend-1)*ni, posy);
    E_Float pt2 = (*f)(jend*ni,posy);
    E_Float pt3 = hl-eh/2.;
    E_Float pt4 = hl+eh/2.;

    E_Float deltal = pt1a - pt1;
    E_Float deltar = pt2 - pt2a;

    if (add > 0 && isVerbose)
    {
      if (eh > deltal)
      {
        MSGY;
        printf("add=" SF_D_ " > 0  while eh=" SF_F_ " ",add,eh);
        printf("is greater than the first step on left side=" SF_F_ ".\n", deltal);
      }
      if (eh > deltar)
      {
        MSGY;
        printf("add=" SF_D_ " > 0  while eh=" SF_F_ " ",add,eh);
        printf("is greater than the first step on right side=" SF_F_ ".\n", deltar);
      }
    }
    if (add < 0 && isVerbose)
    {
      if (eh < deltal)
      {
        MSGY;
        printf("add=" SF_D_ " < 0  while eh=" SF_F_ " ",add,eh);
        printf("is lower than the first step on left side=" SF_F_ ".\n", deltal);
      }
      if (eh < deltar)
      {
        MSGY;
        printf("add=" SF_D_ " < 0 while eh=" SF_F_ " ",add,eh);
        printf("is lower than the last step on right side=" SF_F_ ".\n", deltar);
      }
    }

    // Distribution a gauche
    E_Int npl = suppl + add + 2;
    snl.malloc(npl);
    snl2.malloc(npl);
    k6stretch_(pt1, pt4, snl2.begin(), npl, eh, deltal, 2, 1);
    for (E_Int pp = 0 ; pp < npl ; pp++)
      snl[pp] = -snl2[npl-1-pp]+pt4+pt1;

    // Verification de la distribution : decroissance des tailles de maille a gauche si add > 0
    pb = false;
    if (add > 0)
    {
      for (j = 0 ; j < npl-2 ; j++)
      {
        if (snl[j+2]-snl[j+1] >= snl[j+1]-snl[j])
          pb = true;
      }
      if (pb == true && isVerbose)
      {
        MSGY;
        printf("non-decreasing step-size in left-distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }
    // Verification de la distribution : croissance des tailles de maille a gauche si add < 0
    pb = false;
    if (add < 0)
    {
      for (j = 0 ; j < npl-2 ; j++)
      {
        if (snl[j+2]-snl[j+1] <= snl[j+1]-snl[j])
          pb = true;
      }
      if (pb == true && isVerbose)
      {
        MSGY;
        printf("non-increasing step-size in left-distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }
    pb = false;

    // Distribution a droite
    E_Int npr = suppr + add + 2;
    snr.malloc(npr);
    k6stretch_(pt3, pt2, snr.begin(), npr, eh, deltar, 2, 1);

    // Verification de la distribution : croissance des tailles de maille a droite si add > 0
    pb = false;
    if (add > 0)
    {
      for (j = 0 ; j < npr-2 ; j++)
      {
        if (snr[j+2]-snr[j+1] <= snr[j+1]-snr[j])
          pb = true;
      }
      if (pb == true && isVerbose)
      {
        MSGY;
        printf("non-increasing step-size in right-distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }
    // Verification de la distribution : decroissance des tailles de maille a droite si add < 0
    pb = false;
    if (add < 0)
    {
      for (j = 0 ; j < npr-2 ; j++)
      {
        if (snr[j+2]-snr[j+1] >= snr[j+1]-snr[j])
          pb = true;
      }
      if (pb == true && isVerbose)
      {
        MSGY;
        printf("non-decreasing step-size in right-distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }
    pb = false;

    // Distribution finale
    E_Int np = nj+2*add;
    coord.malloc(ni*np*nk, 3);
    coord.setAllValuesAtNull();

    for (j = 0; j < jstart; j++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = ic + j*ni + k*ni*np;
        coord(ind,1) = (*f)(ind, posx);
        coord(ind,2) = (*f)(ind, posy);
        coord(ind,3) = (*f)(ind, posz);
      }

    for (j = jstart; j < jstart+npl-1; j++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = ic + j*ni + k*ni*np;
        coord(ind,1) = (*f)(ic,posx);
        coord(ind,2) = snl[j-jstart];
        coord(ind,3) = (*f)(ic,posz);
      }

    for (j = jstart+npl-1; j < jstart+npl-1+npr-1; j++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = ic + j*ni + k*ni*np;
        coord(ind,1) = (*f)(ic,posx);
        coord(ind,2) = snr[j-jstart-npl+2];
        coord(ind,3) = (*f)(ic,posz);
      }

    for (j = jstart+npl-1+npr-1; j < nj+2*add; j++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = ic + j*ni + k*ni*np;
        coord(ind,1) = (*f)(ic+(j-2*add)*ni+k*ni*nj, posx);
        coord(ind,2) = (*f)(ic+(j-2*add)*ni+k*ni*nj, posy);
        coord(ind,3) = (*f)(ic+(j-2*add)*ni+k*ni*nj, posz);
      }

    // Verification de la nouvelle distribution globale : croissance des abscisses
    pb = false;
    for (j = 0 ; j < nj+2*add-1 ; j++)
      for (k = 0 ; k < nk ; k++)
      {
        E_Int ind = j*ni + k*ni*np;
        E_Int ind2 = (j+1)*ni + k*ni*np;
        if (coord(ind,2) >= coord(ind2,2))
        {
          if (isVerbose)
          {
            MSGY;
            printf("Y coordinates are not increasing:\n");
            printf("   Y(" SF_D_ ") = " SF_F_ " and Y(" SF_D_ ") = " SF_F_ ".\n",
                   j, coord(ind,2), j+1, coord(ind2,2));
          }
          pb = true;
        }
      }
    if (pb == true && isVerbose)
      printf("      => Please change add=" SF_D_ " or supp=" SF_D_ " or eh=" SF_F_ ".\n",
             add, supp, eh);
    pb = false;

    // Verification de la nouvelle distribution globale : eh
    ehr = coord((jstart+npl-1)*ni,2)-coord((jstart+npl-2)*ni,2)-eh ;
    if ((ehr > eh/1000. || ehr < -eh/1000.) && isVerbose)
    {
      MSGY;
      printf("eh different from the specified value.\n");
      printf("in new distribution, eh = " SF_F_ ".\n", ehr+eh);
    }

    ic++;

    // For debug
    // In reality, we must parametrize the curve
    while (ic < ni)
    {
      for (j = 0; j < np; j++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = ic + j*ni + k*ni*np;
          E_Int ind0 = j*ni + k*ni*np;
          E_Int ind1 = ic + 0 + k*ni*np;
          coord(ind,1) = (*f)(ind1,posx);
          coord(ind,2) = coord(ind0,2);
          coord(ind,3) = (*f)(ind1,posz);
        }
      ic++;
    }

    RELEASESHAREDS(array, f);  
    // Build array
    PyObject* tpl = K_ARRAY::buildArray3(coord, varString, ni, np, nk, api);
    return tpl;
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "enforceY: not used for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforceY: unknown type of array.");
    return NULL;
  }
}

// ============================================================================
/* Enforce a j line in a distribution */
// ============================================================================
PyObject* K_GENERATOR::enforcePlusYMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int supp;    // supressed points
  E_Int add;  // added points
  E_Float eh, ehr;  // enforce height
  E_Int verbose;
  if (!PYPARSETUPLE_(args, O_ R_ TII_ I_,
                    &array, &eh, &supp, &add, &verbose)) return NULL;

  bool isVerbose = (verbose == 1) ? true : false;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  FldArrayF coord;
  FldArrayF sn1;

  E_Int api = f->getApi();

  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "enforcePlusY: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;
    E_Int j,k;
    E_Int jstart, jend;
    E_Bool pb = false;

    // Verification des parametres supp et add et modification si necessaire
    j = 0 ;
    jstart = 0 ;
    if (supp <= 0)
    {
      supp = 1;
      if (isVerbose)
      {
        MSGPLUSY;
        printf("supp <= 0.\n");
        printf("   =>supp set to: " SF_D_ ".\n", supp);
      }
    }
    if (supp > nj-2)
    {
      supp = nj-2;
      if (isVerbose)
      {
        MSGPLUSY;
        printf("cannot suppress " SF_D_ " points.\n", supp);
        printf("   =>supp set to: " SF_D_ ".\n", supp);
      }
    }
    jend = jstart + supp + 1 ;
    if (supp+add < 2) // supp+add < 2 => probleme sur-determine
    {
      if (isVerbose)
      {
        MSGPLUSY;
        printf("supp+add=" SF_D_ " < 2.\n", supp+add);
      }
      while (supp+add < 2) {add++;}
      if (isVerbose) printf("   => add set to: " SF_D_ ".\n", add);
    }
    if (add == 0) // add = 0 => risque de probleme sur la croissance ou decroissance des tailles de maille
    {
      add++;
      if (isVerbose)
      {
        MSGPLUSY;
        printf("add=0 => problems will occur on step-size.\n");
        printf("   => add set to: " SF_D_ ".\n", add);
      }
    }

    // Distribution finale
    E_Int np = nj+add;
    coord.malloc(ni*np*nk, 3);
    coord.setAllValuesAtNull();

    // Nouvelle distribution sur la zone modifiee : jstart <= j <= jend (jstart et jend non modifies)
    E_Int np1 = supp+add+2;
    sn1.malloc(np1);

    // Distributions pour i = 0
    E_Int ic = 0;

    // Verification : comparaison entre eh et delta1
    E_Float delta1 = (*f)(jend*ni,posy)-(*f)((jend-1)*ni,posy);

    if (add > 0)
    {
      if (eh > delta1 && isVerbose)
      {
        MSGPLUSY;
        printf("   add=" SF_D_ " >= 0 while eh=" SF_F_ "", add, eh);
        printf(" is greater than the last step=" SF_F_ ".\n", delta1);
      }
    }
    if (add < 0 && isVerbose)
    {
      if (eh < delta1)
      {
        MSGPLUSY;
        printf("   add=" SF_D_ " <= 0  while eh=" SF_F_ "", add, eh);
        printf(" is lower than the last step=" SF_F_ ".\n",  delta1);
      }
    }

    // Taille de la premiere maille = eh ; taille de la derniere maille inchangee = x(iend) - x(iend-1)
    E_Float pt1 = (*f)(jstart*ni,posy);
    E_Float pt2 = (*f)(jend*ni, posy);
    k6stretch_(pt1, pt2, sn1.begin(), np1, eh, delta1, 2, 1);

    // Verification de la distribution : croissance des points
    for (j = 0 ; j < np1-1 ; j++)
    {
      if (sn1[j+1] <= sn1[j])
      {
        if (isVerbose)
        {
            MSGPLUSY;
            printf("Y coordinates are not increasing:\n");
            printf("   Y(" SF_D_ ") = " SF_F_ "", j, sn1[j]);
            printf("   Y(" SF_D_ ") = " SF_F_ ".\n", j+1, sn1[j+1]);
        }
        pb = true;
      }
    }

    if (pb == true && isVerbose)
    {
      printf("      => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
    }

    // Verification de la distribution : croissance des tailles de maille si add > 0
    pb = false;
    if (add > 0)
    {
      for (j = 0 ; j < np1-2 ; j++)
      {
        if (sn1[j+2]-sn1[j+1] <= sn1[j+1]-sn1[j])
        {
          pb = true;
        }
      }
      if (pb == true && isVerbose)
      {
        MSGPLUSY;
        printf("non-increasing step-size in distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }

    // Verification de la distribution : decroissance des tailles de maille si add < 0
    pb = false;
    if (add < 0)
    {
      for (j = 0 ; j < np1-2 ; j++)
      {
        if (sn1[j+2]-sn1[j+1] >= sn1[j+1]-sn1[j])
          pb = true;
      }
      if (pb == true && isVerbose)
      {
        MSGPLUSY;
        printf("non-decreasing step-size in distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }

    // Creation de la nouvelle distribution
    for (j = 0 ; j < np1 ; j++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = ic + j*ni + k*ni*np;
        coord(ind,1) = (*f)(ic,posx);
        coord(ind,2) = sn1[j];
        coord(ind,3) = (*f)(ic,posz);
      }

    for (j = np1 ; j < np ; j++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = ic + j*ni + k*ni*nj;
        coord(ind,1) = (*f)((j-add)*ni+k*ni*nj,posx);
        coord(ind,2) = (*f)((j-add)*ni+k*ni*nj,posy);
        coord(ind,3) = (*f)((j-add)*ni+k*ni*nj,posz);
      }

    // Verification de la nouvelle distribution globale : croissance des ordonnees
    pb = false;
    for (j = 0 ; j < np-1 ; j++)
      for (k = 0 ; k < nk ; k++)
      {
        E_Int ind = j*ni + k*ni*np;
        E_Int ind2 = (j+1)*ni + k*ni*np;
        if (coord(ind,2) >= coord(ind2,2))
        {
          if (isVerbose)
          {
            MSGPLUSY;
            printf("Y coordinates are not increasing:\n");
            printf("   Y(" SF_D_ ") = " SF_F_ "", j, coord(ind,2));
            printf("   Y(" SF_D_ ") = " SF_F_ ".\n", j+1, coord(ind2,2));
          }
          pb = true;
        }
      }
    if (pb == true && isVerbose)
    {
      printf("      => Please change add=" SF_D_ " or supp=" SF_D_ " or eh=" SF_F_ ".\n", add, supp, eh);
    }
    pb = false;

    // Verification de la nouvelle distribution globale : eh
    ehr = coord((jstart+1)*ni,2)-coord(jstart*ni,2)-eh ;
    if ((ehr > eh/1000. || ehr < -eh/1000.) && isVerbose)
    {
      MSGPLUSY;
      printf("eh different from the specified value.\n");
      printf("   in new distribution, eh=" SF_F_ ".\n", coord(1,2)-coord(0,2));
    }

    ic++;

    // For debug
    // In reality, we must parametrize the curve
    while (ic < ni)
    {
      for (j = 0; j < np; j++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = ic + j*ni + k*ni*np;
          E_Int ind0 = j*ni + k*ni*np;
          E_Int ind1 = ic + 0 + k*ni*nj;
          coord(ind,1) = (*f)(ind1,posx);
          coord(ind,2) = coord(ind0,2);
          coord(ind,3) = (*f)(ind1,posz);
        }
      ic++;
    }
    RELEASESHAREDS(array, f);
    // Build array
    PyObject* tpl = K_ARRAY::buildArray3(coord, varString, ni, np, nk, api);
    return tpl;
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "enforcePlusY: not used for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforcePlusY: invalid type of array.");
    return NULL;
  }
}

// ============================================================================
/* Enforce a j line in a distribution */
// ============================================================================
PyObject* K_GENERATOR::enforceMoinsYMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int supp;    // supressed points
  E_Int add;  // added points
  E_Float eh, ehr;  // enforce height
  E_Int verbose;
  if (!PYPARSETUPLE_(args, O_ R_ TII_ I_,
                    &array, &eh, &supp, &add, &verbose)) return NULL;

  bool isVerbose = (verbose == 1) ? true : false;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  FldArrayF coord;
  FldArrayF sn1, sn2;

  E_Int api = f->getApi();

  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_TypeError,
                        "enforceMoinsY: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int j,k;
    E_Int jstart, jend;
    E_Bool pb = false;

    // Verification des parametres supp et add et modification si necessaire
    jend = nj-1;
    if (supp <= 0)
    {
      supp = 1;
      if (isVerbose)
      {
        MSGMOINSY;
        printf("supp <= 0.\n");
        printf("   =>supp set to: supp=" SF_D_ ".\n", supp);
      }
    }
    if (supp > nj-2)
    {
      supp = nj-2;
      if (isVerbose)
      {
        MSGMOINSY;
        printf("cannot suppress " SF_D_ " points.\n", supp);
        printf("   =>supp set to: supp=" SF_D_ ".\n", supp);
      }
    }
    jstart = jend - supp - 1;
    if (supp+add < 2) // supp+add < 2 => probleme sur-determine
    {
      if (isVerbose)
      {
        MSGMOINSY;
        printf("supp+add=" SF_D_ " < 2.\n", supp+add);
      }
      while (supp+add < 2) {add++;}
      if (isVerbose) printf("   => add set to: add=" SF_D_ ".\n", add);
    }
    if (add == 0) // add = 0 => risque de probleme sur la croissance ou decroissance des tailles de maille
    {
      add++;
      if (isVerbose)
      {
        MSGMOINSY;
        printf("add=0 => problems will occur on step-size.\n");
        printf("   => add set to: add=" SF_D_ ".\n", add);
      }
    }

    // Distribution finale
    E_Int np = nj+add;
    coord.malloc(ni*np*nk, 3);
    coord.setAllValuesAtNull();

    // Nouvelle distribution sur la zone modifiee : jstart <= j <= jend (jstart et jend non modifies)
    E_Int np1 = supp+add+2;
    sn1.malloc(np1);
    sn2.malloc(np1);

    // Distributions pour i = 0
    E_Int ic = 0;

    // Verification : comparaison entre eh et delta1
    E_Float delta1 = (*f)((jstart+1)*ni,posy)-(*f)(jstart*ni,posy);

    if (add > 0 && isVerbose)
    {
      if (eh > delta1)
      {
        MSGMOINSY;
        printf("   add=" SF_D_ " > 0 while eh=" SF_F_ " is greater than the last step=" SF_F_ ".\n", add, eh, delta1);
      }
    }
    if (add < 0 && isVerbose)
    {
      if (eh < delta1)
      {
        MSGMOINSY;
        printf("   add=" SF_D_ " < 0 while eh=" SF_F_ " is lower than the last step=" SF_F_ ".\n", add, eh, delta1);
      }
    }

    // Taille de la premiere maille = delta1 ; taille de la derniere maille inchangee = x(iend) - x(iend-1)
    E_Float pt1 = (*f)(jstart*ni, posy);
    E_Float pt2 = (*f)(jend*ni, posy);
    k6stretch_(pt1, pt2, sn2.begin(), np1, eh, delta1, 2, 1);
    for (j = 0 ; j < np1 ; j++)
      sn1[j] = -sn2[np1-1-j]+pt2+pt1;

    // Verification de la distribution : decroissance des tailles de maille si add > 0
    pb = false;
    if (add > 0)
    {
      for (j = 0 ; j < np1-2 ; j++)
      {
        if (sn1[j+2]-sn1[j+1] >= sn1[j+1]-sn1[j])
        {
          pb = true; break;
        }
      }
      if (pb == true && isVerbose)
      {
        MSGMOINSY;
        printf("Non-decreasing step-size in distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }

    // Verification de la distribution : croissance des tailles de maille si add < 0
    pb = false;
    if (add < 0)
    {
      for (j = 0 ; j < np1-2 ; j++)
      {
        if (sn1[j+2]-sn1[j+1] <= sn1[j+1]-sn1[j])
        {
          pb = true; break;
        }
      }
      if (pb == true && isVerbose)
      {
        MSGMOINSY;
        printf("non-increasing step-size in distribution.\n");
        printf("   => Please change add=" SF_D_ " or eh=" SF_F_ ".\n", add, eh);
      }
    }

    // Creation de la nouvelle distribution
    for (j = 0 ; j < jstart ; j++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = ic + j*ni + k*ni*np;
        E_Int ind0 = ic+j*ni + k*ni*np;
        coord(ind,1) = (*f)(ind0, posx);
        coord(ind,2) = (*f)(ind0, posy);
        coord(ind,3) = (*f)(ind0, posz);
      }

    for (j = jstart ; j < np ; j++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = ic + j*ni + k*ni*nj;
        coord(ind,1) = (*f)(0, posx);
        coord(ind,2) = sn1[j-jstart];
        coord(ind,3) = (*f)(0, posz);
      }

    // Verification de la nouvelle distribution globale : croissance des ordonnees
    pb = false;
    for (j = 0 ; j < np-1 ; j++)
      for (k = 0 ; k < nk ; k++)
      {
        E_Int ind = j*ni + k*ni*np;
        E_Int ind2 = (j+1)*ni + k*ni*np;
        if (coord(ind,2) >= coord(ind2,2))
        {
          if (isVerbose)
          {
            MSGMOINSY;
            printf("Y coordinates are not increasing:\n");
            printf("   Y(" SF_D_ ") = " SF_F_ "", j, coord(ind,2));
            printf("   Y(" SF_D_ ") = " SF_F_ ".\n", j+1, coord(ind2,2));
          }
          pb = true;
        }
      }
    if (pb == true && isVerbose)
    {
      printf("      => Please change add=" SF_D_ ".\n", add);
      printf("      => or supp=" SF_D_ ".\n", supp);
      printf("      => or eh=" SF_F_ ".\n", eh);
    }
    pb = false;

    // Verification de la nouvelle distribution globale : eh
    ehr = coord((np-1)*ni, 2) - coord((np-2)*ni, 2) - eh;
    if (ehr > eh/1000. || ehr < -eh/1000.)
    {
      char msg[256];
      sprintf(msg,
              "enforceMoinsY (eh=" SF_F_ "): eh different from the specified value in new distribution eh=" SF_F_ ".", eh, coord((np-1)*ni,2)-coord((np-2)*ni,2));
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError, msg);
      return NULL;
    }

    ic++;

    // For debug
    // In reality, we must parametrize the curve
    while (ic < ni)
    {
      for (j = 0; j < np; j++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = ic + j*ni + k*ni*np;
          E_Int ind0 = j*ni + k*ni*np;
          E_Int ind1 = ic + 0 + k*ni*nj;
          coord(ind,1) = (*f)(ind1, posx);
          coord(ind,2) = coord(ind0,2);
          coord(ind,3) = (*f)(ind1, posz);
        }
      ic++;
    }
    RELEASESHAREDS(array, f);
    // Build array
    PyObject* tpl = K_ARRAY::buildArray3(coord, varString, ni, np, nk, api);
    return tpl;
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "enforceMoinsY: not used for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforceMoinsY: invalid type of array.");
    return NULL;
  }
}

// ============================================================================
/* Enforce a line in a distribution */
// ============================================================================
PyObject* K_GENERATOR::enforceLineMesh(PyObject* self, PyObject* args)
{
  PyObject* array1; PyObject* array2;
  E_Int supp; // supressed points
  E_Int add;  // added points
  E_Float eh; // enforce heigh
  const E_Float EPS = 1.e-3;
  if (!PYPARSETUPLE_(args, OO_ R_ TII_,
                    &array1, &array2, &eh, &supp, &add)) return NULL;
                    
  // Check array
  E_Int ni, nj, nk, im2, jm2, km2;
  FldArrayF* f1; FldArrayF* f2;
  FldArrayI* cn1; FldArrayI* cn2;
  char *varString1; char *varString2;
  char* eltType1; char* eltType2;
  FldArrayF coord;
  vector<E_Int> pos1; vector<E_Int> pos2;
  FldArrayF sn1, sn2, sn3;

  // Maillage de base a modifier
  E_Int res1 =
    K_ARRAY::getFromArray3(array1, varString1, f1, ni, nj, nk,
                           cn1, eltType1);
  // Courbe a enforcer
  E_Int res2 =
    K_ARRAY::getFromArray3(array2, varString2, f2, im2, jm2, km2,
                           cn2, eltType2);

  E_Int posx1, posy1, posz1;
  E_Int posx2, posy2, posz2;

  E_Int api = f1->getApi();

  if (res1 == 1 && res2 == 1)
  {
    char* varString = new char[strlen(varString1)+strlen(varString2)+4];
    E_Int res0 =
      K_ARRAY::getPosition(varString1, varString2, pos1, pos2, varString);
    if (res0 == -1)
    {
      RELEASESHAREDS(array1, f1);
      RELEASESHAREDS(array2, f2);
      PyErr_SetString(PyExc_TypeError,
                      "enforceLine: common variables list is empty.");
      return NULL;
    }
    else if (res0 == 0 )
    {
      printf("Warning: enforceLine: some variables are different.\n");
      printf("Only common variables are kept.\n");
    }
    posx1 = K_ARRAY::isCoordinateXPresent(varString1);
    posy1 = K_ARRAY::isCoordinateYPresent(varString1);
    posz1 = K_ARRAY::isCoordinateZPresent(varString1);
    posx2 = K_ARRAY::isCoordinateXPresent(varString2);
    posy2 = K_ARRAY::isCoordinateYPresent(varString2);
    posz2 = K_ARRAY::isCoordinateZPresent(varString2);

    if (posx1 == -1 || posy1 == -1 || posz1 == -1 ||
        posx2 == -1 || posy2 == -1 || posz2 == -1)
    {
      RELEASESHAREDS(array1, f1);
      RELEASESHAREDS(array2, f2);
      PyErr_SetString(PyExc_TypeError,
                      "enforceLine: coordinates not found.");
      return NULL;
    }

    posx1++; posy1++; posz1++;
    posx2++; posy2++; posz2++;

    E_Float sl, hl;
    E_Int il,jl,j,k;
    E_Int jstart, jend;
    E_Float ehr = eh;

    // Distribution finale
    E_Int np = nj+2*add;
    coord.malloc(ni*np*nk, 3);
    coord.setAllValuesAtNull();

    // Determination du tube
    E_Int ic = 0;
    E_Int jmin = 1000000;
    E_Int jmax = -1000000;

    E_Float* f1x = f1->begin(posx1);
    E_Float* f1y = f1->begin(posy1);
    E_Float* f1z = f1->begin(posz1);
    E_Float* f2x = f2->begin(posx2);
    E_Float* f2y = f2->begin(posy2);
    //E_Float* f2z = f2->begin(posz2);

    while (ic < ni)
    {
      E_Float sc = f1x[ic];

      // Determine hl : le y de la courbe a enforcer
      for (il = 0; il < im2; il++)
      {
        sl = f2x[il];
        if (sl >= sc - EPS) break;
      }
      if (il == im2)
      {
        RELEASESHAREDS(array1, f1);
        RELEASESHAREDS(array2, f2);
        char msg[256];
        sprintf(msg,
                "enforceLine: cannot find sc=" SF_F_ " in line, stopped.", sc);
        PyErr_SetString(PyExc_TypeError, msg);
        return NULL;
      }
      else if (il == 0) hl = f2y[il];
      else hl = f2y[il-1];

      // to avoid problems
      hl = E_max(hl, ehr);

      // Determination de la cellule de distrib contenant hl
      for (j = 0; j < nj; j++)
      for (k = 0; k < nk; k++)
      {
        E_Int ind = ic + j*ni + k*ni*nj;
        E_Float h = f1y[ind];
        if (h > hl) goto exit;
      }
      exit:
      if (j == nj)
      {
        delete f1; delete f2;
        PyErr_SetString(PyExc_TypeError,
                        "enforceLine: cannot find hl in distrib.");
        return NULL;
      }
      j = E_max(j-1,0);

      jl = j;
      jstart = E_max(jl - supp, 0);
      jend = E_min(jl + supp, nj)+1;
      jmin = E_min(jstart, jmin);
      jmax = E_max(jend, jmax);
      ic++;
    }

    sn1.malloc(add+jmax-jmin+2);
    sn2.malloc(add+jmax-jmin+2);
    sn3.malloc(add+jmax-jmin+2);

    // Distributions
    jstart = jmin;
    jend = jmax;
    jl = E_Int(jmin + (jmax-jmin)*0.5);
    ic = 0;
    while (ic < ni)
    {
      E_Float sc = f1x[ic];

      // Determine la hauteur du resserrement
      for (il = 0; il < im2; il++)
      {
        sl = f2x[il];
        if (sl >= sc - EPS) break;
      }
      if (il == im2)
      {
        RELEASESHAREDS(array1, f1);
        RELEASESHAREDS(array2, f2);
        char msg[256];
        sprintf(msg,
                "enforceLine: cannot find sc=" SF_F_ " in line, stopped.", sc);
        PyErr_SetString(PyExc_TypeError,
                        msg);
        return NULL;
      }
      else if (il == 0) hl = f2y[il];
      else hl = f2y[il-1];
      //printf("index " SF_D_ ": y=" SF_F_ "\n", ic, hl);

      // to avoid problems
      hl = E_max(hl, ehr);

      // Diminish the eh considering delta
      E_Float delta1 =
        f1y[ic+jmax*ni] - f1y[ic+(jmax-1)*ni];
      E_Float delta2 =
        f1y[ic+(jmin+1)*ni] - f1y[ic+jmin*ni];
      eh = E_min(ehr*0.5, 0.5*delta1);
      eh = E_min(eh, 0.5*delta2);

      // premiere distribution
      jstart = jl;
      jend = jmax;
      E_Float pt2a = f1y[ic+(jend-1)*ni];
      E_Float pt2 = f1y[ic+jend*ni];
      E_Float pt1 = hl-eh;
      E_Float pt1a = f1y[ic+(jstart+1)*ni];
      E_Int np1 = add+jend-jstart+1;
      k6stretch_(pt1, pt2, sn1.begin(), np1, 2*eh, pt2-pt2a, 2, 1);

      // deuxieme distribution
      jstart = jmin;
      jend = jl+1;
      pt1 = f1y[ic+jstart*ni];
      pt1a = f1y[ic+(jstart+1)*ni];
      pt2a = f1y[ic+(jend-1)*ni];
      pt2 = hl+eh;
      E_Int np2 = add+jend-jstart+1;
      k6stretch_(pt1, pt2, sn3.begin(), np2, 2*eh, pt1a-pt1, 2, 1);

      for (E_Int i = 0; i < np2; i++) sn2[i] = -sn3[np2-i-1]+pt2+pt1;

      jend = jmax;

      for (j = 0; j < jstart; j++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = ic + j*ni + k*ni*np;
          E_Int ind0 = ic+j*ni+k*ni*nj;
          coord(ind,1) = f1x[ind0];
          coord(ind,2) = f1y[ind0];
          coord(ind,3) = f1z[ind0];
        }

      for (j = jstart; j < jstart+np2; j++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = ic + j*ni + k*ni*np;
          coord(ind,1) = f1x[ic];
          coord(ind,2) = sn2[j-jstart];
          coord(ind,3) = f1z[ic];
        }

      for (j = jstart+np2; j < jstart+np2+np1-2; j++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = ic + j*ni + k*ni*np;
          coord(ind,1) = f1x[ic];
          coord(ind,2) = sn1[j-jstart-np2+2];
          coord(ind,3) = f1z[ic];
        }

      for (j = jstart+np1+np2-2; j < nj+2*add; j++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = ic + j*ni + k*ni*np;
          E_Int ind0 = ic+(j-2*add)*ni+k*ni*nj;
          coord(ind,1) = f1x[ind0];
          coord(ind,2) = f1y[ind0];
          coord(ind,3) = f1z[ind0];
        }
      ic++;
    }

    RELEASESHAREDS(array1, f1);
    RELEASESHAREDS(array2, f2);  
    // Build array
    PyObject* tpl = K_ARRAY::buildArray3(coord, varString, ni, np, nk, api);
    delete [] varString;
    return tpl;
  }
  else if (res1 == 2 || res2 == 2)
  {
    RELEASESHAREDB(res1, array1, f1, cn1);
    RELEASESHAREDB(res2, array2, f2, cn2);      
    PyErr_SetString(PyExc_TypeError,
                    "enforceLine: not used for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforceLine: invalid type of array.");
    return NULL;
  }
}

// ============================================================================
/* Enforce a point in a distribution */
// ============================================================================
PyObject* K_GENERATOR::enforcePoint(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0; // enforced point
  if (!PYPARSETUPLE_(args, O_ R_, &array, &x0)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  FldArrayF coord;

  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);

  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "enforcePoint: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int i, j, k;
    E_Int api = f->getApi();

    // Distribution finale
    coord.malloc(ni*nj*nk, 3);
    coord.setAllValuesAtNull();

    // Determination du point le plus proche de x0
    E_Int jc = 0;
    E_Int kc = 0;
    E_Float dmax = 1.e6;
    E_Float xs = 0.;
    E_Int is = 0;

    for (i = 0; i < ni; i++)
    {
      E_Int ind = i + jc*ni + kc*ni*nj;
      E_Float x = (*f)(ind,posx);
      E_Float d = E_abs(x-x0);
      if (d < dmax)
      {
        dmax = d;
        xs = x;
        is = i;
      }
    }

    E_Float epsilon = x0-xs; // shift vector
    E_Float delta;
    if (is < ni)
      delta = (*f)(is+1,posx)-xs;
    else
      delta = xs - (*f)(is-1,posx);

    if (E_abs(epsilon) < 1.e-12) // nothing to do
    {
      PyObject* tpl = K_ARRAY::buildArray3(*f, varString, ni, nj, nk, api);
      RELEASESHAREDS(array, f);
      return tpl;
    }

    // Determination de sigma
    E_Float sigma = 1./delta*10;

    if (epsilon > 0)
    {
      while (sigma > 1./delta)
      {
        E_Bool cool = true;
        for (i = 0; i < ni-1; i++)
        {
          E_Float x2 = (*f)(i+1,posx);
          E_Float x1 = (*f)(i,posx);
          E_Float alph2 = exp(-sigma*(x2-xs)*(x2-xs));
          E_Float alph1 = exp(-sigma*(x1-xs)*(x1-xs));
          if (alph2 <= alph1 + (x1 - x2)/epsilon )
          {
            cool = false;
            break;
          }
        }
        if (cool == true) goto exit;
        sigma = sigma - 1./delta;
      }
    }
    else
    {
      while (sigma > 1./delta)
      {
        E_Bool cool = true;
        for (i = 0; i < ni-1; i++)
        {
          E_Float x2 = (*f)(i+1,posx);
          E_Float x1 = (*f)(i,posx);
          E_Float alph2 = exp(-sigma*(x2-xs)*(x2-xs));
          E_Float alph1 = exp(-sigma*(x1-xs)*(x1-xs));
          if (alph2 >= alph1 + (x1 - x2)/epsilon )
          {
            cool = false;
            break;
          }
        }
        if (cool == true) goto exit;
        sigma = sigma - 1./delta;
      }
    }
    exit:

    // Shift coordinates
    for (i = 0; i < ni; i++)
      for (j = 0; j < nj; j++)
        for (k = 0; k < nk; k++)
        {
          E_Int ind = i + j*ni + k*ni*nj;
          E_Float x = (*f)(ind,posx);
          E_Float alph = exp(-sigma*(x-xs)*(x-xs));

          coord(ind,1) = (*f)(ind,posx) + alph*epsilon;
          coord(ind,2) = (*f)(ind,posy);
          coord(ind,3) = (*f)(ind,posz);
        }

    // Check
    if (E_abs(coord(0,1)-(*f)(0,posx)) > 1.e-12)
      printf("Warning: enforcePoint: a bound has been moved.");
    if (E_abs(coord(ni-1,1)-(*f)(ni-1,posx)) > 1.e-12)
      printf("Warning: enforcePoint: a bound has been moved.\n");

    printf("Info: enforcePoint: index of enforced point is: " SF_D_ ".\n", is);

    RELEASESHAREDS(array, f);
    // Build array
    PyObject* tpl = K_ARRAY::buildArray3(coord, varString, ni, nj, nk, api);
    return tpl;
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);    
    PyErr_SetString(PyExc_TypeError,
                    "enforcePoint: not used for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforcePoint: invalid type of array.");
    return NULL;
  }
}

// ============================================================================
/* Add a point in in a distribution */
// ============================================================================
PyObject* K_GENERATOR::addPointInDistribution(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int ind;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &ind)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  FldArrayF coord;

  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);

  E_Int api = f->getApi();

  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "addPointInDistribution: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int nip1 = ni+1;
    E_Int i, j, k, l, m;
    E_Float dx;

    if (ind > ni)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "addPointInDistribution: given index is greater than distribution number of points.");
      return NULL;
    }

    // Distribution finale
    coord.malloc(nip1*nj*nk, 3);
    coord.setAllValuesAtNull();

    // Comp of shift
    m = ind;
    if (ind < ni)
      dx = (*f)(m,posx) - (*f)(m-1,posx);
    else
      dx = (*f)(m-1,posx) - (*f)(m-2,posx);

    // Recopie
    for (j = 0; j < nj; j++)
      for (k = 0; k < nk; k++)
      {
        for (i = 0; i < ind; i++)
        {
          l = i + j*nip1 + k*nip1*nj;
          m = i + j*ni + k*ni*nj;
          coord(l, 1) = (*f)(m,posx);
          coord(l, 2) = (*f)(m,posy);
          coord(l, 3) = (*f)(m,posz);
        }
        l = ind + j*nip1 + k*nip1*nj;
        coord(l, 1) = coord(l-1, 1) + dx;
        coord(l, 2) = coord(l-1, 2);
        coord(l, 3) = coord(l-1, 3);

        for (i = ind; i < ni; i++)
        {
          l = i + j*nip1 + k*nip1*nj;
          m = i + j*ni + k*ni*nj;
          coord(l+1, 1) = (*f)(m,posx) + dx;
          coord(l+1, 2) = (*f)(m,posy);
          coord(l+1, 3) = (*f)(m,posz);
        }
      }

    // Contract to stay in [0,1]
    E_Float inv = 1./(1.+dx);
    for (i = 0; i < nip1; i++)
      for (j = 0; j < nj; j++)
        for (k = 0; k < nk; k++)
        {
          l = i + j*nip1 + k*nip1*nj;
          coord(l, 1) = coord(l, 1) * inv;
        }

    RELEASESHAREDS(array, f);
    // Build array
    PyObject* tpl = K_ARRAY::buildArray3(coord, varString, nip1, nj, nk, api);
    return tpl;
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "addPointInDistribution: not used for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "addPointInDistribution: invalid array.");
    return NULL;
  }
}

