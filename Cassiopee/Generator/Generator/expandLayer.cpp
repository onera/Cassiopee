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

# include "generator.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Modify the indicator to expand layers
  checkType = 0 : none
  checkType = 1 : check blanking of neighbouring cells
  checkType = 2 : check spacing of neighbouring cells
  checkType = 3 : raffine si un voisin est masque et la cellule voisine masquee
  est plus fine que la cellule courante.
*/
//=============================================================================
PyObject* K_GENERATOR::modifyIndicToExpandLayer(PyObject* self, PyObject* args)
{
  PyObject *octree, *indicator; 
  E_Int level, corners, checkType;
  if (!PYPARSETUPLE_(args, OO_ III_, &octree, &indicator, &level, &corners,
                    &checkType)) return NULL;
  if (level < 0) {printf("Warning: expandLayer: level is set to 0.\n"); level = 0;}

  // Verif octree HEXA/QUAD
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(octree, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  if (res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "expandLayer: array must be unstructured.");
    RELEASESHAREDB(res, octree, f, cn); return NULL;   
  }
 
  if (strcmp(eltType, "HEXA") != 0 && strcmp(eltType, "QUAD") != 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "expandLayer: array must be HEXA or QUAD.");
    RELEASESHAREDU(octree, f, cn); return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "expandLayer: coordinates not found in array.");
    RELEASESHAREDU(octree, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  // Verif indicator
  E_Int nii, nji, nki;
  FldArrayF* fi; FldArrayI* cni;
  char* varStringi; char* eltTypei;
  E_Int resi = K_ARRAY::getFromArray(indicator, varStringi, fi, nii, nji, nki, 
                                     cni, eltTypei, true);
  if (resi != 1 && resi != 2) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "expandLayer: indic array is not valid.");
    RELEASESHAREDU(octree, f, cn); return NULL;
  }
  E_Int posi = K_ARRAY::isNamePresent("indicator", varStringi);
  if (posi == -1) 
  { 
    RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: expandLayer: no refinement indicator given. Nothing done."); 
    return indicator;
  }
  posi++;
  if (fi->getSize() != cn->getSize())
  {
    RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: expandLayer: refinement indicator size must be equal to the number of elements. Nothing done."); 
    return indicator;
  }
  /*-----------FIN DES VERIFS ------------------*/
  E_Int nelts = cn->getSize(); E_Int npts = f->getSize();
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  // calcul du pas sur chaque grille 
  FldArrayF dht(nelts);
  E_Float* dhtp = dht.begin();
  E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
  E_Float dhmin = K_CONST::E_MAX_FLOAT;

  E_Int nthreads = __NUMTHREADS__;
  E_Float* dhminl = new E_Float [nthreads];
  for (E_Int i = 0; i < nthreads; i++) dhminl[i] = K_CONST::E_MAX_FLOAT;
    
#pragma omp parallel default(shared)
  {
    E_Int ind1, ind2;
    E_Int ithread = __CURRENT_THREAD__;
    
#pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    {
      ind1 = cn1[et]-1; ind2 = cn2[et]-1;
      dhtp[et] = xt[ind2]-xt[ind1];   
      dhminl[ithread] = K_FUNC::E_min(dhminl[ithread], dhtp[et]);
    }
  }
  for (E_Int i = 0; i < nthreads; i++) dhmin = K_FUNC::E_min(dhmin, dhminl[i]);

  delete [] dhminl;

  vector< vector<E_Int> > cEEN(nelts);
  getNeighbourElts(npts, xt, yt, zt, *cn, cEEN, corners, dhmin); 
  E_Float* indict = fi->begin(posi);
  E_Int posc = K_ARRAY::isCellNatureField1Present(varStringi);
  E_Float* cellNp = NULL;
  if (posc > -1)
  {
    posc++;
    cellNp = fi->begin(posc);
  }

  // detection des elements de niveau l 
  E_Float eps = 1.e-10; E_Float dhl = pow(2.,level) * dhmin;
  E_Boolean voisinBlanked;
  E_Int etv; E_Float dhet; E_Float dhetv;
  E_Float dhleps = dhl + eps;

  if (checkType == 1) // check blanking only
  {
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (cellNp[et] == 0.)
      {
        vector<E_Int>& voisins = cEEN[et];
        for (size_t nov = 0; nov < voisins.size(); nov++)
        {
          etv = voisins[nov];
          dhetv = dhtp[etv];
          if (cellNp[etv]==1.)
          {
            if (dhet > dhetv+eps) indict[et]=1.;           
            else if (dhetv > dhet+eps) indict[etv]=1.;  
          }
        }
      }
    }
  }
  else if (checkType == 2)
  {
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (K_FUNC::fEqualZero(dhet-dhl,eps) == true) // niveau l?
      {
        vector<E_Int>& voisins = cEEN[et];
        for (size_t nov = 0; nov < voisins.size(); nov++)
        {
          etv = voisins[nov];
          dhetv = dhtp[etv];
          if (dhetv > dhleps) indict[etv] = 1.; // raffine le voisin si niveau du voisin plus grand
        }
      }
    }
  }
  else if (checkType == 3) // 
  {
    // premiere passe pour les points non masques
    // on les raffine si la cellule est plus grossiere qu'une cellule de corps voisine
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (cellNp[et] > 0.1) // pas masque
      {
        // voisine masquee?
        voisinBlanked = false;
        E_Float voisinStep = K_CONST::E_MAX_FLOAT;
        vector<E_Int>& voisins = cEEN[et];
        for (size_t nov = 0; nov < voisins.size(); nov++)
        {
          etv = voisins[nov];
          if (cellNp[etv] == 0.) 
          {
            voisinBlanked = true;
            dhetv = dhtp[etv];
            if (dhetv < voisinStep) voisinStep = dhetv;
          }
        }
        if (voisinBlanked && voisinStep < dhet-eps) indict[et] = 1.;
      }
    }
    // deuxieme passe pour les points masques
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (K_FUNC::E_abs(cellNp[et]) < eps) // masque
      {
        // voisine non masquee?
        E_Boolean voisinNonBlanked = false;
        E_Float voisinStep = K_CONST::E_MAX_FLOAT;
        vector<E_Int>& voisins = cEEN[et];
        for (size_t nov = 0; nov < voisins.size(); nov++)
        {
          etv = voisins[nov];
          if (cellNp[etv] > 0.1) 
          {
            dhetv = dhtp[etv];
            if (indict[etv] == 1) dhetv = dhetv*0.5;
            voisinNonBlanked = true;
            if (dhetv < voisinStep) voisinStep = dhetv;
          }
        }
        if (voisinNonBlanked && voisinStep < dhet-eps) indict[et] = 1.;
      }
    }
  }
  else if (checkType == 4) // extension au deuxieme niveau
  {
    for (E_Int et = 0; et < nelts; et++) indict[et] = K_CONST::E_MAX_FLOAT;
    
    // propage h dans indict pour les cellules voisines des pts masques
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (cellNp[et] > 0.5) // non masque
      {
        // voisine masquee?
        voisinBlanked = false;
        E_Float voisinStep = K_CONST::E_MAX_FLOAT;
        vector<E_Int>& voisins = cEEN[et];
        for (size_t nov = 0; nov < voisins.size(); nov++)
        {
          etv = voisins[nov];
          if (cellNp[etv] == 0.) 
          {
            dhetv = dhtp[etv];
            voisinBlanked = true;
            if (dhetv < voisinStep) voisinStep = dhetv;
          }
        }
        if (voisinBlanked) indict[et] = min(indict[et],voisinStep);
      }
    }
    // propage voisin du voisin
    E_Float* h2 = new E_Float [nelts];
    for (E_Int et = 0; et < nelts; et++) h2[et] = indict[et];
    
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (cellNp[et] > 0.5) // non masque
      {
        // voisine non masquee?
        voisinBlanked = false;
        E_Float voisinStep = 1.e10;
        vector<E_Int>& voisins = cEEN[et];
        for (size_t nov = 0; nov < voisins.size(); nov++)
        {
          etv = voisins[nov];
          if (cellNp[etv] > 0.1 && indict[etv] < 1.e+10) 
          {
            dhetv = indict[etv];
            if (dhetv < voisinStep) voisinStep = dhetv;
          }
          if (K_FUNC::E_abs(cellNp[etv]) < eps) voisinBlanked = true;
        }
        if (voisinBlanked == false) h2[et] = min(indict[et],voisinStep);
      }
    }
    // Raffine si la cellule est plus grande que la valeur contenue, qui
    // correspond a la taille de la cellule paroi la plus proche
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (h2[et] < 1.e+10 && dhet > h2[et]*1.1) indict[et] = 1.;
      else indict[et] = 0.;
      //indict[et] = 0.;
    }
    delete [] h2;
  }

  else if (checkType == 5) // pure verification que les voisins des pts masque ont bien le meme niveau
  {
    E_Int pb = 0;
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (cellNp[et] > 0.5) // non masque
      {
        voisinBlanked = false;
        vector<E_Int>& voisins = cEEN[et];
        for (size_t nov = 0; nov < voisins.size(); nov++)
        {
          etv = voisins[nov];
          if (K_FUNC::E_abs(cellNp[etv]) < eps) // voisin non masque doit avoir la meme taille de maille 
          {
            dhetv = dhtp[etv];
            if (dhetv < dhet-eps) pb++;
            if (dhetv > dhet+eps) pb++;
            if (dhetv < dhet-eps || dhetv > dhet+eps)
            {
              E_Int ind1 = cn1[et]-1;
              printf("%g %g %g -> point %g %g %g\n", dhet, dhetv, cellNp[etv], xt[ind1], yt[ind1], zt[ind1]);
            }
          }
        }
      }
    }
    printf("==> Pbs: " SF_D_ "\n", pb);
  }

  else if (checkType == 6) // extension au troisieme niveau
  {
    for (E_Int et = 0; et < nelts; et++) indict[et] = K_CONST::E_MAX_FLOAT;
    
    // propage h dans indict pour les cellules voisines des pts masques
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (cellNp[et] > 0.5) // non masque
      {
        // voisine masquee?
        voisinBlanked = false;
        E_Float voisinStep = K_CONST::E_MAX_FLOAT;
        vector<E_Int>& voisins = cEEN[et];
        for (size_t nov = 0; nov < voisins.size(); nov++)
        {
          etv = voisins[nov];
          if (cellNp[etv] == 0.) 
          {
            dhetv = dhtp[etv];
            voisinBlanked = true;
            if (dhetv < voisinStep) voisinStep = dhetv;
          }
        }
        if (voisinBlanked) indict[et] = min(indict[et],voisinStep);
      }
    }
    // propage voisin du voisin
    E_Float* h2 = new E_Float [nelts];
    for (E_Int et = 0; et < nelts; et++) h2[et] = indict[et];
    
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (cellNp[et] > 0.5) // non masque
      {
        // voisine non masquee?
        voisinBlanked = false;
        E_Float voisinStep = 1.e10;
        vector<E_Int>& voisins = cEEN[et];
        for (size_t nov = 0; nov < voisins.size(); nov++)
        {
          etv = voisins[nov];
          if (cellNp[etv] > 0.1 && indict[etv] < 1.e+10) 
          {
            dhetv = indict[etv];
            if (dhetv < voisinStep) voisinStep = dhetv;
          }
          if (K_FUNC::E_abs(cellNp[etv]) < eps) voisinBlanked = true;
        }
        if (voisinBlanked == false) h2[et] = min(indict[et],voisinStep);
      }
    }
    // propage voisin du voisin du voisin (promis apres on arrete)
    E_Float* h3 = new E_Float [nelts];
    for (E_Int et = 0; et < nelts; et++) h3[et] = h2[et];
    
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (cellNp[et] > 0.5) // non masque
      {
        // voisine non masquee?
        voisinBlanked = false;
        E_Float voisinStep = 1.e10;
        vector<E_Int>& voisins = cEEN[et];
        for (size_t nov = 0; nov < voisins.size(); nov++)
        {
          etv = voisins[nov];
          if (cellNp[etv] > 0.1 && h2[etv] < 1.e+10) 
          {
            dhetv = h2[etv];
            if (dhetv < voisinStep) voisinStep = dhetv;
          }
          if (K_FUNC::E_abs(cellNp[etv]) < eps) voisinBlanked = true;
        }
        if (voisinBlanked == false) h3[et] = min(h2[et],voisinStep);
      }
    }
    // Raffine si la cellule est plus grande que la valeur contenue, qui
    // correspond a la taille de la cellule paroi la plus proche
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (h3[et] < 1.e+10 && dhet > h3[et]*1.1) indict[et] = 1.;
      else indict[et] = 0.;
      //indict[et] = 0.;
    }
    delete [] h2;
    delete [] h3;
  }

  /*-----------CONSTRUCTION ARRAY DE SORTIE ------------------*/
  PyObject* tpl;
  if (resi == 1) 
    tpl = K_ARRAY::buildArray(*fi, varStringi, nii, nji, nki);
  else 
    tpl = K_ARRAY::buildArray(*fi, varStringi, *cni, -1, eltTypei, false);
  RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn);
  return tpl;
}
