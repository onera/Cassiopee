/*    
    Copyright 2013-2019 Onera.

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
#ifdef E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOlll", &octree, &indicator, &level, &corners, &checkType)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "OOiii", &octree, &indicator, &level, &corners, &checkType)) return NULL;
#endif
  if (level < 0) {printf("Warning: expandLayer: level is set to 0.\n"); level = 0;}

  printf("expanding=%d\n", checkType);

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
    for (E_Int et = 0; et < nelts; et++)
    {
      dhet = dhtp[et];
      if (cellNp[et] > 0.1) // pas masque
      {
        // voisine masquee?
        E_Boolean voisinBlanked = false;
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
      if (cellNp[et] == 0.) // masque
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

  /*-----------CONSTRUCTION ARRAY DE SORTIE ------------------*/
  PyObject* tpl;
  if (resi == 1) 
    tpl = K_ARRAY::buildArray(*fi, varStringi, nii, nji, nki);
  else 
    tpl = K_ARRAY::buildArray(*fi, varStringi, *cni, -1, eltTypei, false);
  RELEASESHAREDB(resi, indicator, fi, cni);
  return tpl;
}
