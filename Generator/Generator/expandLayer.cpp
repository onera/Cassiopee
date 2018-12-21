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
/* Modify the indicator to expand the layer of level l */
//=============================================================================
PyObject* K_GENERATOR::modifyIndicToExpandLayer(PyObject* self, PyObject* args)
{
  PyObject *octree, *indicator; 
  E_Int level; E_Int corners;
#ifdef E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOll", &octree, &indicator, &level, &corners)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "OOii", &octree, &indicator, &level, &corners)) return NULL;
#endif
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
                    "expandLayer: indic array must be structured.");
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

  // detection des elements de niveau l 
  E_Float eps = 1.e-10; E_Float dhl = pow(2.,level) * dhmin;
  
  E_Int etv;
  E_Float dhleps = dhl + eps;
  for (E_Int et = 0; et < nelts; et++)
  {
    if (K_FUNC::fEqualZero(dhtp[et]-dhl,eps) == true) 
    {
      vector<E_Int>& voisins = cEEN[et];
      for (size_t nov = 0; nov < voisins.size(); nov++)
      {
        etv = voisins[nov];
        if (dhtp[etv] > dhleps) indict[etv] = 1.;
      }
    }
  }

  /*-----------CONSTRUCTION ARRAY DE SORTIE ------------------*/
  PyObject* tpl;
  if (resi == 1) 
    tpl = K_ARRAY::buildArray(*fi, varStringi, nii, nji, nki);
  else 
    tpl = K_ARRAY::buildArray(*fi, varStringi, *cni, -1, eltTypei, 
                              false);
  RELEASESHAREDB(resi, indicator, fi, cni);
  return tpl;
}
