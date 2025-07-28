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

// refine elements of a surface triangular mesh

# include "post.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Raffine un maillage surfacique TRI en fonction d'un critere de raffinement.
   Retourne un maillage surfacique TRI modifie */
//=============================================================================
PyObject* K_POST::refine(PyObject* self, PyObject* args)
{
  // surf: maillage a deraffiner (x,y,z+sol)
  // indic: indicateur de deraffinement: vaut 0 ou 1
  PyObject* surf; PyObject* aindic;
  if (!PyArg_ParseTuple(args, "OO", &surf, &aindic)) return NULL;  

  /*-----------------------------------------------*/
  /* Extraction des donnees du maillage surfacique */ 
  /*-----------------------------------------------*/
  char* varString0; char* eltType0;
  FldArrayF* f; FldArrayI* cn;
  E_Int nil, njl, nkl;
  E_Int res = 
    K_ARRAY::getFromArray3(surf, varString0, f, nil, njl, nkl, cn, eltType0);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "refine: input array is invalid.");
    return NULL;
  }
  if (res != 2 || strcmp(eltType0, "TRI") != 0)
  {
    RELEASESHAREDB(res, surf, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "refine: input array must be TRI.");
    return NULL;
  }
  
  // Check size of array
  E_Int posxu = K_ARRAY::isCoordinateXPresent(varString0);
  E_Int posyu = K_ARRAY::isCoordinateYPresent(varString0);
  E_Int poszu = K_ARRAY::isCoordinateZPresent(varString0);

  if (posxu == -1 || posyu == -1 || poszu == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "refine: coordinates not found in array.");
    RELEASESHAREDU(surf, f, cn);
    return NULL;
  }
  posxu++; posyu++; poszu++;
  
  /*-------------------------------------------*/
  /* Extraction de l'indicateur de raffinement */ 
  /*-------------------------------------------*/
  char* varString1; char* eltType1;
  FldArrayF* findic; FldArrayI* cn1;
  E_Int ni1, nj1, nk1;
  E_Int res1 = K_ARRAY::getFromArray3(aindic, varString1, findic, 
                                      ni1, nj1, nk1, cn1, eltType1);
  E_Int ne = cn->getSize();
  if (res1 == 1)
  {
    if (ni1*nj1*nk1 != ne)
    {
      RELEASESHAREDU(surf, f, cn);
      RELEASESHAREDS(aindic, findic);
      PyErr_SetString(
        PyExc_TypeError,
        "refine: dimension of refinement indicator array must be equal to the number of elements.");
      return NULL;
    }
  }
  else if (res1 == 2)
  {
    if ( findic->getSize() != ne )
    {
      RELEASESHAREDU(surf, f, cn);
      RELEASESHAREDU(aindic, findic, cn1);
      PyErr_SetString(
        PyExc_TypeError,
        "refine: dimension of refinement indicator array must be equal to the number of elements.");
      return NULL;
    }
  }
  else
  {
    RELEASESHAREDU(surf, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "refine: refinement indicator array is invalid.");
    return NULL;
  }
  
  // Passage a un indic en noeuds
  E_Int npts = f->getSize();
  E_Int n1, n2, n3;
  FldArrayIS indic(npts); indic.setAllValuesAtNull();
  E_Int* cnp1 = cn->begin(1);
  E_Int* cnp2 = cn->begin(2);
  E_Int* cnp3 = cn->begin(3);
  for (E_Int i = 0; i < ne; i++)
  {
    if ((*findic)[i] > 0)
    {
      n1 = cnp1[i]-1; n2 = cnp2[i]-1; n3 = cnp3[i]-1;
      indic[n1] = 1; indic[n2] = 1; indic[n3] = 1;
    }
  }
  refineElements(*f, *cn, indic);

  E_Float tolc = 1.e-12; 
  K_CONNECT::cleanConnectivity(posxu, posyu, poszu, tolc, 
                               "TRI", *f, *cn);

  PyObject* t = K_ARRAY::buildArray(*f, varString0, *cn, -1, "TRI");

  RELEASESHAREDB(res1, aindic, findic, cn1);
  RELEASESHAREDU(surf, f, cn);
      
  return t;
}

//=============================================================================
void K_POST::refineElements(FldArrayF& f, FldArrayI& cn, FldArrayIS& indic)
{
  E_Int i, v, n1, n2, n3;
  short i1, i2, i3;
  // Nombre de champs
  E_Int nf = f.getNfld();
  // Nombre d'elements sur le maillage d'origine
  E_Int ne = cn.getSize();
  // Nombre de noeuds sur le maillage d'origine
  E_Int np = f.getSize();

  // Nouvelle connectivite
  FldArrayI co(4*ne, cn.getNfld());

  // Nouveaux coord des points
  FldArrayF fo(6*ne, nf);
  for (v = 1; v <= nf; v++)
  {
    E_Float* fop = fo.begin(v);
    E_Float* fp =f.begin(v);
    for (i = 0; i < np; i++) fop[i] = fp[i];
  }

  // no de l'element courant dans le nouveau maillage
  E_Int nelt = 0;
  // no du point courant dans le nouveau maillage
  E_Int npt = np;
  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2); E_Int* cn3 = cn.begin(3);
  E_Int* co1 = co.begin(1); E_Int* co2 = co.begin(2); E_Int* co3 = co.begin(3);
  for (i = 0; i < ne; i++)
  {
    n1 = cn1[i]-1; n2 = cn2[i]-1; n3 = cn3[i]-1;
    i1 = indic[n1]; i2 = indic[n2]; i3 = indic[n3];
    if (i1 > 0 && i2 > 0 && i3 > 0)
    {
      // Ajoute 3 pts et 4 elements
      for (v = 1; v <= nf; v++)
      {
        fo(npt, v) = 0.5*(f(n1,v) + f(n2,v));
        fo(npt+1, v) = 0.5*(f(n1,v) + f(n3,v));
        fo(npt+2, v) = 0.5*(f(n2,v) + f(n3,v));
      }
      co1[nelt] = n1+1;
      co2[nelt] = npt+1;
      co3[nelt] = npt+2;
      co1[nelt+1] = npt+2;
      co2[nelt+1] = npt+3;
      co3[nelt+1] = n3+1;
      co1[nelt+2] = npt+1;
      co2[nelt+2] = npt+3;
      co3[nelt+2] = npt+2;
      co1[nelt+3] = npt+1;
      co2[nelt+3] = n2+1;
      co3[nelt+3] = npt+3;
      npt += 3;
      nelt += 4;
    }
    else if (i1 > 0 && i2 > 0)
    {
      // Ajoute 1 pt et 2 elements
      for (v = 1; v <= nf; v++) fo(npt, v) = 0.5*(f(n1,v) + f(n2,v));
      co1[nelt] = n1+1;
      co2[nelt] = npt+1;
      co3[nelt] = n3+1;
      co1[nelt+1] = npt+1;
      co2[nelt+1] = n2+1;
      co3[nelt+1] = n3+1;
      npt += 1;
      nelt += 2;
    }
    else if (i1 > 0 && i3 > 0)
    {
      // Ajoute 1 pt et 2 elements
      for (v = 1; v <= nf; v++) fo(npt, v) = 0.5*(f(n1,v) + f(n3,v));
      co1[nelt] = n1+1;
      co2[nelt] = n2+1;
      co3[nelt] = npt+1;
      co1[nelt+1] = npt+1;
      co2[nelt+1] = n2+1;
      co3[nelt+1] = n3+1;
      npt = npt + 1;
      nelt = nelt + 2;
    }
    else if (i2 > 0 && i3 > 0)
    {
      // Ajoute 1 pt et 2 elements
      for (v = 1; v <= nf; v++) fo(npt, v) = 0.5*(f(n2,v) + f(n3,v));
      co1[nelt] = n1+1;
      co2[nelt] = n2+1;
      co3[nelt] = npt+1;
      co1[nelt+1] = n1+1;
      co2[nelt+1] = npt+1;
      co3[nelt+1] = n3+1;
      npt += 1;
      nelt += 2;
    }
    else
    {
      // Ajoute aucun point et 1 element
      co1[nelt] = n1+1;
      co2[nelt] = n2+1;
      co3[nelt] = n3+1;
      nelt += 1;
    }
  }

  // Reallocation
  co.reAllocMat(nelt, 3);
  fo.reAllocMat(npt, nf);
  
  cn = co;
  f = fo;
}
