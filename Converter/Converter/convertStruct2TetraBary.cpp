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

// Convert structured array to tetra "star" array

# include <stdlib.h>
# include "converter.h"
# include "kcore.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Convert structured array to a tetraedrical mesh with addition of points 
   (barycenter of elements and faces) 
   The method deals only with fields located at nodes */
// ============================================================================
PyObject* K_CONVERTER::convertStruct2TetraBary(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, ni, nj, nk, cn, eltType, true);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2TetraBary: array is invalid.");
    return NULL;
  }
  if (res == 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2TetraBary: array must be structured.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if ( posx == 0 || posy == 0 || posz == 0 )
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2TetraBary: coordinates not found.");
    RELEASESHAREDS(array, f); return NULL;
  }

  // dimensions du maillage structure
  E_Int dim0 = 3;
  if (ni == 1)
  {
    if (nj == 1 || nk == 1 ) dim0 = 1;
    else {dim0 = 2; ni = nj; nj = nk; nk = 1;}
  }
  else if (nj == 1)
  {
    if (ni == 1 || nk == 1 ) dim0 = 1;
    else {dim0 = 2; nj = nk; nk = 1;}
  }
  else if (nk == 1)
  {
    if (ni == 1 || nj == 1 ) dim0 = 1;
    else dim0 = 2;
  }
  
  // variables locales
  E_Int ni1 = E_max(1, E_Int(ni)-1);
  E_Int nj1 = E_max(1, E_Int(nj)-1);
  E_Int nk1 = E_max(1, E_Int(nk)-1);
  E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees
  E_Int ninj = ni*nj;
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();
  E_Int nelts = ncells;      // nb d elements pour la nouvelle connectivite TETRA
  E_Int nconnect = 2;        // nb de points par element
  E_Int nptsNewField = npts; // taille des nouveaux champs 
  if ( dim0 == 1)
  {
    nelts = ncells*2;
    nconnect = 2;
    nptsNewField = npts+ncells;
  }
  else if (dim0 == 2) 
  {
    nelts = ncells*4;
    nconnect = 3;
    nptsNewField = npts + ncells; // ajout des barycentres des cellules
  }
  else if (dim0 == 3)
  {
    nelts = ncells*24;
    nconnect = 4;
    nptsNewField = npts + 7*ncells; // ajout des barycentres des cellules et de chaque face
  }  
  // Nouvelle connectivite
  FldArrayI cnnew(nelts, nconnect);
  // pointeurs sur la nouvelle connectivite TETRA
  E_Int* cn1 = cnnew.begin(1);
  E_Int* cn2 = cnnew.begin(2);
  E_Int* cn3 = NULL;
  E_Int* cn4 = NULL;
  if (dim0 != 1) cn3 = cnnew.begin(3);
  if (dim0 == 3) cn4 = cnnew.begin(4);
  // nouveau champ
  FldArrayF fnew(nptsNewField, nfld);
  // pointeurs sur le nouveau champ
  vector<E_Float*> fnewp(nfld);
  for (E_Int p = 0; p < nfld; p++) {fnewp[p] = fnew.begin(p+1);}
  // pointeurs sur l ancien champ
  vector<E_Float*> fp(nfld);
  for (E_Int p = 0; p < nfld; p++) {fp[p] = f->begin(p+1);}
  E_Float* fpp; E_Float* fnpp;

  // Type d'element de la nouvelle connectivite
  char newEltType[256];

  // Construction du maillage non structure
  switch (dim0)
  {
    case 1:
    {
      strcpy(newEltType, "BAR");
      E_Int indp1, indp2;
      E_Int noind = 0; E_Int et = 0;
      E_Int N;
      if (nk1 == 1 && nj1 == 1) N = ni1;
      else if (ni1 == 1 && nj1 == 1) N = nk1;
      else N = nj1;
      
      for (E_Int i = 0; i < N; i++)
      {
        indp1 = i; indp2 = i+1;//point a gauche et droite de la cellule
        cn1[et] = noind+1;//point a gauche
        cn2[et] = noind+2;//barycentre
        cn1[et+1] = noind+2;//barycentre
        cn2[et+1] = noind+3;//point a droite
        
        for (E_Int p = 0; p < nfld; p++) 
        {
          fpp = fp[p]; fnpp = fnewp[p];
          fnpp[noind] = fpp[indp1];
          fnpp[noind+1] = (fpp[indp1]+fpp[indp2])*0.5;
        }
        noind+= 2; et+= 2;
      }
      //dernier point
      for (E_Int p = 0; p < nfld; p++) 
      {
        fpp = fp[p]; fnpp = fnewp[p];
        fnpp[noind] = fpp[N];
      }      
    }
    break;
    case 2:
    {
      strcpy(newEltType, "TRI");
      // tableau des indices des points la cellule
      E_Int ind1, ind2, ind3, ind4;
      E_Int indelt = 0; // indice d element de la nouvelle connectivite TETRA
      E_Int indpt = npts; // indice des nouveaux points
      E_Int indNewElt; // indice du nouvel elmt
      
      for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          ind1 = i + j*ni;      //A(  i,  j,1)
          ind2 = ind1 + 1;      //B(i+1,  j,1)
          ind3 = ind1 + ni;     //C(  i,j+1,1)
          ind4 = ind3 + 1;      //D(i+1,j+1,1)
          // champs pour le barycentre de la cellule
          for (E_Int p = 0; p < nfld; p++) 
          {
            fpp = fp[p]; fnpp = fnewp[p];
            fnpp[indpt] = (fpp[ind1]+fpp[ind2]+fpp[ind3]+fpp[ind4])*0.25;
          }
          indpt++; indNewElt = indpt;
          // connectivite nouveaux elements
          cn1[indelt] = ind1+1;
          cn2[indelt] = ind2+1;
          cn3[indelt] = indNewElt; indelt++;
          
          cn1[indelt] = ind2+1;
          cn2[indelt] = ind4+1;
          cn3[indelt] = indNewElt; indelt++;
            
          cn1[indelt] = ind4+1;
          cn2[indelt] = ind3+1;
          cn3[indelt] = indNewElt; indelt++;
          
          cn1[indelt] = ind3+1;
          cn2[indelt] = ind1+1;
          cn3[indelt] = indNewElt; indelt++;
          // champs du nouvel element tetraedrique pour les neouds [indf1,indf2,indf3,indf4]
          for (E_Int p = 0; p < nfld; p++) 
          {
            fpp = fp[p]; fnpp = fnewp[p];
            fnpp[ind1] = fpp[ind1]; fnpp[ind2] = fpp[ind2]; 
            fnpp[ind3] = fpp[ind3]; fnpp[ind4] = fpp[ind4];
          }
        }
    }
    break;
    case 3:
    { 
      strcpy(newEltType, "TETRA");
      // tableau des indices des points la cellule
      E_Int ind[8];
      // indices des points la face
      E_Int indf1, indf2, indf3, indf4;
      // indir: indirection pour les numerotations des noeuds pour chaque face
      //  ordre des faces dans indir : ABCD, EFGH, ADHE, BCGF, ABFE, DCGH
      E_Int indir[24] = {0,1,2,3,4,5,6,7,0,3,7,4,1,2,6,5,0,1,5,4,3,2,6,7};
      E_Int indelt = 0; // indice d element de la nouvelle connectivite TETRA
      E_Int indpt = npts; // indice des nouveaux points
      E_Int indNewElt, indNewFace; // indices du nouvel elmt et de la nouvelle face
  
      // Parcours des cellules
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind[0] = i + j*ni + k*ninj;        //A(  i,  j,k)
            ind[1] = ind[0] + 1;               //B(i+1,  j,k)
            ind[2] = ind[1] + ni;              //C(i+1,j+1,k)
            ind[3] = ind[2] - 1;               //D(  i,j+1,k)
            ind[4] = ind[0] + ninj;            //E(  i,  j,k+1)
            ind[5] = ind[1] + ninj;            //F(i+1,  j,k+1)
            ind[6] = ind[2] + ninj;            //G(i+1,j+1,k+1)
            ind[7] = ind[3] + ninj;            //H(  i,j+1,k+1) 

            // champs pour le barycentre de la cellule
            for (E_Int p = 0; p < nfld; p++) 
            {
              fnpp = fnewp[p]; fpp = fp[p];
              fnpp[indpt] = (fpp[ind[0]]+fpp[ind[1]]+fpp[ind[2]]+fpp[ind[3]]+fpp[ind[4]]+fpp[ind[5]]+fpp[ind[6]]+fpp[ind[7]])*0.125;
            }
            indpt ++; indNewElt = indpt;
            for (E_Int fa = 0; fa < 6; fa++)
            {
              // champs pour le barycentre de la face
              indf1 = ind[indir[4*fa]]; 
              indf2 = ind[indir[4*fa+1]]; 
              indf3 = ind[indir[4*fa+2]]; 
              indf4 = ind[indir[4*fa+3]]; 
              for (E_Int p = 0; p < nfld; p++) 
              {
                fnpp = fnewp[p]; fpp = fp[p];
                fnpp[indpt] = (fpp[indf1]+fpp[indf2]+fpp[indf3]+fpp[indf4])*0.25;
              }
              indpt ++; indNewFace = indpt;
              // connectivite nouveaux elements lies a la face
              cn1[indelt] = indf1+1;
              cn2[indelt] = indf2+1;
              cn3[indelt] = indNewElt;
              cn4[indelt] = indNewFace; indelt++;
          
              cn1[indelt] = indf2+1;
              cn2[indelt] = indf3+1;
              cn3[indelt] = indNewElt;
              cn4[indelt] = indNewFace; indelt++;
          
              cn1[indelt] = indf3+1;
              cn2[indelt] = indf4+1;
              cn3[indelt] = indNewElt;
              cn4[indelt] = indNewFace; indelt++;
          
              cn1[indelt] = indf4+1;
              cn2[indelt] = indf1+1;
              cn3[indelt] = indNewElt;
              cn4[indelt] = indNewFace; indelt++;
             // champs du nouvel element tetraedrique pour les neouds [indf1,indf2,indf3,indf4]
              for (E_Int p = 0; p < nfld; p++) 
              {
                fnpp = fnewp[p]; fpp = fp[p];
                fnpp[indf1] = fpp[indf1]; fnpp[indf2] = fpp[indf2]; 
                fnpp[indf3] = fpp[indf3]; fnpp[indf4] = fpp[indf4];
              }
            }
          }
    }
    break;
  }

  // Nettoyage de la connectivite TETRA
  K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-12, newEltType, fnew, cnnew);

  // Objet python retourne
  PyObject* tpl = K_ARRAY::buildArray(fnew, varString, cnnew, -1, newEltType);

  // Liberation de la memoire
  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}
// ============================================================================
/* Convert structured array to a tetraedrical mesh with addition of points 
   (barycenter of elements and faces) 
   The method deals only with fields located at nodes */
// ============================================================================
PyObject* K_CONVERTER::convertStruct2TetraBaryBoth(PyObject* self, PyObject* args)
{
  PyObject *array, *arrayc;
  if (!PyArg_ParseTuple(args, "OO", &array, &arrayc)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, ni, nj, nk, cn, eltType, true);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2TetraBary: array is invalid.");
    return NULL;
  }
  if (res == 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2TetraBary: array must be structured.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }
  E_Int nic, njc, nkc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  E_Int resc = K_ARRAY::getFromArray(arrayc, varStringc, 
                                     fc, nic, njc, nkc, cnc, eltTypec, true);
  if (resc != 1 && resc != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2TetraBary: arrayc is invalid.");
    RELEASESHAREDS(array, f);
    return NULL;
  }
  if (resc == 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2TetraBary: arrayc must be structured.");
    RELEASESHAREDS(array, f);
    RELEASESHAREDU(arrayc, fc, cnc);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if ( posx == 0 || posy == 0 || posz == 0 )
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2TetraBary: coordinates not found.");
    RELEASESHAREDS(array, f);  RELEASESHAREDS(arrayc, fc); return NULL;
  }

  // dimensions du maillage structure
  E_Int dim0 = 3;
  if (ni == 1)
  {
    if (nj == 1 || nk == 1 ) dim0 = 1;
    else {dim0 = 2; ni = nj; nj = nk; nk = 1;}
  }
  else if (nj == 1)
  {
    if (ni == 1 || nk == 1 ) dim0 = 1;
    else {dim0 = 2; nj = nk; nk = 1;}
  }
  else if (nk == 1)
  {
    if (ni == 1 || nj == 1 ) dim0 = 1;
    else dim0 = 2;
  }
  
  // variables locales
  E_Int ni1 = E_max(1, E_Int(ni)-1);
  E_Int nj1 = E_max(1, E_Int(nj)-1);
  E_Int nk1 = E_max(1, E_Int(nk)-1);
  E_Int ni1nj1 = ni1*nj1;
  E_Int ncells = ni1nj1*nk1; // nb de cellules structurees
  E_Int ninj = ni*nj;
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();
  E_Int nelts = ncells;      // nb d elements pour la nouvelle connectivite TETRA
  E_Int nconnect = 2;        // nb de points par element
  E_Int nptsNewField = npts; // taille des nouveaux champs 
  if ( dim0 == 1)
  {
    nelts = ncells*2;
    nconnect = 2;
    nptsNewField = npts+ncells;
  }
  else if (dim0 == 2) 
  {
    nelts = ncells*4;
    nconnect = 3;
    nptsNewField = npts + ncells; // ajout des barycentres des cellules
  }
  else if (dim0 == 3)
  {
    nelts = ncells*24;
    nconnect = 4;
    nptsNewField = npts + 7*ncells; // ajout des barycentres des cellules et de chaque face
  }  
  // Nouvelle connectivite
  FldArrayI cnnew(nelts, nconnect);
  E_Int nfldc = fc->getNfld();
  FldArrayF fcnew(nelts, nfldc);
  
  // pointeurs sur la nouvelle connectivite TETRA
  E_Int* cn1 = cnnew.begin(1);
  E_Int* cn2 = cnnew.begin(2);
  E_Int* cn3 = NULL;
  E_Int* cn4 = NULL;
  if (dim0 != 1) cn3 = cnnew.begin(3);
  if (dim0 == 3) cn4 = cnnew.begin(4);
  // nouveau champ
  FldArrayF fnew(nptsNewField, nfld);
  // pointeurs sur le nouveau champ
  vector<E_Float*> fnewp(nfld);
  for (E_Int p = 0; p < nfld; p++) {fnewp[p] = fnew.begin(p+1);}

  vector<E_Float*> fcnewp(nfldc);
  for (E_Int p = 0; p < nfldc; p++) {fcnewp[p] = fcnew.begin(p+1);}

  // pointeurs sur l ancien champ
  vector<E_Float*> fp(nfld);
  for (E_Int p = 0; p < nfld; p++) {fp[p] = f->begin(p+1);}
  vector<E_Float*> fcp(nfldc);
  for (E_Int p = 0; p < nfldc; p++) {fcp[p] = fc->begin(p+1);}

  E_Float* fpp; E_Float* fnpp;

  // Type d'element de la nouvelle connectivite
  char newEltType[256];

  // Construction du maillage non structure
  switch (dim0)
  {
    case 1:
    {
      strcpy(newEltType, "BAR");
      E_Int indp1, indp2;
      E_Int noind = 0; E_Int et = 0;
      E_Int N;
      if (nk1 == 1 && nj1 == 1) N = ni1;
      else if (ni1 == 1 && nj1 == 1) N = nk1;
      else N = nj1;
      
      for (E_Int i = 0; i < N; i++)
      {
        indp1 = i; indp2 = i+1;//point a gauche et droite de la cellule
        cn1[et] = noind+1;//point a gauche
        cn2[et] = noind+2;//barycentre
        cn1[et+1] = noind+2;//barycentre
        cn2[et+1] = noind+3;//point a droite
        for (E_Int p = 0; p < nfldc; p++) 
        {
          fpp = fcp[p]; fnpp = fcnewp[p];
          fnpp[et] = fpp[i];
          fnpp[et+1] = fpp[i];
        }
        for (E_Int p = 0; p < nfld; p++) 
        {
          fpp = fp[p]; fnpp = fnewp[p];
          fnpp[noind] = fpp[indp1];
          fnpp[noind+1] = (fpp[indp1]+fpp[indp2])*0.5;
        }
        noind+= 2; et+= 2;
      }
      //dernier point
      for (E_Int p = 0; p < nfld; p++) 
      {
        fpp = fp[p]; fnpp = fnewp[p];
        fnpp[noind] = fpp[N];
      }      
    }
    break;
    case 2:
    {
      strcpy(newEltType, "TRI");
      // tableau des indices des points la cellule
      E_Int ind1, ind2, ind3, ind4, indcell;
      E_Int indelt = 0; // indice d element de la nouvelle connectivite TETRA
      E_Int indpt = npts; // indice des nouveaux points
      E_Int indNewElt; // indice du nouvel elmt
      
      for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          ind1 = i + j*ni;      //A(  i,  j,1)
          ind2 = ind1 + 1;      //B(i+1,  j,1)
          ind3 = ind1 + ni;     //C(  i,j+1,1)
          ind4 = ind3 + 1;      //D(i+1,j+1,1)
          indcell = i + j*ni1;
          // champs pour le barycentre de la cellule
          for (E_Int p = 0; p < nfld; p++) 
          {
            fpp = fp[p]; fnpp = fnewp[p];
            fnpp[indpt] = (fpp[ind1]+fpp[ind2]+fpp[ind3]+fpp[ind4])*0.25;
          }
          indpt++; indNewElt = indpt;
          for (E_Int p = 0; p < nfldc; p++) 
          {
            fpp = fcp[p]; fnpp = fcnewp[p];        
            fnpp[indelt]   = fpp[indcell];
            fnpp[indelt+1] = fpp[indcell];
            fnpp[indelt+2] = fpp[indcell];
            fnpp[indelt+3] = fpp[indcell];
          }
          // connectivite nouveaux elements
          cn1[indelt] = ind1+1;
          cn2[indelt] = ind2+1;
          cn3[indelt] = indNewElt; indelt++;
          
          cn1[indelt] = ind2+1;
          cn2[indelt] = ind4+1;
          cn3[indelt] = indNewElt; indelt++;
            
          cn1[indelt] = ind4+1;
          cn2[indelt] = ind3+1;
          cn3[indelt] = indNewElt; indelt++;
          
          cn1[indelt] = ind3+1;
          cn2[indelt] = ind1+1;
          cn3[indelt] = indNewElt; indelt++;
          // champs du nouvel element tetraedrique pour les neouds [indf1,indf2,indf3,indf4]
          for (E_Int p = 0; p < nfld; p++) 
          {
            fpp = fp[p]; fnpp = fnewp[p];
            fnpp[ind1] = fpp[ind1]; fnpp[ind2] = fpp[ind2]; 
            fnpp[ind3] = fpp[ind3]; fnpp[ind4] = fpp[ind4];
          }
        }
    }
    break;
    case 3:
    { 
      strcpy(newEltType, "TETRA");
      // tableau des indices des points la cellule
      E_Int ind[8];
      // indices des points la face
      E_Int indf1, indf2, indf3, indf4, indcell;
      // indir: indirection pour les numerotations des noeuds pour chaque face
      //  ordre des faces dans indir : ABCD, EFGH, ADHE, BCGF, ABFE, DCGH
      E_Int indir[24] = {0,1,2,3,4,5,6,7,0,3,7,4,1,2,6,5,0,1,5,4,3,2,6,7};
      E_Int indelt = 0; // indice d element de la nouvelle connectivite TETRA
      E_Int indpt = npts; // indice des nouveaux points
      E_Int indNewElt, indNewFace; // indices du nouvel elmt et de la nouvelle face
  
      // Parcours des cellules
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind[0] = i + j*ni + k*ninj;        //A(  i,  j,k)
            ind[1] = ind[0] + 1;               //B(i+1,  j,k)
            ind[2] = ind[1] + ni;              //C(i+1,j+1,k)
            ind[3] = ind[2] - 1;               //D(  i,j+1,k)
            ind[4] = ind[0] + ninj;            //E(  i,  j,k+1)
            ind[5] = ind[1] + ninj;            //F(i+1,  j,k+1)
            ind[6] = ind[2] + ninj;            //G(i+1,j+1,k+1)
            ind[7] = ind[3] + ninj;            //H(  i,j+1,k+1) 
            indcell = i +j*ni1 + k*ni1nj1;

            // champs pour le barycentre de la cellule
            for (E_Int p = 0; p < nfld; p++) 
            {
              fnpp = fnewp[p]; fpp = fp[p];
              fnpp[indpt] = (fpp[ind[0]]+fpp[ind[1]]+fpp[ind[2]]+fpp[ind[3]]+fpp[ind[4]]+fpp[ind[5]]+fpp[ind[6]]+fpp[ind[7]])*0.125;
            }
            indpt ++; indNewElt = indpt;
            for (E_Int fa = 0; fa < 6; fa++)
            {
              // champs pour le barycentre de la face
              indf1 = ind[indir[4*fa]]; 
              indf2 = ind[indir[4*fa+1]]; 
              indf3 = ind[indir[4*fa+2]]; 
              indf4 = ind[indir[4*fa+3]]; 
              for (E_Int p = 0; p < nfld; p++) 
              {
                fnpp = fnewp[p]; fpp = fp[p];
                fnpp[indpt] = (fpp[indf1]+fpp[indf2]+fpp[indf3]+fpp[indf4])*0.25;
              }
              indpt ++; indNewFace = indpt;
              
              //champs en centres 
              for (E_Int p = 0; p < nfldc; p++) 
              {
                fpp = fcp[p]; fnpp = fcnewp[p];        
                fnpp[indelt]   = fpp[indcell];
                fnpp[indelt+1] = fpp[indcell];
                fnpp[indelt+2] = fpp[indcell];
                fnpp[indelt+3] = fpp[indcell];                
              }

              // connectivite nouveaux elements lies a la face
              cn1[indelt] = indf1+1;
              cn2[indelt] = indf2+1;
              cn3[indelt] = indNewElt;
              cn4[indelt] = indNewFace; indelt++;
          
              cn1[indelt] = indf2+1;
              cn2[indelt] = indf3+1;
              cn3[indelt] = indNewElt;
              cn4[indelt] = indNewFace; indelt++;
          
              cn1[indelt] = indf3+1;
              cn2[indelt] = indf4+1;
              cn3[indelt] = indNewElt;
              cn4[indelt] = indNewFace; indelt++;
          
              cn1[indelt] = indf4+1;
              cn2[indelt] = indf1+1;
              cn3[indelt] = indNewElt;
              cn4[indelt] = indNewFace; indelt++;
             // champs du nouvel element tetraedrique pour les noeuds [indf1,indf2,indf3,indf4]
              for (E_Int p = 0; p < nfld; p++) 
              {
                fnpp = fnewp[p]; fpp = fp[p];
                fnpp[indf1] = fpp[indf1]; fnpp[indf2] = fpp[indf2]; 
                fnpp[indf3] = fpp[indf3]; fnpp[indf4] = fpp[indf4];
              }
            }
          }
    }
    break;
  }

  // Nettoyage de la connectivite TETRA
  K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-12, newEltType, fnew, cnnew);

  // Objet python retourne
  PyObject* l = PyList_New(0);

  PyObject* tpl1 = K_ARRAY::buildArray(fnew, varString, cnnew, -1, newEltType);
  PyList_Append(l, tpl1); Py_DECREF(tpl1);
  PyObject* tpl2 = K_ARRAY::buildArray(fcnew, varStringc, cnnew, -1, newEltType);
  PyList_Append(l, tpl2); Py_DECREF(tpl2);
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDB(resc, arrayc, fc, cnc);

  return l;
}
