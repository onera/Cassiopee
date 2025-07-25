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
# include <vector>
# include <stdlib.h>
# include "converter.h"
# include "kcore.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Convert non-polyhedral 3D array to a tetraedrical mesh with addition of 
   points (barycenter of elements and faces) 
   The method deals only with fields located at nodes */
// ============================================================================
PyObject* K_CONVERTER::convertArray2TetraBary(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString, 
                               f, ni, nj, nk, cn, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2TetraBary: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2TetraBary: array must be unstructured.");
    return NULL;
  }

  if (strcmp(eltType, "NGON") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "convertArray2TetraBary: array must not be NGon.");
    return NULL;
  }

  if ((strcmp(eltType, "NODE") == 0)||
      // (strcmp(eltType, "BAR") == 0) ||
      (strcmp(eltType, "TRI") == 0) ||
      (strcmp(eltType, "TETRA") == 0))
  {
    PyObject* tpl = K_ARRAY::buildArray(*f, varString, *cn, -1, eltType);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if (posx == 0 || posy == 0 || posz == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2TetraBary: coordinates not found.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }

  // nb d'elements
  E_Int nelts = cn->getSize(); 

  // specific datas depending on element type
  E_Int newElts=0; // nb nouveaux elements tetra
  E_Int nfaces=0; // nb de faces par element
  E_Int npoints=0; // nb points par element
  E_Int dim = 3; 
  if (strcmp(eltType, "QUAD") == 0)
  {
    newElts = 4*nelts; nfaces = 0; npoints = 4; // nfaces = 0 au lieu de 4 car on ajoute pas de points sur les faces
    dim = 2;
  }
  else if (strcmp(eltType, "HEXA") == 0)
  {
    newElts = 24*nelts; nfaces = 6; npoints = 8;
    dim = 3;
  }
  else if (strcmp(eltType, "PENTA") == 0)
  {
    newElts = 18*nelts; nfaces = 5; npoints = 6;
    dim = 3;  
  }
  else if (strcmp(eltType, "PYRA") == 0)
  {
    newElts = 16*nelts; nfaces = 5; npoints = 5;
    dim = 3;
  }
  else if (strcmp(eltType, "BAR") == 0)
  {
    newElts = 2*nelts; nfaces = 0; npoints = 2;
    dim = 1;
  }

  // Calcul de la taille du tableau des nouveaux champs
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();
  E_Int sizeFields = npts;

  for (E_Int elt = 0; elt < nelts; elt++)
  {
    for (E_Int fa = 0; fa < nfaces; fa++) {sizeFields += 1;}
    sizeFields += 1;
  }

  // nouveau champ
  FldArrayF fnew(sizeFields, nfld);

  // pointeurs sur le nouveau champ
  vector<E_Float*> fnewp(nfld);
  for (E_Int p=0; p < nfld; p++) {fnewp[p] = fnew.begin(p+1);}
  // pointeurs sur l ancien champ
  vector<E_Float*> fp(nfld);
  for (E_Int p=0; p < nfld; p++) {fp[p] = f->begin(p+1);}  

  // Nouvelle connectivite
  E_Int nconnect = 2; // nb de points par element
  if (dim == 2) nconnect = 3; // elmt triangle
  else if (dim == 1) nconnect = 2; // elt BAR
  else nconnect = 4;          // elmt tetraedre
  FldArrayI newcn(newElts,nconnect);  
  // Pointeurs sur la nouvelle connectivite
  vector<E_Int*> newcnp(nconnect);
  for (E_Int p=0; p < nconnect; p++) {newcnp[p] = newcn.begin(p+1);}
  // Pointeurs sur l ancienne connectivite
  vector<E_Int*> cnp(npoints);
  for (E_Int p = 0; p < npoints; p++) {cnp[p] = cn->begin(p+1);}

  // Type d'element de la nouvelle connectivite
  char newEltType[256];
  if (dim == 1) // BAR
  {
    strcpy(newEltType,"BAR");
    E_Int indp1, indp2; 
    E_Int indpt = npts; //indices des nouveaux points
    E_Int indelt = 0; // indice d elt de la nouvelle connectivite

    for (E_Int et = 0; et < nelts; et++)
    {
      indp1 = cnp[0][et]-1;
      indp2 = cnp[1][et]-1;
    
      // Barycentre de l'element
      for (E_Int p = 0; p < nfld; p++) 
        fnewp[p][indpt] = 0.5*(fp[p][indp1]+fp[p][indp2]);      
      indpt++; 

      // construction des nouveaux elements
      // connectivite
      newcnp[0][indelt] = indp1+1;
      newcnp[1][indelt] = indpt;
      indelt++;

      newcnp[0][indelt] = indpt;
      newcnp[1][indelt] = indp2+1;
      indelt++;
      //coordonnees des 2 nouveaux elements
      for (E_Int p = 0; p < nfld; p++) 
      {
        E_Float* fnewpp = fnewp[p];
        E_Float* fpp = fp[p]; 
        fnewpp[indp1] = fpp[indp1];
        fnewpp[indp2] = fpp[indp2];       
      }
    }
  }
  else if (dim == 2) // Elmt QUAD
  {
    strcpy(newEltType, "TRI");
    
    // Calcul de la nouvelle connectivite et des nouvelles coordonnees
    E_Int ind1, ind2, ind3, ind4; // indices des points la cellule
    E_Int indelt = 0; // indice d element de la nouvelle connectivite
    E_Int indpt = npts; // indice des nouveaux points
    E_Int indNewElt; // indices du nouvel elmt et de la nouvelle face
    for (E_Int elt = 0; elt < nelts; elt++)
    {
      ind1 = cnp[0][elt]-1;
      ind2 = cnp[1][elt]-1;
      ind3 = cnp[2][elt]-1;
      ind4 = cnp[3][elt]-1;

      // calcul du barycentre de l'element
      for (E_Int p = 0; p < nfld; p++) 
      {
        fnewp[p][indpt] = K_CONST::ONE_FOURTH*(fp[p][ind1]+fp[p][ind2]+fp[p][ind3]+fp[p][ind4]);
      }
      indNewElt = indpt; indpt++; 

      // construction des nouveaux elements tetraedriques
      // connectivite des nouveaux elements
      newcnp[0][indelt] = ind1+1;
      newcnp[1][indelt] = ind2+1;
      newcnp[2][indelt] = indNewElt+1; indelt++;
      
      newcnp[0][indelt] = ind2+1;
      newcnp[1][indelt] = ind3+1;
      newcnp[2][indelt] = indNewElt+1; indelt++;
      
      newcnp[0][indelt] = ind3+1;
      newcnp[1][indelt] = ind4+1;
      newcnp[2][indelt] = indNewElt+1; indelt++;
      
      newcnp[0][indelt] = ind4+1;
      newcnp[1][indelt] = ind1+1;
      newcnp[2][indelt] = indNewElt+1; indelt++;
      
     // coordonnees du nouvel element tetraedrique
      for (E_Int p = 0; p < nfld; p++) 
      {
        E_Float* fnewpp = fnewp[p];
        E_Float* fpp = fp[p]; 
        fnewpp[ind1] = fpp[ind1];
        fnewpp[ind2] = fpp[ind2];
        fnewpp[ind3] = fpp[ind3];
        fnewpp[ind4] = fpp[ind4];
      }
    }
  }
  else // Elmts 3D
  {
    strcpy(newEltType, "TETRA");

    vector< vector<E_Int> > faces(nfaces); // vecteur contenant pour un element les indices "locaux" des points de chaque face
    buildFaceIndices(eltType, faces);

    // Calcul de la nouvelle connectivite et des nouvelles coordonnees
    vector<E_Float> fbe(nfld); // champs du barycentre de l element
    vector<E_Float> fbf(nfld); // champs du barycentre de la face
    E_Int ind, ind1, ind2; // indices de noeud 
    E_Int indelt = 0; // indice d element de la nouvelle connectivite
    E_Int indpt = npts; // indice des nouveaux points
    E_Int indNewElt, indNewFace; // indices du nouvel elmt et de la nouvelle face
    E_Int nptsInFace; // nb points pour la face donnee
    for (E_Int elt = 0; elt < nelts; elt++)
    {
      // calcul du barycentre de l'element
      for (E_Int p = 0; p < nfld; p++) 
      {
        fbe[p] = 0.;
        for (E_Int n = 0; n < npoints; n++) 
        {
          ind = cnp[n][elt] -1;
          fbe[p] += fp[p][ind];
        }
        fbe[p] = fbe[p]/npoints;
        
        fnewp[p][indpt] = fbe[p];
      }
      indNewElt = indpt; indpt++; 
      for (E_Int fa = 0; fa < nfaces; fa++)
      {
        nptsInFace = faces[fa].size();
        // calcul du barycentre de la face
        for (E_Int p = 0; p < nfld; p++) 
        {
          fbf[p] = 0.;
          for (E_Int n = 0; n < nptsInFace; n++) 
          {
            ind = cnp[faces[fa][n]][elt] -1;
            fbf[p] += fp[p][ind];
          }
          fbf[p] = fbf[p]/nptsInFace;
          fnewp[p][indpt] = fbf[p];
        }
        indNewFace = indpt; indpt++;
        // construction des nouveaux elements tetraedriques
        for (E_Int n = 0; n < nptsInFace; n++)
        {
          ind1 =  cnp[faces[fa][n]][elt] -1;
          ind2 =  cnp[faces[fa][(n+1)%(nptsInFace)]][elt] -1;

          // connectivite du nouvel element
          newcnp[0][indelt] = ind1+1;
          newcnp[1][indelt] = ind2+1;
          newcnp[2][indelt] = indNewElt+1;
          newcnp[3][indelt] = indNewFace+1;
          indelt++;

          // coordonnees du nouvel element tetraedrique
          for (E_Int p = 0;p < nfld; p++) 
          {
            fnewp[p][ind1] = fp[p][ind1]; // premier point de l arete de l ancien element
            fnewp[p][ind2] = fp[p][ind2]; // second point de l arete de l ancien element
          }
        }
      }
    }
  }

  // Nettoyage de la connectivite creee
  K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-12, newEltType, fnew, newcn);

  // Objet python retourne
  PyObject* tpl = K_ARRAY::buildArray(fnew, varString, newcn, -1, newEltType);

  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}
// ============================================================================
/* Convert  non-polyhedral 3D array to a tetraedrical mesh with addition of 
   points (barycenter of elements and faces) 
   The method deals with fields at nodes AND centers */
// ============================================================================
PyObject* K_CONVERTER::convertArray2TetraBaryBoth(PyObject* self, PyObject* args)
{
  PyObject *array, *arrayc;
  if (!PyArg_ParseTuple(args, "OO", &array, &arrayc)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString, 
                               f, ni, nj, nk, cn, eltType);
  
  if (res != 2 && res != 1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2TetraBary: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2TetraBary: array must be unstructured.");
    return NULL;
  }

  if (strcmp(eltType, "NGON") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "convertArray2TetraBary: array must not be NGon.");
    return NULL;
  }

  if ((strcmp(eltType, "NODE") == 0)||
      // (strcmp(eltType, "BAR") == 0) ||
      (strcmp(eltType, "TRI") == 0) ||
      (strcmp(eltType, "TETRA") == 0))
  {
    PyObject* tpl = K_ARRAY::buildArray(*f, varString, *cn, -1, eltType);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  E_Int nic, njc, nkc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  E_Int resc = K_ARRAY::getFromArray3(arrayc, varStringc, 
                                      fc, nic, njc, nkc, cnc, eltTypec);
  if (resc != 1 && resc != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2TetraBary: arrayc is invalid.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }
  if (resc == 1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2TetraBary: arrayc must be unstructured.");
    RELEASESHAREDU(array, f, cn);
    RELEASESHAREDS(arrayc, fc);
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if (posx == 0 || posy == 0 || posz == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2TetraBary: coordinates not found.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }

  // nb d'elements
  E_Int nelts = cn->getSize(); 

  // specific datas depending on element type
  E_Int newElts=0; // nb nouveaux elements tetra
  E_Int nfaces=0; // nb de faces par element
  E_Int npoints=0; // nb points par element
  E_Int dim = 3; 
  if (strcmp(eltType, "QUAD") == 0)
  {
    newElts = 4*nelts; nfaces = 0; npoints = 4; // nfaces = 0 au lieu de 4 car on ajoute pas de points sur les faces
    dim = 2;
  }
  else if (strcmp(eltType, "HEXA") == 0)
  {
    newElts = 24*nelts; nfaces = 6; npoints = 8;
    dim = 3;
  }
  else if (strcmp(eltType, "PENTA") == 0)
  {
    newElts = 18*nelts; nfaces = 5; npoints = 6;
    dim = 3;  
  }
  else if (strcmp(eltType, "PYRA") == 0)
  {
    newElts = 16*nelts; nfaces = 5; npoints = 5;
    dim = 3;
  }
  else if (strcmp(eltType, "BAR") == 0)
  {
    newElts = 2*nelts; nfaces = 0; npoints = 2;
    dim = 1;
  }

  // Calcul de la taille du tableau des nouveaux champs
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();
  E_Int sizeFields = npts;

  for (E_Int elt = 0; elt < nelts; elt++)
  {
    for (E_Int fa = 0; fa < nfaces; fa++) {sizeFields += 1;}
    sizeFields += 1;
  }

  // nouveau champ
  FldArrayF fnew(sizeFields, nfld);
  E_Int nfldc = fc->getNfld();
  FldArrayF fcnew(newElts, nfldc);

  // pointeurs sur le nouveau champ
  vector<E_Float*> fnewp(nfld);
  for (E_Int p=0;p<nfld;p++) {fnewp[p] = fnew.begin(p+1);}

  vector<E_Float*> fcnewp(nfldc);
  for (E_Int p=0; p<nfldc;p++) {fcnewp[p] = fcnew.begin(p+1);}

  // pointeurs sur l ancien champ
  vector<E_Float*> fp(nfld);
  for (E_Int p=0;p<nfld;p++) {fp[p] = f->begin(p+1);}  
  vector<E_Float*> fcp(nfldc);
  for (E_Int p = 0; p < nfldc; p++) {fcp[p] = fc->begin(p+1);}

  // Nouvelle connectivite
  E_Int nconnect = 2; // nb de points par element
  if (dim == 2) nconnect = 3; // elmt triangle
  else if (dim == 1) nconnect = 2;// elt BAR
  else nconnect = 4;          // elmt tetraedre
  FldArrayI newcn(newElts,nconnect);  
  // Pointeurs sur la nouvelle connectivite
  vector<E_Int*> newcnp(nconnect);
  for (E_Int p=0; p<nconnect; p++) {newcnp[p] = newcn.begin(p+1);}
  // Pointeurs sur l ancienne connectivite
  vector<E_Int*> cnp(npoints);
  for (E_Int p = 0; p < npoints; p++) {cnp[p] = cn->begin(p+1);}

  // Type d'element de la nouvelle connectivite
  char newEltType[256];
  E_Float* fpp; E_Float* fnpp;

  if (dim == 1) // BAR
  {
    strcpy(newEltType,"BAR");
    E_Int indp1, indp2; 
    E_Int indpt = npts; //indices des nouveaux points
    E_Int indelt = 0; // indice d elt de la nouvelle connectivite

    for (E_Int et = 0; et < nelts; et++)
    {
      indp1 = cnp[0][et]-1;
      indp2 = cnp[1][et]-1;
    
      // Barycentre de l'element
      for (E_Int p = 0; p < nfld; p++) 
        fnewp[p][indpt] = 0.5*(fp[p][indp1]+fp[p][indp2]);      
      indpt++; 

      for (E_Int p = 0; p < nfldc; p++) 
      {
        fnpp = fcnewp[p]; fpp = fcp[p]; 
        fnpp[indelt] = fpp[et];
        fnpp[indelt+1] = fpp[et];
      }

      // construction des nouveaux elements
      // connectivite
      newcnp[0][indelt] = indp1+1;
      newcnp[1][indelt] = indpt;
      indelt++;
      
      
      newcnp[0][indelt] = indpt;
      newcnp[1][indelt] = indp2+1;
      indelt++;
      //coordonnees des 2 nouveaux elements
      for (E_Int p = 0; p < nfld; p++) 
      {
        fnpp = fnewp[p]; fpp = fp[p]; 
        fnpp[indp1] = fpp[indp1];
        fnpp[indp2] = fpp[indp2];       
      }
    }
  }
  else if (dim == 2) // Elmt QUAD
  {
    strcpy(newEltType, "TRI");

    // Calcul de la nouvelle connectivite et des nouvelles coordonnees
    E_Int ind1, ind2, ind3, ind4; // indices des points la cellule
    E_Int indelt = 0; // indice d element de la nouvelle connectivite
    E_Int indpt = npts; // indice des nouveaux points
    E_Int indNewElt; // indices du nouvel elmt et de la nouvelle face
    for (E_Int elt = 0; elt < nelts; elt++)
    {
      ind1 = cnp[0][elt]-1;
      ind2 = cnp[1][elt]-1;
      ind3 = cnp[2][elt]-1;
      ind4 = cnp[3][elt]-1;

      // calcul du barycentre de l'element
      for (E_Int p = 0; p < nfld; p++) 
      {
        fnewp[p][indpt] = K_CONST::ONE_FOURTH*(fp[p][ind1]+fp[p][ind2]+fp[p][ind3]+fp[p][ind4]);
      }
      indNewElt = indpt; indpt++; 

      // building new tetra elts 
      // fields located at centers
      for (E_Int p = 0; p < nfldc; p++) 
      {
        fnpp = fcnewp[p]; fpp = fcp[p]; 
        fnpp[indelt] = fpp[elt];
        fnpp[indelt+1] = fpp[elt];
        fnpp[indelt+2] = fpp[elt];
        fnpp[indelt+3] = fpp[elt];
      }  

      // connectivite des nouveaux elements
      newcnp[0][indelt] = ind1+1;
      newcnp[1][indelt] = ind2+1;
      newcnp[2][indelt] = indNewElt+1; indelt++;
      
      newcnp[0][indelt] = ind2+1;
      newcnp[1][indelt] = ind3+1;
      newcnp[2][indelt] = indNewElt+1; indelt++;
      
      newcnp[0][indelt] = ind3+1;
      newcnp[1][indelt] = ind4+1;
      newcnp[2][indelt] = indNewElt+1; indelt++;
      
      newcnp[0][indelt] = ind4+1;
      newcnp[1][indelt] = ind1+1;
      newcnp[2][indelt] = indNewElt+1; indelt++;
      
     // coordonnees du nouvel element tetraedrique
      for (E_Int p = 0; p < nfld; p++) 
      {
        E_Float* fnewpp = fnewp[p];
        E_Float* fpp = fp[p]; 
        fnewpp[ind1] = fpp[ind1];
        fnewpp[ind2] = fpp[ind2];
        fnewpp[ind3] = fpp[ind3];
        fnewpp[ind4] = fpp[ind4];
      }
    }
  }
  else // Elmts 3D
  {
    strcpy(newEltType, "TETRA");

    vector< vector<E_Int> > faces(nfaces); // vecteur contenant pour un element les indices "locaux" des points de chaque face
    buildFaceIndices(eltType, faces);

    // Calcul de la nouvelle connectivite et des nouvelles coordonnees
    vector<E_Float> fbe(nfld); // champs du barycentre de l element
    vector<E_Float> fbf(nfld); // champs du barycentre de la face
    E_Int ind, ind1, ind2; // indices de noeud 
    E_Int indelt = 0; // indice d element de la nouvelle connectivite
    E_Int indpt = npts; // indice des nouveaux points
    E_Int indNewElt, indNewFace; // indices du nouvel elmt et de la nouvelle face
    E_Int nptsInFace; // nb points pour la face donnee
    for (E_Int elt = 0; elt < nelts; elt++)
    {
      // calcul du barycentre de l'element
      for (E_Int p = 0; p < nfld; p++) 
      {
        fbe[p] = 0.;
        for (E_Int n = 0; n < npoints; n++) 
        {
          ind = cnp[n][elt] -1;
          fbe[p] += fp[p][ind];
        }
        fbe[p] = fbe[p]/npoints;
        
        fnewp[p][indpt] = fbe[p];
      }
      indNewElt = indpt; indpt++; 
      for (E_Int fa = 0; fa < nfaces; fa++)
      {
        nptsInFace = faces[fa].size();
        // calcul du barycentre de la face
        for (E_Int p = 0; p < nfld; p++) 
        {
          fbf[p] = 0.;
          for (E_Int n = 0; n < nptsInFace; n++) 
          {
            ind = cnp[faces[fa][n]][elt] -1;
            fbf[p] += fp[p][ind];
          }
          fbf[p] = fbf[p]/nptsInFace;
          fnewp[p][indpt] = fbf[p];
        }
        indNewFace = indpt; indpt ++;
        // construction des nouveaux elements tetraedriques
        for (E_Int n = 0; n < nptsInFace; n++)
        {
          ind1 =  cnp[faces[fa][n]][elt] -1;
          ind2 =  cnp[faces[fa][(n+1)%(nptsInFace)]][elt] -1;

          // connectivite du nouvel element
          newcnp[0][indelt] = ind1+1;
          newcnp[1][indelt] = ind2+1;
          newcnp[2][indelt] = indNewElt+1;
          newcnp[3][indelt] = indNewFace+1;

          // fields located at centers
          for (E_Int p = 0; p < nfldc; p++) 
          {
            fnpp = fcnewp[p]; fpp = fcp[p]; 
            fnpp[indelt] = fpp[elt];
          }
          indelt++;

          // coordonnees du nouvel element tetraedrique
          for (E_Int p = 0;p < nfld; p++) 
          {
            fnewp[p][ind1] = fp[p][ind1]; // premier point de l arete de l ancien element
            fnewp[p][ind2] = fp[p][ind2]; // second point de l arete de l ancien element
          }
        }
      }
    }
  }

  // Nettoyage de la connectivite creee
  K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-12, newEltType, fnew, newcn);
  // Objet python retourne
  PyObject* l = PyList_New(0);

  PyObject* tpl1 = K_ARRAY::buildArray(fnew, varString, newcn, -1, newEltType);
  PyList_Append(l, tpl1); Py_DECREF(tpl1);
  PyObject* tpl2 = K_ARRAY::buildArray(fcnew, varStringc, newcn, -1, newEltType);
  PyList_Append(l, tpl2); Py_DECREF(tpl2);
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDB(resc, arrayc, fc, cnc);
  return l;
}

//=============================================================================
// build local indices for each faces of elements, following CGNS numerotation
//=============================================================================
void K_CONVERTER::buildFaceIndices(char* eltType, vector< vector<E_Int> >& faces)
{ 
  if (strcmp(eltType, "HEXA") == 0)
  {
    faces[0].push_back(0);faces[0].push_back(1);faces[0].push_back(2);faces[0].push_back(3);
    faces[1].push_back(4);faces[1].push_back(5);faces[1].push_back(6);faces[1].push_back(7);
    faces[2].push_back(0);faces[2].push_back(1);faces[2].push_back(5);faces[2].push_back(4);
    faces[3].push_back(3);faces[3].push_back(2);faces[3].push_back(6);faces[3].push_back(7);
    faces[4].push_back(0);faces[4].push_back(3);faces[4].push_back(7);faces[4].push_back(4);
    faces[5].push_back(1);faces[5].push_back(2);faces[5].push_back(6);faces[5].push_back(5);
  }
  else if (strcmp(eltType, "TETRA") == 0)
  {
    faces[0].push_back(0);faces[0].push_back(1);faces[0].push_back(2);
    faces[1].push_back(0);faces[1].push_back(1);faces[1].push_back(3);
    faces[2].push_back(1);faces[2].push_back(2);faces[2].push_back(3);
    faces[3].push_back(2);faces[3].push_back(0);faces[3].push_back(3);
  }
  else if (strcmp(eltType, "PYRA") == 0)
  {
    faces[0].push_back(0);faces[0].push_back(1);faces[0].push_back(2);faces[0].push_back(3);
    faces[1].push_back(0);faces[1].push_back(1);faces[1].push_back(4);
    faces[2].push_back(1);faces[2].push_back(2);faces[2].push_back(4);
    faces[3].push_back(2);faces[3].push_back(3);faces[3].push_back(4);
    faces[4].push_back(3);faces[4].push_back(0);faces[4].push_back(4);
  }
  else if (strcmp(eltType, "PENTA") == 0)
  {
    faces[0].push_back(0);faces[0].push_back(1);faces[0].push_back(2);
    faces[1].push_back(3);faces[1].push_back(4);faces[1].push_back(5);
    faces[2].push_back(0);faces[2].push_back(1);faces[2].push_back(4);faces[2].push_back(3);
    faces[3].push_back(1);faces[3].push_back(2);faces[3].push_back(5);faces[3].push_back(4);
    faces[4].push_back(0);faces[4].push_back(2);faces[4].push_back(5);faces[4].push_back(3);
  }
}
