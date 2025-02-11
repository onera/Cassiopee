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

// Convert polyhedral array to tetra array
# include <vector>
# include <stdlib.h>
# include "converter.h"
# include "kcore.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Convert  polyhedral array to a tetraedrical mesh */
// ============================================================================
PyObject* K_CONVERTER::convertNGon2TetraBary(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  res = K_ARRAY::getFromArray(array, varString, 
                              f, ni, nj, nk, cn, eltType, true);
  if (res == 1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: array must be unstructured NGon.");
    return NULL; 
  }
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: array is invalid.");
    return NULL;
  }
  if (strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "convertNGon2TetraBary: array must be NGon.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if (posx == 0 || posy == 0 || posz == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: coordinates not found.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }

  // Pointeurs sur la connectivite
  // cn1: connectivite face/noeuds
  // cn2: connectivite elmt/faces
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin();
  E_Int sizeFN = cn1[1]; // taille de la connectivite Face/Noeuds
  E_Int nelts = cn2[sizeFN+2]; // nombre d'elements pour la connectivite cn
  cn2 += sizeFN+4;
  
  if ( sizeFN == 0 || nelts == 0 ) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: empty connectivity.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  // Dimensions des champs
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();

  // Tableau de position des faces dans la connectivite
  FldArrayI posFace;
  K_CONNECT::getPosFaces(*cn, posFace);

  // tableau donnant pour chaque element sa dimension
  FldArrayI dimElts;
  K_CONNECT::getDimElts(*cn, posFace, dimElts);
  // dimension des elements (dimension unique dans l'array)
  E_Int dim = dimElts[0];
  // Type d'element de la nouvelle connectivite
  char newEltType[256];
  E_Int nconnect; // nb de points par element
  if (dim == 1) {nconnect = 2; strcpy(newEltType,"BAR");}
  else if (dim == 2) {nconnect = 3; strcpy(newEltType,"TRI");}
  else if (dim == 3) {nconnect = 4; strcpy(newEltType,"TETRA");}
  else 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: wrong dimension for elements.");
    RELEASESHAREDU(array, f, cn); return NULL;    
  }
 
  // Connectivite Element/Noeuds
  vector< vector<E_Int> > cnEV(nelts);
  K_CONNECT::connectNG2EV(*cn, cnEV);

  // Calcul de la taille du tableau des nouveaux champs
  E_Int sizeConnect = 0; E_Int sizeFields = npts;
  E_Int nbFaces, nbNodes, numface, pos;
  switch (dim)
  {
    case 1:
    {
      sizeConnect = nelts;
    }
    break;
    case 2:
    {
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        nbFaces = cn2[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          numface = cn2[fa+1]-1;
          pos = posFace[numface];
          nbNodes = cn1[pos];
          for (E_Int n=0; n < nbNodes; n++) {sizeConnect += 1;}
        }
        sizeFields += 1;
        cn2 += nbFaces+1;
      }
    }
    break;
    case 3:
    {
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        nbFaces = cn2[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          numface = cn2[fa+1]-1;
          pos = posFace[numface];
          nbNodes = cn1[pos];
          for (E_Int n=0; n < nbNodes; n++) {sizeConnect += 1;}
          sizeFields += 1;
        }
        sizeFields += 1;
        cn2 += nbFaces+1;
      }
    }
    break;
  }

  // Nouvelle connectivite
  FldArrayI newcn(sizeConnect,nconnect); // tetraedre
  // Pointeurs sur la nouvelle connectivite
  vector<E_Int*> newcnp(nconnect);
  for (E_Int p = 0; p < nconnect; p++) {newcnp[p] = newcn.begin(p+1);}
  // Nouveau champ
  FldArrayF fnew(sizeFields, nfld);
  // Pointeurs sur les nouveaux champs
  vector<E_Float*> fnewp(nfld);
  for (E_Int p=0; p<nfld; p++) {fnewp[p] = fnew.begin(p+1);}
  // pointeurs sur l'ancien champ
  vector<E_Float*> fp(nfld);
  for (E_Int p=0; p < nfld; p++) {fp[p] = f->begin(p+1);}  
  
  // Calcul de la nouvelle connectivite et des nouveaux champs
  vector<E_Float> fbe(nfld); // champs du barycentre de l element
  vector<E_Float> fbf(nfld); // champs du barycentre de la face
  E_Int ind, ind1, ind2; // indices de noeud 
  E_Int indelt = 0; // indice d element de la nouvelle connectivite
  E_Int indpt = npts; // indice des nouveaux points
  E_Int npoints; // nb points par element ou par face selon le contexte
  E_Int indNewElt, indNewFace; // indices du nouvel elmt et de la nouvelle face
  cn2 = cn->begin(); cn2 += sizeFN+4;

  switch (dim)
  {
    case 1:
    {
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        nbFaces = cn2[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          numface = cn2[fa+1]-1;
          pos = posFace[numface];
          npoints = cn1[pos]; pos++;
          ind = cn1[pos]-1;
          // connectivite du nouvel element BAR
          newcnp[fa][indelt] = ind+1;
            // champs du nouvel element BAR
          for (E_Int p=0;p<nfld;p++) 
          {
            fnewp[p][ind] = fp[p][ind]; // premier point de l'arete de l ancien element
          }
        }
        indelt++;
        cn2 += nbFaces+1;
      }
    }
    break;
    case 2:
    {
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        vector<E_Int>& vertices = cnEV[elt];// sommets associes a l elt
        npoints = vertices.size();
        for (E_Int p=0;p<nfld;p++) 
        {
          fbe[p] = 0.;
          for (E_Int n=0; n<npoints; n++)
          {
            ind = vertices[n]-1;
            fbe[p] += fp[p][ind];
          }
          fbe[p] = fbe[p]/npoints;
          fnewp[p][indpt] = fbe[p];
        }
        indNewElt = indpt; indpt ++; 
        nbFaces = cn2[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          numface = cn2[fa+1]-1;
          pos = posFace[numface];
          npoints = cn1[pos]; pos++;
          // construction des nouveaux elements triedriques
          for (E_Int n=0; n < npoints; n++)
          {
            ind1 = cn1[pos+n]-1;
            ind2 = cn1[pos+(n+1)%(npoints)]-1;

            // connectivite du nouvel element
            newcnp[0][indelt] = ind1+1;
            newcnp[1][indelt] = ind2+1;
            newcnp[2][indelt] = indNewElt+1;
            indelt++;

            // champs du nouvel element triedrique
            for (E_Int p=0;p<nfld;p++) 
            {
              fnewp[p][ind1] = fp[p][ind1]; // premier point de l arete de l ancien element
              fnewp[p][ind2] = fp[p][ind2]; // second point de l arete de l ancien element
            }
          }
        }
        cn2 += nbFaces+1;
      }
    }
    break;
    case 3:
    {
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        vector<E_Int>& vertices = cnEV[elt];// sommets associes a l'elt
        npoints = vertices.size();
        for (E_Int p=0;p<nfld;p++) 
        {
          fbe[p] = 0.;
          for (E_Int n=0; n<npoints; n++)
          {
            ind = vertices[n]-1;
            fbe[p] += fp[p][ind];
          }
          fbe[p] = fbe[p]/npoints;
          fnewp[p][indpt] = fbe[p];
        }
        indNewElt = indpt; indpt++; 
        nbFaces = cn2[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          numface = cn2[fa+1]-1;
          pos = posFace[numface];
          npoints = cn1[pos]; pos++;
          // calcul du barycentre de la face
          for (E_Int p=0;p<nfld;p++) 
          {
            fbf[p] = 0.;
            for (E_Int n=0; n < npoints; n++) 
            {
              ind = cn1[pos+n]-1;
              fbf[p] += fp[p][ind];
            }
            fbf[p] = fbf[p]/npoints;
 
            fnewp[p][indpt] = fbf[p];
          }
          indNewFace = indpt; indpt++;
          // construction des nouveaux elements tetraedriques
          for (E_Int n=0; n < npoints; n++)
          {
            ind1 = cn1[pos+n]-1;
            ind2 = cn1[pos+(n+1)%(npoints)]-1;

            // connectivite du nouvel element
            newcnp[0][indelt] = ind1+1;
            newcnp[1][indelt] = ind2+1;
            newcnp[2][indelt] = indNewElt+1;
            newcnp[3][indelt] = indNewFace+1;
            indelt++;

            // champs du nouvel element tetraedrique
            for (E_Int p=0;p<nfld;p++) 
            {
              fnewp[p][ind1] = fp[p][ind1]; // premier point de l arete de l ancien element
              fnewp[p][ind2] = fp[p][ind2]; // second point de l arete de l ancien element
            }
          }
        }
        cn2 += nbFaces+1;
      }
    }
    break;
  }
  
  K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-12, newEltType, fnew, newcn);

  // Objet python retourne
  PyObject* tpl = K_ARRAY::buildArray(fnew, varString, newcn, -1, newEltType);

  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}
// ============================================================================
/* Convert  polyhedral array to a tetraedrical mesh */
// ============================================================================
PyObject* K_CONVERTER::convertNGon2TetraBaryBoth(PyObject* self, PyObject* args)
{
  PyObject *array, *arrayc;
  if (!PyArg_ParseTuple(args, "OO", &array, &arrayc)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  res = K_ARRAY::getFromArray(array, varString, 
                              f, ni, nj, nk, cn, eltType, true);
  if (res == 1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: array must be unstructured NGon.");
    return NULL; 
  }
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: array is invalid.");
    return NULL;
  }
  if (strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "convertNGon2TetraBary: array must be NGon.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if (posx == 0 || posy == 0 || posz == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: coordinates not found.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  E_Int nic, njc, nkc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  E_Int resc = K_ARRAY::getFromArray(arrayc, varStringc, 
                                     fc, nic, njc, nkc, cnc, eltTypec, true);
  if (resc != 1 && resc != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: arrayc is invalid.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }
  if (resc == 1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: arrayc must be unstructured.");
    RELEASESHAREDU(array, f, cn);
    RELEASESHAREDS(arrayc, fc);
    return NULL;
  }
  if (strcmp(eltTypec, "NGON*") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    RELEASESHAREDU(arrayc, fc, cnc);
    PyErr_SetString(PyExc_ValueError,
                    "convertNGon2TetraBary: arrayc must be NGon.");
    return NULL;
  }

  // Pointeurs sur la connectivite
  // cn1: connectivite face/noeuds
  // cn2: connectivite elmt/faces
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin();
  E_Int sizeFN = cn1[1]; // taille de la connectivite Face/Noeuds
  E_Int nelts = cn2[sizeFN+2]; // nombre d'elements pour la connectivite cn
  cn2 += sizeFN+4;
  
  if ( sizeFN == 0 || nelts == 0 ) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: empty connectivity.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  // Dimensions des champs
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();
  E_Int nfldc = fc->getNfld();

  // Tableau de position des faces dans la connectivite
  FldArrayI posFace;
  K_CONNECT::getPosFaces(*cn, posFace);

  // tableau donnant pour chaque element sa dimension
  FldArrayI dimElts;
  K_CONNECT::getDimElts(*cn, dimElts);
  // dimension des elements (dimension unique dans l array)
  E_Int dim = dimElts[0];
  // Type d element de la nouvelle connectivite
  char newEltType[256];
  E_Int nconnect; // nb de points par element
  if (dim == 1) {nconnect = 2; strcpy(newEltType,"BAR");}
  else if (dim == 2) {nconnect = 3; strcpy(newEltType,"TRI");}
  else if (dim == 3) {nconnect = 4; strcpy(newEltType,"TETRA");}
  else 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertNGon2TetraBary: wrong dimension for elements.");
    RELEASESHAREDU(array, f, cn); return NULL;    
  }
 
  // Connectivite Element/Noeuds
  vector< vector<E_Int> > cnEV(nelts);
  K_CONNECT::connectNG2EV(*cn, cnEV);

  // Calcul de la taille du tableau des nouveaux champs
  E_Int sizeConnect = 0; E_Int sizeFields = npts;
  E_Int nbFaces, nbNodes, numface, pos;
  switch (dim)
  {
    case 1:
    {
      sizeConnect = nelts;
    }
    break;
    case 2:
    {
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        nbFaces = cn2[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          numface = cn2[fa+1]-1;
          pos = posFace[numface];
          nbNodes = cn1[pos];
          for (E_Int n=0; n < nbNodes; n++) {sizeConnect += 1;}
        }
        sizeFields += 1;
        cn2 += nbFaces+1;
      }
    }
    break;
    case 3:
    {
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        nbFaces = cn2[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          numface = cn2[fa+1]-1;
          pos = posFace[numface];
          nbNodes = cn1[pos];
          for (E_Int n=0; n < nbNodes; n++) {sizeConnect += 1;}
          sizeFields += 1;
        }
        sizeFields += 1;
        cn2 += nbFaces+1;
      }
    }
    break;
  }

  // Nouvelle connectivite
  FldArrayI newcn(sizeConnect,nconnect); // tetraedre
  // Pointeurs sur la nouvelle connectivite
  vector<E_Int*> newcnp(nconnect);
  for (E_Int p=0; p<nconnect; p++) {newcnp[p] = newcn.begin(p+1);}
  // Nouveau champ
  FldArrayF fnew(sizeFields, nfld);
  FldArrayF fcnew(sizeConnect,nfldc);

  // Pointeurs sur les nouveaux champs
  vector<E_Float*> fnewp(nfld);
  for (E_Int p=0; p<nfld; p++) {fnewp[p] = fnew.begin(p+1);}
  vector<E_Float*> fcnewp(nfld);
  for (E_Int p=0; p<nfldc; p++) {fcnewp[p] = fcnew.begin(p+1);}

  // pointeurs sur l ancien champ
  vector<E_Float*> fp(nfld);
  for (E_Int p=0;p<nfld;p++) {fp[p] = f->begin(p+1);}  
  vector<E_Float*> fcp(nfldc);
  for (E_Int p=0;p<nfldc;p++) {fcp[p] = fc->begin(p+1);}  
  
  // Calcul de la nouvelle connectivite et des nouveaux champs
  vector<E_Float> fbe(nfld); // champs du barycentre de l element
  vector<E_Float> fbf(nfld); // champs du barycentre de la face
  E_Int ind, ind1, ind2; // indices de noeud 
  E_Int indelt = 0; // indice d element de la nouvelle connectivite
  E_Int indpt = npts; // indice des nouveaux points
  E_Int npoints; // nb points par element ou par face selon le contexte
  E_Int indNewElt, indNewFace; // indices du nouvel elmt et de la nouvelle face
  cn2 = cn->begin(); cn2 += sizeFN+4;

  switch (dim)
  {
    case 1:
    {
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        nbFaces = cn2[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          numface = cn2[fa+1]-1;
          pos = posFace[numface];
          npoints = cn1[pos]; pos++;
          ind = cn1[pos]-1;
          // connectivite du nouvel element BAR
          newcnp[fa][indelt] = ind+1;
                    
          // champs en centres du nouvel element BAR
          for (E_Int p=0;p<nfldc;p++) 
            fcnewp[p][indelt] = fcp[p][elt]; 
          
          // champs du nouvel element BAR
          for (E_Int p=0;p<nfld;p++) 
          {
            fnewp[p][ind] = fp[p][ind]; // premier point de l arete de l ancien element
          }
        }
        indelt++;
        cn2 += nbFaces+1;
      }
    }
    break;
    case 2:
    {
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        vector<E_Int>& vertices = cnEV[elt];// sommets associes a l elt
        npoints = vertices.size();
        for (E_Int p=0;p<nfld;p++) 
        {
          fbe[p] = 0.;
          for (E_Int n=0; n<npoints; n++)
          {
            ind = vertices[n]-1;
            fbe[p] += fp[p][ind];
          }
          fbe[p] = fbe[p]/npoints;
          fnewp[p][indpt] = fbe[p];
        }
        indNewElt = indpt; indpt ++; 
        nbFaces = cn2[0];

        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          numface = cn2[fa+1]-1;
          pos = posFace[numface];
          npoints = cn1[pos]; pos++;
          // construction des nouveaux elements triedriques
          for (E_Int n=0; n < npoints; n++)
          {
            ind1 = cn1[pos+n]-1;
            ind2 = cn1[pos+(n+1)%(npoints)]-1;
            
            //champs en centres
            for (E_Int p=0;p<nfldc;p++) 
              fcnewp[p][indelt] = fcp[p][elt];

            // connectivite du nouvel element
            newcnp[0][indelt] = ind1+1;
            newcnp[1][indelt] = ind2+1;
            newcnp[2][indelt] = indNewElt+1;
            indelt++;

            // champs du nouvel element TRI
            for (E_Int p=0;p<nfld;p++) 
            {
              fnewp[p][ind1] = fp[p][ind1]; // premier point de l arete de l ancien element
              fnewp[p][ind2] = fp[p][ind2]; // second point de l arete de l ancien element
            }
          }
        }
        cn2 += nbFaces+1;
      }
    }
    break;
    case 3:
    {
      for (E_Int elt = 0; elt < nelts; elt++)
      {
        vector<E_Int>& vertices = cnEV[elt];// sommets associes a l elt
        npoints = vertices.size();
        for (E_Int p=0;p<nfld;p++) 
        {
          fbe[p] = 0.;
          for (E_Int n=0; n<npoints; n++)
          {
            ind = vertices[n]-1;
            fbe[p] += fp[p][ind];
          }
          fbe[p] = fbe[p]/npoints;
          fnewp[p][indpt] = fbe[p];
        }
        indNewElt = indpt; indpt ++; 
        nbFaces = cn2[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          numface = cn2[fa+1]-1;
          pos = posFace[numface];
          npoints = cn1[pos]; pos++;
          // calcul du barycentre de la face
          for (E_Int p=0;p<nfld;p++) 
          {
            fbf[p] = 0.;
            for (E_Int n=0; n < npoints; n++) 
            {
              ind = cn1[pos+n]-1;
              fbf[p] += fp[p][ind];
            }
            fbf[p] = fbf[p]/npoints;
 
            fnewp[p][indpt] = fbf[p];
          }
          indNewFace = indpt; indpt ++;
          // construction des nouveaux elements tetraedriques
          for (E_Int n=0; n < npoints; n++)
          {
            ind1 = cn1[pos+n]-1;
            ind2 = cn1[pos+(n+1)%(npoints)]-1;
             
            // champs en centres du nouvel element TETRA
            for (E_Int p=0;p<nfldc;p++) 
              fcnewp[p][indelt] = fcp[p][elt];

            // connectivite du nouvel element
            newcnp[0][indelt] = ind1+1;
            newcnp[1][indelt] = ind2+1;
            newcnp[2][indelt] = indNewElt+1;
            newcnp[3][indelt] = indNewFace+1;
            indelt++;

            // champs du nouvel element tetraedrique
            for (E_Int p=0;p<nfld;p++) 
            {
              fnewp[p][ind1] = fp[p][ind1]; // premier point de l arete de l ancien element
              fnewp[p][ind2] = fp[p][ind2]; // second point de l arete de l ancien element
            }
          }
        }
        cn2 += nbFaces+1;
      }
    }
    break;
  }
  K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-12, newEltType, fnew, newcn);
  PyObject* l = PyList_New(0);

  PyObject* tpl1 = K_ARRAY::buildArray(fnew, varString, newcn, -1, newEltType);
  PyList_Append(l, tpl1); Py_DECREF(tpl1);
  PyObject* tpl2 = K_ARRAY::buildArray(fcnew, varStringc, newcn, -1, newEltType);
  PyList_Append(l, tpl2); Py_DECREF(tpl2);
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDB(resc, arrayc, fc, cnc);
  return l;
 
}
