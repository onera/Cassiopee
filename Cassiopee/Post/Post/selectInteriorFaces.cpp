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

// selectInteriorFaces

# include "post.h"
# include <stdio.h>
# include <string.h>

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
/* Selectionne les facettes interieures d'un array.
   Si l'argument strict vaut 1, les facettes ayant uniquement des noeuds 
   interieurs sont prises en compte.
*/
// ============================================================================
PyObject* K_POST::selectInteriorFaces(PyObject* self, PyObject* args)
{
  PyObject* array; E_Int strict;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &strict))
  {
    return NULL;
  }
  // Extract array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cnp;
  E_Int res;
  E_Int ni, nj, nk;
  res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cnp, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectInteriorFaces: array is invalid.");
    return NULL;
  }

  // Verifications
  E_Int nfaces = 0; // nombre de facettes par element
  if (strcmp(eltType, "TRI") == 0)
  {
    nfaces = 3;
  }
  else if (strcmp(eltType, "QUAD") == 0) 
  {
    nfaces = 4;
  }
  else
  {
    RELEASESHAREDB(res, array, f, cnp);
    PyErr_SetString(PyExc_TypeError,
                    "selectInteriorFaces: only for TRI and QUAD.");
    return NULL;
  }
  
  if (strict != 0 && strict != 1)
  {
    printf("Warning: selectInteriorFaces: strict must be equal to 0 or 1. Set to default value 0.\n");
    strict = 0;          
  }

  // Calcul les facettes interieures
  E_Int node1, node2;
  E_Int nb = f->getSize(); // nbre de noeuds
  E_Int nfld = f->getNfld();// nb de champs
  FldArrayI& cn = *cnp;
  E_Int ne = cn.getSize(); // nbre d'elements
  E_Int nt = cn.getNfld(); // taille des elements
  FldArrayIS faces(nb*nb); faces.setAllValuesAtNull();
  
  for (E_Int e = 0; e < ne; e++)
  {
    for (E_Int n = 1; n < nt; n++)
    {
      node1 = cn(e, n);
      node2 = cn(e, n+1);
      faces[(node1-1)+(node2-1)*nb]++; 
    }
    node1 = cn(e, nt);
    node2 = cn(e, 1);
    faces[(node1-1)+(node2-1)*nb]++; 
    }

  // Liste des facettes interieures
  FldArrayF* fnodes = new FldArrayF(ne*2*nfaces, nfld);
  FldArrayI* connect = new FldArrayI(ne*nfaces, 2);
  E_Int cf = 0; //compteur facettes
  E_Int cc = 0; //compteur connect
  switch (strict) 
  {
    case 0:
      for (E_Int f1 = 0; f1 < nb; f1++)
      {
        for (E_Int f2 = 0; f2 < nb; f2++)
        {
          if (faces[f1+f2*nb] + faces[f2+f1*nb] == 2) // facette non exterieure
          {
            faces[f1+f2*nb] = 0;
            faces[f2+f1*nb] = 0;
            for (E_Int n = 1; n <= nfld; n++)
              (*fnodes)(cf, n) = (*f)(f1, n);
            for (E_Int n = 1; n <= nfld; n++)
              (*fnodes)(cf+1, n) = (*f)(f2, n);
            (*connect)(cc, 1) = cf+1;
            (*connect)(cc, 2) = cf+2;
            cf = cf+2;
            cc = cc+1;
          }
        }
      }
      break;
  
    case 1: //noeuds interieurs
    
      // tag des noeuds exterieurs
      FldArrayIS tag(nb); tag.setAllValuesAtNull(); //tag vaut 1 si noeud ext
      for (E_Int f1 = 0; f1 < nb; f1++)
      {
        for (E_Int f2 = 0; f2 < nb; f2++)
        {
          if (f1 != f2) 
          {
            if (faces[f1+f2*nb] + faces[f2+f1*nb] == 1) // tag des noeuds exterieurs
            {
              tag[f1] = 1; tag[f2] = 1;
            }
          }
        }
      }

      for (E_Int f1 = 0; f1 < nb; f1++)
      {
        for (E_Int f2 = 0; f2 < nb; f2++)
        {
          if (tag[f1] == 0 && tag[f2] == 0 && 
              faces[f1+f2*nb] + faces[f2+f1*nb] == 2) 
          {
            faces[f1+f2*nb] = 0;
            faces[f2+f1*nb] = 0;
            for (E_Int n = 1; n <= nfld; n++)
              (*fnodes)(cf, n) = (*f)(f1, n);
            for (E_Int n = 1; n <= nfld; n++)
              (*fnodes)(cf+1, n) = (*f)(f2, n);
            (*connect)(cc, 1) = cf+1;
            (*connect)(cc, 2) = cf+2;
            cf = cf+2;
            cc = cc+1;
          }
        }
      }
      break;
  }
  fnodes->reAllocMat(cf, nfld);
  connect->reAllocMat(cc, 2);

  // Build array
  PyObject* tpl;
  tpl = K_ARRAY::buildArray(*fnodes, varString, *connect, -1, "BAR");
  delete fnodes; delete connect;
  RELEASESHAREDB(res, array, f, cnp);

  return tpl;
}
