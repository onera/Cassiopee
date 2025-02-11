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

// Selectionne les faces qui sont a la frontiere d'un tag=0 et tag=1

# include <stdio.h>
# include <string.h>
# include "post.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
/* 
   Selectionne les faces a la frontiere d'un tag=0 et tag=1.
*/
// ============================================================================
PyObject* K_POST::frontFaces(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* tagArray;
  if (!PyArg_ParseTuple(args, "OO", &array, &tagArray)) return NULL;
  
  // Extraction array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;
  E_Int res;
  res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, 
                              cn, eltType, true);

  PyObject* tpl = NULL;
  if (res == 2)
  {
    // Extraction tag
    char* varStringt; char* eltTypet;
    FldArrayF* tag; FldArrayI* cnt;
    E_Int nit, njt, nkt;
    E_Int rest;
    rest = K_ARRAY::getFromArray(tagArray, varStringt, tag, nit, njt, nkt, 
                                 cnt, eltTypet, true);
    if (rest != 1 && rest != 2)
    {
      RELEASESHAREDU(array, f, cn); 
      PyErr_SetString(PyExc_TypeError,
                      "frontFaces: tag array is invalid.");
      return NULL;
    }
    if (tag->getSize() != f->getSize())
    {
      RELEASESHAREDB(rest, tagArray, tag, cnt); 
      RELEASESHAREDU(array, f, cn); 
      PyErr_SetString(PyExc_TypeError,
                      "frontFaces: tag must be located on the same grid.");
      return NULL;
    }
    
    tpl = frontFacesUnstructured(varString, *f, *cn, eltType, *tag);
    RELEASESHAREDU(array, f, cn);
    RELEASESHAREDB(rest, tagArray, tag, cnt); 
    return tpl;
  }
  else if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "frontFaces: only for unstructured arays.");
    RELEASESHAREDS(array, f); 
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "frontFaces: invalid array.");
    return NULL;
  }
}

//=============================================================================
// frontFaces pour non-structure
//=============================================================================
PyObject* K_POST::frontFacesUnstructured(char* varString, FldArrayF& f, 
                                         FldArrayI& cn, char* eltType,
                                         FldArrayF& tag)
{
  E_Int nelts = cn.getSize(); // nombre d'elements
  E_Int nvpe = cn.getNfld(); // nbre de noeuds par elements
  E_Int nfld = f.getNfld(); // nbre de champs dans f

  E_Int nfaces = 0; // nb de faces par element
  E_Int nvertex = 0; // nb de sommets par face
  char eltTypeOut[256]; // type des faces de sortie
  E_Float* tagp = tag.begin();
  if (strcmp(eltType, "TRI") == 0)
  {
    nfaces = 3; nvertex = 2;
    strcpy(eltTypeOut, "BAR");
  }
  else if (strcmp(eltType, "QUAD") == 0) 
  {
    nfaces = 4; nvertex = 2;
    strcpy(eltTypeOut, "BAR");
  }
  else if (strcmp(eltType, "TETRA") == 0) 
  {
    nfaces = 4; nvertex = 3;
    strcpy(eltTypeOut, "TRI");
  }
  else if (strcmp(eltType, "HEXA") == 0) 
  {
    nfaces = 6; nvertex = 4;
    strcpy(eltTypeOut, "QUAD");
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "frontFaces: not implemented for this type of elements.");
    return NULL;
  }

  // Sortie
  E_Int tzero, tone, ind, indv;
  E_Int size0 = 0;
  for (E_Int et = 0; et < nelts; et++)
  {
    // Compte les tags=0
    tzero = 0;
    for (E_Int nv = 1; nv <= nvpe; nv++)
    {
      ind = cn(et,nv)-1; if (tagp[ind] == 0) tzero++;
    }
    if ( tzero != 0 && tzero != nvpe) size0 += 1;
  }
  size0 *= nfaces;
  FldArrayI* connect = new FldArrayI(size0, nvertex);
  FldArrayI& connectp = *connect;
  FldArrayF* faces = new FldArrayF(nvertex*size0, nfld);
  FldArrayF& facesp = *faces;
  E_Int cfaces = 0; E_Int celt = 0;
  FldArrayI fi(nfaces, nvertex);

  for (E_Int et = 0; et < nelts; et++)
  {
    // Compte les tags=0
    tzero = 0;
    for (E_Int nv = 1; nv <= nvpe; nv++)
    {
      ind = cn(et, nv)-1; if (tagp[ind] == 0) tzero++;
    }
    if (tzero != 0 && tzero != nvpe)
    {
      buildFaceInfo(et, cn, fi);
      for (E_Int nf = 0; nf < nfaces; nf++) // pour chaque face
      {
        tone = 0;
        for (E_Int n = 1; n <= nvertex; n++)
        {
          indv = fi(nf, n)-1; if (tagp[indv] == 1) tone++;
        }
        if (tone == nvertex)
        {
          // On ajoute cette face
          for (E_Int n = 0; n < nvertex; n++)
          {
            indv = fi(nf, n+1)-1;
            for (E_Int nv = 1; nv <= nfld; nv++)
              facesp(cfaces+n, nv) = f(indv, nv);
          }
          for (E_Int nv = 1; nv <= nvertex; nv++)
            connectp(celt, nv) = cfaces+nv;
          cfaces = cfaces + nvertex; celt++;
        }
      }
    }
  }
  faces->reAllocMat(cfaces, nfld);
  connect->reAllocMat(celt, nvertex);
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString)+1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString)+1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString)+1;
  if (posx != 0 && posy != 0 && posz != 0)
    K_CONNECT::cleanConnectivity(posx, posy, posz, 
                                 1.e-12, eltTypeOut, 
                                 *faces, *connect);
  PyObject* tpl = K_ARRAY::buildArray(*faces, varString, 
                                      *connect, -1, eltTypeOut);
  delete faces; delete connect;
  return tpl;
}
