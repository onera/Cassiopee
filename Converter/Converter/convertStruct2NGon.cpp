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

// Convert structured array to NGON array

# include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Convert a structured array to a NGON array */
// ============================================================================
PyObject* K_CONVERTER::convertStruct2NGon(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, ni, nj, nk, cnl, eltType, true);

  if (res == 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2NGon: array must be structured.");
    RELEASESHAREDU(array, f, cnl); return NULL; 
  }
  else if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2NGon: array is invalid.");
    return NULL;
  }
  // Check ni, nj, nk
  if (ni < 1 || nj < 1 || nk < 1)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "convertStruct2NGon: ni, nj, nk must be >=1.");
    RELEASESHAREDS(array, f); return NULL;
  }
  // 1D, 2D or 3D ?
  E_Int dim0 = 3;
  if (ni == 1)
  {
    if (nj == 1 || nk == 1) dim0 = 1;
    else dim0 = 2;
  }
  else if (nj == 1)
  {
    if (ni == 1 || nk == 1) dim0 = 1;
    else dim0 = 2;
  }
  else if (nk == 1)
  {
    if (ni == 1 || nj == 1) dim0 = 1;
    else dim0 = 2;
  }
  
  E_Int nfld = f->getNfld();
  E_Int ninj = ni*nj; E_Int npts = ninj*nk;
  E_Int ni1 = E_max(1, E_Int(ni)-1);
  E_Int nj1 = E_max(1, E_Int(nj)-1);
  E_Int nk1 = E_max(1, E_Int(nk)-1);
  E_Int ni1nj1 = ni1*nj1;
  E_Int ncells = ni1nj1*nk1; // nb de cellules structurees
  E_Int nfaces = 0; E_Int sizeFN = 0; E_Int sizeEF = 0;
  if (dim0 == 1) 
  {nfaces = npts; sizeFN = nfaces*(1+1); sizeEF = ncells*(2+1);}
  else if (dim0 == 2) 
  {
    if (ni == 1) nfaces = nj*nk1+nj1*nk;
    else if (nj == 1) nfaces = ni*nk1+ni1*nk;
    else nfaces = ni*nj1+ni1*nj;
    sizeFN = (2+1)*nfaces; sizeEF = (4+1)*ncells;
  }
  else 
  { 
    nfaces = ni*nj1*nk1 + ni1*nj*nk1+ ni1nj1*nk;
    sizeFN = (4+1)*nfaces; sizeEF = (6+1)*ncells;
  }
  E_Int csize = sizeFN + sizeEF + 4;
  //PyObject* tpl = K_ARRAY::buildArray(nfld, varString, npts, ncells, -1, 
  //                                    "NGON", false, csize);
  //E_Float* fieldp = K_ARRAY::getFieldPtr(tpl);
  //FldArrayF field(npts, nfld, fieldp, true); field = *f;
  //E_Int* cnp = K_ARRAY::getConnectPtr(tpl);
  //FldArrayI cn(csize, 1, cnp, true);
  FldArrayF field(npts, nfld); field = *f;
  FldArrayI cn(csize, 1);
  E_Int* cnp = cn.begin(1);

  // Build the NGon connectivity
  cnp[0] = nfaces;
  cnp[1] = sizeFN;
  cnp[sizeFN+2] = ncells;
  cnp[sizeFN+3] = sizeEF;
  E_Int* cFN = cnp+2;
  E_Int* cEF = cnp+sizeFN+4;
  // Connectivite FN
  E_Int ind1, ind2, ind3, ind4, ind5, ind6;
  E_Int ninti = ni*nj1*nk1;
  E_Int nintj = ni1*nj*nk1;

  E_Int c = 0;
  E_Int nidim;
  switch (dim0)
  {
    case 1:
      // connectivite FN
      if (nk == 1 && nj == 1) nidim = ni;
      else if (ni == 1 && nk == 1) nidim = nj;
      else nidim = nk; 
      for (E_Int i = 0; i < nidim; i++)
      {
        cFN[c] = 1; // 1 noeud par face
        cFN[c+1] = i+1;
        c += 2;
      }
      // connectivite EF
      c = 0;
      nidim = E_max(nidim-1,1);
      for (E_Int i = 0; i < nidim; i++)
      {
        cEF[c] = 2;
        cEF[c+1] = i+1;
        cEF[c+2] = i+2;
        c += 3;
      }
      break;
    
    case 2:
      if (nk == 1)
      {
        // Faces en i
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni; i++)
          {
            cFN[c] = 2;
            ind1 = i+j*ni; cFN[c+1] = ind1+1;
            ind2 = ind1+ni; cFN[c+2] = ind2+1;
            c += 3;
          }
        // Faces en j
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cFN[c] = 2;
            ind1 = i+j*ni; cFN[c+1] = ind1+1;
            ind2 = ind1+1; cFN[c+2] = ind2+1;
            c += 3;
          }
        // Connectivite EF
        c = 0;
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cEF[c] = 4;
            ind1 = i+j*ni; // faces en i
            ind2 = ind1+1;
            ind3 = ninti+i+j*ni1; // faces en j
            ind4 = ind3+ni1;
            cEF[c+1] = ind1+1;
            cEF[c+2] = ind3+1;
            cEF[c+3] = ind2+1;
            cEF[c+4] = ind4+1;
            c += 5;
          }
      }// fin nk = 1
      else if (nj == 1) 
      {
        // Faces en i
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni; i++)
          {
            cFN[c] = 2;
            ind1 = i+k*ni; cFN[c+1] = ind1+1;
            ind2 = ind1+ni; cFN[c+2] = ind2+1;
            c += 3;
          }
        // Faces en k
        for (E_Int k = 0; k < nk; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cFN[c] = 2;
            ind1 = i+k*ni; cFN[c+1] = ind1+1;
            ind2 = ind1+1; cFN[c+2] = ind2+1;
            c += 3;
          }
        // Connectivite EF
        c = 0;
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cEF[c] = 4;
            ind1 = i+k*ni; // faces en i
            ind2 = ind1+1;
            ind3 = ninti+i+k*ni1; // faces en k
            ind4 = ind3+ni1;

            cEF[c+1] = ind1+1;
            cEF[c+2] = ind3+1;
            cEF[c+3] = ind2+1;
            cEF[c+4] = ind4+1;
            c += 5;
          }
      }
      else // ni = 1
      {
        // Faces en j
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int j = 0; j < nj; j++)
          {
            cFN[c] = 2;
            ind1 = j+k*nj; cFN[c+1] = ind1+1;
            ind2 = ind1+nj; cFN[c+2] = ind2+1;
            c += 3;
          }
  
      // Faces en k
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj1; j++)
          {
            cFN[c] = 2;
            ind1 = j+k*nj; cFN[c+1] = ind1+1;
            ind2 = ind1+1; cFN[c+2] = ind2+1;
            c += 3;
          }

      // Connectivite EF
      c = 0;
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          {
            cEF[c] = 4;
            ind1 = j+k*nj; // faces en j
            ind2 = ind1+1;
            ind3 = nintj+j+k*nj1; // faces en k
            ind4 = ind3+nj1;
            cEF[c+1] = ind1+1;
            cEF[c+2] = ind3+1;
            cEF[c+3] = ind2+1;
            cEF[c+4] = ind4+1;
            c += 5;
          }
      }
      break;
    
    default:
      // Faces en i
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni; i++)
          {
            cFN[c] = 4;
            ind1 = i+j*ni+k*ninj;
            ind2 = ind1+ni;
            ind3 = ind2+ninj;
            ind4 = ind1+ninj;
            cFN[c+1] = ind1+1;
            cFN[c+2] = ind2+1;
            cFN[c+3] = ind3+1;
            cFN[c+4] = ind4+1;
            c += 5;
          }
  
      // Faces en j
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cFN[c] = 4;
            ind1 = i+j*ni+k*ninj;
            ind2 = ind1+1;
            ind3 = ind2+ninj;
            ind4 = ind1+ninj;
            cFN[c+1] = ind1+1;
            cFN[c+2] = ind2+1;
            cFN[c+3] = ind3+1;
            cFN[c+4] = ind4+1;
            c += 5;
          }
  
      // Faces en k
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cFN[c] = 4;
            ind1 = i+j*ni+k*ninj;
            ind2 = ind1+1;
            ind3 = ind2+ni;
            ind4 = ind1+ni;
            cFN[c+1] = ind1+1;
            cFN[c+2] = ind2+1;
            cFN[c+3] = ind3+1;
            cFN[c+4] = ind4+1;
            c += 5;
          }

      // Connectivite EF
      c = 0;
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cEF[c] = 6;
            ind1 = i+j*ni+k*ni*nj1; // faces en i
            ind2 = ind1+1;
            ind3 = ninti+i+j*ni1+k*ni1*nj; // faces en j
            ind4 = ind3+ni1;
            ind5 = ninti+nintj+i+j*ni1+k*ni1nj1; // faces en k
            ind6 = ind5+ni1nj1;
            cEF[c+1] = ind1+1;
            cEF[c+2] = ind2+1;
            cEF[c+3] = ind3+1;
            cEF[c+4] = ind4+1;
            cEF[c+5] = ind5+1;
            cEF[c+6] = ind6+1;
            c += 7;
          }
      break;
  }
  RELEASESHAREDS(array, f);

  /* clean connectivity */
  /* Est-il necessaire? */
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString)+1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString)+1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString)+1;
  E_Float tol = 1.e-12;
  if (posx > 0 && posy > 0 && posz > 0)
    K_CONNECT::cleanConnectivityNGon(posx, posy, posz, tol, field, cn);

  PyObject* tpl = K_ARRAY::buildArray(field, varString, 
                                      cn, -1, "NGON");
  return tpl;
}
