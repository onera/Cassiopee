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
  E_Int api = 1;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &api)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* dummy;
  res = K_ARRAY::getFromArray3(array, varString, 
                               f, ni, nj, nk, cnl, dummy);
  E_Int shift = 1; if (api == 3) shift = 0;

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
  {nfaces = npts; sizeFN = nfaces*(1+shift); sizeEF = ncells*(2+shift);}
  else if (dim0 == 2)
  {
    if (ni == 1) nfaces = nj*nk1+nj1*nk;
    else if (nj == 1) nfaces = ni*nk1+ni1*nk;
    else nfaces = ni*nj1+ni1*nj;
    sizeFN = (2+shift)*nfaces; sizeEF = (4+shift)*ncells;
  }
  else 
  { 
    nfaces = ni*nj1*nk1 + ni1*nj*nk1+ ni1nj1*nk;
    sizeFN = (4+shift)*nfaces; sizeEF = (6+shift)*ncells;
  }

  // Build an empty NGON array
  E_Int ngonType = 1; // CGNSv3 compact array1
  if (api == 2) ngonType = 2; // CGNSv3, array2
  else if (api == 3) ngonType = 3; // force CGNSv4, array3 
  E_Boolean center = false;
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, ncells, nfaces, 
                                       "NGON", sizeFN, sizeEF, ngonType,
                                       center, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

  // Acces non universel sur les nouveaux ptrs
  E_Int* ngon2 = cn2->getNGon();
  E_Int* nface2 = cn2->getNFace();
  E_Int *indPG2 = NULL, *indPH2 = NULL; 
  if (api == 2 || api == 3) // array2 ou array3
  {
    indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
  }
  
  E_Int ninti = ni*nj1*nk1;
  E_Int nintj = ni1*nj*nk1;

#pragma omp parallel
  {
    E_Int c;
    if (dim0 == 1)
    {
      E_Int nidim, nidim2;
      if (nk == 1 && nj == 1) nidim = ni;
      else if (ni == 1 && nk == 1) nidim = nj;
      else nidim = nk;
      nidim2 = E_max(nidim-1,1);

      // connectivite FN
#pragma omp for nowait
      for (E_Int i = 0; i < nidim; i++)
      {
        c = (1+shift)*i;
        ngon2[c] = 1;
        ngon2[c+shift] = i+1;
      }
      // connectivite EF
#pragma omp for nowait
      for (E_Int i = 0; i < nidim2; i++)
      {
        c = (2+shift)*i;
        nface2[c] = 2;
        nface2[c+shift] = i+1;
        nface2[c+shift+1] = i+2;
      }
    }
    else if (dim0 == 2)
    {
      E_Int ind1, ind2, ind3, ind4;
      if (nk == 1)
      {
        // connectivite FN
        // Faces en i
#pragma omp for nowait
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni; i++)
          {
            c = (2+shift)*(j*ni + i);
            ngon2[c] = 2;
            ind1 = i+j*ni; ngon2[c+shift] = ind1+1;
            ind2 = ind1+ni; ngon2[c+shift+1] = ind2+1;
          }
        // Faces en j
#pragma omp for nowait
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            c = (2+shift)*(nj1*ni + j*ni1 + i);
            ngon2[c] = 2;
            ind1 = i+j*ni; ngon2[c+shift] = ind1+1;
            ind2 = ind1+1; ngon2[c+shift+1] = ind2+1;
          }
        // Connectivite EF
#pragma omp for nowait
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            c = (4+shift)*(j*ni1 + i);
            ind1 = i+j*ni; // faces en i
            ind2 = ind1+1;
            ind3 = ninti+i+j*ni1; // faces en j
            ind4 = ind3+ni1;
            nface2[c] = 4;
            nface2[c+shift] = ind1+1;
            nface2[c+shift+1] = ind3+1;
            nface2[c+shift+2] = ind2+1;
            nface2[c+shift+3] = ind4+1;
          }
      }// fin nk = 1
      else if (nj == 1) 
      {
        // connectivite FN
        // Faces en i
#pragma omp for nowait
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni; i++)
          {
            c = (2+shift)*(k*ni + i);
            ngon2[c] = 2;
            ind1 = i+k*ni; ngon2[c+shift] = ind1+1;
            ind2 = ind1+ni; ngon2[c+shift+1] = ind2+1;
          }
        // Faces en k
#pragma omp for nowait
        for (E_Int k = 0; k < nk; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            c = (2+shift)*(nk1*ni + k*ni1 + i);
            ngon2[c] = 2;
            ind1 = i+k*ni; ngon2[c+shift] = ind1+1;
            ind2 = ind1+1; ngon2[c+shift+1] = ind2+1;
          }
        // Connectivite EF
#pragma omp for nowait
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            c = (4+shift)*(k*ni1 + i);
            ind1 = i+k*ni; // faces en i
            ind2 = ind1+1;
            ind3 = ninti+i+k*ni1; // faces en k
            ind4 = ind3+ni1;
            nface2[c] = 4;
            nface2[c+shift] = ind1+1;
            nface2[c+shift+1] = ind3+1;
            nface2[c+shift+2] = ind2+1;
            nface2[c+shift+3] = ind4+1;
          }
      }
      else // ni = 1
      {
        // connectivite FN
        // Faces en j
#pragma omp for nowait
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int j = 0; j < nj; j++)
          {
            c = (2+shift)*(k*nj + j);
            ngon2[c] = 2;
            ind1 = j+k*nj; ngon2[c+shift] = ind1+1;
            ind2 = ind1+nj; ngon2[c+shift+1] = ind2+1;
          }
  
        // Faces en k
#pragma omp for nowait
        for (E_Int k = 0; k < nk; k++)
          for (E_Int j = 0; j < nj1; j++)
          {
            c = (2+shift)*(nk1*nj + k*nj1 + j);
            ngon2[c] = 2;
            ind1 = j+k*nj; ngon2[c+shift] = ind1+1;
            ind2 = ind1+1; ngon2[c+shift+1] = ind2+1;
          }

        // Connectivite EF
#pragma omp for nowait
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int j = 0; j < nj1; j++)
            {
              c = (4+shift)*(k*nj1 + j);
              ind1 = j+k*nj; // faces en j
              ind2 = ind1+1;
              ind3 = nintj+j+k*nj1; // faces en k
              ind4 = ind3+nj1;
              nface2[c] = 4;
              nface2[c+shift] = ind1+1;
              nface2[c+shift+1] = ind3+1;
              nface2[c+shift+2] = ind2+1;
              nface2[c+shift+3] = ind4+1;
            }
      }
    }
    else
    {
      E_Int ind1, ind2, ind3, ind4, ind5, ind6;
      // connectivite FN
      // Faces en i
#pragma omp for nowait
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni; i++)
          {
            c = (4+shift)*(k*nj1*ni + j*ni + i);
            ind1 = i+j*ni+k*ninj;
            ind2 = ind1+ni;
            ind3 = ind2+ninj;
            ind4 = ind1+ninj;
            ngon2[c] = 4;
            ngon2[c+shift] = ind1+1;
            ngon2[c+shift+1] = ind2+1;
            ngon2[c+shift+2] = ind3+1;
            ngon2[c+shift+3] = ind4+1;
          }
  
      // Faces en j
#pragma omp for nowait
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            c = (4+shift)*(nk1*nj1*ni + k*nj*ni1 + j*ni1 + i);
            ind1 = i+j*ni+k*ninj;
            ind2 = ind1+1;
            ind3 = ind2+ninj;
            ind4 = ind1+ninj;
            ngon2[c] = 4;
            ngon2[c+shift] = ind1+1;
            ngon2[c+shift+1] = ind2+1;
            ngon2[c+shift+2] = ind3+1;
            ngon2[c+shift+3] = ind4+1;
          }
  
      // Faces en k
#pragma omp for nowait
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            c = (4+shift)*(nk1*(nj1*ni + nj*ni1) + k*nj1*ni1 + j*ni1 + i);
            ind1 = i+j*ni+k*ninj;
            ind2 = ind1+1;
            ind3 = ind2+ni;
            ind4 = ind1+ni;
            ngon2[c] = 4;
            ngon2[c+shift] = ind1+1;
            ngon2[c+shift+1] = ind2+1;
            ngon2[c+shift+2] = ind3+1;
            ngon2[c+shift+3] = ind4+1;
          }

      // Connectivite EF
#pragma omp for nowait
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            c = (6+shift)*(k*nj1*ni1 + j*ni1 + i);
            ind1 = i+j*ni+k*ni*nj1; // faces en i
            ind2 = ind1+1;
            ind3 = ninti+i+j*ni1+k*ni1*nj; // faces en j
            ind4 = ind3+ni1;
            ind5 = ninti+nintj+i+j*ni1+k*ni1nj1; // faces en k
            ind6 = ind5+ni1nj1;
            nface2[c] = 6;
            nface2[c+shift] = ind1+1;
            nface2[c+shift+1] = ind2+1;
            nface2[c+shift+2] = ind3+1;
            nface2[c+shift+3] = ind4+1;
            nface2[c+shift+4] = ind5+1;
            nface2[c+shift+5] = ind6+1;
          }
    }

    // Start offset indices
    if (api == 2 || api == 3) // array2 ou array3
    {
#pragma omp for nowait
      for (E_Int i = 0; i < nfaces; i++) indPG2[i] = (pow(2,dim0-1)+shift)*i;
#pragma omp for nowait
      for (E_Int i = 0; i < ncells; i++) indPH2[i] = (2*dim0+shift)*i;
    }

    // Copy fields to f2
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* f2p = f2->begin(n);
#pragma omp for
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }

  RELEASESHAREDS(array, f);

  // Clean connectivity
  E_Float tol = 1.e-12;
  E_Bool rmOverlappingPts=true; E_Bool rmOrphanPts=false;
  E_Bool rmDuplicatedFaces=true; E_Bool rmDuplicatedElts=false;
  E_Bool rmDegeneratedFaces=true; E_Bool rmDegeneratedElts=false;
  PyObject* tplClean = K_CONNECT::V_cleanConnectivity(
      varString, *f2, *cn2, "NGON", tol,
      rmOverlappingPts, rmOrphanPts,
      rmDuplicatedFaces, rmDuplicatedElts,
      rmDegeneratedFaces, rmDegeneratedElts
  );
  
  RELEASESHAREDU(tpl, f2, cn2);
  if (tplClean == NULL) return tpl;
  else { Py_DECREF(tpl); return tplClean; }
}
