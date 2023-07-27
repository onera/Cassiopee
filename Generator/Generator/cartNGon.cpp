/*    
    Copyright 2013-2023 Onera.

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

// Generateur non structure NGON a partir d'une grille cartesienne

# include "generator.h"

using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
/* Generateur de grille non structuree NGon
   IN: x0, y0, z0: origine de la grille
   IN: hi, hj, hk: pas de la grille 
   IN: ni, nj, nk: nombre de points
   OUT: array definissant le maillage cree. */
//=============================================================================
PyObject* K_GENERATOR::cartNGon(PyObject* self, PyObject* args)
{
  E_Int ni, nj, nk;
  E_Float xo, yo, zo;
  E_Float hi, hj, hk;
  E_Int api = 1;
  if (!PYPARSETUPLE(args, 
                    "(ddd)(ddd)(lll)l", "(ddd)(ddd)(iii)i",
                    "(fff)(fff)(lll)l", "(fff)(fff)(lll)i",
                    &xo, &yo, &zo, &hi, &hj, &hk, &ni, &nj, &nk, &api))
  {
    return NULL;
  }
  
  // Check ni, nj, nk
  if (ni < 1 || nj < 1 || nk < 1)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cartNGon: ni, nj, nk must be >=1.");
    return NULL;
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

  // Create cartesian mesh
  E_Int ind;
  E_Int ninj = ni*nj; E_Int npts = ninj*nk;
  E_Int ni1 = E_max(1, E_Int(ni)-1);
  E_Int nj1 = E_max(1, E_Int(nj)-1);
  E_Int nk1 = E_max(1, E_Int(nk)-1);  
  E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees

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
    nfaces = ni*nj1*nk1 + ni1*nj*nk1+ ni1*nj1*nk;
    sizeFN = (4+1)*nfaces; sizeEF = (6+1)*ncells;
  }

  E_Int shift = 1;
  if (api == 3) { sizeFN -= nfaces; sizeEF -= ncells; shift = 0; }

  //PyObject* tpl = K_ARRAY::buildArray2(3, "x,y,z", npts, ncells, -1, 
  //                                     "NGON", false, sizeFN, sizeEF, nfaces, api);

  PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", npts, ncells, nfaces, "NGON", 
                                       sizeFN, sizeEF, false, api);

  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  //char* varString; char* eltType;
  //K_ARRAY::getFromArray2(tpl, varString, f, ni, nj, nk, cn, eltType);
  K_ARRAY::getFromArray3(tpl, f, cn);
  
  E_Int* cFN = cn->getNGon();
  E_Int* cEF = cn->getNFace();

  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);

  E_Int nij = ni*nj;
  E_Int nijk = ni*nj*nk;
  E_Int i, j, k;

#pragma omp parallel for default(shared) private(k,j,i,ind)
  for (ind = 0; ind < nijk; ind++)
  {
    k = ind/nij;
    j = (ind-k*nij)/ni;
    i = ind-j*ni-k*nij;
    xt[ind] = xo + i * hi;
    yt[ind] = yo + j * hj;
    zt[ind] = zo + k * hk;
  }

  // Build the NGon connectivity

  // Connectivite FN
  E_Int ind1, ind2, ind3, ind4, ind5, ind6;
  E_Int ninti = ni*nj1*nk1;
  E_Int nintj = ni1*nj*nk1;
  //E_Int nintk = ni1*nj1*nk;

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
        cFN[c+shift] = i+1;
        c += 1+shift;
      }
      // connectivite EF
      c = 0;
      nidim = E_max(nidim-1,1);
      for (E_Int i = 0; i < nidim; i++)
      {
        cEF[c] = 2;
        cEF[c+shift] = i+1;
        cEF[c+shift+1] = i+2;
        c += 2+shift;
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
            ind1 = i+j*ni; cFN[c+shift] = ind1+1;
            ind2 = ind1+ni; cFN[c+shift+1] = ind2+1;
            c += 2+shift;
          }
        // Faces en j
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cFN[c] = 2;
            ind1 = i+j*ni; cFN[c+shift] = ind1+1;
            ind2 = ind1+1; cFN[c+shift+1] = ind2+1;
            c += 2+shift;
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
            cEF[c+shift] = ind1+1; // en 2D, la connectivite EF doit faire
            cEF[c+shift+1] = ind3+1; // un cycle de faces
            cEF[c+shift+2] = ind2+1;
            cEF[c+shift+3] = ind4+1;
            c += 4+shift;
          }
      }// fin nk = 1
      else if (nj == 1) 
      {
        // Faces en i
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni; i++)
          {
            cFN[c] = 2;
            ind1 = i+k*ni; cFN[c+shift] = ind1+1;
            ind2 = ind1+ni; cFN[c+shift+1] = ind2+1;
            c += 2+shift;
          }
        // Faces en k
        for (E_Int k = 0; k < nk; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cFN[c] = 2;
            ind1 = i+k*ni; cFN[c+shift] = ind1+1;
            ind2 = ind1+1; cFN[c+shift+1] = ind2+1;
            c += 2+shift;
          }
        // Connectivite EF
        c = 0;
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cEF[c] = 3+shift;
            ind1 = i+k*ni; // faces en i
            ind2 = ind1+1;
            ind3 = ninti+i+k*ni1; // faces en k
            ind4 = ind3+ni1;

            cEF[c+shift] = ind1+1;
            cEF[c+shift+1] = ind3+1;
            cEF[c+shift+2] = ind2+1;
            cEF[c+shift+3] = ind4+1;
            c += 4+shift;
          }
      }
      else // ni = 1
      {
        // Faces en j
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int j = 0; j < nj; j++)
          {
            cFN[c] = 2;
            ind1 = j+k*nj;cFN[c+shift] = ind1+1;
            ind2 = ind1+nj;cFN[c+shift+1] = ind2+1;
            c += 2+shift;
          }
  
      // Faces en k
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj1; j++)
          {
            cFN[c] = 2;
            ind1 = j+k*nj; cFN[c+shift] = ind1+1;
            ind2 = ind1+1; cFN[c+shift+1] = ind2+1;
            c += 2+shift;
          }

      // Connectivite EF
      c = 0;
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          {
            cEF[c] = 3+shift;
            ind1 = j+k*nj; // faces en j
            ind2 = ind1+1;
            ind3 = nintj+j+k*nj1; // faces en k
            ind4 = ind3+nj1;
            cEF[c+shift] = ind1+1;
            cEF[c+shift+1] = ind3+1;
            cEF[c+shift+2] = ind2+1;
            cEF[c+shift+3] = ind4+1;
            c += 4+shift;
          }
      }
      break;
    
    default:
      // Faces en i
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni; i++)
          {
            cFN[c] = 3+shift;
            ind1 = i+j*ni+k*ninj;
            ind2 = ind1+ni;
            ind3 = ind2+ninj;
            ind4 = ind1+ninj;
            cFN[c+shift] = ind1+1;
            cFN[c+shift+1] = ind2+1;
            cFN[c+shift+2] = ind3+1;
            cFN[c+shift+3] = ind4+1;
            c += 4+shift;
          }
  
      // Faces en j
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cFN[c] = 3+shift;
            ind1 = i+j*ni+k*ninj;
            ind2 = ind1+1;
            ind3 = ind2+ninj;
            ind4 = ind1+ninj;
            cFN[c+shift] = ind1+1;
            cFN[c+shift+1] = ind2+1;
            cFN[c+shift+2] = ind3+1;
            cFN[c+shift+3] = ind4+1;
            c += 4+shift;
          }
  
      // Faces en k
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cFN[c] = 3+shift;
            ind1 = i+j*ni+k*ninj;
            ind2 = ind1+1;
            ind3 = ind2+ni;
            ind4 = ind1+ni;
            cFN[c+shift] = ind1+1;
            cFN[c+shift+1] = ind2+1;
            cFN[c+shift+2] = ind3+1;
            cFN[c+shift+3] = ind4+1;
            c += 4+shift;
          }

      // Connectivite EF
      c = 0;
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            cEF[c] = 5+shift;
            ind1 = i+j*ni+k*ni*nj1; // faces en i
            ind2 = ind1+1;
            ind3 = ninti+i+j*ni1+k*ni1*nj; // faces en j
            ind4 = ind3+ni1;
            ind5 = ninti+nintj+i+j*ni1+k*ni1*nj1; // faces en k
            ind6 = ind5+ni1*nj1;
            cEF[c+shift] = ind1+1;
            cEF[c+shift+1] = ind2+1;
            cEF[c+shift+2] = ind3+1;
            cEF[c+shift+3] = ind4+1;
            cEF[c+shift+4] = ind5+1;
            cEF[c+shift+5] = ind6+1;
            c += 6+shift;
          }
      break;
  }

  // index
  if (api == 2 || api == 3) // array2 ou array3
  {
   E_Int* indPG = cn->getIndPG();
   E_Int* indPH = cn->getIndPH();
#pragma omp parallel
  {
    #pragma omp for
    for (E_Int i = 0; i < nfaces; i++) indPG[i] = (pow(2,dim0-1)+shift)*i;
    #pragma omp for
    for (E_Int i = 0; i < ncells; i++) indPH[i] = (2*dim0+shift)*i;  
  }
  if (api == 3) { indPG[nfaces] = pow(2,dim0-1)*nfaces; indPH[ncells] = 2*dim0*ncells; }
  }
  //RELEASESHAREDU(tpl, f, cn);
  delete f; delete cn;
  Py_DECREF(PyList_GetItem(tpl,1)); Py_DECREF(PyList_GetItem(tpl,2));
  return tpl;
}
