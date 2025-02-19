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

// Convert an unstructured array to NGON array

# include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Convert an unstructured array to a NGON array */
// ============================================================================
PyObject* K_CONVERTER::convertUnstruct2NGon(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int api = 1;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &api)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString, 
                               f, ni, nj, nk, cnl, eltType);

  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertUnstruct2NGon: array is invalid.");
    return NULL;
  }
 
  if (strcmp(eltType, "NGON") == 0)
  { RELEASESHAREDU(array, f, cnl); return array; }

  // Acces universel sur BE/ME
  E_Int nc = cnl->getNConnect();
  E_Int elOffset = 0, fcOffset = 0; // element and face offsets for ME
  E_Int sizeFNOffset = 0, sizeEFOffset = 0; // connectivity offsets
  // Acces universel aux eltTypes
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  // Number of elements, faces, and vertices per connectivity
  std::vector<E_Int> nepc(nc), nfpc(nc), nvpc(nc);
  // Number of faces per element and vertices per face for each connectivity
  std::vector<E_Int> nf(nc), nv(nc);
  // Total number of elements and faces for all connectivities
  E_Int ntotelts = 0, ntotfaces = 0;
  E_Int sizeFN = 0, sizeEF = 0;

  E_Int shift = 1; if (api == 3) shift = 0;

  // Boucle sur toutes les connectivites une premiere fois pour savoir si elles
  // sont valides et calcule les quantites totales (elements, faces, etc)
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cnl->getConnect(ic));
    char* eltTypConn = eltTypes[ic];
    // Check that this connectivity is valid
    if (not(strcmp(eltTypConn,"BAR")==0 || strcmp(eltTypConn,"TRI")==0 || 
            strcmp(eltTypConn,"QUAD")==0 || strcmp(eltTypConn,"TETRA")==0 || 
            strcmp(eltTypConn,"HEXA")==0 || strcmp(eltTypConn,"PENTA")==0 ||
            strcmp(eltTypConn,"PYRA")==0))
    {
      PyErr_SetString(PyExc_TypeError, 
                    "convertUnstruct2NGon: elt type of array (%d) is invalid.");
      for (size_t v = 0; v < eltTypes.size(); v++) delete [] eltTypes[v];
      RELEASESHAREDU(array, f, cnl); return NULL;  
    }

    nepc[ic] = cm.getSize();

    if (strcmp(eltTypConn, "BAR") == 0) // peut avoir des T-Branches
    {
      // 1 sommet par face, 2 faces par elt
      nv[ic] = 1; nf[ic] = 2;
      nfpc[ic] = nf[ic]*nepc[ic];
      nvpc[ic] = nv[ic]*nfpc[ic];
    }
    else if (strcmp(eltTypConn, "TRI") == 0)
    {
      // 2 sommets par face, 3 faces par elt
      nv[ic] = 2; nf[ic] = 3;
      nfpc[ic] = nf[ic]*nepc[ic];
      nvpc[ic] = nv[ic]*nfpc[ic];
    }
    else if (strcmp(eltTypConn, "QUAD") == 0)
    {
      // 2 sommets par face, 4 faces par elt
      nv[ic] = 2; nf[ic] = 4;
      nfpc[ic] = nf[ic]*nepc[ic];
      nvpc[ic] = nv[ic]*nfpc[ic];
    }
    else if (strcmp(eltTypConn, "TETRA") == 0)
    {
      // 3 sommets par face, 4 faces par elt
      nv[ic] = 3; nf[ic] = 4;
      nfpc[ic] = nf[ic]*nepc[ic];
      nvpc[ic] = nv[ic]*nfpc[ic];
    }
    else if (strcmp(eltTypConn, "HEXA") == 0)
    {
      // 4 sommets par face, 6 faces par elt
      nv[ic] = 4; nf[ic] = 6;
      nfpc[ic] = nf[ic]*nepc[ic];
      nvpc[ic] = nv[ic]*nfpc[ic];
    }
    else if (strcmp(eltTypConn, "PENTA") == 0)
    {
      // 3 quad et 2 tri par elt, 5 faces par elt
      nv[ic] = -1; nf[ic] = 5;
      nfpc[ic] = nf[ic]*nepc[ic];
      nvpc[ic] = 18*nepc[ic];
    }
    else if (strcmp(eltTypConn, "PYRA") == 0)
    {
      // 1 quad et 4 tri par elt, 5 faces par elt
      nv[ic] = -1; nf[ic] = 5;
      nfpc[ic] = nf[ic]*nepc[ic]; 
      nvpc[ic] = 16*nepc[ic];
    }

    ntotfaces += nfpc[ic];
    ntotelts += nepc[ic];
    sizeFN += nvpc[ic] + shift*nfpc[ic];
    sizeEF += nfpc[ic] + shift*nepc[ic];
  }

  // Build an empty NGON connectivity
  E_Int npts = f->getSize(); E_Int nfld = f->getNfld();
  E_Int ngonType = 1;
  if (api == 1) ngonType = 1; // CGNSv3 compact array1
  else if (api == 3) ngonType = 3; // force CGNSv4, array3
  else if (api == 2) ngonType = 2; // CGNSv3, array2
  E_Boolean center = false;

  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, ntotelts,
                                       ntotfaces, "NGON", sizeFN, sizeEF,
                                       ngonType, center, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

  // Acces non universel sur les ptrs NGON
  E_Int* ngon2 = cn2->getNGon();
  E_Int* nface2 = cn2->getNFace();
  E_Int *indPG2 = NULL, *indPH2 = NULL; 
  if (api == 2 || api == 3) // array2 ou array3
  {
    indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
  }
  
  // Loop over all ME connectivities to fill the new NGON connectivity
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cnl->getConnect(ic));
    char* eltTypConn = eltTypes[ic];
    
#pragma omp parallel
    {
      E_Int c1, c2, nof;
      // Connectivites FN & EF
      if (strcmp(eltTypConn, "BAR") == 0) 
      {
        E_Int v1, v2;
#pragma omp for
        for (E_Int et = 0; et < nepc[ic]; et++)
        {
          v1 = cm(et,1); v2 = cm(et,2); 
          // connectivite face/noeuds
          c1 = sizeFNOffset + (2+2*shift)*et;
          ngon2[c1] = 1; ngon2[c1+shift] = v1; // face 1
          ngon2[c1+shift+1] = 1; ngon2[c1+2*shift+1] = v2; // face 2
          // connectivite elt/faces
          c2 = sizeEFOffset + (shift+2)*et; nof = shift + fcOffset + nf[ic]*et;
          nface2[c2] = 2; nface2[c2+shift] = nof+1-shift; nface2[c2+shift+1] = nof+2-shift;
        }
      }
      else if (strcmp(eltTypConn, "TRI") == 0) 
      {
        E_Int v1, v2, v3;
#pragma omp for
        for (E_Int et = 0; et < nepc[ic]; et++)
        {
          v1 = cm(et,1); v2 = cm(et,2); v3 = cm(et,3);
          // connectivite face/noeuds
          c1 = sizeFNOffset + (3*shift+6)*et;
          ngon2[c1] = 2; ngon2[c1+shift] = v1; ngon2[c1+shift+1] = v2; // face 1
          ngon2[c1+shift+2] = 2; ngon2[c1+2*shift+2] = v2; ngon2[c1+2*shift+3] = v3; // face 2
          ngon2[c1+2*shift+4] = 2; ngon2[c1+3*shift+4] = v3; ngon2[c1+3*shift+5] = v1; // face 3
          // connectivite elt/faces
          c2 = sizeEFOffset + (shift+3)*et; nof = shift + fcOffset + nf[ic]*et;
          nface2[c2] = 3; nface2[c2+shift] = nof+1-shift; nface2[c2+shift+1] = nof+2-shift; nface2[c2+shift+2] = nof+3-shift;   
        }
      }
      else if (strcmp(eltTypConn, "QUAD") == 0) 
      {
        E_Int v1, v2, v3, v4;
#pragma omp for
        for (E_Int et = 0; et < nepc[ic]; et++)
        {
          v1 = cm(et,1); v2 = cm(et,2); v3 = cm(et,3); v4 = cm(et,4);
          // connectivite face/noeuds
          c1 = sizeFNOffset + (4*shift+8)*et;
          ngon2[c1] = 2; ngon2[c1+shift] = v1; ngon2[c1+shift+1] = v2; // face 1
          ngon2[c1+shift+2] = 2; ngon2[c1+2*shift+2] = v2; ngon2[c1+2*shift+3] = v3; // face 2
          ngon2[c1+2*shift+4] = 2; ngon2[c1+3*shift+4] = v3; ngon2[c1+3*shift+5] = v4; // face 3
          ngon2[c1+3*shift+6] = 2; ngon2[c1+4*shift+6] = v4; ngon2[c1+4*shift+7] = v1; // face 4
          // connectivite elt/faces
          c2 = sizeEFOffset + (shift+4)*et; nof = shift + fcOffset + nf[ic]*et;
          nface2[c2] = 4; nface2[c2+shift] = nof+1-shift; nface2[c2+shift+1] = nof+2-shift;
          nface2[c2+shift+2] = nof+3-shift; nface2[c2+shift+3] = nof+4-shift;      
        }
      }
      else if (strcmp(eltTypConn, "TETRA") == 0) 
      {
        E_Int v1, v2, v3, v4;
#pragma omp for
        for (E_Int et = 0; et < nepc[ic]; et++)
        {
          v1 = cm(et,1); v2 = cm(et,2); v3 = cm(et,3); v4 = cm(et,4);
          // connectivite face/noeuds
          c1 = sizeFNOffset + (4*shift+12)*et;  
          ngon2[c1] = 3; ngon2[c1+shift] = v1; ngon2[c1+shift+1] = v2; ngon2[c1+shift+2] = v3;// face 1
          ngon2[c1+shift+3] = 3; ngon2[c1+2*shift+3] = v1; ngon2[c1+2*shift+4] = v2; ngon2[c1+2*shift+5] = v4;// face 2
          ngon2[c1+2*shift+6] = 3; ngon2[c1+3*shift+6] = v2; ngon2[c1+3*shift+7] = v3; ngon2[c1+3*shift+8] = v4;// face 3
          ngon2[c1+3*shift+9] = 3; ngon2[c1+4*shift+9] = v3; ngon2[c1+4*shift+10] = v1; ngon2[c1+4*shift+11] = v4;// face 4
          // connectivite elt/faces
          c2 = sizeEFOffset + (shift+4)*et; nof = shift + fcOffset + nf[ic]*et;
          nface2[c2] = 4; nface2[c2+shift] = nof+1-shift; nface2[c2+shift+1] = nof+2-shift;
          nface2[c2+shift+2] = nof+3-shift; nface2[c2+shift+3] = nof+4-shift;      
        }
      }
      else if (strcmp(eltTypConn, "HEXA") == 0) 
      {
        E_Int v1, v2, v3, v4, v5, v6, v7, v8;
#pragma omp for
        for (E_Int et = 0; et < nepc[ic]; et++)
        {
          v1 = cm(et,1); v2 = cm(et,2); v3 = cm(et,3); v4 = cm(et,4);
          v5 = cm(et,5); v6 = cm(et,6); v7 = cm(et,7); v8 = cm(et,8);
          // connectivite face/noeuds
          c1 = sizeFNOffset + (6*shift+24)*et; 
          ngon2[c1] = 4; ngon2[c1+shift] = v1; ngon2[c1+shift+1] = v2; ngon2[c1+shift+2] = v3; ngon2[c1+shift+3] = v4;// face 1
          ngon2[c1+shift+4] = 4; ngon2[c1+2*shift+4] = v5; ngon2[c1+2*shift+5] = v6; ngon2[c1+2*shift+6] = v7; ngon2[c1+2*shift+7] = v8;// face 2
          ngon2[c1+2*shift+8] = 4; ngon2[c1+3*shift+8] = v1; ngon2[c1+3*shift+9] = v2; ngon2[c1+3*shift+10] = v6; ngon2[c1+3*shift+11] = v5;// face 3
          ngon2[c1+3*shift+12] = 4; ngon2[c1+4*shift+12] = v4; ngon2[c1+4*shift+13] = v3; ngon2[c1+4*shift+14] = v7; ngon2[c1+4*shift+15] = v8;// face 4
          ngon2[c1+4*shift+16] = 4; ngon2[c1+5*shift+16] = v1; ngon2[c1+5*shift+17] = v4; ngon2[c1+5*shift+18] = v8; ngon2[c1+5*shift+19] = v5;// face 5
          ngon2[c1+5*shift+20] = 4; ngon2[c1+6*shift+20] = v2; ngon2[c1+6*shift+21] = v3; ngon2[c1+6*shift+22] = v7; ngon2[c1+6*shift+23] = v6;// face 6
          // connectivite elt/faces
          c2 = sizeEFOffset + (shift+6)*et; nof = shift + fcOffset + nf[ic]*et;
          nface2[c2] = 6; nface2[c2+shift] = nof+1-shift; nface2[c2+shift+1] = nof+2-shift; nface2[c2+shift+2] = nof+3-shift;
          nface2[c2+shift+3] = nof+4-shift; nface2[c2+shift+4] = nof+5-shift; nface2[c2+shift+5] = nof+6-shift;
        }
      }
      else if (strcmp(eltTypConn, "PENTA") == 0) 
      {
        E_Int v1, v2, v3, v4, v5, v6;
#pragma omp for
        for (E_Int et = 0; et < nepc[ic]; et++)
        {
          v1 = cm(et,1); v2 = cm(et,2); v3 = cm(et,3);
          v4 = cm(et,4); v5 = cm(et,5); v6 = cm(et,6);
          // connectivite face/noeuds
          c1 = sizeFNOffset + (5*shift+18)*et;
          ngon2[c1] = 3; ngon2[c1+shift] = v1; ngon2[c1+shift+1] = v2; ngon2[c1+shift+2] = v3;// face 1 TRI
          ngon2[c1+shift+3] = 3; ngon2[c1+2*shift+3] = v4; ngon2[c1+2*shift+4] = v5; ngon2[c1+2*shift+5] = v6;// face 2 TRI
          ngon2[c1+2*shift+6] = 4; ngon2[c1+3*shift+6] = v1; ngon2[c1+3*shift+7] = v2; ngon2[c1+3*shift+8] = v5; ngon2[c1+3*shift+9] = v4;// face 3 : QUAD
          ngon2[c1+3*shift+10] = 4; ngon2[c1+4*shift+10] = v2; ngon2[c1+4*shift+11] = v3; ngon2[c1+4*shift+12]= v6; ngon2[c1+4*shift+13] = v5;// face 4 : QUAD
          ngon2[c1+4*shift+14] = 4; ngon2[c1+5*shift+14] = v3; ngon2[c1+5*shift+15] = v1; ngon2[c1+5*shift+16] = v4; ngon2[c1+5*shift+17] = v6;// face 5 : QUAD
          // connectivite elt/faces
          c2 = sizeEFOffset + (shift+5)*et; nof = shift + fcOffset + nf[ic]*et;
          nface2[c2] = 5; nface2[c2+shift] = nof+1-shift; nface2[c2+shift+1] = nof+2-shift; nface2[c2+shift+2] = nof+3-shift;
          nface2[c2+shift+3] = nof+4-shift; nface2[c2+shift+4] = nof+5-shift;      
        }
      }
      else if (strcmp(eltTypConn, "PYRA") == 0) 
      {
        E_Int v1, v2, v3, v4, v5;
#pragma omp for
        for (E_Int et = 0; et < nepc[ic]; et++)
        {
          v1 = cm(et,1); v2 = cm(et,2); v3 = cm(et,3);
          v4 = cm(et,4); v5 = cm(et,5);
          // connectivite face/noeuds
          c1 = sizeFNOffset + (5*shift+16)*et;
          ngon2[c1] = 4; ngon2[c1+shift] = v1; ngon2[c1+shift+1] = v2; ngon2[c1+shift+2] = v3; ngon2[c1+shift+3] = v4; // face 1: QUAD
          ngon2[c1+shift+4] = 3; ngon2[c1+2*shift+4] = v1; ngon2[c1+2*shift+5] = v2; ngon2[c1+2*shift+6] = v5;// face 2: TRI
          ngon2[c1+2*shift+7] = 3; ngon2[c1+3*shift+7] = v2; ngon2[c1+3*shift+8] = v3; ngon2[c1+3*shift+9] = v5;// face 3: TRI
          ngon2[c1+3*shift+10] = 3; ngon2[c1+4*shift+10] = v3; ngon2[c1+4*shift+11] = v4; ngon2[c1+4*shift+12] = v5;// face 4: TRI
          ngon2[c1+4*shift+13] = 3; ngon2[c1+5*shift+13] = v4; ngon2[c1+5*shift+14] = v1; ngon2[c1+5*shift+15] = v5;// face 5: TRI
          // connectivite elt/faces
          c2 = sizeEFOffset + (shift+5)*et; nof = shift + fcOffset + nf[ic]*et;
          nface2[c2] = 5; nface2[c2+shift] = nof+1-shift; nface2[c2+shift+1] = nof+2-shift; nface2[c2+shift+2] = nof+3-shift;
          nface2[c2+shift+3] = nof+4-shift; nface2[c2+shift+4] = nof+5-shift;      
        }
      }

      // Start offset indices
      if (api == 2 || api == 3) // array2 ou array3
      {
        E_Int c = 0, d = 0;
        if (strcmp(eltTypConn, "PENTA") == 0)
        {
#pragma omp for
          for (E_Int i = 0; i < nepc[ic]; i++)
          {
            c = fcOffset + i*5;
            d = sizeFNOffset + i*(18+5*shift);
            indPG2[c] = d;
            indPG2[c+1] = d + (3+shift);
            indPG2[c+2] = d + (6+2*shift);
            indPG2[c+3] = d + (10+3*shift);
            indPG2[c+4] = d + (14+4*shift);
          }
        }
        else if (strcmp(eltTypConn, "PYRA") == 0)
        {
#pragma omp for
          for (E_Int i = 0; i < nepc[ic]; i++)
          {
            c = fcOffset + i*5;
            d = sizeFNOffset + i*(16+5*shift);
            indPG2[c] = d;
            indPG2[c+1] = d + (4+shift);
            indPG2[c+2] = d + (7+2*shift);
            indPG2[c+3] = d + (10+3*shift);
            indPG2[c+4] = d + (13+4*shift);
          }
        }
        else
        {
#pragma omp for
          for (E_Int i = 0; i < nfpc[ic]; i++)
          {
            c = fcOffset + i;
            indPG2[c] = sizeFNOffset + (nv[ic]+shift)*i;
          }
        }

#pragma omp for
        for (E_Int i = 0; i < nepc[ic]; i++)
        {
          c = elOffset + i;
          indPH2[c] = sizeEFOffset + (nf[ic]+shift)*i;
        } 
      }
    }

    // Increment offsets
    fcOffset += nfpc[ic];
    elOffset += nepc[ic];
    sizeFNOffset += nvpc[ic] + shift*nfpc[ic];
    sizeEFOffset += nfpc[ic] + shift*nepc[ic];
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

#pragma omp parallel
  {
    // Copy fields to f2
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* f2p = f2->begin(n);
#pragma omp for
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }

  RELEASESHAREDU(array, f, cnl);

  // Clean connectivity
  E_Float tol = 1.e-12;
  E_Bool rmOverlappingPts=true; E_Bool rmOrphanPts=false;
  E_Bool rmDuplicatedFaces=true; E_Bool rmDuplicatedElts=false;
  E_Bool rmDegeneratedFaces=false; E_Bool rmDegeneratedElts=false;
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
