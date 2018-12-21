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
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, ni, nj, nk, cnl, eltType, true);

  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertUnstruct2NGon: array is invalid.");
    return NULL;
  }
 
  if (strcmp(eltType, "NGON") == 0)
  { RELEASESHAREDU(array, f, cnl); return array; }

  if (strcmp(eltType, "TRI")   != 0 && strcmp(eltType, "QUAD") != 0 &&  
      strcmp(eltType, "TETRA") != 0 && strcmp(eltType, "HEXA") != 0 &&
      strcmp(eltType, "PENTA") != 0 && strcmp(eltType, "BAR")  != 0 && 
      strcmp(eltType, "PYRA")  != 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertUnstruct2NGon: elt type of array (%d) is invalid.");
    RELEASESHAREDU(array, f, cnl); return NULL;    
  }
  E_Int npts = f->getSize(); E_Int nfld = f->getNfld();
  E_Int nelts = cnl->getSize();
  E_Int sizeFN = 0; E_Int sizeEF = 0;
  E_Int nfaces = 0;
  if (strcmp(eltType, "BAR") == 0) // peut avoir des T-Branches
  {
    nfaces = 2*nelts;  
    sizeFN = nfaces*(1+1);// 1 sommet par face
    sizeEF = nelts*(2+1);// 2 faces par elt
  }
  else if (strcmp(eltType, "TRI") == 0)
  {
    nfaces = 3*nelts;
    sizeFN = nfaces*(2+1);// 2 sommets par face
    sizeEF = nelts*(3+1); // 3 faces par elt
  }
  else if (strcmp(eltType, "QUAD") == 0)
  {
    nfaces = 4*nelts;
    sizeFN = nfaces*(2+1);// 2 sommets par face
    sizeEF = nelts*(4+1); // 4 faces par elt
  }
  else if (strcmp(eltType, "TETRA") == 0)
  {
    nfaces = 4*nelts;
    sizeFN = nfaces*(3+1);// 3 sommets par face
    sizeEF = nelts*(4+1); // 4 faces par elt
  }
  else if (strcmp(eltType, "HEXA") == 0)
  {
    nfaces = 6*nelts;
    sizeFN = nfaces*(4+1);// 4 sommets par face
    sizeEF = nelts*(6+1); // 6 faces par elt
  }
  else if (strcmp(eltType, "PENTA") == 0)
  {
    nfaces = 5*nelts;
    sizeFN = nelts*(3*(4+1)+2*(3+1)); // 3 quad et 2 tri par elt
    sizeEF = nelts*(5+1); // 5 faces par elt
  }
  else if (strcmp(eltType, "PYRA") == 0)
  {
    nfaces = 5*nelts;
    sizeFN = nelts*(1*(4+1)+4*(3+1));// 1 quad et 4 tri par elt
    sizeEF = nelts*(5+1); // 5 faces par elt
  }
  E_Int csize = sizeFN+sizeEF+4;

  // PyObject* tpl = K_ARRAY::buildArray(nfld, varString, npts, nelts, -1, 
  //                                     "NGON", false, csize);
  // E_Float* fieldp = K_ARRAY::getFieldPtr(tpl);
  // FldArrayF field(npts, nfld, fieldp, true); field = *f;
  // E_Int* cnp = K_ARRAY::getConnectPtr(tpl);
  // FldArrayI cn(csize, 1, cnp, true);

  FldArrayF field(npts, nfld); field = *f;
  FldArrayI cn(csize, 1);
  E_Int* cnp = cn.begin(1);

  // Build the NGon connectivity
  cnp[0] = nfaces;
  cnp[1] = sizeFN;
  cnp[sizeFN+2] = nelts; 
  cnp[sizeFN+3] = sizeEF;
  E_Int* cFN = cnp+2;
  E_Int* cEF = cnp+sizeFN+4;
  E_Int c1 = 0; E_Int c2 = 0; E_Int nof = 1;
  // Connectivite FN
  if (strcmp(eltType, "BAR") == 0) 
  {
    E_Int* cn1 = cnl->begin(1);
    E_Int* cn2 = cnl->begin(2);
    for (E_Int et = 0; et < nelts; et++)
    {
      E_Int v1 = cn1[et]; E_Int v2 = cn2[et]; 
      // connectivite face/noeuds
      cFN[c1]   = 1; cFN[c1+1] = v1;// face 1
      cFN[c1+2] = 1; cFN[c1+3] = v2;// face 2
      c1 += 4;
      // connectivite elt/faces
      cEF[c2] = 2; cEF[c2+1] = nof; cEF[c2+2] = nof+1;
      c2 += 3; nof += 2;
    }
  }
  else if (strcmp(eltType, "TRI") == 0) 
  {
    E_Int* cn1 = cnl->begin(1);
    E_Int* cn2 = cnl->begin(2);
    E_Int* cn3 = cnl->begin(3);
    for (E_Int et = 0; et < nelts; et++)
    {
      E_Int v1 = cn1[et]; E_Int v2 = cn2[et]; E_Int v3 = cn3[et]; 
      // connectivite face/noeuds
      cFN[c1]   = 2; cFN[c1+1] = v1; cFN[c1+2] = v2; // face 1
      cFN[c1+3] = 2; cFN[c1+4] = v2; cFN[c1+5] = v3; // face 2
      cFN[c1+6] = 2; cFN[c1+7] = v3; cFN[c1+8] = v1; // face 3
      c1 += 9;
      // connectivite elt/faces
      cEF[c2] = 3; cEF[c2+1] = nof; cEF[c2+2] = nof+1; cEF[c2+3] = nof+2;   
      c2 += 4; nof += 3;
    }
  }
  else if (strcmp(eltType, "QUAD") == 0) 
  {
    E_Int* cn1 = cnl->begin(1);
    E_Int* cn2 = cnl->begin(2);
    E_Int* cn3 = cnl->begin(3);
    E_Int* cn4 = cnl->begin(4);
    for (E_Int et = 0; et < nelts; et++)
    {
      E_Int v1 = cn1[et]; E_Int v2 = cn2[et]; E_Int v3 = cn3[et]; E_Int v4 = cn4[et]; 
      // connectivite face/noeuds
      cFN[c1]   = 2; cFN[c1+1] = v1;  cFN[c1+2] = v2;// face 1
      cFN[c1+3] = 2; cFN[c1+4] = v2;  cFN[c1+5] = v3;// face 2
      cFN[c1+6] = 2; cFN[c1+7] = v3;  cFN[c1+8] = v4;// face 3
      cFN[c1+9] = 2; cFN[c1+10] = v4; cFN[c1+11] = v1;// face 4
      c1 += 12;
      // connectivite elt/faces
      cEF[c2] = 4; cEF[c2+1] = nof; cEF[c2+2] = nof+1; cEF[c2+3] = nof+2; cEF[c2+4] = nof+3;      
      c2 += 5;nof += 4;
    }
  }
  else if (strcmp(eltType, "TETRA") == 0) 
  {
    E_Int* cn1 = cnl->begin(1);
    E_Int* cn2 = cnl->begin(2);
    E_Int* cn3 = cnl->begin(3);
    E_Int* cn4 = cnl->begin(4);
    E_Int v1, v2, v3, v4;
    for (E_Int et = 0; et < nelts; et++)
    {
      v1 = cn1[et]; v2 = cn2[et]; v3 = cn3[et]; v4 = cn4[et]; 
      // connectivite face/noeuds
      cFN[c1]   = 3; cFN[c1+1] = v1; cFN[c1+2] = v2; cFN[c1+3] = v3;// face 1
      cFN[c1+4] = 3; cFN[c1+5] = v1; cFN[c1+6] = v2; cFN[c1+7] = v4;// face 2
      cFN[c1+8] = 3; cFN[c1+9] = v2; cFN[c1+10] = v3; cFN[c1+11] = v4;// face 3
      cFN[c1+12] =3; cFN[c1+13] = v3; cFN[c1+14] = v1; cFN[c1+15] = v4;// face 4
      c1 += 16;      
      // connectivite elt/faces
      cEF[c2] = 4; cEF[c2+1] = nof; cEF[c2+2] = nof+1; cEF[c2+3] = nof+2; cEF[c2+4] = nof+3;      
      c2 += 5; nof += 4;
    }
  }
  else if (strcmp(eltType, "HEXA") == 0) 
  {
    E_Int* cn1 = cnl->begin(1);
    E_Int* cn2 = cnl->begin(2);
    E_Int* cn3 = cnl->begin(3);
    E_Int* cn4 = cnl->begin(4);
    E_Int* cn5 = cnl->begin(5);
    E_Int* cn6 = cnl->begin(6);
    E_Int* cn7 = cnl->begin(7);
    E_Int* cn8 = cnl->begin(8);
    E_Int v1, v2, v3, v4, v5, v6, v7, v8;
    for (E_Int et = 0; et < nelts; et++)
    {
      v1 = cn1[et]; v2 = cn2[et]; v3 = cn3[et]; v4 = cn4[et]; 
      v5 = cn5[et]; v6 = cn6[et]; v7 = cn7[et]; v8 = cn8[et]; 
      // connectivite face/noeuds
      cFN[c1]   = 4; cFN[c1+1] = v1; cFN[c1+2] = v2; cFN[c1+3] = v3; cFN[c1+4] = v4;// face 1
      cFN[c1+5] = 4; cFN[c1+6] = v5; cFN[c1+7] = v6; cFN[c1+8] = v7; cFN[c1+9] = v8;// face 2
      cFN[c1+10]= 4; cFN[c1+11] = v1; cFN[c1+12] = v2; cFN[c1+13] = v6; cFN[c1+14] = v5;// face 3
      cFN[c1+15]= 4; cFN[c1+16] = v4; cFN[c1+17] = v3; cFN[c1+18] = v7; cFN[c1+19] = v8;// face 4
      cFN[c1+20]= 4; cFN[c1+21] = v1; cFN[c1+22] = v4; cFN[c1+23] = v8; cFN[c1+24] = v5;// face 5
      cFN[c1+25]= 4; cFN[c1+26] = v2; cFN[c1+27] = v3; cFN[c1+28] = v7; cFN[c1+29] = v6;// face 6
      c1+=30; 
      // connectivite elt/faces
      cEF[c2] = 6; cEF[c2+1] = nof; cEF[c2+2] = nof+1; cEF[c2+3] = nof+2; cEF[c2+4] = nof+3; cEF[c2+5] = nof+4; cEF[c2+6] = nof+5;
      c2 += 7; nof += 6;
    }
  }
  else if (strcmp(eltType, "PENTA") == 0) 
  {
    E_Int* cn1 = cnl->begin(1);
    E_Int* cn2 = cnl->begin(2);
    E_Int* cn3 = cnl->begin(3);
    E_Int* cn4 = cnl->begin(4);
    E_Int* cn5 = cnl->begin(5);
    E_Int* cn6 = cnl->begin(6);
    E_Int v1, v2, v3, v4, v5, v6;
    for (E_Int et = 0; et < nelts; et++)
    {
      v1 = cn1[et]; v2 = cn2[et]; v3 = cn3[et]; 
      v4 = cn4[et]; v5 = cn5[et]; v6 = cn6[et];  
      // connectivite face/noeuds
      cFN[c1]   = 3; cFN[c1+1] = v1; cFN[c1+2] = v2; cFN[c1+3] = v3;// face 1 TRI
      cFN[c1+4] = 3; cFN[c1+5] = v4; cFN[c1+6] = v5; cFN[c1+7] = v6;// face 2 TRI
      cFN[c1+8]  = 4; cFN[c1+9] = v1; cFN[c1+10] = v2; cFN[c1+11] = v5; cFN[c1+12] = v4;// face 3 : QUAD
      cFN[c1+13] = 4; cFN[c1+14] = v2; cFN[c1+15] = v3; cFN[c1+16] = v6; cFN[c1+17] = v5;// face 4 : QUAD
      cFN[c1+18] = 4; cFN[c1+19] = v3; cFN[c1+20] = v1; cFN[c1+21] = v4; cFN[c1+22] = v6;// face 5 : QUAD
      c1 += 23;
      // connectivite elt/faces
      cEF[c2] = 5; cEF[c2+1] = nof; cEF[c2+2] = nof+1; cEF[c2+3] = nof+2; cEF[c2+4] = nof+3; cEF[c2+5] = nof+4;      
      c2 += 6; nof += 5;
    }
  }
  else if (strcmp(eltType, "PYRA") == 0) 
  {
    E_Int* cn1 = cnl->begin(1);
    E_Int* cn2 = cnl->begin(2);
    E_Int* cn3 = cnl->begin(3);
    E_Int* cn4 = cnl->begin(4);
    E_Int* cn5 = cnl->begin(5);
    E_Int v1, v2, v3, v4, v5;
    for (E_Int et = 0; et < nelts; et++)
    {
      v1 = cn1[et]; v2 = cn2[et]; v3 = cn3[et]; 
      v4 = cn4[et]; v5 = cn5[et]; 
      // connectivite face/noeuds
      cFN[c1]   = 4; cFN[c1+1] = v1; cFN[c1+2] = v2; cFN[c1+3] = v3; cFN[c1+4] = v4; // face 1: QUAD
      cFN[c1+5] = 3; cFN[c1+6] = v1; cFN[c1+7] = v2; cFN[c1+8] = v5;// face 2: TRI
      cFN[c1+9] = 3; cFN[c1+10] = v2; cFN[c1+11] = v3; cFN[c1+12] = v5;// face 3: TRI
      cFN[c1+13]= 3; cFN[c1+14] = v3; cFN[c1+15] = v4; cFN[c1+16] = v5;// face 4: TRI
      cFN[c1+17]= 3; cFN[c1+18] = v4; cFN[c1+19] = v1; cFN[c1+20] = v5;// face 5: TRI
      c1+=21;
      // connectivite elt/faces
      cEF[c2] = 5; cEF[c2+1] = nof; cEF[c2+2] = nof+1; cEF[c2+3] = nof+2; cEF[c2+4] = nof+3; cEF[c2+5] = nof+4;      
      c2 += 6; nof += 5;
    }
  }
  RELEASESHAREDU(array, f, cnl); 
  
  /* clean connectivity */
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
