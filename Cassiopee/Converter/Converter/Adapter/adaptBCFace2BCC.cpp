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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Convert BCFace numpy to a BCC numpy 
   BCFace: list of indices of faces (global indices)
   BCC: connectivity of boundary */
//=============================================================================
PyObject* K_CONVERTER::adaptBCFace2BCC(PyObject* self, PyObject* args)
{
  PyObject* BCFaceO, *cnO;
  char* eltType;
  if (!PYPARSETUPLE_(args, OO_ S_, &BCFaceO, &cnO, &eltType)) return NULL;

  // Check numpy (BCFace)
  FldArrayI* BCFace;
  E_Int res = K_NUMPY::getFromPointList(BCFaceO, BCFace);
  E_Int* faces = BCFace->begin();
  E_Int nint = BCFace->getSize();
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "adaptBCFace2BCC: BCFace is invalid.");
    return NULL;
  }

  // Check numpy (cn)
  FldArrayI* cn;
  res = K_NUMPY::getFromNumpyArray(cnO, cn);
  if (res == 0)
  {
    RELEASESHAREDN(BCFaceO, BCFace);
    PyErr_SetString(PyExc_TypeError, 
                    "adaptBCFace2BCC: cn is invalid.");
    return NULL;
  }
  E_Int* cnp = cn->begin();
  
  E_Int face[6][4];
  E_Int nfaces = 0; /*E_Int nof = 0;*/ E_Int type = 0;
  /*E_Int np = 0;*/
  if (strcmp(eltType, "BAR") == 0)
  { 
    nfaces = 2; /*nof = 1;*/ type = 0; /*np = 2;*/
    face[0][0] = 1; face[1][0] = 2;
  }
  else if (strcmp(eltType, "QUAD") == 0) 
  {
    nfaces = 4; /*nof = 2;*/ type = 1; /*np = 4;*/
    face[0][0] = 1; face[0][1] = 2;
    face[1][0] = 2; face[1][1] = 3;
    face[2][0] = 3; face[2][1] = 4;
    face[3][0] = 4; face[3][1] = 1;
  }
  else if (strcmp(eltType, "TRI") == 0) 
  {
    nfaces = 3; /*nof = 2;*/ type = 1; /*np = 3;*/
    face[0][0] = 1; face[0][1] = 2;
    face[1][0] = 2; face[1][1] = 3;
    face[2][0] = 3; face[2][1] = 1;
  }
  else if (strcmp(eltType, "HEXA") == 0) 
  {
    nfaces = 6; /*nof = 4;*/ type = 3; /*np = 8;*/
    face[0][0] = 1; face[0][1] = 4; face[0][2] = 3; face[0][3] = 2;
    face[1][0] = 1; face[1][1] = 2; face[1][2] = 6; face[1][3] = 5;
    face[2][0] = 2; face[2][1] = 3; face[2][2] = 7; face[2][3] = 6;
    face[3][0] = 3; face[3][1] = 4; face[3][2] = 8; face[3][3] = 7;
    face[4][0] = 1; face[4][1] = 5; face[4][2] = 8; face[4][3] = 4;
    face[5][0] = 5; face[5][1] = 6; face[5][2] = 7; face[5][3] = 8;
  }
  else if (strcmp(eltType, "TETRA") == 0) 
  {
    nfaces = 4; /*nof = 3;*/ type = 2; /*np = 4;*/
    face[0][0] = 1; face[0][1] = 3; face[0][2] = 2;
    face[1][0] = 1; face[1][1] = 2; face[1][2] = 4;
    face[2][0] = 2; face[2][1] = 3; face[2][2] = 4;
    face[3][0] = 3; face[3][1] = 1; face[3][2] = 4;
  }
  else if (strcmp(eltType, "PYRA") == 0) 
  {
    nfaces = 5; /*nof = 3; np = 5;*/
    type = -1; // In this case, we must output 2 connectivities
    face[0][0] = 1; face[0][1] = 4; face[0][2] = 3; // T
    face[1][0] = 3; face[1][1] = 2; face[1][2] = 1; // T
    face[2][0] = 1; face[2][1] = 2; face[2][2] = 5; // T
    face[3][0] = 2; face[3][1] = 3; face[3][2] = 5; // T
    face[4][0] = 3; face[4][1] = 4; face[4][2] = 5; face[4][3] = 1; // Q
  }
  else if (strcmp(eltType, "PENTA") == 0) 
  {
    nfaces = 5; /*nof = 4; np = 6;*/
    type = -2; // In this case, we must output 2 connectivities
    face[0][0] = 1; face[0][1] = 2; face[0][2] = 5; face[0][3] = 4; // Q
    face[1][0] = 2; face[1][1] = 3; face[1][2] = 6; face[1][3] = 5; // Q
    face[2][0] = 3; face[2][1] = 1; face[2][2] = 4; face[2][3] = 6; // Q
    face[3][0] = 1; face[3][1] = 3; face[3][2] = 2; // T
    face[4][0] = 4; face[4][1] = 5; face[4][2] = 6; // T
  }
  
  E_Int indFace, inde, indf;

  // Compte les triangles et les quads
  E_Int nBARS = 0; E_Int nTRIS = 0; E_Int nQUADS = 0;
  for (E_Int i = 0; i < nint; i++)
  {
    indFace = faces[i]-1;
    inde = indFace / nfaces;
    indf = indFace - inde*nfaces;
    if (type == 1) nBARS++;
    else if (type == 2) nTRIS++;
    else if (type == 3) nQUADS++;
    else if (type == -1) { if (indf>3) nQUADS++; else nTRIS++; }
    else if (type == -2) { if (indf>2) nTRIS++; else nQUADS++; }
  }
  printf("count " SF_D3_ "\n", nBARS, nTRIS, nQUADS);
  PyObject* tplB = K_NUMPY::buildNumpyArray(nBARS*2, 1, 1, 1);
  PyObject* tplT = K_NUMPY::buildNumpyArray(nTRIS*3, 1, 1, 1);
  PyObject* tplQ = K_NUMPY::buildNumpyArray(nQUADS*4, 1, 1, 1);

  E_Int* cB = K_NUMPY::getNumpyPtrI(tplB);
  E_Int* cT = K_NUMPY::getNumpyPtrI(tplT);
  E_Int* cQ = K_NUMPY::getNumpyPtrI(tplQ);
 
  E_Int cBAR = 0; E_Int cTRI = 0; E_Int cQUAD = 0;
  for (E_Int i = 0; i < nint; i++)
  {
    indFace = faces[i]-1;
    inde = indFace / nfaces;
    indf = indFace - inde*nfaces;
    printf("inde " SF_D2_ "\n", inde, indf);
    if (type == 1) // BARS
    {
      for (E_Int f = 0; f < 2; f++) 
      { 
        cB[cBAR+f] = cnp[2*inde+face[indf][f]-1];
      }
      cBAR += 2;
    }
    else if (type == 2) // TRI
    {
      for (E_Int f = 0; f < 3; f++) 
      { 
        cT[cTRI+f] = cnp[3*inde+face[indf][f]-1]; 
      }
      cTRI += 3;
    }
    else if (type == 3) // QUAD
    {
      for (E_Int f = 0; f < 4; f++) 
      { 
        cQ[cQUAD+f] = cnp[inde*4+face[indf][f]-1]; 
      }
      cQUAD += 4;
    }
    else if (type == -1) 
    {
      if (indf>3) 
      { 
        for (E_Int f = 0; f < 4; f++) 
        { 
          cQ[cQUAD+f] = cnp[inde*5+face[indf][f]-1]; 
        }
        cQUAD += 4;
      }
      else 
      { 
        for (E_Int f = 0; f < 3; f++) 
        { 
          cT[cTRI+f] = cnp[inde*5+face[indf][f]-1]; 
        }
        cTRI += 3;
      }
    }
    else if (type == -2)
    {
      if (indf>2) 
      {
        for (E_Int f = 0; f < 3; f++) 
        { 
          cT[cTRI+f] = cnp[inde*6+face[indf][f]-1]; 
        }
        cTRI += 3; 
      }
      else 
      { 
        for (E_Int f = 0; f < 4; f++) 
        { 
          cQ[cQUAD+f] = cnp[inde*6+face[indf][f]-1]; 
        }
        cQUAD += 4;
      }
    }
  }
  //for (E_Int i = 0; i < 4*nQUADS; i++) printf("%d - %d\n", i, cQ[i]);

  RELEASESHAREDN(BCFaceO, BCFace);
  RELEASESHAREDN(cnO, cn);

  return Py_BuildValue("(OOO)", tplB, tplT, tplQ);
}
