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
/* Convert BCC numpy to a BCFace numpy
   BCC: connectivity of boundary 
   BCFace: list of indices of faces (global indices)
   To be finished.
*/
//=============================================================================
PyObject* K_CONVERTER::adaptBCC2BCFace(PyObject* self, PyObject* args)
{
  PyObject* BCCO, *cnO;
  char* eltType;
  if (!PyArg_ParseTuple(args, "OOs", &BCCO, &cnO, &eltType)) return NULL;

  // Check numpy (BCFace)
  FldArrayI* BCC;
  E_Int res = K_NUMPY::getFromNumpyArray(BCCO, BCC, true);
  //E_Int* bcc = BCC->begin();
  //E_Int nint = BCC->getSize();
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "adaptBCC2BCFace: BCC is invalid.");
    return NULL;
  }

  // Check numpy (cn - connectivite volumique)
  FldArrayI* cn;
  res = K_NUMPY::getFromNumpyArray(cnO, cn, true);
  if (res == 0)
  {
    RELEASESHAREDN(BCCO, BCC);
    PyErr_SetString(PyExc_TypeError, 
                    "adaptBCC2BCFace: cn is invalid.");
    return NULL;
  }
  //E_Int* cnp = cn->begin();
  
  printf("eltType %s\n", eltType);
  E_Int face[6][4];
  E_Int nfaces = 0; E_Int nof = 0; E_Int type = 0;
  E_Int np = 0;
  if (strcmp(eltType, "BAR") == 0)
  { 
    nfaces = 2; nof = 1; type = 0; np = 2;
    face[0][0] = 1; face[1][0] = 2;
  }
  else if (strcmp(eltType, "QUAD") == 0) 
  {
    nfaces = 4; nof = 2; type = 1; np = 4;
    face[0][0] = 1; face[0][1] = 2;
    face[1][0] = 2; face[1][1] = 3;
    face[2][0] = 3; face[2][1] = 4;
    face[3][0] = 4; face[3][1] = 1;
  }
  else if (strcmp(eltType, "TRI") == 0) 
  {
    nfaces = 3; nof = 2; type = 1; np = 3;
    face[0][0] = 1; face[0][1] = 2;
    face[1][0] = 2; face[1][1] = 3;
    face[2][0] = 3; face[2][1] = 1;
  }
  else if (strcmp(eltType, "HEXA") == 0) 
  {
    nfaces = 6; nof = 4; type = 3; np = 8;
    face[0][0] = 1; face[0][1] = 4; face[0][2] = 3; face[0][3] = 2;
    face[1][0] = 1; face[1][1] = 2; face[1][2] = 6; face[1][3] = 5;
    face[2][0] = 2; face[2][1] = 3; face[2][2] = 7; face[2][3] = 6;
    face[3][0] = 3; face[3][1] = 4; face[3][2] = 8; face[3][3] = 7;
    face[4][0] = 1; face[4][1] = 5; face[4][2] = 8; face[4][3] = 4;
    face[5][0] = 5; face[5][1] = 6; face[5][2] = 7; face[5][3] = 8;
  }
  else if (strcmp(eltType, "TETRA") == 0) 
  {
    nfaces = 4; nof = 3; type = 2; np = 4;
    face[0][0] = 1; face[0][1] = 3; face[0][2] = 2;
    face[1][0] = 1; face[1][1] = 2; face[1][2] = 4;
    face[2][0] = 2; face[2][1] = 3; face[2][2] = 4;
    face[3][0] = 3; face[3][1] = 1; face[3][2] = 4;
  }
  else if (strcmp(eltType, "PYRA") == 0) 
  {
    nfaces = 5; nof = 3;  np = 5;
    type = -1; // In this case, we must output 2 connectivities
    face[0][0] = 1; face[0][1] = 4; face[0][2] = 3; // T
    face[1][0] = 3; face[1][1] = 2; face[1][2] = 1; // T
    face[2][0] = 1; face[2][1] = 2; face[2][2] = 5; // T
    face[3][0] = 2; face[3][1] = 3; face[3][2] = 5; // T
    face[4][0] = 3; face[4][1] = 4; face[4][2] = 5; face[4][3] = 1; // Q
  }
  else if (strcmp(eltType, "PENTA") == 0) 
  {
    nfaces = 5; nof = 4; np = 6;
    type = -2; // In this case, we must output 2 connectivities
    face[0][0] = 1; face[0][1] = 2; face[0][2] = 5; face[0][3] = 4; // Q
    face[1][0] = 2; face[1][1] = 3; face[1][2] = 6; face[1][3] = 5; // Q
    face[2][0] = 3; face[2][1] = 1; face[2][2] = 4; face[2][3] = 6; // Q
    face[3][0] = 1; face[3][1] = 3; face[3][2] = 2; // T
    face[4][0] = 4; face[4][1] = 5; face[4][2] = 6; // T
  }

  /* ... */

  RELEASESHAREDN(BCCO, BCC);
  RELEASESHAREDN(cnO, cn);

  return Py_BuildValue("O", Py_None);
}
