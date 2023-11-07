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
#include "proto.h"

PyObject *K_XCORE::extractLeafMesh(PyObject *self, PyObject *args)
{
  PyObject *o;

  if (!PYPARSETUPLE_(args, O_, &o)) {
    PyErr_SetString(PyExc_ValueError, "adaptMesh(): wrong input.");
    return NULL;
  }

  // Unpack mesh
  if (!PyCapsule_IsValid(o, "adaptMesh")) {
    PyErr_SetString(PyExc_TypeError, "extractLeafMesh: bad mesh hook.");
    return NULL;
  }
  mesh *M = (mesh *)PyCapsule_GetPointer(o, "adaptMesh");
  
  tree &CT = M->ctree;
  tree &FT = M->ftree;

  auto &ecells = CT.l2g;
  ecells.clear();
  for (E_Int i = 0; i < M->ncells; i++) {
    if (CT.enabled[i]) ecells.push_back(i);
  }
  E_Int necells = ecells.size();

  std::unordered_map<E_Int, E_Int> efaces;
  std::unordered_map<E_Int, E_Int> epoints;
  E_Int nefaces = 0, nepoints = 0;

  std::vector<E_Int> NFACE;
  std::vector<E_Int> NGON;
  std::vector<E_Int> XCELLS(1, 0);
  std::vector<E_Int> XFACES(1, 0);

  // Make NFACE
  for (E_Int i = 0; i < necells; i++) {
    E_Int cell = ecells[i];
    for (E_Int j = M->xcells[cell]; j < M->xcells[cell+1]; j++) {
      E_Int face = M->NFACE[j];
      NFACE.push_back(face);
    }
    XCELLS.push_back(M->xcells[cell+1]-M->xcells[cell]);
  }
  
  for (E_Int i = 0; i < necells; i++)
    XCELLS[i+1] += XCELLS[i];

  // Make NGON
  FT.l2g.clear();
  for (E_Int i = 0; i < necells; i++) {
    for (E_Int j = XCELLS[i]; j < XCELLS[i+1]; j++) {
      E_Int face = NFACE[j];
      if (efaces.find(face) == efaces.end()) {
        FT.l2g.push_back(face);
        efaces[face] = nefaces++;

        for (E_Int k = M->xfaces[face]; k < M->xfaces[face+1]; k++) {
          E_Int point = M->NGON[k];
          NGON.push_back(point);

          if (epoints.find(point) == epoints.end()) {
            M->pointMap[nepoints] = point;
            epoints[point] = nepoints++;
          }
        }

        XFACES.push_back(M->xfaces[face+1]-M->xfaces[face]);
      }
    }
  }

  for (E_Int i = 0; i < nefaces; i++)
    XFACES[i+1] += XFACES[i];

  // Build array
  const char *varStringOut = "CoordinateX,CoordinateY,CoordinateZ";

  PyObject* m = K_ARRAY::buildArray3(3, varStringOut, nepoints, necells,
    nefaces, "NGON", XFACES[nefaces], XCELLS[necells], 3, false, 3);

  FldArrayF *fo;
  FldArrayI *cno;
  K_ARRAY::getFromArray3(m, fo, cno);

  E_Float *pxo = fo->begin(1);
  E_Float *pyo = fo->begin(2);
  E_Float *pzo = fo->begin(3);

  // Fill out coordinates
  for (E_Int i = 0; i < nepoints; i++) {
    E_Int point = M->pointMap[i];
    E_Float *ptr = &M->xyz[3*point]; 
    pxo[i] = ptr[0];
    pyo[i] = ptr[1];
    pzo[i] = ptr[2];
  }

  // Fill out indPG and NGON
  E_Int *indPGo = cno->getIndPG();
  for (E_Int i = 0; i < nefaces+1; i++) indPGo[i] = XFACES[i];

  E_Int *ngono = cno->getNGon();
  E_Int *ptr = ngono;
  for (E_Int i = 0; i < nefaces; i++) {
    E_Int start = XFACES[i];
    E_Int end = XFACES[i+1];
    
    for (E_Int j = start; j < end; j++)
      *ptr++ = epoints[NGON[j]]+1;
  }

  // Fill out indPH and NFACE
  E_Int *indPHo = cno->getIndPH();
  for (E_Int i = 0; i < necells+1; i++) indPHo[i] = XCELLS[i];

  E_Int *nfaceo = cno->getNFace();
  ptr = nfaceo;
  for (E_Int i = 0; i < necells; i++) {
    E_Int start = XCELLS[i];
    E_Int end = XCELLS[i+1];
    for (E_Int j = start; j < end; j++)
      *ptr++ = efaces[NFACE[j]]+1;
  }

  return m;
}
