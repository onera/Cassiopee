/*    
    Copyright 2013-2026 ONERA.

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

#include "occ.h"

#include "TopoDS.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "Bnd_Box.hxx"
#include "BRepBndLib.hxx"
#include "BRep_Builder.hxx"

// ============================================================================
/* Return the bounding box of given face */
// ============================================================================
PyObject* K_OCC::getBoundingBox(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listFaces;
  if (!PYPARSETUPLE_(args, O_ O_, &hook, &listFaces)) return NULL;

  GETSHAPE;
  GETMAPSURFACES;

  E_Float xmin = E_MAXFLOAT;
  E_Float ymin = E_MAXFLOAT;
  E_Float zmin = E_MAXFLOAT;
  E_Float xmax = -E_MAXFLOAT;
  E_Float ymax = -E_MAXFLOAT;
  E_Float zmax = -E_MAXFLOAT;

  // By face
  if (listFaces == Py_None) // all faces of topshape
  {
    Bnd_Box box; 
    BRepBndLib::Add(*shape, box); 
    // compute bounding box 
    box.SetGap(0.0);
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  }
  else
  {
    // make compound
    E_Int nfaces = PyList_Size(listFaces);
    BRep_Builder builder;
    TopoDS_Compound shc;
    builder.MakeCompound(shc);
    
    for (E_Int no = 0; no < nfaces; no++)
    {
      PyObject* noFaceO = PyList_GetItem(listFaces, no);
      E_Int noFace = PyInt_AsLong(noFaceO);
      const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
      builder.Add(shc, F);
    }
    Bnd_Box box; 
    BRepBndLib::Add(shc, box); 
    // compute bounding box 
    box.SetGap(0.0);
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  }
  return Py_BuildValue("dddddd", xmin,ymin,zmin,xmax,ymax,zmax);
} 