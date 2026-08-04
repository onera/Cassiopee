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

#include "TColStd_HSequenceOfTransient.hxx"
#include "TopoDS.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "ShapeAnalysis.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"

// ============================================================================
/* Return the face area of given face */
// ============================================================================
PyObject* K_OCC::getFaceArea(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listFaces;
  if (!PYPARSETUPLE_(args, O_ O_, &hook, &listFaces)) return NULL;

  GETPACKET;
  GETMAPSURFACES;

  E_Float area = 0.;
  if (listFaces == Py_None) // all faces of topshape
  {
    for (E_Int i=1; i <= surfaces.Extent(); i++)
    {
      GProp_GProps props;
      BRepGProp::SurfaceProperties(surfaces(i), props);
      area += props.Mass();
    }
  }
  else
  {
    E_Int nfaces = PyList_Size(listFaces);
    for (E_Int no = 0; no < nfaces; no++)
    {
      PyObject* noFaceO = PyList_GetItem(listFaces, no);
      E_Int noFace = PyInt_AsLong(noFaceO);

      GProp_GProps props;
      BRepGProp::SurfaceProperties(surfaces(noFace), props);
      area += props.Mass();
    }
  }
  return Py_BuildValue("d", area);
} 