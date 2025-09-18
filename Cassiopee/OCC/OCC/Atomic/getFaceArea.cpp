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

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  //TopoDS_Shape* shp = (TopoDS_Shape*) packet[0];
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  //TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  TopExp_Explorer expl;

  //const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));

  // By wire
  /*
  for (expl.Init(surfaces(noFace), TopAbs_WIRE); expl.More(); expl.Next())
  {
    const TopoDS_Wire& W = TopoDS::Wire(expl.Current());
    E_Float surface = ShapeAnalysis::ContourArea(W);
    //printf("wire surface=%f\n", surface);	
    area = std::max(area, surface);
  }
  printf("area1=%g\n", area);
  */

  // By face
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