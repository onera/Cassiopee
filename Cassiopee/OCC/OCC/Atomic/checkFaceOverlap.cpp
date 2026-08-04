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
// Check if two faces overlap

#include "occ.h"

#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS.hxx"
#include "BRepAlgoAPI_Common.hxx"
#include "BRepGProp.hxx"
#include "GProp_GProps.hxx"

// ============================================================================
/* Return the face area of given face */
// ============================================================================
PyObject* K_OCC::checkFaceOverlap(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int face1; E_Int face2;
  E_Float tol;
  if (!PYPARSETUPLE_(args, O_ II_ R_, &hook, &face1, &face2, &tol)) return NULL;

  GETPACKET;
  GETMAPSURFACES;
  const TopoDS_Face& F1 = TopoDS::Face(surfaces(face1));
  const TopoDS_Face& F2 = TopoDS::Face(surfaces(face2));
  
  // Find common part
  BRepAlgoAPI_Common common(F1, F2);

  // Explore result in faces
  TopExp_Explorer expl1;
  TopTools_IndexedMapOfShape sf;
  TopExp::MapShapes(common, TopAbs_FACE, sf);
  E_Int nof = sf.Extent(); 
  //if (nof == 0) printf("no overlap.\n");
  //else printf("two faces overlap.\n");
  
  E_Float returnValue = 0.0;

  // compute commmon area
  GProp_GProps props;
  BRepGProp::SurfaceProperties(F1, props);
  E_Float area1 = props.Mass();
  BRepGProp::SurfaceProperties(F2, props);
  //E_Float area2 = props.Mass();
  if (nof > 0)
  {
    BRepGProp::SurfaceProperties(sf(1), props);
    E_Float area = props.Mass();
    if (area > area1 - tol) 
    { 
      //printf("full overlap.\n");
      returnValue = 1.0;
    }
    else 
    {
      //printf("partial overlap %g.\n", area/area1); 
      returnValue = area/area1;
    }
  }

  // Explore result in edges
  TopTools_IndexedMapOfShape se;
  TopExp::MapShapes(common, TopAbs_EDGE, se);
  E_Int noe = se.Extent(); 
  if (noe == 0) 
  {
    //printf("no overlap.\n");
  }
  else if (noe == 1) 
  { 
    //printf("overlap in one edge/contact.\n");
    returnValue = -1.0;
  }
  else 
  {
    //printf("overlap in multiple edges.\n");
    returnValue = -noe;
  }

  // return the percentage of overlap (0 to 1)
  return Py_BuildValue("d", returnValue);
}
