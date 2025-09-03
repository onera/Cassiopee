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
#ifndef _OCC_OCC_H_
#define _OCC_OCC_H_

#include "kcore.h"

namespace K_OCC
{
  PyObject* convertCAD2Arrays0(PyObject* self, PyObject* args); // with OCC internal
  PyObject* convertCAD2Arrays1(PyObject* self, PyObject* args); // with T3Mesher
  PyObject* convertCAD2Arrays2(PyObject* self, PyObject* args); // with T3Mesher
  
  PyObject* readCAD(PyObject* self, PyObject* args);
  PyObject* writeCAD(PyObject* self, PyObject* args);
  PyObject* createEmptyCAD(PyObject* self, PyObject* args);
  PyObject* mergeCAD(PyObject* self, PyObject* args);
  PyObject* freeHook(PyObject* self, PyObject* args);

  PyObject* printOCAF(PyObject* self, PyObject* args);
  PyObject* printShapeOCAF(PyObject* self, PyObject* args);
  PyObject* getFaceNameInOCAF(PyObject* self, PyObject* args);
  PyObject* getFaceNameInOCAF2(PyObject* self, PyObject* args);
  PyObject* getEdgeNameInOCAF2(PyObject* self, PyObject* args);

  PyObject* bottle(PyObject* self, PyObject* args);
  PyObject* addSphere(PyObject* self, PyObject* args);
  PyObject* addCylinder(PyObject* self, PyObject* args);
  PyObject* addBox(PyObject* self, PyObject* args);
  PyObject* addSquare(PyObject* self, PyObject* args);
  PyObject* addLine(PyObject* self, PyObject* args);
  PyObject* addCircle(PyObject* self, PyObject* args);
  PyObject* addSpline(PyObject* self, PyObject* args);
  
  PyObject* getNbFaces(PyObject* self, PyObject* args);
  PyObject* getNbEdges(PyObject* self, PyObject* args);
  PyObject* getFileAndFormat(PyObject* self, PyObject* args);

  PyObject* meshGlobalEdges1(PyObject* self, PyObject* args);
  PyObject* meshGlobalEdges2(PyObject* self, PyObject* args);
  PyObject* meshGlobalEdges3(PyObject* self, PyObject* args);
  PyObject* meshGlobalEdges4(PyObject* self, PyObject* args);
  PyObject* meshEdgesByFace(PyObject* self, PyObject* args);
  PyObject* meshEdgesByFace2(PyObject* self, PyObject* args);
  PyObject* meshEdgesByFace3(PyObject* self, PyObject* args);
  PyObject* getEdgeNoByFace(PyObject* self, PyObject* args);
  PyObject* identifyLoopsInEdges(PyObject* self, PyObject* args);
  PyObject* evalEdge(PyObject* self, PyObject* args);
  PyObject* evalFace(PyObject* self, PyObject* args);
  PyObject* projectOnFaces(PyObject* self, PyObject* args);
  PyObject* projectOnEdges(PyObject* self, PyObject* args);
  PyObject* linkNodes2CAD(PyObject* self, PyObject* args);
  PyObject* updateFcadidFromNcadid(PyObject* self, PyObject* args);
  PyObject* updateNcadidFromFcadid(PyObject* self, PyObject* args);
  PyObject* getNodalParameters(PyObject* self, PyObject* args);
  PyObject* trimesh(PyObject* self, PyObject* args);

  PyObject* meshOneEdge(PyObject* self, PyObject* args);
  PyObject* meshEdgesOfFace(PyObject* self, PyObject* args);

  PyObject* analyseEdges(PyObject* self, PyObject* args);
  PyObject* getFaceArea(PyObject* self, PyObject* args);
  PyObject* getFaceOrientation(PyObject* self, PyObject* args);
  PyObject* areEdgeIdentical(PyObject* self, PyObject* args);

  PyObject* splitFaces(PyObject* self, PyObject* args);
  PyObject* fixShape(PyObject* self, PyObject* args);
  PyObject* trimFaces(PyObject* self, PyObject* args);
  PyObject* sewing(PyObject* self, PyObject* args);
  PyObject* removeFaces(PyObject* self, PyObject* args);
  PyObject* fillHole(PyObject* self, PyObject* args);
  PyObject* addFillet(PyObject* self, PyObject* args);
  PyObject* mergeFaces(PyObject* self, PyObject* args);

  PyObject* translate(PyObject* self, PyObject* args);
  PyObject* scale(PyObject* self, PyObject* args);
  PyObject* rotate(PyObject* self, PyObject* args);

  PyObject* getOppData(PyObject* self, PyObject* args);
  PyObject* identifyTags(PyObject* self, PyObject* args);
}

#endif
