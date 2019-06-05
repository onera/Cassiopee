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

#ifndef _INTERSECTOR_INTERSECTOR_H_
#define _INTERSECTOR_INTERSECTOR_H_

# include "kcore.h"

namespace K_INTERSECTOR
{
  PyObject* conformUnstr(PyObject* self, PyObject* args);
  
  PyObject* booleanIntersection(PyObject* self, PyObject* args);
  PyObject* booleanUnion(PyObject* self, PyObject* args);
  PyObject* booleanMinus(PyObject* self, PyObject* args);
  PyObject* booleanIntersectionBorder(PyObject* self, PyObject* args);
  PyObject* booleanModifiedSolid(PyObject* self, PyObject* args);
  PyObject* DiffSurf(PyObject* self, PyObject* args);
  PyObject* XcellN(PyObject* self, PyObject* args);
  PyObject* P1ConservativeChimeraCoeffs(PyObject* self, PyObject* args);

  PyObject* selfX(PyObject* self, PyObject* args);
  PyObject* triangulateExteriorFaces(PyObject* self, PyObject* args);
  PyObject* triangulateSpecifiedFaces(PyObject* self, PyObject* args);
  PyObject* convexifyFaces(PyObject* self, PyObject* args);
  PyObject* closeOctalCells(PyObject* self, PyObject* args);
  
  PyObject* collapseUncomputableFaces(PyObject* self, PyObject* args);
  PyObject* removeNonManifoldExternalCells(PyObject* self, PyObject* args);
  
  PyObject* prepareCellsSplit(PyObject* self, PyObject* args);
  PyObject* splitNonStarCells(PyObject* self, PyObject* args);
  PyObject* simplifyCells(PyObject* self, PyObject* args);

  PyObject* agglomerateSmallCells(PyObject* self, PyObject* args);
  PyObject* agglomerateNonStarCells(PyObject* self, PyObject* args);
  //PyObject* agglomerateUncomputableCells(PyObject* self, PyObject* args);
  PyObject* agglomerateCellsWithSpecifiedFaces(PyObject* self, PyObject* args);

  PyObject* getOverlappingFaces(PyObject* self, PyObject* args);

  PyObject* adaptCells(PyObject* self, PyObject* args);
  PyObject* adaptBox(PyObject* self, PyObject* args);
  
  PyObject* extractUncomputables(PyObject* self, PyObject* args);
  PyObject* extractPathologicalCells(PyObject* self, PyObject* args);
  PyObject* extractOuterLayers(PyObject* self, PyObject* args);
  PyObject* extractNthCell(PyObject* self, PyObject* args);
  PyObject* extractNthFace(PyObject* self, PyObject* args);
  PyObject* removeNthCell(PyObject* self, PyObject* args);

  PyObject* statsUncomputableFaces(PyObject* self, PyObject* args);
  PyObject* statsSize(PyObject* self, PyObject* args);
  
  PyObject* computeAspectRatio(PyObject* self, PyObject* args);

  PyObject* centroids(PyObject* self, PyObject* args);

  PyObject* diffMesh(PyObject* self, PyObject* args);

  PyObject* checkCellsClosure(PyObject* self, PyObject* args);
  PyObject* checkForDegenCells(PyObject* self, PyObject* args);
  PyObject* edgeLengthExtrema(PyObject* self, PyObject* args);
  PyObject* reorientExternalFaces(PyObject* self, PyObject* args);
  PyObject* removeBaffles(PyObject* self, PyObject* args);

  PyObject* convert2Polyhedron(PyObject* self, PyObject* args);
  PyObject* oneZonePerCell(PyObject* self, PyObject* args);

  PyObject* extrudeUserDefinedBC(PyObject* self, PyObject* args);

  PyObject* convertNGON2DToNGON3D(PyObject* self, PyObject* args);  

  /////////// syncronizing the tree ///////////
  PyObject* updatePointLists(PyObject* self, PyObject* args);
  /////////////////////////////////////////////
  
  E_Int check_is_NGON(PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType);
  
  /*PyObject* conservative_transfer(PyObject* self, PyObject* args);
  PyObject* total_mass(PyObject* self, PyObject* args);
  PyObject* deltaMass(PyObject* self, PyObject* args);
  PyObject* normL1(PyObject* self, PyObject* args);*/
}
#endif
