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

#ifndef _INTERSECTOR_INTERSECTOR_H_
#define _INTERSECTOR_INTERSECTOR_H_

# include "kcore.h"

namespace K_INTERSECTOR
{
  PyObject* conformUnstr(PyObject* self, PyObject* args);
  
  PyObject* booleanIntersection(PyObject* self, PyObject* args);
  PyObject* booleanUnion(PyObject* self, PyObject* args);
  PyObject* booleanUnionMZ(PyObject* self, PyObject* args);
  PyObject* booleanMinus(PyObject* self, PyObject* args);
  PyObject* booleanIntersectionBorder(PyObject* self, PyObject* args);
  PyObject* booleanModifiedSolid(PyObject* self, PyObject* args);
  PyObject* DiffSurf(PyObject* self, PyObject* args);

  PyObject* XcellN(PyObject* self, PyObject* args);
  PyObject* superMesh(PyObject* self, PyObject* args);
  PyObject* superMeshCompSurf(PyObject* self, PyObject* args);
  PyObject* computeTNCFields(PyObject* self, PyObject* args);
  
  PyObject* P1ConservativeInterpolation(PyObject* self, PyObject* args);
  PyObject* P1ConservativeChimeraCoeffs(PyObject* self, PyObject* args);

  PyObject* selfX(PyObject* self, PyObject* args);
  PyObject* triangulateExteriorFaces(PyObject* self, PyObject* args);
  PyObject* triangulateSpecifiedFaces(PyObject* self, PyObject* args);
  PyObject* triangulateNFaces(PyObject* self, PyObject* args);
  PyObject* convexifyFaces(PyObject* self, PyObject* args);
  
  PyObject* closeCells(PyObject* self, PyObject* args);
  PyObject* closeCells_mpi(PyObject* self, PyObject* args);
  
  PyObject* replaceFaces(PyObject* self, PyObject* args);
  
  PyObject* collapseUncomputableFaces(PyObject* self, PyObject* args);
  PyObject* collapseSmallCells(PyObject* self, PyObject* args);
  PyObject* collapseSmallEdges(PyObject* self, PyObject* args);

  PyObject* removeNonManifoldExternalCells(PyObject* self, PyObject* args);
  
  PyObject* prepareCellsSplit(PyObject* self, PyObject* args);
  PyObject* splitNonStarCells(PyObject* self, PyObject* args);
  PyObject* simplifyCells(PyObject* self, PyObject* args);
  PyObject* simplifySurf(PyObject* self, PyObject* args);
  PyObject* simplifyFaces(PyObject* self, PyObject* args);
  PyObject* syncMacthPeriodicFaces(PyObject* self, PyObject* args);

  PyObject* agglomerateSmallCells(PyObject* self, PyObject* args);
  PyObject* shellAgglomerateSmallCells(PyObject* self, PyObject* args);
  PyObject* agglomerateNonStarCells(PyObject* self, PyObject* args);
  //PyObject* agglomerateUncomputableCells(PyObject* self, PyObject* args);
  PyObject* agglomerateCellsWithSpecifiedFaces(PyObject* self, PyObject* args);

  PyObject* immerseNodes(PyObject* self, PyObject* args);

  PyObject* getOverlappingFaces(PyObject* self, PyObject* args);
  PyObject* getCollidingTopFaces(PyObject* self, PyObject* args);
  PyObject* getCollidingCells(PyObject* self, PyObject* args);
  PyObject* getAnisoInnerFaces(PyObject* self, PyObject* args);

  PyObject* getFaceIdsWithCentroids(PyObject* self, PyObject* args);
  PyObject* getFaceIdsCollidingVertex(PyObject* self, PyObject* args);
  PyObject* getCells(PyObject* self, PyObject* args);
  PyObject* getFaces(PyObject* self, PyObject* args);

  PyObject* getNthNeighborhood(PyObject* self, PyObject* args);

  PyObject* estimateAdapReq(PyObject* self, PyObject* args);

  PyObject* adaptCells(PyObject* self, PyObject* args);
  PyObject* adaptCells_mpi(PyObject* self, PyObject* args);
  
  PyObject* adaptBox(PyObject* self, PyObject* args);

  PyObject* initForAdaptCells(PyObject* self, PyObject* args);
  PyObject* createHMesh(PyObject* self, PyObject* args);
  PyObject* deleteHMesh(PyObject* self, PyObject* args);
  
  PyObject* createSensor(PyObject* self, PyObject* args);
  PyObject* assignData2Sensor(PyObject* self, PyObject* args);
  PyObject* deleteSensor(PyObject* self, PyObject* args);
  
  PyObject* conformizeHMesh(PyObject* self, PyObject* args);
  PyObject* interpolateHMeshNodalField(PyObject* self, PyObject* args);
  
  PyObject* extractUncomputables(PyObject* self, PyObject* args);
  PyObject* extractPathologicalCells(PyObject* self, PyObject* args);
  PyObject* extractOuterLayers(PyObject* self, PyObject* args);
  PyObject* extractNthCell(PyObject* self, PyObject* args);
  PyObject* extractNthFace(PyObject* self, PyObject* args);
  PyObject* extractBiggestCell(PyObject* self, PyObject* args);
  PyObject* removeNthCell(PyObject* self, PyObject* args);
  PyObject* removeNthFace(PyObject* self, PyObject* args);
  PyObject* extractBadVolCells(PyObject* self, PyObject* args);
  PyObject* extractOverConnectedCells(PyObject* self, PyObject* args);

  PyObject* statsUncomputableFaces(PyObject* self, PyObject* args);
  PyObject* statsSize(PyObject* self, PyObject* args);
  
  PyObject* computeGrowthRatio(PyObject* self, PyObject* args);

  PyObject* centroids(PyObject* self, PyObject* args);
  PyObject* volumes(PyObject* self, PyObject* args);
  PyObject* volume(PyObject* self, PyObject* args);

  PyObject* diffMesh(PyObject* self, PyObject* args);

  PyObject* checkCellsClosure(PyObject* self, PyObject* args);
  PyObject* checkForDegenCells(PyObject* self, PyObject* args);
  PyObject* checkForBigCells(PyObject* self, PyObject* args);
  PyObject* checkCellsFlux(PyObject* self, PyObject* args);
  PyObject* checkCellsVolume(PyObject* self, PyObject* args);
  PyObject* checkCellsVolumeAndGrowthRatio(PyObject* self, PyObject* args);
  PyObject* checkAngularExtrema(PyObject* self, PyObject* args);
  PyObject* edgeLengthExtrema(PyObject* self, PyObject* args);
  PyObject* edgeLengthMax(PyObject* self, PyObject* args);
  
  PyObject* detectIdenticalCells(PyObject* self, PyObject* args);
  PyObject* detectOverConnectedFaces(PyObject* self, PyObject* args);
  
  PyObject* externalFaces(PyObject* self, PyObject* args);
  PyObject* reorient(PyObject* self, PyObject* args);
  PyObject* reorientSpecifiedFaces(PyObject* self, PyObject* args);

  PyObject* removeBaffles(PyObject* self, PyObject* args);

  PyObject* convert2Polyhedron(PyObject* self, PyObject* args);
  PyObject* oneZonePerCell(PyObject* self, PyObject* args);
  PyObject* oneZonePerFace(PyObject* self, PyObject* args);

  PyObject* extrudeBC(PyObject* self, PyObject* args);
  PyObject* extrudeSurf(PyObject* self, PyObject* args);
  PyObject* extrudeRevolSurf(PyObject* self, PyObject* args);

  PyObject* convertNGON2DToNGON3D(PyObject* self, PyObject* args);
  PyObject* convertBasic2NGONFaces(PyObject* self, PyObject* args);
  PyObject* oneph(PyObject* self, PyObject* args);  
  PyObject* drawOrientation(PyObject* self, PyObject* args);

  /////////// syncronizing the tree ///////////
  PyObject* updatePointLists(PyObject* self, PyObject* args);
  PyObject* exchangePointLists(PyObject* self, PyObject* args);
  PyObject* transposePointLists(PyObject* self, PyObject* args);
  /////////////////////////////////////////////

  PyObject* merge(PyObject* self, PyObject* args);
  PyObject* concatenate(PyObject* self, PyObject* args);
  
  E_Int check_is_of_type(const std::vector<std::string>& types, PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType);
  E_Int check_is_NGON(PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType);
  E_Int check_is_BAR(PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType);
  E_Int check_is_BASICF(PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType);
  
  enum eType { UNKN = -1, TETRA, PYRA, PRISM3, HEXA, /*PRISMN,*/ LAYER, BASIC };
  
  eType check_has_NGON_BASIC_ELEMENT(const K_FLD::IntArray & cnt);

  E_Int get_of_type(const std::vector<std::string>& types, PyObject* arr, K_FLD::FloatArray& f1, bool only_coords, K_FLD::IntArray& cn1, char*& varString, char*& eltType);
  E_Int getFromNGON(PyObject* arr, K_FLD::FloatArray& f1, bool only_coords, K_FLD::IntArray& cn1, char*& varString, char*& eltType);
  E_Int getFromBAR(PyObject* arr, K_FLD::FloatArray& f1, bool only_coords, K_FLD::IntArray& cn1, char*& varString, char*& eltType);

  PyObject* testmain(PyObject* self, PyObject* args);

  /*PyObject* conservative_transfer(PyObject* self, PyObject* args);
  PyObject* total_mass(PyObject* self, PyObject* args);
  PyObject* deltaMass(PyObject* self, PyObject* args);
  PyObject* normL1(PyObject* self, PyObject* args);*/

  PyObject* deleteCOM(PyObject*, PyObject*);

}
#endif
