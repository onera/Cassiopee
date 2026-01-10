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
#define K_ARRAY_UNIQUE_SYMBOL
#include "intersector.h"
#include <sstream>

#include "Nuga/include/Hexahedron.h"
#include "Nuga/include/Prism.h"
#include "Nuga/include/Pyramid.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyintersector [] =
{
  
  {"conformUnstr", K_INTERSECTOR::conformUnstr, METH_VARARGS},
  
  {"booleanIntersection", K_INTERSECTOR::booleanIntersection, METH_VARARGS},
  {"booleanUnion", K_INTERSECTOR::booleanUnion, METH_VARARGS},
  {"booleanUnionMZ", K_INTERSECTOR::booleanUnionMZ, METH_VARARGS},
  {"booleanMinus", K_INTERSECTOR::booleanMinus, METH_VARARGS},
  {"booleanIntersectionBorder", K_INTERSECTOR::booleanIntersectionBorder, METH_VARARGS},
  {"booleanModifiedSolid", K_INTERSECTOR::booleanModifiedSolid, METH_VARARGS},
  {"DiffSurf", K_INTERSECTOR::DiffSurf, METH_VARARGS},

  {"XcellN", K_INTERSECTOR::XcellN, METH_VARARGS},
  {"superMesh", K_INTERSECTOR::superMesh, METH_VARARGS},
  {"superMeshCompSurf", K_INTERSECTOR::superMeshCompSurf, METH_VARARGS},
  {"computeTNCFields", K_INTERSECTOR::computeTNCFields, METH_VARARGS},
  
  {"P1ConservativeInterpolation", K_INTERSECTOR::P1ConservativeInterpolation, METH_VARARGS},
  {"P1ConservativeChimeraCoeffs", K_INTERSECTOR::P1ConservativeChimeraCoeffs, METH_VARARGS},
  
  {"selfX", K_INTERSECTOR::selfX, METH_VARARGS},
  {"triangulateExteriorFaces", K_INTERSECTOR::triangulateExteriorFaces, METH_VARARGS},
  {"triangulateSpecifiedFaces", K_INTERSECTOR::triangulateSpecifiedFaces, METH_VARARGS},
  {"triangulateNFaces", K_INTERSECTOR::triangulateNFaces, METH_VARARGS},
  {"convexifyFaces", K_INTERSECTOR::convexifyFaces, METH_VARARGS},
  {"prepareCellsSplit", K_INTERSECTOR::prepareCellsSplit, METH_VARARGS},
  {"simplifyCells", K_INTERSECTOR::simplifyCells, METH_VARARGS},
  {"simplifySurf", K_INTERSECTOR::simplifySurf, METH_VARARGS},
  {"simplifyFaces", K_INTERSECTOR::simplifyFaces, METH_VARARGS},
  {"syncMacthPeriodicFaces", K_INTERSECTOR::syncMacthPeriodicFaces, METH_VARARGS},
  {"splitNonStarCells", K_INTERSECTOR::splitNonStarCells, METH_VARARGS},
  {"collapseUncomputableFaces", K_INTERSECTOR::collapseUncomputableFaces, METH_VARARGS},
  {"collapseSmallCells", K_INTERSECTOR::collapseSmallCells, METH_VARARGS},
  {"collapseSmallEdges", K_INTERSECTOR::collapseSmallEdges, METH_VARARGS },
  {"removeNonManifoldExternalCells", K_INTERSECTOR::removeNonManifoldExternalCells, METH_VARARGS},
  {"agglomerateSmallCells", K_INTERSECTOR::agglomerateSmallCells, METH_VARARGS},
  {"shellAgglomerateSmallCells", K_INTERSECTOR::shellAgglomerateSmallCells, METH_VARARGS},
  {"agglomerateNonStarCells", K_INTERSECTOR::agglomerateNonStarCells, METH_VARARGS},
  //{"agglomerateUncomputableCells", K_INTERSECTOR::agglomerateUncomputableCells, METH_VARARGS},
  {"immerseNodes", K_INTERSECTOR::immerseNodes, METH_VARARGS},
  {"agglomerateCellsWithSpecifiedFaces", K_INTERSECTOR::agglomerateCellsWithSpecifiedFaces, METH_VARARGS},
  
  {"adaptCells", K_INTERSECTOR::adaptCells, METH_VARARGS},
  {"adaptCells_mpi", K_INTERSECTOR::adaptCells_mpi, METH_VARARGS},
  
  {"adaptBox", K_INTERSECTOR::adaptBox, METH_VARARGS},
  
  {"initForAdaptCells", K_INTERSECTOR::initForAdaptCells, METH_VARARGS},
  {"createHMesh", K_INTERSECTOR::createHMesh, METH_VARARGS},
  {"deleteHMesh", K_INTERSECTOR::deleteHMesh, METH_VARARGS},
  {"conformizeHMesh", K_INTERSECTOR::conformizeHMesh, METH_VARARGS},
  
  {"createSensor", K_INTERSECTOR::createSensor, METH_VARARGS},
  {"deleteSensor", K_INTERSECTOR::deleteSensor, METH_VARARGS},
  {"assignData2Sensor", K_INTERSECTOR::assignData2Sensor, METH_VARARGS},
  
  {"interpolateHMeshNodalField", K_INTERSECTOR::interpolateHMeshNodalField, METH_VARARGS},
  
  {"closeCells", K_INTERSECTOR::closeCells, METH_VARARGS},
  {"closeCells_mpi", K_INTERSECTOR::closeCells_mpi, METH_VARARGS},

  {"replaceFaces", K_INTERSECTOR::replaceFaces, METH_VARARGS},
  
  {"extractUncomputables", K_INTERSECTOR::extractUncomputables, METH_VARARGS},
  {"extractPathologicalCells", K_INTERSECTOR::extractPathologicalCells, METH_VARARGS},
  {"extractOuterLayers", K_INTERSECTOR::extractOuterLayers, METH_VARARGS},
  {"extractNthCell", K_INTERSECTOR::extractNthCell, METH_VARARGS},
  {"extractNthFace", K_INTERSECTOR::extractNthFace, METH_VARARGS},
  {"extractBiggestCell", K_INTERSECTOR::extractBiggestCell, METH_VARARGS},
  {"removeNthCell", K_INTERSECTOR::removeNthCell, METH_VARARGS},
  {"removeNthFace", K_INTERSECTOR::removeNthFace, METH_VARARGS},
  {"extractBadVolCells", K_INTERSECTOR::extractBadVolCells, METH_VARARGS},
  {"extractOverConnectedCells", K_INTERSECTOR::extractOverConnectedCells, METH_VARARGS},

  {"getNthNeighborhood", K_INTERSECTOR::getNthNeighborhood, METH_VARARGS},

  {"getOverlappingFaces", K_INTERSECTOR::getOverlappingFaces, METH_VARARGS},
  {"getCollidingTopFaces", K_INTERSECTOR::getCollidingTopFaces, METH_VARARGS},
  {"getCollidingCells", K_INTERSECTOR::getCollidingCells, METH_VARARGS},
  {"getAnisoInnerFaces", K_INTERSECTOR::getAnisoInnerFaces, METH_VARARGS},
  {"estimateAdapReq", K_INTERSECTOR::estimateAdapReq, METH_VARARGS},
  {"getFaceIdsWithCentroids", K_INTERSECTOR::getFaceIdsWithCentroids, METH_VARARGS},
  {"getFaceIdsCollidingVertex", K_INTERSECTOR::getFaceIdsCollidingVertex, METH_VARARGS},
  {"getCells", K_INTERSECTOR::getCells, METH_VARARGS},
  {"getFaces", K_INTERSECTOR::getFaces, METH_VARARGS},

  {"statsUncomputableFaces", K_INTERSECTOR::statsUncomputableFaces, METH_VARARGS},
  {"statsSize", K_INTERSECTOR::statsSize, METH_VARARGS},
  
  {"computeGrowthRatio", K_INTERSECTOR::computeGrowthRatio, METH_VARARGS},
  {"centroids", K_INTERSECTOR::centroids, METH_VARARGS},
  {"volumes", K_INTERSECTOR::volumes, METH_VARARGS},
  {"volume", K_INTERSECTOR::volume, METH_VARARGS},
  
  {"diffMesh", K_INTERSECTOR::diffMesh, METH_VARARGS},

  { "checkCellsClosure", K_INTERSECTOR::checkCellsClosure, METH_VARARGS },
  { "checkForDegenCells", K_INTERSECTOR::checkForDegenCells, METH_VARARGS },
  { "checkForBigCells", K_INTERSECTOR::checkForBigCells, METH_VARARGS },
  { "checkCellsFlux", K_INTERSECTOR::checkCellsFlux, METH_VARARGS },
  { "checkCellsVolume", K_INTERSECTOR::checkCellsVolume, METH_VARARGS },
  { "checkCellsVolumeAndGrowthRatio", K_INTERSECTOR::checkCellsVolumeAndGrowthRatio, METH_VARARGS },
  { "checkAngularExtrema", K_INTERSECTOR::checkAngularExtrema, METH_VARARGS },
  { "detectIdenticalCells", K_INTERSECTOR::detectIdenticalCells, METH_VARARGS },
  { "detectOverConnectedFaces", K_INTERSECTOR::detectOverConnectedFaces, METH_VARARGS },
  { "edgeLengthExtrema", K_INTERSECTOR::edgeLengthExtrema, METH_VARARGS },
  { "edgeLengthMax", K_INTERSECTOR::edgeLengthMax, METH_VARARGS },
  { "removeBaffles", K_INTERSECTOR::removeBaffles, METH_VARARGS },
  { "convert2Polyhedron", K_INTERSECTOR::convert2Polyhedron, METH_VARARGS },
  { "oneZonePerCell", K_INTERSECTOR::oneZonePerCell, METH_VARARGS },
  { "oneZonePerFace", K_INTERSECTOR::oneZonePerFace, METH_VARARGS },
  
  { "extrudeBC", K_INTERSECTOR::extrudeBC, METH_VARARGS },
  { "extrudeSurf", K_INTERSECTOR::extrudeSurf, METH_VARARGS },
  { "extrudeRevolSurf", K_INTERSECTOR::extrudeRevolSurf, METH_VARARGS },

  { "externalFaces", K_INTERSECTOR::externalFaces, METH_VARARGS },
  { "reorient", K_INTERSECTOR::reorient, METH_VARARGS },
  { "reorientSpecifiedFaces", K_INTERSECTOR::reorientSpecifiedFaces, METH_VARARGS },

  { "convertNGON2DToNGON3D", K_INTERSECTOR::convertNGON2DToNGON3D, METH_VARARGS },
  { "convertBasic2NGONFaces", K_INTERSECTOR::convertBasic2NGONFaces, METH_VARARGS },
  { "oneph", K_INTERSECTOR::oneph, METH_VARARGS },
  { "drawOrientation", K_INTERSECTOR::drawOrientation, METH_VARARGS },


  /////////// syncronizing the tree ///////////
  { "updatePointLists", K_INTERSECTOR::updatePointLists, METH_VARARGS },
  { "exchangePointLists", K_INTERSECTOR::exchangePointLists, METH_VARARGS },
  { "transposePointLists", K_INTERSECTOR::transposePointLists, METH_VARARGS },
  /////////////////////////////////////////////
  { "merge", K_INTERSECTOR::merge, METH_VARARGS },
  { "concatenate", K_INTERSECTOR::concatenate, METH_VARARGS },

  { "testmain", K_INTERSECTOR::testmain, METH_VARARGS },

  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "intersector",
        NULL,
        -1,
        Pyintersector
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_intersector();
  PyMODINIT_FUNC PyInit_intersector()
#else
  PyMODINIT_FUNC initintersector();
  PyMODINIT_FUNC initintersector()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("intersector", Pyintersector);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

E_Int K_INTERSECTOR::check_is_of_type(const std::vector<std::string>& types, PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType)
{
  E_Int ni, nj, nk;
  
  E_Int res = K_ARRAY::getFromArray(arr, varString, f1, ni, nj, nk,
                                    cn1, eltType);

  //std::cout << "eltType ???????" << eltType << std::endl;
     
  bool err = (res !=2);

  for (size_t i=0; (i < types.size()) && !err; ++i)
    err &= (strcmp(eltType, types[i].c_str()) != 0);

  if (err)
  {
    std::stringstream o;
    o << "input error : " << eltType << " is an invalid array, must be a ";
    for (size_t i=0; i < types.size()-1; ++i){
      o << types[i];
      if (i < types.size()-2) o << ", ";
    }
    o << " or " << types[types.size()-1] << " array." ;
    PyErr_SetString(PyExc_TypeError, o.str().c_str());//fixme triangulateExteriorFaces : PASS A STRING AS INPUT
    //delete f1; delete cn1;
    //f1 = nullptr; cn1 = nullptr;
    return 1;
  }

  // Check coordinates.
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if ((posx == -1) || (posy == -1) || (posz == -1))
  {
    PyErr_SetString(PyExc_TypeError, "input error : can't find coordinates in array.");//fixme  conformUnstr
    //delete f1; delete cn1;
    //f1 = nullptr; cn1 = nullptr;
    return 1;
  }
  
  return 0;
}

E_Int K_INTERSECTOR::get_of_type
(const std::vector<std::string>& types, PyObject* arr, K_FLD::FloatArray& f1, bool only_coords, K_FLD::IntArray& cn1, char*& varString, char*& eltType)
{
  E_Int ni, nj, nk;
  
  E_Int res = K_ARRAY::getFromArray(arr, varString, f1, ni, nj, nk, cn1, eltType);

  //std::cout << "eltType ???????" << eltType << std::endl;
     
  bool err = (res != 2);

  if (!err)
  {
    std::set<std::string> stypes(types.begin(), types.end());
    std::string selt(eltType);
    err = (stypes.find(selt) == stypes.end()); // eltType is not in the input list
  }

  if (err)
  {
    std::stringstream o;
    o << "input error : invalid array, must be a ";
    for (size_t i=0; i < types.size()-1; ++i){
      o << types[i];
      if (i < types.size() - 2) o << ", ";
    }
    if (types.size() > 1) o << " or ";
    o << types[types.size()-1] << " array." ;
    PyErr_SetString(PyExc_TypeError, o.str().c_str());//fixme triangulateExteriorFaces : PASS A STRING AS INPUT
    //delete f1; delete cn1;
    //f1 = nullptr; cn1 = nullptr;
    return 1;
  }

  // Check coordinates.
  E_Int pos[3];
  pos[0] = K_ARRAY::isCoordinateXPresent(varString);
  pos[1] = K_ARRAY::isCoordinateYPresent(varString);
  pos[2] = K_ARRAY::isCoordinateZPresent(varString);

  if ((pos[0] == -1) || (pos[1] == -1) || (pos[2] == -1))
  {
    PyErr_SetString(PyExc_TypeError, "input error : can't find coordinates in array.");//fixme  conformUnstr
    //delete f1; delete cn1;
    //f1 = nullptr; cn1 = nullptr;
    return 1;
  }

  if (pos[0] == 0 && pos[1] == 1 && pos[2] == 2) only_coords = false; //nothing to do

  if (only_coords)
  {
    E_Int npts = f1.cols();
    K_FLD::FloatArray crd(3, npts);
    for (E_Int i = 0; i < npts; ++i)
    {
      for (size_t k=0; k < 3; ++k)
        crd(k, i) = f1(pos[k], i);
    }
    f1 = std::move(crd);
    //std::cout << "COORDS : " << f1.rows() << std::endl;
  }
  //else
  	//std::cout << " U COORDS : " << f1.rows() << std::endl;
  
  return 0;
}

E_Int K_INTERSECTOR::check_is_NGON(PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType)
{
  std::vector<std::string> types;
  types.push_back(std::string("NGON"));
  return check_is_of_type(types, arr, f1, cn1, varString, eltType);
}

E_Int K_INTERSECTOR::getFromNGON(PyObject* arr, K_FLD::FloatArray& f1, bool only_coords, K_FLD::IntArray& cn1, char*& varString, char*& eltType)
{
  std::vector<std::string> types;
  types.push_back(std::string("NGON"));
  return get_of_type(types, arr, f1, only_coords, cn1, varString, eltType);
}

E_Int K_INTERSECTOR::check_is_BAR(PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType)
{
  std::vector<std::string> types;
  types.push_back(std::string("BAR"));
  return check_is_of_type(types, arr, f1, cn1, varString, eltType);
}

E_Int K_INTERSECTOR::getFromBAR(PyObject* arr, K_FLD::FloatArray& f1, bool only_coords, K_FLD::IntArray& cn1, char*& varString, char*& eltType)
{
  std::vector<std::string> types;
  types.push_back(std::string("BAR"));
  return get_of_type(types, arr, f1, only_coords, cn1, varString, eltType);
}

E_Int K_INTERSECTOR::check_is_BASICF(PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType)
{
  std::vector<std::string> types;
  types.push_back(std::string("TRI"));
  types.push_back(std::string("QUAD"));

  return check_is_of_type(types, arr, f1, cn1, varString, eltType);
}


K_INTERSECTOR::eType K_INTERSECTOR::check_has_NGON_BASIC_ELEMENT(const K_FLD::IntArray & cnt)
{
  using ngon_type = ngon_t<K_FLD::IntArray>;
  ngon_type ng(cnt); //fixme: temporary hack
  E_Int s1(0), s2(0), s3(0), s4(0);  
  //E_Int err = 0;
  for (E_Int i = 0; (i < ng.PHs.size()); ++i){
    if (K_MESH::Hexahedron::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s1;
    else if (K_MESH::Tetrahedron::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s2;
    else if (K_MESH::Pyramid::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s3;
    else if (K_MESH::Prism::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s4;    
  }

  if (ng.PHs.size()==s1)      return eType::HEXA; 
  else if (ng.PHs.size()==s2) return eType::TETRA;
  else if (ng.PHs.size()==s4) return eType::PRISM3;
  else if (s1+s2+s3+s4 > 0)   return eType::BASIC;
  else return eType::UNKN;
}
