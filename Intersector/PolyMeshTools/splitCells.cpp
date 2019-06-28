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


//#define FLAG_STEP
//#define DEBUG_2019
#ifdef DEBUG_2019
#include <chrono>
#include <ctime>
#endif

# include <string>
# include <sstream> 
# include "intersector.h"
# include "Fld/ngon_t.hxx"
# include "Nuga/Delaunay/Triangulator.h"
# include "Nuga/Boolean/Splitter.h"

#include "Nuga/include/tree.hxx"
#include "Nuga/include/geom_sensor.hxx"
#include "Nuga/include/xsensor.hxx"
#include "Nuga/include/adaptor.hxx"
#include "Nuga/include/hierarchical_mesh.hxx"

#include "Search/BbTree.h"

//#include <iostream>

using namespace std;
using namespace NUGA;

#ifdef FLAG_STEP
E_Int chrono::verbose = 0;
#endif

#define HX8 8
#define TH4 4
#define PYRA5 5
#define PENTA6 6
#define MIXED 9  // means that the mesh contains at least one basic-type element

typedef ngon_t<K_FLD::IntArray> ngon_type;

E_Int check_has_NGON_BASIC_ELEMENT(ngon_type& ng)
{
  E_Int s1(0), s2(0), s3(0), s4(0);  
  E_Int err = 0;
  for (E_Int i = 0; (i < ng.PHs.size()); ++i){
    if (K_MESH::Polyhedron<0>::is_HX8(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s1;
    else if (K_MESH::Polyhedron<0>::is_TH4(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s2;
    else if (K_MESH::Polyhedron<0>::is_PY5(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s3;
    else if (K_MESH::Polyhedron<0>::is_PR6(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s4;    
  }
#ifdef DEBUG_2019   
//  std::cout << "ng.PHs.size()= " << ng.PHs.size() << std::endl;  
//  std::cout << "s1= " << s1 << std::endl;
//  std::cout << "s2= " << s2 << std::endl;
//  std::cout << "s3= " << s3 << std::endl;
//  std::cout << "s4= " << s4 << std::endl;
#endif

  if (ng.PHs.size()==s1) return HX8;           // pure HX8 
  else if (ng.PHs.size()==s2) return TH4;      // pure TH4
  else if (ng.PHs.size()==s4) return PENTA6; // pure PENTA
  else if (s1+s2+s3+s4 > 0)                    // mixed basic (or Pyra which is mixed after first split in iso mode)
    return MIXED;
  else err=-1;

  if (err==-1)
  {
    //std::cout << "input error : err => " << err << std::endl;
    //std::cout << "input error : eltType => " << eltType << std::endl;
    PyErr_SetString(PyExc_TypeError, "input error : invalid array, must contain some basic elements (TETRA and HEXA currently) in NGON format.");
  }

  return err;
}


//=============================================================================
/* Agglomerate superfuous faces (overdefined polyhedra) */
//=============================================================================
PyObject* K_INTERSECTOR::splitNonStarCells(PyObject* self, PyObject* args)
{
  PyObject *arr(0);
  E_Float PH_conc_threshold(1./3.);
  E_Float PH_cvx_threshold(0.05);
  E_Float PG_cvx_threshold(1.e-8);

  if (!PYPARSETUPLEF(args, "Oddd", "Offf", &arr, &PH_conc_threshold, &PH_cvx_threshold, &PG_cvx_threshold)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  
  ngon_type ngi(cnt), ngo;
  
  NUGA::Splitter::split_non_star_phs<DELAUNAY::Triangulator>(crd, ngi, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold, ngo);

  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;

  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* Agglomerate superfuous faces (overdefined polyhedra) */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCells(PyObject* self, PyObject* args)
{
  PyObject *arr(nullptr), *arrS(nullptr);
  E_Int sensor_type(0), itermax(0), force_basic(0);
  

  if (!PYPARSETUPLEI(args, "OOlll", "OOiii", &arr, &arrS, &sensor_type, &itermax, &force_basic)) return NULL;

  K_FLD::FloatArray* f(nullptr), *fS(nullptr);
  K_FLD::IntArray* cn(nullptr), *cnS(nullptr);
  char* varString, *eltType, *varString2, *eltType2;
  // Check the mesh (NGON)
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err)
  {
    delete f; delete cn;
    return NULL;
  }

    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  E_Int ni, nj, nk;
  E_Int res2 = K_ARRAY::getFromArray(arrS, varString2, fS, ni, nj, nk,
                                     cnS, eltType2);
  
  K_FLD::FloatArray & crdS = *fS;
  const K_FLD::IntArray & cntS = *cnS;

  if (sensor_type == 1 /*xensor*/)
  {
    if ( (res2 != 2) || ((res2 == 2) && (strcmp(eltType2, "HEXA") != 0) ) )
    {
     PyErr_SetString(PyExc_ValueError,
       "adaptCells: xsensor currently support only HEXA mesh as source mesh.");
     delete f; delete cn;
     delete fS; delete cnS;
     return NULL;
    }
  }
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  PyObject *l(PyList_New(0));

  E_Int elt_type = check_has_NGON_BASIC_ELEMENT(ngi);
  if (elt_type==-1)
  {
    PyErr_SetString(PyExc_ValueError,
       "adaptCells: input mesh to adapt must have basic elements (Only Tets and Hexas are currently adapted) and must be in NGON format.");
    delete f; delete cn;
    delete fS; delete cnS;
    return NULL;
  }
  else if (elt_type==HX8 && !force_basic)
  {
    using mesh_type = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;

    mesh_type hmesh(crd, ngi);
  
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t sensor(hmesh, 1/*max_pts per cell*/, itermax);
      NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }
    else if (sensor_type == 1) //xsensor
    {
      using sensor_t = NUGA::xsensor<K_MESH::Hexahedron, mesh_type>;
      sensor_t sensor(hmesh, cntS, itermax);
      NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }
    else
    {
      using sensor_t = NUGA::geom_sensor2<mesh_type>;
      sensor_t sensor(hmesh, 1/*max_pts per cell*/, itermax);
      NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }
  
    std::vector<E_Int> oids;
    hmesh.conformize(oids);

    K_FLD::IntArray cnto;
    hmesh._ng.export_to_array(cnto);
  
    // pushing out the mesh
    PyObject *tpl = K_ARRAY::buildArray(hmesh._crd, varString, cnto, -1, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // pushing out PG history
    tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
 
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

  }
  else if (elt_type==TH4 && !force_basic)
  {
    using mesh_type = NUGA::hierarchical_mesh<K_MESH::Tetrahedron, NUGA::ISO>;
    mesh_type hmesh(crd, ngi);

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported
  
    //std::cout << "sensor_type : " << sensor_type << std::endl;
  
    if (sensor_type == 0)
    {
      //std::cout << "adapting..." << std::endl;
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t sensor(hmesh, 1/*max_pts per cell*/, itermax);
      NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }
    else if (sensor_type == 1) //xsensor
    {
      //using sensor_t = NUGA::xsensor<K_MESH::Tetrahedron, mesh_type>;
      //sensor_t sensor(hmesh, cntS, itermax);
      //NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }
    else
    {
      //std::cout << "adapting..." << std::endl;
      using sensor_t = NUGA::geom_sensor2<mesh_type>;
      sensor_t sensor(hmesh, 1/*max_pts per cell*/, itermax);
      NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }   

    std::vector<E_Int> oids;
    hmesh.conformize(oids);
    K_FLD::IntArray cnto;
    hmesh._ng.export_to_array(cnto);
      
    // pushing out the mesh
    PyObject *tpl = K_ARRAY::buildArray(hmesh._crd, varString, cnto, -1, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // pushing out PG history
    tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
 
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  else if (elt_type==PENTA6 && !force_basic)
  {
    using mesh_type = NUGA::hierarchical_mesh<K_MESH::Prism, NUGA::ISO>;
    mesh_type hmesh(crd, ngi);

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported
  
    //std::cout << "sensor_type : " << sensor_type << std::endl;
  
    if (sensor_type == 0)
    {
      std::cout << "force_basic" << std::endl;
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t sensor(hmesh, 1/*max_pts per cell*/, itermax);
      NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }
    else if (sensor_type == 1) //xsensor
    {
      //using sensor_t = NUGA::xsensor<K_MESH::Tetrahedron, mesh_type>;
      //sensor_t sensor(hmesh, cntS, itermax);
      //NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }
    else
    {
      //std::cout << "adapting..." << std::endl;
      using sensor_t = NUGA::geom_sensor2<mesh_type>;
      sensor_t sensor(hmesh, 1/*max_pts per cell*/, itermax);
      NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }

    std::vector<E_Int> oids;
    hmesh.conformize(oids);

    K_FLD::IntArray cnto;
    hmesh._ng.export_to_array(cnto);
  
    // pushing out the mesh
    PyObject *tpl = K_ARRAY::buildArray(hmesh._crd, varString, cnto, -1, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // pushing out PG history
    tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
 
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  else if (elt_type==MIXED || force_basic)
  {    
    using mesh_type = NUGA::hierarchical_mesh<K_MESH::Basic, NUGA::ISO>;
    mesh_type hmesh(crd, ngi);

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported
  
    //std::cout << "sensor_type : " << sensor_type << std::endl;
  
    if (sensor_type == 0)
    {
      //std::cout << "adapting..." << std::endl;
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t sensor(hmesh, 1/*max_pts per cell*/, itermax);
      NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    } 
    else if (sensor_type == 1) //xsensor
    {
      //using sensor_t = NUGA::xsensor<K_MESH::Basic, mesh_type>;
      //sensor_t sensor(hmesh, cntS, itermax);
      //NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }
    else
    {
      //std::cout << "adapting..." << std::endl;
      using sensor_t = NUGA::geom_sensor2<mesh_type>;
      sensor_t sensor(hmesh, 1/*max_pts per cell*/, itermax);
      NUGA::adaptor<mesh_type, sensor_t>::run(hmesh, sensor, crdS);
    }   

    std::vector<E_Int> oids;
    hmesh.conformize(oids);
    K_FLD::IntArray cnto;
    hmesh._ng.export_to_array(cnto);
  
    // pushing out the mesh
    PyObject *tpl = K_ARRAY::buildArray(hmesh._crd, varString, cnto, -1, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);

    // pushing out PG history
    tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
 
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  delete f; delete cn;
  delete fS; delete cnS;
  return l;
}

//=============================================================================
/* todo */
//=============================================================================
PyObject* K_INTERSECTOR::adaptBox(PyObject* self, PyObject* args)
{
  PyObject *arrS(nullptr);
  E_Float bratio(10.);
  E_Int itermax(0);
  

  if (!PYPARSETUPLE(args, "Odl", "Odi", "Ofl", "Ofi", &arrS, &bratio, &itermax)) return NULL;

  if (bratio < 1.)bratio = 1.;

  //std::cout << "in K_INTERSECTOR::adaptBox" << std::endl;


  K_FLD::FloatArray* f(nullptr);
  K_FLD::IntArray* cn(nullptr);
  char* varString, *eltType;

  E_Int ni, nj, nk;
  /*E_Int res2 = */K_ARRAY::getFromArray(arrS, varString, f, ni, nj, nk,
                                     cn, eltType);

  K_FLD::FloatArray & crdS = *f;
  //std::cout << "crd : " << crdS.cols() << "/" << crdS.rows() << std::endl;

  //std::cout << "compute box..." << std::endl;

  
  // Create the box
  K_SEARCH::BBox3D box;
  K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd(crdS);
  box.compute(acrd);

  // Box center and deltas
  E_Float C[3], Delta[3], DeltaMax=0.;
  for (int i=0; i < 3; ++i){
  	C[i] = 0.5*(box.minB[i] + box.maxB[i]);
  	Delta[i] = box.maxB[i] - box.minB[i];
  	DeltaMax = std::max(DeltaMax, Delta[i]);
  }
  // Enlarge it 
  for (int i=0; i < 3; ++i)
  {
  	box.minB[i] = C[i] - (bratio*DeltaMax/*Delta[i]*/);
  	box.maxB[i] = C[i] + (bratio*DeltaMax/*Delta[i]*/);
  }

  //std::cout << "convert box..." << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi;
  K_FLD::FloatArray crd;
  K_FLD::IntArray cnt;
  box.convert2NG<ngon_type>(crd, ngi);

  //std::cout << "adapt box..." << std::endl;

  using mesh_type = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;
  using sensor_type = NUGA::geom_sensor/*geom_static_sensor*/<mesh_type>;
  
  mesh_type hmesh(crd, ngi);
  sensor_type sensor(hmesh, 1/*max pts per cell*/, itermax);
  
  NUGA::adaptor<mesh_type, sensor_type>::run(hmesh, sensor, crdS);

  //std::cout << "output leaves..." << std::endl;
  std::vector<E_Int> oids;
  hmesh.conformize(oids);

  K_FLD::IntArray cnto;
  hmesh._ng.export_to_array(cnto);

  //std::cout << "output ..." << std::endl;
  
  PyObject* tpl = K_ARRAY::buildArray(hmesh._crd, varString, cnto, -1, "NGON", false);

  delete f; delete cn;
  return tpl;
}


//=======================  Intersector/PolyMeshTools/splitCells.cpp ====================
