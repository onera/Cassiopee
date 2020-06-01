/*    
    Copyright 2013-2020 Onera.

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
#include "Nuga/include/nodal_sensor.hxx"
#include "Nuga/include/adaptor.hxx"
#include "Nuga/include/hierarchical_mesh.hxx"

#include "Search/BbTree.h"
#include <memory>
//#include <iostream>

using namespace std;
using namespace NUGA;

#ifdef FLAG_STEP
E_Int chrono::verbose = 0;
#endif

#define HMESH_HOOK_ID 77
#define SENSOR_HOOK_ID 78

using ngon_type = ngon_t<K_FLD::IntArray>;
using elt_t = K_MESH::Polyhedron<0>::eType;
using subdiv_t = NUGA::eSUBDIV_TYPE;

elt_t check_has_NGON_BASIC_ELEMENT(const K_FLD::IntArray & cnt)
{
  ngon_type ng(cnt); //fixme: temporary hack
  E_Int s1(0), s2(0), s3(0), s4(0);  
  E_Int err = 0;
  for (E_Int i = 0; (i < ng.PHs.size()); ++i){
    if (K_MESH::Hexahedron::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s1;
    else if (K_MESH::Tetrahedron::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s2;
    else if (K_MESH::Pyramid::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s3;
    else if (K_MESH::Prism::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s4;    
  }
#ifdef DEBUG_2019   
//  std::cout << "ng.PHs.size()= " << ng.PHs.size() << std::endl;  
//  std::cout << "s1= " << s1 << std::endl;
//  std::cout << "s2= " << s2 << std::endl;
//  std::cout << "s3= " << s3 << std::endl;
//  std::cout << "s4= " << s4 << std::endl;
#endif

  if (ng.PHs.size()==s1)      return elt_t::HEXA; 
  else if (ng.PHs.size()==s2) return elt_t::TETRA;
  else if (ng.PHs.size()==s4) return elt_t::PRISM3;
  else if (s1+s2+s3+s4 > 0)   return elt_t::BASIC;
  else return elt_t::UNKN;
}

template <NUGA::eSUBDIV_TYPE STYPE>
void* __createHM(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt);

// ISO strategy
template<>
void* __createHM<NUGA::ISO>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt)
{
  if (typ == elt_t::UNKN)
  {
    PyErr_SetString(PyExc_ValueError,
      "createHMesh: input mesh to adapt must have basic elements and must be in NGON format.");
    return nullptr;
  }
  else if (typ == elt_t::HEXA)
  {
    using mesh_type = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;
    return new mesh_type(crd, cnt);
  }
  else if (typ == elt_t::TETRA)
  {
    using mesh_type = NUGA::hierarchical_mesh<K_MESH::Tetrahedron, NUGA::ISO>;
    return new mesh_type(crd, cnt);
  }
  else if (typ == elt_t::PRISM3)
  {
    using mesh_type = NUGA::hierarchical_mesh<K_MESH::Prism, NUGA::ISO>;
    return new mesh_type(crd, cnt);
  }
  else if (typ == elt_t::BASIC)
  {
    using mesh_type = NUGA::hierarchical_mesh<K_MESH::Basic, NUGA::ISO>;
    return new mesh_type(crd, cnt);
  }
  return NULL;
}

// ISO_HEX strategy
template<>
void* __createHM<NUGA::ISO_HEX>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt)
{
  using mesh_type = NUGA::hierarchical_mesh<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>;
  return new mesh_type(crd, cnt);
}

//============================================================================
/* Create a hmesh and returns a hook */
//============================================================================
PyObject* K_INTERSECTOR::createHMesh(PyObject* self, PyObject* args)
{
  PyObject* hook;
  void** packet = new void*[5];  // hook_ID, hmesh ptr, subdiv policy, elt type, varString

  E_Int* hookid = new E_Int;  packet[0] = hookid;
  //void* hmesh_ptr = nullptr;  packet[1] = hmesh_ptr;// templated hmesh type to build 
  E_Int* subtype = new E_Int; packet[2] = subtype;
  elt_t* etyp = new elt_t;    packet[3] = etyp;
  std::string* vString = new std::string; packet[4] = vString;

  *hookid = HMESH_HOOK_ID;

  PyObject *arr;
  E_Int SUBDIV_TYPE(0);/*ISO*/
  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arr, &SUBDIV_TYPE)) return nullptr;

  *subtype = SUBDIV_TYPE; // 0 : ISO, 1 :ISO_HEX ...

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char *eltType, *varString;
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  *vString = varString;

  std::unique_ptr<K_FLD::FloatArray> afmesh(f);  // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<K_FLD::IntArray> acmesh(cn); // to avoid to call explicit delete at several places in the code.
  if (err) return nullptr;
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  *etyp = elt_t::UNKN;// Polyhedron
  if (SUBDIV_TYPE == 0) //ISO
    *etyp = check_has_NGON_BASIC_ELEMENT(cnt);

  if (SUBDIV_TYPE == 0) //ISO
    packet[1] = __createHM<NUGA::ISO>(*etyp, crd, cnt);
  else if (SUBDIV_TYPE == 1) //ISO_HEX
    packet[1] = __createHM<NUGA::ISO_HEX>(*etyp, crd, cnt);

  if (packet[1] == nullptr) return nullptr;

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif
  
  return hook;
}

//=============================================================================
/* 
   Confomize a hmesh
 */
//=============================================================================
template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void __conformizeHM(const void* hmesh_ptr, K_FLD::FloatArray*& crdo, K_FLD::IntArray& cnto, std::vector<E_Int>& oids)
{
  using mesh_type = NUGA::hierarchical_mesh<ELT_t, STYPE>;

  mesh_type* hmesh = (mesh_type*)hmesh_ptr;

  ngon_type ngo;
  hmesh->conformize(ngo, oids);
  ngo.export_to_array(cnto);

  crdo = &hmesh->_crd;
}

template <NUGA::eSUBDIV_TYPE STYPE>
void __conformizeHM(E_Int etype, const void* hmesh_ptr, K_FLD::FloatArray*& crdo, K_FLD::IntArray& cnto, std::vector<E_Int>& oids)
{
  if (etype == elt_t::HEXA)
    __conformizeHM<K_MESH::Hexahedron, STYPE>(hmesh_ptr, crdo, cnto, oids);
  else if (etype == (E_Int)elt_t::TETRA)
    __conformizeHM<K_MESH::Tetrahedron, STYPE>(hmesh_ptr, crdo, cnto, oids);
  else if (etype == (E_Int)elt_t::PRISM3)
    __conformizeHM<K_MESH::Prism, STYPE>(hmesh_ptr, crdo, cnto, oids);
  else if (etype == (E_Int)elt_t::BASIC)
    __conformizeHM<K_MESH::Basic, STYPE>(hmesh_ptr, crdo, cnto, oids);
}

template <>
void __conformizeHM<NUGA::ISO_HEX>(E_Int etype/*dummy*/, const void* hmesh_ptr, K_FLD::FloatArray*& crdo, K_FLD::IntArray& cnto, std::vector<E_Int>& oids)
{
  __conformizeHM<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh_ptr, crdo, cnto, oids);
}

//
PyObject* K_INTERSECTOR::conformizeHMesh(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PyArg_ParseTuple(args, "O", &hook)) return nullptr;

  //todo CW : utiliser unpackHMesh

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0]; // type of hook
  E_Int* sub_type = (E_Int*)packet[2]; // subdivision type ISO, ISO_HEX, DIR...
  E_Int* elt_type = (E_Int*)packet[3]; // type of elements in hmesh
  std::string* vString = (std::string*)packet[4];
  assert (*type == HMESH_HOOK_ID);

  
  PyObject *l(PyList_New(0));
  K_FLD::IntArray cnto;
  std::vector<E_Int> oids;
  K_FLD::FloatArray* crd(nullptr);

  if (*sub_type == 0) // ISO
    __conformizeHM<NUGA::ISO>(*elt_type, packet[1], crd, cnto, oids);
  else if (*sub_type == 1) // ISO_HEX
    __conformizeHM<NUGA::ISO_HEX>(*elt_type, packet[1], crd, cnto, oids);

  // pushing out the mesh
  PyObject *tpl = K_ARRAY::buildArray(*crd, vString->c_str(), cnto, -1, "NGON", false);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // pushing out PG history
  tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  return l;
}

//============================================================================
/* Deletes a hmesh */
//============================================================================
template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void __deleteHM(const void* hmesh_ptrs)
{
  using mesh_type = NUGA::hierarchical_mesh<ELT_t, STYPE>;
  mesh_type* hmesh = (mesh_type*)hmesh_ptrs;
  delete hmesh;
}

template <NUGA::eSUBDIV_TYPE STYPE>
void __deleteHM(E_Int etype, const void* hmesh_ptr)
{
  if (etype == elt_t::HEXA)
    __deleteHM<K_MESH::Hexahedron, STYPE>(hmesh_ptr);
  else if (etype == (E_Int)elt_t::TETRA)
    __deleteHM<K_MESH::Tetrahedron, STYPE>(hmesh_ptr);
  else if (etype == (E_Int)elt_t::PRISM3)
    __deleteHM<K_MESH::Prism, STYPE>(hmesh_ptr);
  else if (etype == (E_Int)elt_t::BASIC)
    __deleteHM<K_MESH::Basic, STYPE>(hmesh_ptr);
}

template <>
void __deleteHM<NUGA::ISO_HEX>(E_Int etype/*dummy*/, const void* hmesh_ptr)
{
  __deleteHM<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh_ptr);
}

PyObject* K_INTERSECTOR::deleteHMesh(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PyArg_ParseTuple(args, "O", &hook))
  {
      return NULL;
  }

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0]; // type of hook
  E_Int* sub_type = (E_Int*)packet[2]; // subdivision type ISO or DIR
  E_Int* elt_type = (E_Int*)packet[3]; // type of elements in hmesh
  std::string* vString = (std::string*)packet[4];
  assert (*type == HMESH_HOOK_ID);

  if (*sub_type == 0) // ISO
    __deleteHM<NUGA::ISO>(*elt_type, packet[1]);
  else if (*sub_type == 1) // ISO_HEX
    __deleteHM<NUGA::ISO_HEX>(*elt_type, packet[1]);
  
  delete type; delete vString;
  delete sub_type; delete elt_type;
  delete [] packet;
  Py_INCREF(Py_None);
  return Py_None;
}

//============================================================================
/* Create a Geom sensor and returns a hook */
//============================================================================

template<typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void* __createGeomSensor(void* hmesh, E_Int smoothing_type, E_Int itermax)
{
  using hmesh_t = NUGA::hierarchical_mesh<ELT_t, STYPE>;
  using sensor_t = NUGA::geom_sensor<hmesh_t>;

  return new sensor_t(*(hmesh_t*)hmesh, NUGA::eSmoother(smoothing_type), 1/*max_pts_per_cell*/, itermax);
}


 template<typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
 void* __createXSensor(void* hmesh, E_Int smoothing_type, E_Int itermax)
 {
   //todo CW
 }

template<NUGA::eSUBDIV_TYPE STYPE>
void* __createGeomSensor(E_Int etype, void* hmesh, E_Int sensor_type, E_Int smoothing_type, E_Int itermax);

template<>
void* __createGeomSensor<NUGA::ISO>(E_Int elt_type, void* hmesh, E_Int sensor_type, E_Int smoothing_type, E_Int itermax)
{
  if (elt_type == elt_t::HEXA)
  {
    if (sensor_type == 1)
      return __createXSensor<K_MESH::Hexahedron, NUGA::ISO>(hmesh, smoothing_type, itermax);
    else
      return __createGeomSensor<K_MESH::Hexahedron, NUGA::ISO>(hmesh, smoothing_type, itermax);
  }
  else if (elt_type == elt_t::TETRA)
    return __createGeomSensor<K_MESH::Tetrahedron, NUGA::ISO>(hmesh, smoothing_type, itermax);
  else if (elt_type == elt_t::PRISM3)
    return __createGeomSensor<K_MESH::Prism, NUGA::ISO>(hmesh, smoothing_type, itermax);
  else if (elt_type == elt_t::BASIC)
    return __createGeomSensor<K_MESH::Basic, NUGA::ISO>(hmesh, smoothing_type, itermax);
  else 
  {
    PyErr_SetString(PyExc_ValueError, "createSensor: wrong element type in hmesh for ISO strategy.");
    return nullptr;
  }
}

//
template<>
void* __createGeomSensor<NUGA::ISO_HEX>(E_Int etype, void* hmesh, E_Int sensor_type, E_Int smoothing_type, E_Int itermax)
{
  return __createGeomSensor<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh, smoothing_type, itermax);
}

//============================================================================
/* Creates a Geom Sensor */
//============================================================================
PyObject* K_INTERSECTOR::createGeomSensor(PyObject* self, PyObject* args)
{
  PyObject *hook_sensor;

  PyObject *hook_hmesh(nullptr);
  E_Int smoothing_type(0), sensor_type(0), itermax(0); // sensor_type = 0(geom_sensor) or 1 (xsensor) 

  if (!PYPARSETUPLEI(args, "Olll", "Oiii", &hook_hmesh, &sensor_type, &smoothing_type, &itermax)) return NULL;


  // Unpack hmesh hook
  // ==================
  void** packet_hm = NULL;
  void * hmesh = nullptr;
  elt_t elt_type(elt_t::UNKN);

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet_hm = (void**)PyCObject_AsVoidPtr(hook_hmesh);
#else
  packet_hm = (void**)PyCapsule_GetPointer(hook_hmesh, NULL);
#endif

  E_Int              hook_hm_id = *((E_Int*)(packet_hm[0]));
  hmesh = (void*)packet_hm[1];
  NUGA::eSUBDIV_TYPE subdiv_type = *((NUGA::eSUBDIV_TYPE*)(packet_hm[2]));
  elt_type = *((elt_t*)(packet_hm[3]));

  if (hook_hm_id != HMESH_HOOK_ID)
  {
    PyErr_SetString(PyExc_TypeError,
      "createSensor: this function requires a identify hmesh hook.");
    return NULL;
  }

  // Create packet for sensor hook 
  // ==============================
  void** packet_ss = new void*[4]; // hook id, sensor type, hook, smoothing

  E_Int* hook_sensor_id = new E_Int;   packet_ss[0] = hook_sensor_id;
  E_Int* psens_type = new E_Int;       packet_ss[1] = psens_type;
                                       packet_ss[2] = nullptr;
  E_Int* psmooth_type = new E_Int;     packet_ss[3] = psmooth_type;
  
  *hook_sensor_id = SENSOR_HOOK_ID;
  *psens_type = sensor_type;
  *psmooth_type = smoothing_type;

  if (subdiv_type == 0) // ISO
    packet_ss[2] = __createGeomSensor<NUGA::ISO>(elt_type, hmesh, sensor_type, smoothing_type, itermax);
  else if (subdiv_type == 1) // ISO_HEX
    packet_ss[2] = __createGeomSensor<NUGA::ISO_HEX>(elt_type, hmesh, sensor_type, smoothing_type, itermax);
  
  // Create sensor hook 
  // ===================
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook_sensor = PyCObject_FromVoidPtr(packet_ss, NULL);
#else
  hook_sensor = PyCapsule_New(packet_ss, NULL, NULL);
#endif

  return hook_sensor;
}


// //============================================================================
// /* Deletes a hsensor */
// //============================================================================
// PyObject* K_INTERSECTOR::deleteSensor(PyObject* self, PyObject* args)
// {
//   PyObject* hook_sensor;
//   if (!PyArg_ParseTuple(args, "O", &hook_sensor))
//   {
//       return NULL;
//   }

//   // recupere le hook
//   void** packet_ss = NULL;

//   #if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
//     packet_ss = (void**) PyCObject_AsVoidPtr(hook_sensor);
//   #else
//     packet_ss = (void**) PyCapsule_GetPointer(hook_sensor, NULL);
//   #endif

//   E_Int hook_ss_id  = *((E_Int*)(packet_ss[0]));

//   if (hook_ss_id != SENSOR_HOOK_ID)
//   {
//     PyErr_SetString(PyExc_TypeError, 
//                     "unpackSensor: this function requires a identify sensor hook.");
//     return NULL;
//   }   

//   sensor_type = *((E_Int*)(packet_ss[1]));


//   // delete type; delete varString;
//   // delete sub_type; delete elt_type;
//   // delete [] packet;
//   // Py_INCREF(Py_None);
//   // return Py_None;
// }

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

template <typename MESH_t, typename sensor_t>
E_Int __adapt
(MESH_t* hmesh, bool ishook,
 sensor_t& sensor, typename sensor_t::sensor_data_t& data,
 const char* varString, PyObject *out)
{
  hmesh->init();

  E_Int err(0);
  err = sensor.assign_data(data);

  NUGA::adaptor<MESH_t, sensor_t>::run(*hmesh, sensor);

  std::vector<E_Int> oids;
  K_FLD::IntArray cnto;

  if (!ishook)
  {
    ngon_type ngo;
    hmesh->conformize(ngo, oids);
    ngo.export_to_array(cnto);
  }
  else
    hmesh->_ng.export_to_array(cnto);

  // pushing out the mesh
  PyObject *tpl = K_ARRAY::buildArray(hmesh->_crd, varString, cnto, -1, "NGON", false);
  PyList_Append(out, tpl);
  Py_DECREF(tpl);

  // pushing out PG history
  if (!ishook)
  {
    tpl = K_NUMPY::buildNumpyArray(&oids[0], oids.size(), 1, 0);
    PyList_Append(out, tpl);
    Py_DECREF(tpl);
  }

  if (!ishook) delete hmesh;
  return 0;
}

template <typename MESH_t>
MESH_t* get_hmesh(K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, void* hook_hmesh)
{
  MESH_t * hmesh = nullptr;
  if (hook_hmesh != nullptr)
    hmesh = (MESH_t*)hook_hmesh;
  else
    hmesh = new MESH_t(crd, cnt);
  return hmesh;
}

template <subdiv_t STYPE>
E_Int __adapt_wrapper
(K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int elt_type, E_Int sensor_type, E_Int smoothing_type,
  K_FLD::FloatArray& gdata, std::vector<E_Int>& ndata,
  E_Int itermax, const char* varString, PyObject *out, void* hookhm, const K_FLD::IntArray* cntS = nullptr);

///
template <>
E_Int __adapt_wrapper<NUGA::ISO>
(K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int elt_type, E_Int sensor_type, E_Int smoothing_type, 
 K_FLD::FloatArray& gdata, std::vector<E_Int>& ndata, 
 E_Int itermax, const char* varString, PyObject *out, void* hookhm, const K_FLD::IntArray* cntS)
{
  //std::cout << "sensor_type : " << sensor_type << std::endl;

  bool is_hook = (hookhm != nullptr);

  NUGA::eSmoother smoo_typ = (smoothing_type == 0) ? NUGA::eSmoother::V1_NEIGH : NUGA::eSmoother::SHELL;

  E_Int err(0);
  if (elt_type==elt_t::HEXA)
  {
    using ELT_type = K_MESH::Hexahedron;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, NUGA::ISO>;

    mesh_type* hmesh = get_hmesh<mesh_type>(crd, cnt, hookhm); //sync the hook or build new hmesh
  
    if (sensor_type == 0) // geom sensor
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;

      sensor_t sensor(*hmesh, smoo_typ, 1/*max_pts per cell*/, itermax);

      err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, gdata, varString, out);
    }
    else if (sensor_type == 1) // xsensor
    {
      using sensor_t = NUGA::xsensor<ELT_type, mesh_type>;
      sensor_t sensor(*hmesh, smoo_typ, 1/*max_pts per cell*/, itermax, cntS);

      err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, gdata, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;
      sensor_t sensor(*hmesh);

      err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, ndata, varString, out);
    }
  }
  else if (elt_type==elt_t::TETRA)
  {
    using ELT_type = K_MESH::Tetrahedron;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, NUGA::ISO>;

    mesh_type* hmesh = get_hmesh<mesh_type>(crd, cnt, hookhm); //sync the hook or build new one

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported
    
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t sensor(*hmesh, smoo_typ, 1/*max_pts per cell*/, itermax);

      err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, gdata, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;
      sensor_t sensor(*hmesh);

      err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, ndata, varString, out);
    }
  }
  else if (elt_type==elt_t::PRISM3)
  {
    using ELT_type = K_MESH::Prism;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, NUGA::ISO>;

    mesh_type* hmesh = get_hmesh<mesh_type>(crd, cnt, hookhm); //sync the hook or build new one

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported
    
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t sensor(*hmesh, smoo_typ, 1/*max_pts per cell*/, itermax);

      err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, gdata, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;
      sensor_t sensor(*hmesh);

      err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, ndata, varString, out);
    }
  }
  else if (elt_type==elt_t::BASIC)
  {    
    using ELT_type = K_MESH::Basic;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, NUGA::ISO>;

    mesh_type* hmesh = get_hmesh<mesh_type>(crd, cnt, hookhm); //sync the hook or build new one

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported

    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t sensor(*hmesh, smoo_typ, 1/*max_pts per cell*/, itermax);

      err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, gdata, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;
      sensor_t sensor(*hmesh);

      err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, ndata, varString, out);
    }
  }
  return err;
}

///
template <>
E_Int __adapt_wrapper<NUGA::ISO_HEX>
(K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int elt_type, E_Int sensor_type, E_Int smoothing_type,
  K_FLD::FloatArray& gdata, std::vector<E_Int>& ndata,
  E_Int itermax, const char* varString, PyObject *out, void* hookhm, const K_FLD::IntArray* cntS)
{
  //std::cout << "sensor_type : " << sensor_type << std::endl;

  bool is_hook = (hookhm != nullptr);

  NUGA::eSmoother smoo_typ = (smoothing_type == 0) ? NUGA::eSmoother::V1_NEIGH : NUGA::eSmoother::SHELL;

  E_Int err(0);
  using ELT_type = K_MESH::Polyhedron<0>;
  using mesh_type = NUGA::hierarchical_mesh<ELT_type, NUGA::ISO_HEX>;

  mesh_type* hmesh = get_hmesh<mesh_type>(crd, cnt, hookhm); //sync the hook or build new hmesh

  if (sensor_type == 0) // geom sensor
  {
    using sensor_t = NUGA::geom_sensor<mesh_type>;

    sensor_t sensor(*hmesh, smoo_typ, 1/*max_pts per cell*/, itermax);

    err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, gdata, varString, out);
  }
  //else if (sensor_type == 1) // xsensor
  //{
  //  using sensor_t = NUGA::xsensor<ELT_type, mesh_type>;
  //  sensor_t sensor(*hmesh, smoo_typ, 1/*max_pts per cell*/, itermax, cntS);

  //  err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, gdata, varString, out);
  //}
  else if (sensor_type == 2) //nodal sensor
  {
    using sensor_t = NUGA::nodal_sensor<mesh_type>;
    sensor_t sensor(*hmesh);

    err = __adapt<mesh_type, sensor_t>(hmesh, is_hook, sensor, ndata, varString, out);
  }
  
  return err;
}

//=============================================================================
/* Agglomerate superfuous faces (overdefined polyhedra) */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCells(PyObject* self, PyObject* args)
{
  PyObject *arr(nullptr), *arrS(nullptr), *hook_hmesh(nullptr);
  E_Int sensor_type(0), smoothing_type(0/*face*/), itermax(0), subdiv_type(0);

  //
  if (!PYPARSETUPLEI(args, "OOllllO", "OOiiiiO", &arr, &arrS, &sensor_type, &smoothing_type, &itermax, &subdiv_type, &hook_hmesh)) return NULL;

  NUGA::eSUBDIV_TYPE sub_type = subdiv_t(subdiv_type);

  K_FLD::FloatArray* f(nullptr), *fS(nullptr);
  K_FLD::IntArray* cn(nullptr), *cnS(nullptr);
  char* varString, *eltType, *varString2, *eltType2;
  // Check the mesh (NGON)
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);

  std::unique_ptr<K_FLD::FloatArray> afmesh(f);  // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<K_FLD::IntArray> acmesh(cn); // to avoid to call explicit delete at several places in the code.

  if (err) return nullptr;

  //
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  E_Int ni, nj, nk;
  E_Int res2 = K_ARRAY::getFromArray(arrS, varString2, fS, ni, nj, nk,
                                     cnS, eltType2);
  
  std::unique_ptr<K_FLD::FloatArray> afmesh2(fS);  // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<K_FLD::IntArray> acmesh2(cnS);

  K_FLD::FloatArray & crdS = *fS;

  if (sensor_type == 1 /*xensor*/)
  {
    if ( (res2 != 2) || ((res2 == 2) && (strcmp(eltType2, "HEXA") != 0) ) )
    {
     PyErr_SetString(PyExc_ValueError,
       "adaptCells: xsensor currently support only HEXA mesh as source mesh.");
     return nullptr;
    }
  }

  bool update_hook = (hook_hmesh != Py_None);
  elt_t elt_type (elt_t::UNKN);
  void * hookhm = nullptr;
  if (update_hook)
  {
    void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    packet = (void**) PyCObject_AsVoidPtr(hook_hmesh);
#else
    packet = (void**) PyCapsule_GetPointer(hook_hmesh, NULL);
#endif
    E_Int type = *((E_Int*)(packet[0]));
    sub_type = *((NUGA::eSUBDIV_TYPE*)(packet[2]));
    elt_type = *((elt_t*)(packet[3]));

    hookhm = (void*)packet[1];

    if (type != HMESH_HOOK_ID)
    {
      PyErr_SetString(PyExc_TypeError, 
                    "adaptCells: this function requires a identify hmesh hook.");
      return nullptr;
    }
  }
  else
    elt_type = check_has_NGON_BASIC_ELEMENT(cnt);

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  PyObject *l(PyList_New(0));

  if (elt_type==-1)
  {
    PyErr_SetString(PyExc_ValueError,
       "adaptCells: input mesh to adapt must have basic elements and must be in NGON format.");
    return nullptr;
  }

  std::vector<E_Int> dummy;
  if (sub_type == NUGA::ISO)
    err = __adapt_wrapper<ISO>(crd, cnt, elt_type, sensor_type, smoothing_type, crdS, dummy, itermax, varString, l, hookhm, cnS);
  else if (sub_type == NUGA::ISO_HEX)
    err = __adapt_wrapper<ISO_HEX>(crd, cnt, elt_type, sensor_type, smoothing_type, crdS, dummy, itermax, varString, l, hookhm, cnS);

  return (err) ? nullptr : l;
}


//=============================================================================
/* get hmesh hook  */
//=============================================================================
void* unpackHMesh(PyObject* hook_hmesh, subdiv_t &subdiv_type, elt_t &elt_type, std::string& vString)
{
  //todo CW : retourner des pointeurs plutot que des valeurs, pour pouvoir nettoyer dans deleteHMesh

  void** packet = NULL;

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**)PyCObject_AsVoidPtr(hook_hmesh);
#else
  packet = (void**)PyCapsule_GetPointer(hook_hmesh, NULL);
#endif

  E_Int* hook_hm_id    = (E_Int*)packet[0];               // type of hook
  void* hmesh          = packet[1];                       // uuntyped hmesh ptr
  subdiv_type          = subdiv_t (*((E_Int*)packet[2])); // subdivision type ISO, ISO_HEX, DIR...
  elt_type             = elt_t(*((E_Int*)packet[3]));     // type of elements in hmesh
  vString              = *((std::string*)packet[4]);      // for buildArray
  
  if (*hook_hm_id != HMESH_HOOK_ID)
  {
    PyErr_SetString(PyExc_TypeError,
      "unpackHMesh: this function requires a identify hmesh hook.");
    return nullptr;
  }
  return hmesh;
}

//=============================================================================
/* get sensor hook  */
//=============================================================================
E_Int unpackSensor(PyObject* hook_sensor, void* &sensor, E_Int &sensor_type, E_Int& smoothing_type)
{
  void** packet_ss = NULL;

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet_ss = (void**)PyCObject_AsVoidPtr(hook_sensor);
#else
  packet_ss = (void**)PyCapsule_GetPointer(hook_sensor, NULL);
#endif

  E_Int hook_ss_id = *((E_Int*)(packet_ss[0]));

  if (hook_ss_id != SENSOR_HOOK_ID)
  {
    PyErr_SetString(PyExc_TypeError,
      "unpackSensor: this function requires a identify sensor hook.");
    return NULL;
  }

  sensor_type = *((E_Int*)(packet_ss[1]));
  sensor = (void*)packet_ss[2];
  smoothing_type = *((E_Int*)(packet_ss[3]));//todo CW : veriier que ca marche
  return 1;
}

//=============================================================================
/* Dynamic cells adaptation */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCellsDyn(PyObject* self, PyObject* args)
{

  PyObject *arrS(nullptr), *hook_hmesh(nullptr), *hook_sensor(nullptr);
  if (!PYPARSETUPLEI(args, "OOO", "OOO", &arrS, &hook_hmesh, &hook_sensor)) return NULL;


  // Unpack hmesh hook
  // ==================
  void * hmesh = nullptr;
  elt_t elt_type;
  NUGA::eSUBDIV_TYPE subdiv_type;
  std::string vString;

  std::cout << " unpack hmesh" << std::endl;
  hmesh = unpackHMesh(hook_hmesh, subdiv_type, elt_type, vString);
  if (hmesh == nullptr) return nullptr;
  std::cout << " unpack hmesh OK " << std::endl;

  // Unpack sensor hook
  // ==================
  void *sensor = nullptr;
  E_Int sensor_type;
  E_Int smoothing_type;
  std::cout << " unpack sensor" << std::endl;
  unpackSensor(hook_sensor, sensor, sensor_type, smoothing_type);
  std::cout << " unpack sensor OK" << std::endl;

  //todo CW

  PyObject* out(PyList_New(0));
  return out; // TMP (juste pour retourner qqchose)
}

// ###################################################################################


//=============================================================================
/* Adapt cells with respect to the nodal subdivisions query */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCellsNodal(PyObject* self, PyObject* args)
{
  PyObject *arr(nullptr), *nodal_vals(nullptr), *hook_hmesh(nullptr);

  if (!PyArg_ParseTuple(args, "OOO", &arr, &nodal_vals, &hook_hmesh)) return NULL;

  K_FLD::FloatArray* f(nullptr), *fS(nullptr);
  K_FLD::IntArray* cn(nullptr), *cnS(nullptr);
  char* varString, *eltType, *varString2, *eltType2;
  // Check the mesh (NGON)
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  //// Gestion interne HEXA/TETRA format : Commente car resultat different lorsque fait dans le python
  // E_Int ni, nj, nk;
  
  // E_Int res = K_ARRAY::getFromArray(arr, varString, f, ni, nj, nk,
  //                                   cn, eltType);
     
  // bool err = (res !=2);
  // bool is_TH4 = (strcmp(eltType, "TETRA") == 0);
  // bool is_HX8 = (strcmp(eltType, "HEXA") == 0);
  // err |= (strcmp(eltType, "NGON") != 0) && !is_TH4 && !is_HX8;
  // if (err)
  // { 
  //   PyErr_SetString(PyExc_TypeError, "input error : invalid array, must be a unstructured NGON, HEXA or TETRA array.");//fixme triangulateExteriorFaces : PASS A STRING AS INPUT
  //   delete f; delete cn;
  //   return NULL;
  // }

  // // Check coordinates.
  // E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  // E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  // E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  // if ((posx == -1) || (posy == -1) || (posz == -1))
  // {
  //   PyErr_SetString(PyExc_TypeError, "input error : can't find coordinates in array.");//fixme  conformUnstr
  //   delete f; delete cn;
  //   return NULL;
  // }

  std::unique_ptr<K_FLD::FloatArray> afmesh(f);  // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<K_FLD::IntArray> acmesh(cn); // to avoid to call explicit delete at several places in the code.
  if (err)
    return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  E_Int *nodv, size, nfld;
  E_Int res2 = K_NUMPY::getFromNumpyArray(nodal_vals, nodv, size, nfld, true/* shared*/, false/* inverse*/);

  bool update_hook = (hook_hmesh != Py_None);
  elt_t elt_type (elt_t::UNKN);
  void * hookhm = nullptr;
  if (update_hook)
  {
    //todo CW : utiliser unpackHMesh
    void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    packet = (void**) PyCObject_AsVoidPtr(hook_hmesh);
#else
    packet = (void**) PyCapsule_GetPointer(hook_hmesh, NULL);
#endif
    E_Int type = *((E_Int*)(packet[0]));
    NUGA::eSUBDIV_TYPE sub_type = *((NUGA::eSUBDIV_TYPE*)(packet[2]));
    elt_type = *((elt_t*)(packet[3]));

    hookhm = (void*)packet[1];

    if (type != HMESH_HOOK_ID)
    {
      PyErr_SetString(PyExc_TypeError, 
                    "adaptCellsNodal: this function requires a identify hmesh hook.");
      return NULL;
    }
  }
  else
    elt_type = check_has_NGON_BASIC_ELEMENT(cnt);

  ///
  std::vector<E_Int> nodal_data(crd.cols(), 0);
  E_Int imax = std::min(size, crd.cols());
  
  for (E_Int i=0; i < imax; ++i)
  {
      E_Int v = nodv[i];
      v = std::max((E_Int)0,v); //disable agglo currently
      nodal_data[i] = v;
  }
  
  PyObject *l(PyList_New(0));

  if (elt_type == elt_t::UNKN)
  {
    PyErr_SetString(PyExc_ValueError,
       "adaptCells: input mesh to adapt must have basic elements (Only Tets and Hexas are currently adapted) and must be in NGON format.");
    return NULL;
  }

  K_FLD::FloatArray dummy;
  //if (sub_type == NUGA::ISO)
    err = __adapt_wrapper<NUGA::ISO>(crd, cnt, elt_type, 2/*nodal sensor*/, -1, dummy, nodal_data, -1/*itermax*/, varString, l, hookhm);
  //else if (sub_type == NUGA::ISO_HEX)
    //err = __adapt_wrapper<NUGA::ISO_HEX>(crd, cnt, elt_type, 2/*nodal sensor*/, -1, dummy, nodal_data, -1/*itermax*/, varString, l, hookhm);

  return (err) ? nullptr : l;
}

//=============================================================================
/* todo */
//=============================================================================
PyObject* K_INTERSECTOR::adaptBox(PyObject* self, PyObject* args)
{
  PyObject *arrS(nullptr);
  E_Float bratio(10.);
  E_Int itermax(0), smoothing_type(0);

  if (!PYPARSETUPLE(args, "Odll", "Odii", "Ofll", "Ofii", &arrS, &bratio, &smoothing_type, &itermax)) return NULL;

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
  box.convert2NG<ngon_type>(crd, ngi);
  K_FLD::IntArray cnt;
  ngi.export_to_array(cnt);//fixme : convoluted

  //std::cout << "adapt box..." << std::endl;

  using mesh_type = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;
  using sensor_type = NUGA::geom_sensor/*geom_static_sensor*/<mesh_type>;
  
  mesh_type hmesh(crd, cnt);
  hmesh.init();
  sensor_type sensor(hmesh, NUGA::eSmoother(smoothing_type), 1/*max pts per cell*/, itermax);
  E_Int err(0);
  err = sensor.assign_data(crdS);
  
  NUGA::adaptor<mesh_type, sensor_type>::run(hmesh, sensor);

  //std::cout << "output leaves..." << std::endl;
  std::vector<E_Int> oids;
  ngon_type ngo;
  hmesh.conformize(ngo, oids);

  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);

  //std::cout << "output ..." << std::endl;
  
  PyObject* tpl = K_ARRAY::buildArray(hmesh._crd, varString, cnto, -1, "NGON", false);

  delete f; delete cn;
  return tpl;
}


//=======================  Intersector/PolyMeshTools/splitCells.cpp ====================
