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
# include "Nuga/include/ngon_t.hxx"
# include "Nuga/include/Triangulator.h"
# include "Nuga/include/Splitter.h"

#include "Nuga/include/tree.hxx"
#include "Nuga/include/geom_sensor.hxx"
#include "Nuga/include/xsensor.hxx"
#include "Nuga/include/nodal_sensor.hxx"
#include "Nuga/include/adaptor.hxx"
#include "Nuga/include/hierarchical_mesh.hxx"

#include "Nuga/include/BbTree.h"
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
/* get hmesh hook  */
//=============================================================================
void* unpackHMesh(PyObject* hook_hmesh, E_Int *&hook_hm_id, E_Int *&subdiv_type, E_Int *&elt_type, std::string *&vString, void **&packet)
{

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**)PyCObject_AsVoidPtr(hook_hmesh);
#else
  packet = (void**)PyCapsule_GetPointer(hook_hmesh, NULL);
#endif

  hook_hm_id           = (E_Int*)packet[0];        // type of hook
 
  if (*hook_hm_id != HMESH_HOOK_ID)
  {
    PyErr_SetString(PyExc_TypeError,
      "unpackHMesh: this function requires a identify hmesh hook.");
    return nullptr;
  }
  
  void* hmesh          = packet[1];                // untyped hmesh ptr
  subdiv_type          = (E_Int*)packet[2];        // subdivision type ISO, ISO_HEX, DIR...  
  elt_type             = (E_Int*)packet[3];        // type of elements in hmesh
  vString              = (std::string*)packet[4];  // for buildArray

  return hmesh;
}

//=============================================================================
/* get sensor hook  */
//=============================================================================
void* unpackSensor(PyObject* hook_sensor, E_Int *&hook_ss_id, E_Int *&sensor_type, E_Int *&smoothing_type, E_Int *&subdiv_type, E_Int *&elt_type, void **&packet_ss)
{

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet_ss = (void**)PyCObject_AsVoidPtr(hook_sensor);
#else
  packet_ss = (void**)PyCapsule_GetPointer(hook_sensor, NULL);
#endif

  hook_ss_id = (E_Int*)packet_ss[0];

  if (*hook_ss_id != SENSOR_HOOK_ID)
  {
    PyErr_SetString(PyExc_TypeError,
      "unpackSensor: this function requires a identify sensor hook.");
    return nullptr;
  }

  sensor_type    = (E_Int*)packet_ss[1];
  smoothing_type = (E_Int*)packet_ss[3];
  elt_type       = (E_Int*)packet_ss[4];
  subdiv_type    = (E_Int*)packet_ss[5];
  void* sensor   = packet_ss[2];

  return sensor;
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

  void * hmesh = nullptr;
  E_Int* sub_type, *elt_type, *hook_id;
  std::string* vString;
  void** packet = NULL;
  hmesh = unpackHMesh(hook, hook_id, sub_type, elt_type, vString, packet);

  PyObject *l(PyList_New(0));
  K_FLD::IntArray cnto;
  std::vector<E_Int> oids;
  K_FLD::FloatArray* crd(nullptr);

  if (*sub_type == 0) // ISO
    __conformizeHM<NUGA::ISO>(*elt_type, hmesh, crd, cnto, oids);
  else if (*sub_type == 1) // ISO_HEX
    __conformizeHM<NUGA::ISO_HEX>(*elt_type, hmesh, crd, cnto, oids);

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
  void * hmesh = nullptr;
  E_Int* sub_type, *elt_type, *hook_id;
  std::string* vString;
  void** packet = NULL ; 
  hmesh = unpackHMesh(hook, hook_id, sub_type, elt_type, vString, packet);
  
  if (*sub_type == 0) // ISO
    __deleteHM<NUGA::ISO>(*elt_type, hmesh);
  else if (*sub_type == 1) // ISO_HEX
    __deleteHM<NUGA::ISO_HEX>(*elt_type, hmesh);

  delete hook_id;
  delete vString;
  delete sub_type;
  delete elt_type;
  delete [] packet;
  Py_INCREF(Py_None);
  return Py_None;
}

//============================================================================
/* Create a Geom sensor and returns a hook */
//============================================================================

template<typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void* __createSensor(void* hmesh, E_Int smoothing_type, E_Int itermax, E_Int sensor_type)
{
  if (sensor_type == 0)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t = NUGA::geom_sensor<hmesh_t>;
    return new sensor_t(*(hmesh_t*)hmesh, NUGA::eSmoother(smoothing_type), 1/*max_pts_per_cell*/, itermax);
  }
  /*else if (sensor_type == 1)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t = NUGA::xsensor<hmesh_t, ELT_t>;
    return new sensor_t(*(hmesh_t*)hmesh, NUGA::eSmoother(smoothing_type), 1, itermax);
  }*/
  else if (sensor_type == 2)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t = NUGA::nodal_sensor<hmesh_t>;
    return new sensor_t(*(hmesh_t*)hmesh);
  }
}

template<>
void* __createSensor<K_MESH::Hexahedron, NUGA::ISO>(void* hmesh, E_Int smoothing_type, E_Int itermax, E_Int sensor_type)
{
  if (sensor_type == 0)
  {
    using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;
    using sensor_t = NUGA::geom_sensor<hmesh_t>;
    return new sensor_t(*(hmesh_t*)hmesh, NUGA::eSmoother(smoothing_type), 1/*max_pts_per_cell*/, itermax);
  }
  else if (sensor_type == 1)
  {
    using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;
    using sensor_t = NUGA::xsensor<K_MESH::Hexahedron, hmesh_t>;
    return new sensor_t(*(hmesh_t*)hmesh, NUGA::eSmoother(smoothing_type), 1, itermax);
  }
  else if (sensor_type == 2)
  {
    using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;
    using sensor_t = NUGA::nodal_sensor<hmesh_t>;
    return new sensor_t(*(hmesh_t*)hmesh);
  }
}

template<>
void* __createSensor<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(void* hmesh, E_Int smoothing_type, E_Int itermax, E_Int sensor_type)
{
  if (sensor_type == 0)
  {
    using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>;
    using sensor_t = NUGA::geom_sensor<hmesh_t>;
    return new sensor_t(*(hmesh_t*)hmesh, NUGA::eSmoother(smoothing_type), 1/*max_pts_per_cell*/, itermax);
  }
  else if (sensor_type == 2)
  {
    using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>;
    using sensor_t = NUGA::nodal_sensor<hmesh_t>;
    return new sensor_t(*(hmesh_t*)hmesh);
  }
}

template<NUGA::eSUBDIV_TYPE STYPE>
void* __createSensor(E_Int etype, void* hmesh, E_Int sensor_type, E_Int smoothing_type, E_Int itermax);

template<>
void* __createSensor<NUGA::ISO>(E_Int elt_type, void* hmesh, E_Int sensor_type, E_Int smoothing_type, E_Int itermax)
{
  if (elt_type == elt_t::HEXA)
    return __createSensor<K_MESH::Hexahedron, NUGA::ISO>(hmesh, smoothing_type, itermax, sensor_type);
  else if (elt_type == elt_t::TETRA)
    return __createSensor<K_MESH::Tetrahedron, NUGA::ISO>(hmesh, smoothing_type, itermax, sensor_type);
  else if (elt_type == elt_t::PRISM3)
    return __createSensor<K_MESH::Prism, NUGA::ISO>(hmesh, smoothing_type, itermax, sensor_type);
  else if (elt_type == elt_t::BASIC)
    return __createSensor<K_MESH::Basic, NUGA::ISO>(hmesh, smoothing_type, itermax, sensor_type);
  else 
  {
    PyErr_SetString(PyExc_ValueError, "createSensor: wrong element type in hmesh for ISO strategy.");
    return nullptr;
  }
}

//
template<>
void* __createSensor<NUGA::ISO_HEX>(E_Int etype, void* hmesh, E_Int sensor_type, E_Int smoothing_type, E_Int itermax)
{
  return __createSensor<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh, smoothing_type, itermax, sensor_type);
}

//============================================================================
/* Creates a Sensor */
//============================================================================
PyObject* K_INTERSECTOR::createSensor(PyObject* self, PyObject* args)
{
  PyObject *hook_sensor(nullptr);

  PyObject *hook_hmesh(nullptr);
  E_Int smoothing_type(0), sensor_type(0), itermax(0); // sensor_type = 0(geom_sensor) or 1 (xsensor) 

  if (!PYPARSETUPLEI(args, "Olll", "Oiii", &hook_hmesh, &sensor_type, &smoothing_type, &itermax)) return NULL;


  // // Unpack hmesh hook
  // // ==================
  E_Int* subtype_hm, *elttype_hm, *hook_id;
  std::string* vString;
  void ** packet_h = NULL;
  void* hmesh = unpackHMesh(hook_hmesh, hook_id, subtype_hm, elttype_hm, vString, packet_h);

  // // Create packet for sensor hook 
  // // ==============================
  void** packet_ss = new void*[6];  // hook ID, sensor type, hmesh ptr, smoothing type, elt type, subdiv type

  packet_ss[2] = nullptr;
  packet_ss[4] = nullptr;

  // HOOK ID
  E_Int* hook_sensor_id = new E_Int;
  packet_ss[0]          = hook_sensor_id;
  *hook_sensor_id       = SENSOR_HOOK_ID;

  // SENSOR TYPE
  E_Int* sens_type      = new E_Int;
  packet_ss[1]          = sens_type;
  *sens_type            = sensor_type;

  // SMOOTHING TYPE
  E_Int* smooth_type    = new E_Int;
  packet_ss[3]          = smooth_type;
  *smooth_type          = smoothing_type;

  // ELT TYPE
  E_Int* elt_type       = new E_Int;
  packet_ss[4]          = elt_type;
  *elt_type             = *elttype_hm;

  // SUBDIV TYPE
  E_Int* sub_type       = new E_Int;
  packet_ss[5]          = sub_type;
  *sub_type             = *subtype_hm;

  if (sensor_type == 2 && *sub_type != NUGA::ISO)
  {
    PyErr_SetString(PyExc_ValueError,
       "adaptCells: nodal sensor only supports ISO subdivision currently.");
    return nullptr;
  }

  // HMESH PTR
  if (*subtype_hm == 0) // ISO
    packet_ss[2] = __createSensor<NUGA::ISO>(*elt_type, hmesh, sensor_type, smoothing_type, itermax);
  else if (*subtype_hm == 1) // ISO_HEX
    packet_ss[2] = __createSensor<NUGA::ISO_HEX>(*elt_type, hmesh, sensor_type, smoothing_type, itermax);

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
// ============================================================================
template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void __deleteSensor(E_Int sensor_type, const void* sensor_ptrs)
{
  if (sensor_type == 0)
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_type   = NUGA::geom_sensor<mesh_type>;
    sensor_type* sensor = (sensor_type*)sensor_ptrs;
    delete sensor;    
  }
  else if (sensor_type == 1)
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_type   = NUGA::xsensor<ELT_t, mesh_type>;
    sensor_type* sensor = (sensor_type*)sensor_ptrs;
    delete sensor;    
  }
  else if (sensor_type == 2)
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_type   = NUGA::nodal_sensor<mesh_type>;
    sensor_type* sensor = (sensor_type*)sensor_ptrs;
    delete sensor;    
  }
}

template <NUGA::eSUBDIV_TYPE STYPE>
void __deleteSensor(E_Int etype, E_Int sensor_type, const void* sensor_ptr)
{
  if (etype == elt_t::HEXA)
    __deleteSensor<K_MESH::Hexahedron, STYPE>(sensor_type, sensor_ptr);
  else if (etype == (E_Int)elt_t::TETRA)
    __deleteSensor<K_MESH::Tetrahedron, STYPE>(sensor_type, sensor_ptr);
  else if (etype == (E_Int)elt_t::PRISM3)
    __deleteSensor<K_MESH::Prism, STYPE>(sensor_type, sensor_ptr);
  else if (etype == (E_Int)elt_t::BASIC)
    __deleteSensor<K_MESH::Basic, STYPE>(sensor_type, sensor_ptr);
}

template <>
void __deleteSensor<NUGA::ISO_HEX>(E_Int etype/*dummy*/, E_Int sensor_type, const void* sensor_ptr)
{
  __deleteSensor<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(sensor_type, sensor_ptr);
}


PyObject* K_INTERSECTOR::deleteSensor(PyObject* self, PyObject* args)
{
  PyObject* hook_sensor;
  if (!PyArg_ParseTuple(args, "O", &hook_sensor))
  {
      return NULL;
  }

  // recupere le hook
  void *sensor = nullptr;
  E_Int *hook_ss_id, *sensor_type,*smoothing_type;
  E_Int *subdiv_type, *elt_type;
  void** packet;
  sensor = unpackSensor(hook_sensor, hook_ss_id, sensor_type, smoothing_type, subdiv_type, elt_type, packet);  

  if (*subdiv_type == 0) // ISO
    __deleteSensor<NUGA::ISO>(*elt_type, *sensor_type, sensor);
  else if (*subdiv_type == 1) // ISO_HEX
    __deleteSensor<NUGA::ISO_HEX>(*elt_type, *sensor_type, sensor);

  delete hook_ss_id, sensor_type, smoothing_type, subdiv_type, elt_type;
  delete [] packet;
  Py_INCREF(Py_None);
  return Py_None;
}

// //============================================================================
// /* assign data to a sensor*/
// ============================================================================
template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void __assign_sensor_data
(E_Int sensor_type, void* psensor, K_FLD::FloatArray& crdS, K_FLD::IntArray& cntS, std::vector<E_Int>& nodal_data)
{
  if (sensor_type == 0)
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t   = NUGA::geom_sensor<mesh_type>;
    sensor_t* sensor = (sensor_t*)psensor;

    sensor->assign_data(crdS);
  }
  else if (sensor_type == 2)
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t   = NUGA::nodal_sensor<mesh_type>;
    sensor_t* sensor = (sensor_t*)psensor;

    sensor->assign_data(nodal_data);
  }
}

// xsensor is only available for HEXA/ISO 
template <>
void __assign_sensor_data<K_MESH::Hexahedron, NUGA::ISO>
(E_Int sensor_type, void* psensor, K_FLD::FloatArray& crdS, K_FLD::IntArray& cntS, std::vector<E_Int>& nodal_data)
{
  if (sensor_type == 0)
  {
    using mesh_type     = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;
    using sensor_t   = NUGA::geom_sensor<mesh_type>;
    sensor_t* sensor = (sensor_t*)psensor;

    sensor->assign_data(crdS);
  }
  else if (sensor_type == 1)
  {
    using mesh_type     = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;
    using sensor_t   = NUGA::xsensor<K_MESH::Hexahedron, mesh_type>;
    sensor_t* sensor = (sensor_t*)psensor;

    NUGA::mesh_s data;
    data.crd = crdS;
    data.cnt = cntS;
    sensor->assign_data(data);
  }
  else if (sensor_type == 2)
  {
    using mesh_type     = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;
    using sensor_t   = NUGA::nodal_sensor<mesh_type>;
    sensor_t* sensor = (sensor_t*)psensor;

    sensor->assign_data(nodal_data);
  }
}


template <NUGA::eSUBDIV_TYPE STYPE>
void __assign_sensor_data
(E_Int etype, E_Int sensor_type, void* sensor_ptr, K_FLD::FloatArray& crdS, K_FLD::IntArray& cntS, std::vector<E_Int>& nodal_data)
{
  if (etype == elt_t::HEXA)
    __assign_sensor_data<K_MESH::Hexahedron, STYPE>(sensor_type, sensor_ptr, crdS, cntS, nodal_data);
  else if (etype == (E_Int)elt_t::TETRA)
    __assign_sensor_data<K_MESH::Tetrahedron, STYPE>(sensor_type, sensor_ptr, crdS, cntS, nodal_data);
  else if (etype == (E_Int)elt_t::PRISM3)
    __assign_sensor_data<K_MESH::Prism, STYPE>(sensor_type, sensor_ptr, crdS, cntS, nodal_data);
  else if (etype == (E_Int)elt_t::BASIC)
    __assign_sensor_data<K_MESH::Basic, STYPE>(sensor_type, sensor_ptr, crdS, cntS, nodal_data);
}

template <>
void __assign_sensor_data<NUGA::ISO_HEX>
(E_Int etype, E_Int sensor_type, void* sensor_ptr, K_FLD::FloatArray& crdS, K_FLD::IntArray& cntS, std::vector<E_Int>& nodal_data)
{
  __assign_sensor_data<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(sensor_type, sensor_ptr, crdS, cntS, nodal_data);
}


PyObject* K_INTERSECTOR::assignData2Sensor(PyObject* self, PyObject* args)
{
  PyObject* hook_sensor, *dataSensor(nullptr);
  if (!PyArg_ParseTuple(args, "OO", &hook_sensor, &dataSensor))
  {
      return NULL;
  }

  // recupere le hook
  void *sensor = nullptr;
  E_Int *hook_ss_id, *sensor_type,*smoothing_type;
  E_Int *subdiv_type, *elt_type;
  void** packet;
  sensor = unpackSensor(hook_sensor, hook_ss_id, sensor_type, smoothing_type, subdiv_type, elt_type, packet);  

  
  //geom/xsensor or nodal_sensor data ?
  std::vector<E_Int> nodal_data;
  K_FLD::FloatArray fS;
  K_FLD::IntArray cnS;

  if (PyList_Check(dataSensor)) // Array (mesh or coordinates)
  {
    E_Int ni, nj, nk;
    char* varString, *eltType;
    E_Int res = K_ARRAY::getFromArray(dataSensor, varString, fS, ni, nj, nk, cnS, eltType);
    //std::cout << "res/eltType/fs sz/cns sz : " << res << "/" << eltType << "/" << fS.cols() << "/" << cnS.cols() << std::endl;
  }
  else // assuming numpy for nodal
  {
     E_Int *nodv, size, nfld;
     E_Int res2 = K_NUMPY::getFromNumpyArray(dataSensor, nodv, size, nfld, true/* shared*/, false/* inverse*/);
     nodal_data.resize(size);  
     for (E_Int i=0; i < size; ++i)
     {
       E_Int v = nodv[i];
       v = std::max((E_Int)0,v); //disable agglo currently
       nodal_data[i] = v;
     }
  }

  if (*subdiv_type == 0) // ISO
    __assign_sensor_data<NUGA::ISO>(*elt_type, *sensor_type, sensor, fS, cnS, nodal_data);
  else if (*subdiv_type == 1) // ISO_HEX
    __assign_sensor_data<NUGA::ISO_HEX>(*elt_type, *sensor_type, sensor, fS, cnS, nodal_data);

  Py_INCREF(Py_None);
  return Py_None;
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


template <typename MESH_t, typename sensor_t>
E_Int __adapt
(MESH_t* hmesh, sensor_t* sensor, const char* varString, PyObject *out)
{
  hmesh->init();
  NUGA::adaptor<MESH_t, sensor_t>::run(*hmesh, *sensor);

  std::vector<E_Int> oids;
  K_FLD::IntArray cnto;

  hmesh->_ng.export_to_array(cnto);

  // pushing out the mesh
  PyObject *tpl = K_ARRAY::buildArray(hmesh->_crd, varString, cnto, -1, "NGON", false);
  PyList_Append(out, tpl);
  Py_DECREF(tpl);

  return 0;
}

///
template <subdiv_t STYPE>
E_Int __adapt_wrapper
(E_Int elt_type, E_Int sensor_type,
const char* varString, PyObject *out,
void* hookhm, void* hooksens)
{
  E_Int err(0);
  if (elt_type==elt_t::HEXA)
  {
    using ELT_type = K_MESH::Hexahedron;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, STYPE>;

    mesh_type* hmesh = (mesh_type*)hookhm;

    if (sensor_type == 0) // geom sensor
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t* sensor = (sensor_t*)hooksens;

      err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
    }
    else if (sensor_type == 1) // xsensor
    {
      using sensor_t = NUGA::xsensor<ELT_type, mesh_type>;
      sensor_t* sensor = (sensor_t*)hooksens;
      
      err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;
      sensor_t* sensor = (sensor_t*)hooksens;

      err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
    }
  }
  else if (elt_type==elt_t::TETRA)
  {
    using ELT_type = K_MESH::Tetrahedron;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, STYPE>;

    mesh_type* hmesh = (mesh_type*)hookhm;

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported
    
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t* sensor = (sensor_t*)hooksens;

      err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;
      sensor_t* sensor = (sensor_t*)hooksens;

      err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
    }
  }
  else if (elt_type==elt_t::PRISM3)
  {
    using ELT_type = K_MESH::Prism;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, STYPE>;

    mesh_type* hmesh = (mesh_type*)hookhm;

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported
    
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t* sensor = (sensor_t*)hooksens;

      err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;
      sensor_t* sensor = (sensor_t*)hooksens;

      err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
    }
  }
  else if (elt_type==elt_t::BASIC)
  {    
    using ELT_type = K_MESH::Basic;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, STYPE>;

    mesh_type* hmesh = (mesh_type*)hookhm;

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported

    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;
      sensor_t* sensor = (sensor_t*)hooksens;

      err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;
      sensor_t* sensor = (sensor_t*)hooksens;

      err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
    }
  }
  return err;
}

template <>
E_Int __adapt_wrapper<NUGA::ISO_HEX>
(E_Int elt_type/*dummy*/, E_Int sensor_type,
const char* varString, PyObject *out,
void* hookhm, void* hooksens)
{
  E_Int err(0);
  
  using ELT_type = K_MESH::Polyhedron<0>;
  using mesh_type = NUGA::hierarchical_mesh<ELT_type, NUGA::ISO_HEX>;

  mesh_type* hmesh = (mesh_type*)hookhm;

  if (sensor_type == 0) // geom sensor
  {
    using sensor_t = NUGA::geom_sensor<mesh_type>;
    sensor_t* sensor = (sensor_t*)hooksens;

    err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
  }
  // else if (sensor_type == 1) // xsensor
  // {
  //   using sensor_t = NUGA::xsensor<ELT_type, mesh_type>;
  //   sensor_t* sensor = (sensor_t*)hooksens;
    
  //   err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
  // }
  // else if (sensor_type == 2) //nodal sensor
  // {
  //   using sensor_t = NUGA::nodal_sensor<mesh_type>;
  //   sensor_t* sensor = (sensor_t*)hooksens;

  //   err = __adapt<mesh_type, sensor_t>(hmesh, sensor, varString, out);
  // }

  return err;
}

//=============================================================================
/* Hierarchical Mesh Adaptation */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCells(PyObject* self, PyObject* args)
{

  PyObject *hook_hmesh(nullptr), *hook_sensor(nullptr);

  if (!PyArg_ParseTuple(args, "OO", &hook_hmesh, &hook_sensor)) return NULL;

  // Unpack hmesh hook
  // ==================
  void * hmesh = nullptr;
  E_Int* elt_type, *subdiv_type, *hook_id;
  std::string* vString;
  void** packet = NULL ; 

  // std::cout << " unpack hmesh" << std::endl;
  hmesh = unpackHMesh(hook_hmesh, hook_id, subdiv_type, elt_type, vString, packet);
  if (hmesh == nullptr) return nullptr;
  // std::cout << " unpack hmesh OK " << std::endl;

  // Unpack sensor hook
  // ==================
  void *sensor = nullptr;
  E_Int *hook_ss_id, *sensor_type, *smoothing_type;
  E_Int *subdiv_type_ss, *elt_type_ss;
  void** packet_ss;
  // std::cout << " unpack sensor" << std::endl;
  sensor = unpackSensor(hook_sensor, hook_ss_id, sensor_type, smoothing_type, subdiv_type_ss, elt_type_ss, packet_ss);  
  // std::cout << " unpack sensor OK" << std::endl;

  // Check basic elements presence     # to do in createHMesh step? 
  // =============================
  if (*elt_type==-1)
  {
    PyErr_SetString(PyExc_ValueError,
       "adaptCells: input mesh to adapt must have basic elements and must be in NGON format.");
    return nullptr;
  }

  if (*elt_type_ss != *elt_type)
  {
    PyErr_SetString(PyExc_ValueError,
      "adaptCells: type of elements : inconsistency between sensor and hmesh.");
    //std::cout << "hm/sens : " << *elt_type << " vs " << *elt_type_ss << std::endl;
    return nullptr;
  }
  
  if (*sensor_type == 1 /*xsensor*/)
  {
    if (*elt_type_ss != elt_t::HEXA)
    {
     PyErr_SetString(PyExc_ValueError,
       "adaptCells: xsensor currently support only HEXA mesh as source mesh.");
     return nullptr;
    }
  }

  // Adaptation
  // ==========
  PyObject *l(PyList_New(0));
  std::vector<E_Int> dummy;

  E_Int err(0);
  if (*subdiv_type == NUGA::ISO)
    err = __adapt_wrapper<ISO>(*elt_type, *sensor_type, vString->c_str(), l, hmesh, sensor);
  else if (*subdiv_type == NUGA::ISO_HEX)
    err = __adapt_wrapper<ISO_HEX>(*elt_type, *sensor_type, vString->c_str(), l, hmesh, sensor);

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
