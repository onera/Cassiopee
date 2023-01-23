/*    
    Copyright 2013-2023 Onera.

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

# include <string>
# include <sstream> 
# include "intersector.h"
# include "Nuga/include/ngon_t.hxx"
# include "Nuga/include/Triangulator.h"
# include "Nuga/include/Splitter.h"

#include "Nuga/include/tree.hxx"
#include "Nuga/include/geom_sensor.hxx"
#include "Nuga/include/xsensor2.hxx"
#include "Nuga/include/nodal_sensor.hxx"
#include "Nuga/include/cell_sensor.hxx"
#include "Nuga/include/adaptor_mpi.hxx"
#include "Nuga/include/hierarchical_mesh.hxx"
#include "Nuga/include/smoother.hxx"

#include "Nuga/include/BbTree.h"
#include <memory>
//#include <iostream>

#include "adaptCells_hook.h"


using namespace std;
using namespace NUGA;

#ifdef FLAG_STEP
E_Int chrono::verbose = 0;
#endif

#define HMESH_HOOK_ID 77
#define SENSOR_HOOK_ID 78
#define COM_HOOK_ID 79
#define PACK_SIZE 5

using ngon_type = ngon_t<K_FLD::IntArray>;
using subdiv_t = NUGA::eSUBDIV_TYPE;
using elt_t = K_INTERSECTOR::eType;

//// DICO / MAP utils //////////////////////////
 
void convert_dico_to_map___int_int_vecint
(
  PyObject *py_zone_to_zone_to_list_owned,
  std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_zone_to_list_owned)
{
  if (PyDict_Check(py_zone_to_zone_to_list_owned))
  {
    //E_Int nzid = PyDict_Size(py_zone_to_zone_to_list_owned);

    PyObject *py_zid/*key*/, *py_zone_to_list_owned /*value : map jzid to ptlist*/;
    Py_ssize_t pos = 0;

    while (PyDict_Next(py_zone_to_zone_to_list_owned, &pos, &py_zid, &py_zone_to_list_owned))
    {
      int zid = (int) PyInt_AsLong(py_zid);

      assert (PyDict_Check(py_zone_to_list_owned) == 1); // it s a map

      PyObject *py_jzid/*key*/, *py_ptlist_owned /*value : ptlist*/;
      Py_ssize_t pos1 = 0;

      while (PyDict_Next(py_zone_to_list_owned, &pos1, &py_jzid, &py_ptlist_owned))
      {
        int jzid = (int) PyInt_AsLong(py_jzid);

        assert (PyArray_Check(py_ptlist_owned) == 1) ; // it s a numpy
        
        PyArrayObject* pyarr = reinterpret_cast<PyArrayObject*>(py_ptlist_owned);

        long ndims = PyArray_NDIM(pyarr);
        assert (ndims == 1); // vector
        npy_intp* dims = PyArray_SHAPE(pyarr);

        E_Int ptl_sz = dims[0];
        
        //long* dataPtr = static_cast<long*>(PyArray_DATA(pyarr));
        E_Int* dataPtr = (E_Int*)PyArray_DATA(pyarr);

        std::vector<E_Int> ptl(ptl_sz);
        for (size_t u=0; u < ptl_sz; ++u) ptl[u] = dataPtr[u];

        //std::cout << "max in C is : " << *std::max_element(ALL(ptl)) << std::endl;

        zone_to_zone_to_list_owned[zid][jzid]=ptl;

      }
    }
  }
}

void convert_dico_to_map__int_pairint
(
  PyObject *py_rid_to_zones,
  std::map<int, std::pair<int,int>>& rid_to_zones)
{
  if (PyDict_Check(py_rid_to_zones))
  {
    // E_Int nzid = PyDict_Size(py_rid_to_zones);

    PyObject *py_rid/*key*/, *py_pair_owned /*value : map zid to ptlist*/;
    Py_ssize_t pos = 0;

    while (PyDict_Next(py_rid_to_zones, &pos, &py_rid, &py_pair_owned))
    {
      int rid = (int) PyInt_AsLong(py_rid);

      assert (PyTuple_Check(py_pair_owned) == 1); // is it a tuple ?

      PyTupleObject* pytup = reinterpret_cast<PyTupleObject*>(py_pair_owned);    
      Py_ssize_t nb = PyTuple_GET_SIZE(pytup);

      // -----

      std::pair<int,int> pair_zid;

      PyObject * z1 PyTuple_GET_ITEM(pytup, 0);
      PyObject * z2 PyTuple_GET_ITEM(pytup, 1);

      pair_zid.first  = (double) PyFloat_AsDouble(z1);
      pair_zid.second = (double) PyFloat_AsDouble(z2);

      rid_to_zones[rid] = pair_zid;
    }
  }
}

struct transf_t {
  double t[6];
  bool operator==(const transf_t& r) const
  {
    for (size_t k=0; k < 6; ++k)
      if (t[k] != r.t[k]) return false;
    return true;
  }
  bool operator<(const transf_t& r) const
  {
    if (*this == r) return false;
    for (size_t k=0; k < 6; ++k)
      if (t[k] <r.t[k]) return true;
    return false;
  }
};

int convert_dico_to_map___transfo_to_vecint
(
  PyObject *py_transfo_to_list,
  std::map<transf_t, std::vector<int>>& transfo_to_list
)
{
  if (PyDict_Check(py_transfo_to_list) == 0) return 1;
  
  //E_Int nzid = PyDict_Size(transfo_to_list);

  PyObject *py_transfo/*key*/, *py_vecint /*value : vector<int>*/;
  Py_ssize_t pos = 0;

  transf_t t;

  while (PyDict_Next(py_transfo_to_list, &pos, &py_transfo, &py_vecint))
  {
    // key
    assert (PyTuple_Check(py_transfo) == 1) ; // it s a tuple (Xx, Yc, Zc, R)
    PyTupleObject* pytup = reinterpret_cast<PyTupleObject*>(py_transfo);
    Py_ssize_t nb = PyTuple_GET_SIZE(pytup);

    assert (nb == 6);
    for (size_t i=0; i < 6; ++i)
    {
      PyObject * p PyTuple_GET_ITEM(pytup, i);
      t.t[i] = (double) PyFloat_AsDouble(p);
      //std::cout << "transfo " << i << " : " << t.t[i] << std::endl;
    }

    // val
    assert (PyArray_Check(py_vecint) == 1) ; // it s a numpy
    PyArrayObject* pyarr = reinterpret_cast<PyArrayObject*>(py_vecint);

    long ndims = PyArray_NDIM(pyarr);
    assert (ndims == 1); // vector
    npy_intp* dims = PyArray_SHAPE(pyarr);

    E_Int sz = dims[0];
    
    //long* dataPtr = static_cast<long*>(PyArray_DATA(pyarr));
    E_Int* dataPtr = (E_Int*)PyArray_DATA(pyarr);

    transfo_to_list[t].resize(sz);
    for (size_t u=0; u < sz; ++u) transfo_to_list[t][u] = dataPtr[u];
  }
}

////////////////////////////////////////////////


///
template <NUGA::eSUBDIV_TYPE STYPE>
void* __createHM(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int zid);

// ISO strategy
template<>
void* __createHM<NUGA::ISO>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int zid)
{
  bool reorient =false;
  bool sync_match=false;

  if (typ == elt_t::UNKN)
  {
    PyErr_WarnEx(PyExc_Warning,
      "createHMesh: input mesh to adapt must have basic elements and must be in NGON format : no adaptation for this mesh.", 1);
    return nullptr;
  }
  else if (typ == elt_t::HEXA)
  {
    using elt_type = K_MESH::Hexahedron;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;

    hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt), reorient, sync_match);
    hm->zid = zid;
    return hm;
  }
  else if (typ == elt_t::TETRA)
  {
    using elt_type = K_MESH::Tetrahedron;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;

    hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt), reorient, sync_match);
    hm->zid = zid;
    return hm;
  }
  else if (typ == elt_t::PRISM3)
  {
    using elt_type = K_MESH::Prism;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;

    hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt), reorient, sync_match);
    hm->zid = zid;
    return hm;
  }
  else if (typ == elt_t::BASIC)
  {
    using elt_type = K_MESH::Basic;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;

    hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt), reorient, sync_match);
    hm->zid = zid;
    return hm;
  }
  return NULL;
}

// ISO_HEX strategy
template<>
void* __createHM<NUGA::ISO_HEX>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int zid)
{
  bool reorient =false;
  bool sync_match=false;

  using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>;

  hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt), reorient, sync_match);
  hm->zid = zid;
  return hm;
}

// DIR strategy
template<>
void* __createHM<NUGA::DIR>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int zid)
{
  bool reorient =false;
  bool sync_match=false;

  if (typ != elt_t::HEXA)
  {
    PyErr_WarnEx(PyExc_Warning,
      "createHMesh: directionnal policy is only supported with Hexahedral mesh currently.", 1);
    return nullptr;
  }

  using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::DIR>;
  
  hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt), reorient, sync_match);
  hm->zid = zid;
  return hm;
}

//============================================================================
/* Create a (list of) hmesh/hzone(s) and returns a (list of) hook(s) */
//============================================================================
PyObject* K_INTERSECTOR::createHMesh2(PyObject* self, PyObject* args)
{
  PyObject* hook;
  void** packet = new void*[HMESH_PACK_SIZE];  // hook_ID, hmesh ptr, subdiv policy, elt type, varString

  E_Int* hookid = new E_Int;  packet[0] = hookid;
  //void* hmesh_ptr = nullptr;  packet[1] = hmesh_ptr;// templated hmesh type to build 
  E_Int* subtype = new E_Int; packet[2] = subtype;
  elt_t* etyp = new elt_t;    packet[3] = etyp;
  std::string* vString = new std::string; packet[4] = vString;
  
  E_Int*zid = new E_Int; packet[5] = zid;

  *hookid = HMESH_HOOK_ID;

  PyObject *arr;

  if (!PYPARSETUPLEI(args, "Oll", "Oii", &arr, subtype, zid)) return nullptr;

  // mesh
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

  // ISO_HEX +> Polyhedron
  *etyp = elt_t::UNKN;// Polyhedron
  if (*subtype == NUGA::ISO || *subtype == NUGA::DIR)
    *etyp = check_has_NGON_BASIC_ELEMENT(cnt);

  if (*subtype == NUGA::ISO)
    packet[1] = __createHM<NUGA::ISO>(*etyp, crd, cnt, *zid);
  else if (*subtype == NUGA::ISO_HEX)
    packet[1] = __createHM<NUGA::ISO_HEX>(*etyp, crd, cnt, *zid);
  else if (*subtype == NUGA::DIR)
    packet[1] = __createHM<NUGA::DIR>(*etyp, crd, cnt, *zid);

  if (packet[1] == nullptr) return Py_None;// the input mesh does not have basic elts

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif
  
  return hook;
}

template <typename hmesh_t, typename sensor_t>
E_Int __adapt
(std::vector<hmesh_t*>& hmeshes, std::vector<sensor_t*>& sensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  if (hmeshes.empty()) return 0;
  if (sensors.empty()) return 0;

  size_t nb_meshes = hmeshes.size();
  if (nb_meshes != sensors.size()) return 1;

  using adaptor_t = NUGA::adaptor_mpi<hmesh_t, sensor_t>;

  adaptor_t a;
  a.sensors = sensors;
  a.do_agglo = false;
  std::vector<int> zids(hmeshes.size());
  for (size_t i=0; i < hmeshes.size(); ++i)zids[i]=hmeshes[i]->zid;

  a.run(hmeshes, zids, zone_to_rid_to_list_owned, rid_to_zones, zonerank, COM);

  for (size_t i=0; i < nb_meshes; ++i)
  {
    std::vector<E_Int> oids;
    K_FLD::IntArray cnto;

    hmeshes[i]->_ng.export_to_array(cnto);

    // pushing out the mesh
    PyObject *tpl = K_ARRAY::buildArray(hmeshes[i]->_crd, varString, cnto, -1, "NGON", false);
    PyList_Append(out, tpl);
    Py_DECREF(tpl);
  }

  //std::cout << "__adapt : DONE" << std::endl;

  return 0;
}

///
template <subdiv_t STYPE>
int __adapt_wrapper
(int elt_type, int sensor_type,
std::vector<void*>&hookhmes, std::vector<void*>&hooksensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_zone_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  int err(0);
  size_t nb_meshes = hookhmes.size();

  assert (nb_meshes == hooksensors.size());


  if (elt_type==elt_t::HEXA)
  {
    using ELT_type = K_MESH::Hexahedron;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, STYPE>;

    std::vector<mesh_type*> hmeshes(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) hmeshes[i] = (mesh_type*)hookhmes[i];
   
    if (sensor_type == 0) // geom sensor
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      //std::cout << "before __adapt " << std::endl;

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);

      //std::cout << "after __adapt " << std::endl;
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      //std::cout << "before __adapt " << std::endl;

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);

      //std::cout << "after __adapt " << std::endl;
    }
    else if (sensor_type == 1 || sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
  }
  else if (elt_type==elt_t::TETRA)
  {
    using ELT_type = K_MESH::Tetrahedron;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, STYPE>;

    std::vector<mesh_type*> hmeshes(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) hmeshes[i] = (mesh_type*)hookhmes[i];
    
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      //std::cout << "before __adapt " << std::endl;

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);

      //std::cout << "after __adapt " << std::endl;
    }
    else if (sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
  }
  else if (elt_type==elt_t::PRISM3)
  {
    using ELT_type = K_MESH::Prism;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, STYPE>;

    std::vector<mesh_type*> hmeshes(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) hmeshes[i] = (mesh_type*)hookhmes[i];
    
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      //std::cout << "before __adapt " << std::endl;

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);

      //std::cout << "after __adapt " << std::endl;
    }
    else if (sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
  }
  else if (elt_type==elt_t::BASIC)
  {
    using ELT_type = K_MESH::Basic;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, STYPE>;

    std::vector<mesh_type*> hmeshes(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) hmeshes[i] = (mesh_type*)hookhmes[i];

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported

    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      //std::cout << "before __adapt " << std::endl;

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);

      //std::cout << "after __adapt " << std::endl;
    }
    else if (sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
  }
  return err;
}

template <>
int __adapt_wrapper<NUGA::ISO_HEX>
(int elt_type/*dummy*/, int sensor_type,
std::vector<void*>& hookhmes, std::vector<void*>& hooksensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_zone_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  int err(0);
  /*size_t nb_meshes = hookhmes.size();
  
  using ELT_type = K_MESH::Polyhedron<0>;
  using mesh_type = NUGA::hierarchical_mesh<ELT_type, NUGA::ISO_HEX>;

  std::vector<mesh_type*> hmeshes(nb_meshes);
  for (size_t i=0; i < nb_meshes; ++i) hmeshes[i] = (mesh_type*)hookhmes[i];

  if (sensor_type == 0) // geom sensor
  {
    using sensor_t = NUGA::geom_sensor<mesh_type>;

    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

    err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }
  else if (sensor_type == 1) // xsensor
  {
    using sensor_t = NUGA::xsensor<ELT_type, mesh_type>;
    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];
  
     err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }
  else if (sensor_type == 2) //nodal sensor
  {
    using sensor_t = NUGA::nodal_sensor<mesh_type>;
   std::vector<sensor_t*> sensors(nb_meshes);
   for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

   err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }
  else if (sensor_type == 3) //cell sensor
  {
    using sensor_t = NUGA::cell_sensor<mesh_type>;

    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

    //std::cout << "before __adapt " << std::endl;

    err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);

    //std::cout << "after __adapt " << std::endl;
  }
  else if (sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
*/
  return err;
}

template <>
int __adapt_wrapper<NUGA::DIR>
(int elt_type/*dummy*/, int sensor_type,
std::vector<void*>& hookhmes, std::vector<void*>& hooksensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_zone_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  int err(0);
  /*size_t nb_meshes = hookhmes.size();
  
  using ELT_type = K_MESH::Hexahedron;
  using mesh_type = NUGA::hierarchical_mesh<ELT_type, NUGA::DIR>;

  std::vector<mesh_type*> hmeshes(nb_meshes);
  for (size_t i=0; i < nb_meshes; ++i)
    hmeshes[i] = (mesh_type*)(hookhmes[i]);

  if (sensor_type == 0) // geom sensor
  {
    using sensor_t = NUGA::geom_sensor<mesh_type>;

    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i)
      sensors[i] = (sensor_t*)(hooksensors[i]);

    err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }
  else if (sensor_type == 1) // xsensor
  {
    using sensor_t = NUGA::xsensor<ELT_type, mesh_type>;
    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i)
      sensors[i] = (sensor_t*)(hooksensors[i]);
  
     err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }
  else if (sensor_type == 2) //nodal sensor
  {
    using sensor_t = NUGA::nodal_sensor<mesh_type>;
   std::vector<sensor_t*> sensors(nb_meshes);
   for (size_t i=0; i < nb_meshes; ++i)
     sensors[i] = (sensor_t*)(hooksensors[i]);

   err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }
  else if (sensor_type == 3) //cell sensor
  {
    using sensor_t = NUGA::cell_sensor<mesh_type>;

    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i)
      sensors[i] = (sensor_t*)(hooksensors[i]);

    //std::cout << "before __adapt " << std::endl;

    err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);

    //std::cout << "after __adapt " << std::endl;
  }
  else if (sensor_type == 4) // xsensor2
  {
    using sensor_t = NUGA::xsensor2<mesh_type>;

    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

    err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }
*/
  return err;
}

//=============================================================================
/* Initialize the mesh (shift_gem, roerient..) */
//=============================================================================
PyObject* K_INTERSECTOR::initForAdaptCells(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_dict_transfo_to_list;

  if (!PyArg_ParseTuple(args, "OO", &arr, &py_dict_transfo_to_list)) return nullptr;

  // 1. Get mesh and check is NGON
  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check mesh is NGON
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return nullptr;
  
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  // std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  // std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  //2. dico to map
  std::map<transf_t, std::vector<int>> transfo_to_list;
  err = convert_dico_to_map___transfo_to_vecint(py_dict_transfo_to_list, transfo_to_list);
  if (err)
  {
    std::cout << "adaptCells_mpi : input is not a dictionary" << std::endl;
    return nullptr;
  }


  // We reorient the PG of our NGON
  ngi.flag_externals(1);

  DELAUNAY::Triangulator dt;
  bool has_been_reversed;
  err = ngon_type::reorient_skins(dt, crd, ngi, has_been_reversed); //orientate normal outwards
  if (err) return nullptr;

  for (auto& transfo_map : transfo_to_list) //loop over each transformation (match, transla, rota, etc)
  {
    auto & transfo = transfo_map.first.t;
    auto & ptlist  = transfo_map.second;

    size_t sz = ptlist.size();

    E_Float center_rota[3];
    E_Float axis_rota[3];
    for (E_Int i = 0; i < 3; ++i)
    {
      center_rota[i] = transfo[i];
      axis_rota[i] = transfo[i+3];
    }

    E_Float angle = NUGA::normalize<3>(axis_rota);
    angle *= NUGA::PI / 180; //degree --> radian

    if (angle != 0.) //positive rotation periodicity
    {
      auto crd_tmp = crd;

      NUGA::axial_rotate(crd_tmp, center_rota, axis_rota, -angle);

      for (size_t i=0; i < sz; i++)
      {
        E_Int pti = ptlist[i] -1;

        E_Int* nodes = ngi.PGs.get_facets_ptr(pti);
        int nnodes = ngi.PGs.stride(pti);

        K_MESH::Polygon::shift_geom(crd_tmp, nodes, nnodes, 1);
      }
    }

    else //translation, match, negative rotation periodicity
    {
      for (size_t i=0; i < sz; i++)
      {
        E_Int pti = ptlist[i] -1;

        E_Int* nodes = ngi.PGs.get_facets_ptr(pti);
        int nnodes = ngi.PGs.stride(pti);

        K_MESH::Polygon::shift_geom(crd, nodes, nnodes, 1);
      }
    }
  }

  K_FLD::IntArray ng_arr;
  ngi.export_to_array(ng_arr);
  
  PyObject* m = K_ARRAY::buildArray(crd, varString, ng_arr, 8, "NGON", false);

  return m;
}


//=============================================================================
/* Hierarchical Mesh Adaptation : MPI version (has MPI calls) */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCells_mpi(PyObject* self, PyObject* args)
{

  // std::cout << "adaptCells : begin" << std::endl;
  PyObject *hook_hmeshes(nullptr), *hook_sensors(nullptr), *py_zone_to_rid_to_list_owned(nullptr);
  PyObject *py_zonerank(nullptr), *py_rid_to_zones(nullptr);
  MPI_Comm COM = MPI_COMM_WORLD;

  if (!PyArg_ParseTuple(args, "OOOOO", &hook_hmeshes, &hook_sensors, &py_zone_to_rid_to_list_owned, &py_zonerank, &py_rid_to_zones)) return NULL;
  //std::cout << "adaptCells : after parse tuple" << std::endl;

  // 1. GET MESHES AND SENSORS

  int nb_meshes{1};
  int nb_sensors{1};
  bool input_is_list{false};
  if (PyList_Check(hook_hmeshes))
  {
    nb_meshes = PyList_Size(hook_hmeshes);
    nb_sensors = PyList_Size(hook_sensors);
    input_is_list=true;
    //std::cout << "is list (m/s) : " << nb_meshes << "/" << nb_sensors << std::endl;
  }

  if (nb_meshes != nb_sensors)
  {
    PyErr_SetString(PyExc_ValueError,
       "adaptCells: nb of sensors and meshes does not match.");
    return nullptr;
  }

  std::vector<void*> hmeshes, sensors;
  //for unpacking hmeshes
  int* elt_type{ nullptr }, *subdiv_type{ nullptr }, *hook_id{ nullptr };
  std::string* vString{ nullptr };
  //for unpacking sensors
  int *hook_ss_id{ nullptr }, *sensor_type{ nullptr }, *smoothing_type{ nullptr }, *zid{nullptr};
  int *subdiv_type_ss{ nullptr }, *elt_type_ss{ nullptr };

  //std::cout << "adaptCells : before loop" << std::endl;
  std::vector<E_Int> zids;
  for (int m = 0; m < nb_meshes; ++m)
  {
    PyObject* hook_hm = nullptr;
    if (input_is_list)
      hook_hm = PyList_GetItem(hook_hmeshes, m);
    else hook_hm = hook_hmeshes;

    // Unpack hmesh hook
    // ==================
    void** packet{ nullptr };

    //std::cout << " unpack hmesh : " <<  hook_hm << std::endl;
    void * hmesh = unpackHMesh(hook_hm, hook_id, subdiv_type, elt_type, zid, vString, packet);
    if (hmesh == nullptr) return nullptr;
    //std::cout << " unpack hmesh OK " << std::endl;

    hmeshes.push_back(hmesh);
    zids.push_back(*zid);

    PyObject* hook_ss = nullptr;
    if (input_is_list)
      hook_ss = PyList_GetItem(hook_sensors, m);
    else hook_ss = hook_sensors;

    // Unpack sensor hook
    // ==================
    void** packet_ss{ nullptr };
    //std::cout << " unpack sensor" << std::endl;
    void* sensor = unpackSensor(hook_ss, hook_ss_id, sensor_type, smoothing_type, subdiv_type_ss, elt_type_ss, packet_ss);
    //std::cout << " unpackSensor OK " << std::endl;

    //std::cout << "sensor ptr retrieved : " << sensor << std::endl;
    sensors.push_back(sensor);

    // Check basic elements presence     # to do in createHMesh step? 
    // =============================
    if ( (*subdiv_type != NUGA::ISO_HEX) && *elt_type==-1) // K_MESH::Polyhedron<0>::eType::UNKN
    {
      PyErr_SetString(PyExc_ValueError,
         "adaptCells: input mesh to adapt must have basic elements for that subdivision policy (ISO/DIR).");
      return nullptr;
    }

    if (*elt_type_ss != *elt_type)
    {
      PyErr_SetString(PyExc_ValueError,
        "adaptCells: type of elements : inconsistency between sensor and hmesh.");
      std::cout << "hm/sens : " << *elt_type << " vs " << *elt_type_ss << std::endl;
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
  }

  // 2. GET ZONERANK

  std::vector<int> zonerank;
  if (PyDict_Check(py_zonerank))
  {
    E_Int nranks = PyDict_Size(py_zonerank);
    for (E_Int r = 0; r < nranks; ++r)
    {
      PyObject* key = PyInt_FromLong ((long) r);
      PyObject* py_rank = PyDict_GetItem(py_zonerank,key);
      assert (py_rank);
      int rank = (int) PyInt_AsLong(py_rank);
      //std::cout << "key/item : " << r << "/" << my_val <<std::endl;//<< *item << std::endl;
      zonerank.push_back(rank);
    }
  }

  // 3. GET zone_to_rid_to_list_owned
  std::map<int, std::map<int, std::vector<E_Int>>> zone_to_rid_to_list_owned;
  convert_dico_to_map___int_int_vecint(py_zone_to_rid_to_list_owned, zone_to_rid_to_list_owned);
  //assert (zone_to_zone_to_list_owned == nb_meshes);

  /*std::cout << "adaptCells : before __adapt_wrapper" << std::endl;
  std::cout << "sub type : " << *subdiv_type << std::endl;
  std::cout << "elt_type : " << *elt_type << std::endl;
  std::cout << "sensor_type : " << *sensor_type << std::endl;
  std::cout << "vString : " << vString << std::endl;
  std::cout << "hmeshes : " << hmeshes.size() << std::endl;
  std::cout << "sensors : " << sensors.size() << std::endl;*/

  // 4. py_rid_to_zones => rid_to_zones
  std::map<int, std::pair<int,int>> rid_to_zones;
  convert_dico_to_map__int_pairint(py_rid_to_zones, rid_to_zones);
  // assert (zone_to_zone_to_list_owned == nb_meshes);

  // Adaptation
  // ==========
  PyObject *l(PyList_New(0));
  std::vector<E_Int> dummy;

  E_Int err(0);
  if (*subdiv_type == NUGA::ISO)
    err = __adapt_wrapper<ISO>(*elt_type, *sensor_type, hmeshes, sensors, rid_to_zones, zonerank, zone_to_rid_to_list_owned, COM, vString->c_str(), l);
  /*else if (*subdiv_type == NUGA::ISO_HEX)
    err = __adapt_wrapper<ISO_HEX>(*elt_type, *sensor_type, hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, vString->c_str(), l);
  else if (*subdiv_type == NUGA::DIR)
    err = __adapt_wrapper<DIR>(*elt_type, *sensor_type, hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, vString->c_str(), l);
  */

  return (err) ? nullptr : l;
}

//=============================================================================
/* Conformize and dump the enabled cells in a hierarchical mesh : (has MPI calls) */
//=============================================================================
template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void __conformizeHM(const void* hmesh_ptr, K_FLD::FloatArray& crdo, K_FLD::IntArray& cnto, 
                    std::map<E_Int, std::vector<E_Int>>& rid_to_ptlist,
                    std::map<int, std::pair<int,int>>& rid_to_zones,
                    std::vector<std::vector<E_Int>>&     bcptlists,
                    std::vector<std::vector<E_Float>>& fieldsC,
                    std::vector<std::vector<E_Float>>& fieldsN,
                    std::vector<std::vector<E_Float>>& fieldsF)
{
  using mesh_type = NUGA::hierarchical_mesh<ELT_t, STYPE>;

  mesh_type* hmesh = (mesh_type*)hmesh_ptr;

  // update the bc pointlists in hmesh
  for (size_t b=0; b < bcptlists.size(); ++b)
  {
    hmesh->update_pointlist(bcptlists[b]);
  }

  int zid = hmesh->zid;

  for (auto& r_to_ptl: rid_to_ptlist)
  {
    int rid = r_to_ptl.first;
    auto itr2z =  rid_to_zones.find(rid);
    assert (itr2z != rid_to_zones.end());
    int jzid = (itr2z->second.first == zid) ? itr2z->second.second : itr2z->second.first; // get opposed zone

    bool reverse = (jzid < zid);
    
    if (jzid == zid) // self join : the first half is not reverse, the second half is
    {
      std::vector<E_Int> ptlL, ptlR;
      int sz = r_to_ptl.second.size()/2;
      ptlL.insert(ptlL.end(), r_to_ptl.second.begin(), r_to_ptl.second.begin() + sz);
      ptlR.insert(ptlR.end(), r_to_ptl.second.begin()+sz, r_to_ptl.second.end());

      hmesh->update_pointlist(ptlL, false);

      hmesh->update_pointlist(ptlR, true);

      r_to_ptl.second = ptlL;
      r_to_ptl.second.insert(r_to_ptl.second.end(), ALL(ptlR));
    }
    else
      hmesh->update_pointlist(r_to_ptl.second, reverse);
  }
  
  if (hmesh->pghids0.empty()) //the 3 are empty/full at the same time
  {
    K_CONNECT::IdTool::init_inc(hmesh->pthids0, hmesh->_nb_pts0);
    K_CONNECT::IdTool::init_inc(hmesh->pghids0, hmesh->_nb_pgs0);
    K_CONNECT::IdTool::init_inc(hmesh->phhids0, hmesh->_nb_phs0);
  }

  ngon_type ngo;
  std::vector<E_Int> pghids1, phhids1, hmpgid_to_confpgid;
  hmesh->conformize(ngo, hmpgid_to_confpgid, pghids1, phhids1);

  // HISTORY BETWEEN PREVIOUS OUTPUT AND CUURENT ONE
  if (hmpgid_to_confpgid.empty()) // hmesh==exported <=> no adaptation
    K_CONNECT::IdTool::init_inc(hmpgid_to_confpgid, hmesh->_ng.PGs.size());

  NUGA::history_t histo;
  hmesh->build_histo_between_two_enabled_states(hmesh->pghids0, hmesh->phhids0, pghids1, phhids1, hmpgid_to_confpgid, histo);
  
  // EXPORT MESH

  crdo = hmesh->_crd;
  //std::cout << "sz bfore : " << crdo.cols() << std::endl;
  std::vector<E_Int> ptnids1;
  ngo.compact_to_used_nodes(ngo.PGs, crdo, ptnids1);
  //std::cout << "sz after : " << crdo.cols() << std::endl;

  ngo.export_to_array(cnto);

   // JOINS and BC POINTLIST updates
  
  // convert ids to conformed view
  if (!hmpgid_to_confpgid.empty())
  {
    for (size_t b=0; b < bcptlists.size(); ++b)
      for (size_t j=0; j < bcptlists[b].size(); ++j)
        bcptlists[b][j] = hmpgid_to_confpgid[bcptlists[b][j]-1]+1;

    for (auto& r_to_ptl: rid_to_ptlist)
    {
      auto & ptl = r_to_ptl.second;
      for (size_t j=0; j < ptl.size(); ++j)
        ptl[j] = hmpgid_to_confpgid[ptl[j]-1]+1;
    }
  }

  // TRANSFER CENTER SOLUTION FIELDS
  hmesh->project_cell_center_sol_order1(hmesh->phhids0, fieldsC);

  // TRANSFER FACE SOLUTION FIELDS
  // todo : something very similar to above : intensive quantities mgt
  //hmesh->project_face_center_sol_order1(hmesh->pghids0, ...);

  // TRANSFER FACE FLAGS (e.g 'face CAD id')
  // rule : agglo set to -1 (if different) , subdiv inherits parent value
  std::vector<E_Int> src_ids;
  K_CONNECT::IdTool::init_inc(src_ids, hmesh->pghids0.size()); //previous state is always a compressed view of the hmesh

  std::vector<std::vector<E_Float>> new_fieldsF(fieldsF.size());
  for (size_t f = 0; f < fieldsF.size(); ++f)
  {
    histo.transfer_pg_colors(src_ids, fieldsF[f], new_fieldsF[f]);
    new_fieldsF[f].resize(ngo.PGs.size(), -1);
  }
  fieldsF = new_fieldsF;

  // TRANSFER NODE SOLUTION FIELDS == NODE FLAGS
  std::vector<std::vector<E_Float>> new_fieldsN(fieldsN.size());
  for (size_t f = 0; f < fieldsN.size(); ++f)
  {
    new_fieldsN[f].resize(crdo.cols(), NUGA::FLOAT_MAX);

    for (size_t i=0; i < hmesh->pthids0.size(); ++i)
    {
      E_Int tgtid = ptnids1[hmesh->pthids0[i]];
      if (tgtid == IDX_NONE) continue;
      new_fieldsN[f][tgtid]=fieldsN[f][i];
    }

  }
  fieldsN = new_fieldsN;

  // update of the current enabling state
  hmesh->pghids0 = pghids1;
  hmesh->phhids0 = phhids1;
  K_CONNECT::IdTool::reverse_indirection(ptnids1, hmesh->pthids0);

}

template <NUGA::eSUBDIV_TYPE STYPE>
void __conformizeHM(E_Int etype, const void* hmesh_ptr, K_FLD::FloatArray& crdo, K_FLD::IntArray& cnto,
                    std::map<E_Int, std::vector<E_Int>>& rid_to_ptlist,
                    std::map<int, std::pair<int,int>>& rid_to_zones,
                    std::vector<std::vector<E_Int>>&     bcptlists,
                    std::vector<std::vector<E_Float>>&   fieldsC,
                    std::vector<std::vector<E_Float>>&   fieldsN,
                    std::vector<std::vector<E_Float>>&   fieldsF)
{
  if (etype == elt_t::HEXA)
    __conformizeHM<K_MESH::Hexahedron, STYPE>(hmesh_ptr, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF);
  else if (etype == (E_Int)elt_t::TETRA)
    __conformizeHM<K_MESH::Tetrahedron, STYPE>(hmesh_ptr, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF);
  else if (etype == (E_Int)elt_t::PRISM3)
    __conformizeHM<K_MESH::Prism, STYPE>(hmesh_ptr, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF);
  else if (etype == (E_Int)elt_t::BASIC)
    __conformizeHM<K_MESH::Basic, STYPE>(hmesh_ptr, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF);
}

template <>
void __conformizeHM<NUGA::ISO_HEX>(E_Int etype, const void* hmesh_ptr, K_FLD::FloatArray& crdo, K_FLD::IntArray& cnto,
                                   std::map<E_Int, std::vector<E_Int>>& rid_to_ptlist,
                                   std::map<int, std::pair<int,int>>& rid_to_zones,
                                   std::vector<std::vector<E_Int>>&     bcptlists,
                                   std::vector<std::vector<E_Float>>&   fieldsC,
                                   std::vector<std::vector<E_Float>>&   fieldsN,
                                   std::vector<std::vector<E_Float>>&   fieldsF)
{
  __conformizeHM<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh_ptr, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF);
}

template <>
void __conformizeHM<NUGA::DIR>(E_Int etype, const void* hmesh_ptr, K_FLD::FloatArray& crdo, K_FLD::IntArray& cnto,
                                   std::map<E_Int, std::vector<E_Int>>& rid_to_ptlist,
                                   std::map<int, std::pair<int,int>>& rid_to_zones,
                                   std::vector<std::vector<E_Int>>&     bcptlists,
                                   std::vector<std::vector<E_Float>>&   fieldsC,
                                   std::vector<std::vector<E_Float>>&   fieldsN,
                                   std::vector<std::vector<E_Float>>&   fieldsF)
{
  __conformizeHM<K_MESH::Hexahedron, NUGA::DIR>(hmesh_ptr, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF);
}

///
PyObject* K_INTERSECTOR::conformizeHMesh2(PyObject* self, PyObject* args)
{ 
  //hooks[i], bcptlists, jzone_to_ptlist, fieldsC, fieldsN, fieldsF
  PyObject* hook, *py_bcptlists, *py_rid_to_ptlist, *py_rid_to_zones;
  PyObject* pyfieldsC, *pyfieldsN, *pyfieldsF;
  
  if (!PyArg_ParseTuple(args, "OOOOOOO", &hook, &py_bcptlists, &py_rid_to_ptlist, &py_rid_to_zones, &pyfieldsC, &pyfieldsN, &pyfieldsF)) return nullptr;

  int* sub_type{ nullptr }, *elt_type{ nullptr }, *hook_id{ nullptr }, *zid(nullptr);
  std::string* vString{ nullptr };
  void** packet{ nullptr };
  void* hmesh = unpackHMesh(hook, hook_id, sub_type, elt_type, zid, vString, packet);

  // GET RID_TO_ZONES MAP 
  std::map<int, std::pair<int,int>> rid_to_zones;
  convert_dico_to_map__int_pairint(py_rid_to_zones, rid_to_zones);

  // BCs
  std::vector<std::vector<E_Int>> bcptlists;
  if (py_bcptlists != Py_None)
  {
    E_Int nb_ptl = PyList_Size(py_bcptlists);
    bcptlists.resize(nb_ptl);

    for (E_Int i=0; i < nb_ptl; ++i)
    {
      PyObject * pyBCptList = PyList_GetItem(py_bcptlists, i);
      E_Int *ptL, size, nfld;
      /*E_Int res2 = */K_NUMPY::getFromPointList(pyBCptList, ptL, size, nfld, true/* shared*/);
      //std::cout << "res2/size/nfld : " << res2 << "/" << size << "/" << nfld << std::endl;

      std::vector<E_Int> vPtL(ptL, ptL+size);
      bcptlists[i] = vPtL;
    }
  }

  // Joins
  std::map<E_Int, std::vector<E_Int>> rid_to_ptlist;
  if (py_rid_to_ptlist != Py_None && PyDict_Check(py_rid_to_ptlist))
  {
    PyObject *py_rid/*key*/, *py_ptlist /*value : ptlist*/;
    Py_ssize_t pos = 0;

    while (PyDict_Next(py_rid_to_ptlist, &pos, &py_rid, &py_ptlist))
    {
      int rid = (int) PyInt_AsLong(py_rid);

      assert (PyArray_Check(py_ptlist) == 1) ; // it s a numpy
      //std::cout << "est ce un numpy ??? " << isnumpy << std::endl;

      PyArrayObject* pyarr = reinterpret_cast<PyArrayObject*>(py_ptlist);

      long ndims = PyArray_NDIM(pyarr);
      assert (ndims == 1); // vector
      npy_intp* dims = PyArray_SHAPE(pyarr);

      E_Int ptl_sz = dims[0];
      
      //long* dataPtr = static_cast<long*>(PyArray_DATA(pyarr));
      E_Int* dataPtr = (E_Int*)PyArray_DATA(pyarr);

      std::vector<E_Int> ptl(ptl_sz);
      for (size_t u=0; u < ptl_sz; ++u) ptl[u] = dataPtr[u];

      //std::cout << "max in C is : " << *std::max_element(ALL(ptl)) << std::endl;

      rid_to_ptlist[rid]=ptl;
    }
  }

  // FIELDS
  bool has_fieldsC = (pyfieldsC != Py_None);
  bool has_fieldsN = (pyfieldsN != Py_None);
  bool has_fieldsF = (pyfieldsF != Py_None);
  
  char* fvarStringsC, *fvarStringsN, *feltType;
  std::vector<std::vector<E_Float>> fieldsC, fieldsN, fieldsF;
  E_Int nfields{0};

  if (has_fieldsC)
  {
    
    E_Int ni, nj, nk;
    K_FLD::FloatArray fldsC;
    K_FLD::IntArray cn;

    E_Int res = 
      K_ARRAY::getFromArray(pyfieldsC, fvarStringsC, fldsC, ni, nj, nk, cn, feltType);

    /*std::cout << "res : " << res << std::endl;
    std::cout << "var : " << fvarStrings[0] << std::endl;
    std::cout << "field C : " << fldsC.rows() << "/" << fldsC.cols() << std::endl;
    std::cout << "cn : " << cn.rows() << "/" << cn.cols() << std::endl;*/

    nfields = fldsC.rows();
    E_Int nvals = fldsC.cols();
    fieldsC.resize(nfields);
    for (E_Int j = 0; j < nfields; ++j)
      fieldsC[j].resize(nvals);

    for (E_Int i = 0; i < nvals; ++i)
      for (E_Int j = 0; j < nfields; ++j)
        fieldsC[j][i] = fldsC(j, i);
  }

  if (has_fieldsN)
  {
    
    E_Int ni, nj, nk;
    K_FLD::FloatArray fldsN;
    K_FLD::IntArray cn;

    E_Int res = 
      K_ARRAY::getFromArray(pyfieldsN, fvarStringsN, fldsN, ni, nj, nk, cn, feltType);

    /*td::cout << "res : " << res << std::endl;
    std::cout << "var : " << fvarStrings[0] << std::endl;
    std::cout << "field N : " << fldsN.rows() << "/" << fldsN.cols() << std::endl;
    std::cout << "cn : " << cn.rows() << "/" << cn.cols() << std::endl;*/

    nfields = fldsN.rows();
    E_Int nvals = fldsN.cols();
    fieldsN.resize(nfields);
    for (E_Int j = 0; j < nfields; ++j)
      fieldsN[j].resize(nvals);

    for (E_Int i = 0; i < nvals; ++i)
      for (E_Int j = 0; j < nfields; ++j)
        fieldsN[j][i] = fldsN(j, i);
  }

  if (has_fieldsF)
  {

    E_Int nfieldsF = PyList_Size(pyfieldsF);

    std::vector<E_Float*> pFid(nfieldsF);

    fieldsF.resize(nfieldsF);

    for (E_Int j = 0; j < nfieldsF; ++j)
    {
      PyObject* fieldFi = PyList_GetItem(pyfieldsF, j);

      FldArrayF* Fid;
      E_Int res = K_NUMPY::getFromNumpyArray(fieldFi, Fid, true);
      pFid[j] = Fid->begin(1);
      fieldsF[j].resize(Fid->getSize());
    }

    for (E_Int j = 0; j < nfieldsF; ++j)
      for (E_Int i = 0; i < fieldsF[j].size(); ++i)
      {
        fieldsF[j][i] = pFid[j][i];
        //std::cout << "pFid[j][i] : " << pFid[j][i] << std::endl;
      }
  }

  // CONFORMIZE, UPDATE POINTLISTS and TRANSFER FIELDS

  PyObject *l(PyList_New(0));
  K_FLD::IntArray cnto;
  std::vector<E_Int> dummy_oids, jzids;
  
  K_FLD::FloatArray crdo;

  if (*sub_type == NUGA::ISO)
    __conformizeHM<NUGA::ISO>(*elt_type, hmesh, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF);
  else if (*sub_type == NUGA::ISO_HEX)
    __conformizeHM<NUGA::ISO_HEX>(*elt_type, hmesh, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF);
  else if (*sub_type == NUGA::DIR)
    __conformizeHM<NUGA::DIR>(*elt_type, hmesh, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF);

  //  0 : pushing out the mesh
  PyObject *tpl = K_ARRAY::buildArray(crdo, vString->c_str(), cnto, -1, "NGON", false);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // 1 : pushing out the pc pointlist List
  PyObject* bc_ptlist_List = PyList_New(0);
  for (size_t b = 0; b < bcptlists.size(); ++b)
  {
    PyObject* np = K_NUMPY::buildNumpyArray(&bcptlists[b][0], bcptlists[b].size(), 1, 0);
    PyList_Append(bc_ptlist_List, np);
    Py_DECREF(np);
  }
  PyList_Append(l, bc_ptlist_List);
  Py_DECREF(bc_ptlist_List);

  // 2 : pushing out joins pointlist map : jzid to ptlist
  PyObject * rid_to_ptlist_dict = PyDict_New();
  for (auto& r_to_ptl : rid_to_ptlist)
  {
    E_Int rid = r_to_ptl.first;
    auto& ptl = r_to_ptl.second;

    PyObject* key = Py_BuildValue("i", rid);
    PyObject* np = K_NUMPY::buildNumpyArray(&ptl[0], ptl.size(), 1, 0);

    PyDict_SetItem(rid_to_ptlist_dict, key, np);
    Py_DECREF(np);
  }
  PyList_Append(l, rid_to_ptlist_dict);
  Py_DECREF(rid_to_ptlist_dict);

  // FIELDS

  // center fields
  {
    PyObject* tpl = Py_None;
    K_FLD::FloatArray farr;
    if (has_fieldsC)
    {
      farr.resize(nfields, fieldsC[0].size());
      for (size_t i=0; i < fieldsC.size(); ++i)
      {
        std::vector<double>& fld = fieldsC[i];
        for (size_t j = 0; j < fld.size(); ++j)farr(i, j) = fld[j];
      }
      tpl = K_ARRAY::buildArray(farr, fvarStringsC, cnto, -1, feltType, false);
    }
    
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  // node fields
  {
    PyObject* tpl = Py_None;
    K_FLD::FloatArray farr;
    if (has_fieldsN)
    {
      farr.resize(nfields, fieldsN[0].size());
      for (size_t i=0; i < fieldsN.size(); ++i)
      {
        std::vector<double>& fld = fieldsN[i];
        for (size_t j = 0; j < fld.size(); ++j)farr(i, j) = fld[j];
      }
      tpl = K_ARRAY::buildArray(farr, fvarStringsN, cnto, -1, feltType, false);
    }

    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  // face fields (currently only fcadid)
  {
    PyObject* tpl = Py_None;
    if (has_fieldsF) // fcadid only for now
      tpl = K_NUMPY::buildNumpyArray(&fieldsF[0][0], fieldsF[0].size(), 1, 0);

    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  return l;
}


//=============================================================================
/* Exchange the owned PointLists : to update PointListDonor (has MPI calls) */
//=============================================================================
PyObject* K_INTERSECTOR::exchangePointLists(PyObject* self, PyObject* args)
{
  // zonerank, Cmpi.rank, Cmpi.size, zone_to_zone_to_list_owned
  PyObject *py_rid_to_zones(nullptr), *py_zonerank(nullptr), *py_zone_to_rid_to_list_owned(nullptr);
  MPI_Comm COM = MPI_COMM_WORLD;
  E_Int rank{0}, nranks{1};

  if (!PYPARSETUPLEI(args, "OOllO", "OOiiO", &py_rid_to_zones, &py_zonerank, &rank, &nranks, &py_zone_to_rid_to_list_owned)) return nullptr;

  // 1. GET ZONERANK
  std::vector<int> zonerank;
  if (PyDict_Check(py_zonerank))
  {
    E_Int nranks = PyDict_Size(py_zonerank);
    for (E_Int r = 0; r < nranks; ++r)
    {
      PyObject* key = PyInt_FromLong ((long) r);
      PyObject* py_rank = PyDict_GetItem(py_zonerank,key);
      assert (py_rank);
      int rank = (int) PyInt_AsLong(py_rank);
      //std::cout << "key/item : " << r << "/" << my_val <<std::endl;//<< *item << std::endl;
      zonerank.push_back(rank);
    }
  }

  // 2. GET POINTLISTS MAP 
  std::map<int, std::map<int, std::vector<E_Int>>> zone_to_rid_to_list_owned;
  convert_dico_to_map___int_int_vecint(py_zone_to_rid_to_list_owned, zone_to_rid_to_list_owned);
  //assert (zone_to_zone_to_list_owned.size() == nb_meshes);

  // 3. GET RID_TO_ZONES MAP 
  std::map<int, std::pair<int,int>> rid_to_zones;
  convert_dico_to_map__int_pairint(py_rid_to_zones, rid_to_zones);

  // 3. EXCHANGE
  std::map<int, std::map<int, std::vector<E_Int>>> zone_to_rid_to_list_opp;
  NUGA::pointlist_msg_type::exchange_pointlists(rid_to_zones, zonerank, COM, rank, nranks, zone_to_rid_to_list_owned, zone_to_rid_to_list_opp);

  // 4. pushing out joins pointlist map : 'zid to jzid to ptlist'
  PyObject * zone_to_rid_to_list_opp_dict = PyDict_New();
  for (auto& z_to_rid_to_ptl : zone_to_rid_to_list_opp)
  {
    E_Int zid = z_to_rid_to_ptl.first;
    auto& rid_to_list_opp = z_to_rid_to_ptl.second;

    PyObject* key_zid = Py_BuildValue("i", zid);

    PyObject* rid_to_list_opp_dict = PyDict_New();

    for (auto& rid_to_ptl : rid_to_list_opp)
    {
      E_Int rid = rid_to_ptl.first;
      auto & ptl = rid_to_ptl.second;
      PyObject* key_rid = Py_BuildValue("i", rid);

      PyObject* np = K_NUMPY::buildNumpyArray(&ptl[0], ptl.size(), 1, 0);
      PyDict_SetItem(rid_to_list_opp_dict, key_rid, np);
      Py_DECREF(np);
    }

    PyDict_SetItem(zone_to_rid_to_list_opp_dict, key_zid, rid_to_list_opp_dict);
    Py_DECREF(rid_to_list_opp_dict);
  }

  return zone_to_rid_to_list_opp_dict;

}
