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
#include "Nuga/include/adaptor_para.hxx"
#include "Nuga/include/hierarchical_mesh.hxx"
#include "Nuga/include/smoother.hxx"

#include "Nuga/include/BbTree.h"
#include <memory>
//#include <iostream>
#include "dico_to_stl.h"
#include "adaptCells_hook.h"

using namespace std;
using namespace NUGA;

#ifdef FLAG_STEP
E_Int chrono::verbose = 0;
#endif

#define HMESH_HOOK_ID 77
#define SENSOR_HOOK_ID 78
#define PACK_SIZE 5

using ngon_type = ngon_t<K_FLD::IntArray>;
using subdiv_t = NUGA::eSUBDIV_TYPE;
using elt_t = K_INTERSECTOR::eType;

# include "mpi.h"
# include "mpi4py/mpi4py.h"

//=============================================================================
/* generic adatptation function */
//=============================================================================
template <typename hmesh_t, typename sensor_t>
E_Int __adapt_lvl0
(std::vector<void*>& hookhmeshes, std::vector<void*>& hooksensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  if (hookhmeshes.empty()) return 0;
  if (hooksensors.empty()) return 0;

  size_t nb_meshes = hookhmeshes.size();
  if (nb_meshes != hooksensors.size()) return 1;

  using para_algo_t = hybrid_para_algo<hmesh_t, E_Int>; //hybrid MPI-OMP

  using adaptor_t = adaptor_para<para_algo_t, hmesh_t, sensor_t>;

  // static cast from void* to known types
  std::vector<hmesh_t*> hmeshes(nb_meshes);
  for (size_t i=0; i < nb_meshes; ++i) hmeshes[i] = (hmesh_t*)hookhmeshes[i];
  std::vector<sensor_t*> sensors(nb_meshes);
  for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

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
std::vector<void*>&hmeshes, std::vector<void*>&sensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_zone_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  int err(0);

  if (elt_type==elt_t::HEXA)
  {
    using mesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, STYPE>;
   
    if (sensor_type == 0) // geom sensor
    {
      using sensor_t = NUGA::geom_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 1 || sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
  }
  else if (elt_type==elt_t::TETRA)
  {
    using mesh_t = NUGA::hierarchical_mesh<K_MESH::Tetrahedron, STYPE>;
    
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 1 || sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
  }
  else if (elt_type==elt_t::PRISM3)
  {
    using mesh_t = NUGA::hierarchical_mesh<K_MESH::Prism, STYPE>;

    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 1 || sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
  }
  else if (elt_type==elt_t::BASIC)
  {
    using mesh_t = NUGA::hierarchical_mesh<K_MESH::Basic, STYPE>;

    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    else if (sensor_type == 1 || sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
  }
  return err;
}

///
template <>
int __adapt_wrapper<NUGA::DIR_PROTO>
(int elt_type, int sensor_type,
std::vector<void*>&hmeshes, std::vector<void*>&sensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_zone_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  E_Int err(0);

  using mesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::DIR_PROTO>;

  if (sensor_type == 0) // geom sensor
  {
    using sensor_t = NUGA::geom_sensor<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }
  else if (sensor_type == 2) //nodal sensor
  {
    using sensor_t = NUGA::nodal_sensor<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }
  else if (sensor_type == 3) //cell sensor
  {
    using sensor_t = NUGA::cell_sensor<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }
  else if (sensor_type == 1 || sensor_type == 4) // xsensor2
  {
    using sensor_t = NUGA::xsensor2<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
  }

  return err;
}

///
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
  assert(false);//todo
  return 1;
}

///
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
  assert(false);//todo
  return 1;
}

//=============================================================================
/* Hierarchical Mesh Adaptation : MPI version (has MPI calls) */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCells_mpi(PyObject* self, PyObject* args)
{

  // std::cout << "adaptCells : begin" << std::endl;
  PyObject *hook_hmeshes(nullptr), *hook_sensors(nullptr), *py_zone_to_rid_to_list_owned(nullptr);
  PyObject *py_zonerank(nullptr), *py_rid_to_zones(nullptr), *mpi4pyCom(nullptr);

  if (!PYPARSETUPLE_(args, OOOO_ OO_, &hook_hmeshes, &hook_sensors, &py_zone_to_rid_to_list_owned, &py_zonerank, &py_rid_to_zones, &mpi4pyCom)) return NULL;
  //std::cout << "adaptCells : after parse tuple" << std::endl;

  void* pt_comm = (void*)&(((PyMPICommObject*)mpi4pyCom)->ob_mpi);
  MPI_Comm COM = *((MPI_Comm*) pt_comm);

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
/* Exchange the owned PointLists : to update PointListDonor (has MPI calls) */
//=============================================================================
PyObject* K_INTERSECTOR::exchangePointLists(PyObject* self, PyObject* args)
{
  // zonerank, Cmpi.rank, Cmpi.size, zone_to_zone_to_list_owned
  PyObject *py_rid_to_zones(nullptr), *py_zonerank(nullptr), *py_zone_to_rid_to_list_owned(nullptr), *mpi4pyCom(nullptr);
  E_Int rank{0}, nranks{1};

  if (!PYPARSETUPLE_(args, OO_ II_ OO_, &py_rid_to_zones, &py_zonerank, &rank, &nranks, &py_zone_to_rid_to_list_owned, &mpi4pyCom)) return nullptr;

  void* pt_comm = (void*)&(((PyMPICommObject*)mpi4pyCom)->ob_mpi);
  MPI_Comm COM = *((MPI_Comm*) pt_comm);

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
