/*    
    Copyright 2013-2022 Onera.

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
#include "Nuga/include/xsensor.hxx"
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


///
template <NUGA::eSUBDIV_TYPE STYPE>
void* __createHM(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int zid);

// ISO strategy
template<>
void* __createHM<NUGA::ISO>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int zid)
{
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

    hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt));
    hm->zid = zid;
    return hm;
  }
  else if (typ == elt_t::TETRA)
  {
    using elt_type = K_MESH::Tetrahedron;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;

    hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt));
    hm->zid = zid;
    return hm;
  }
  else if (typ == elt_t::PRISM3)
  {
    using elt_type = K_MESH::Prism;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;

    hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt));
    hm->zid = zid;
    return hm;
  }
  else if (typ == elt_t::BASIC)
  {
    using elt_type = K_MESH::Basic;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;

    hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt));
    hm->zid = zid;
    return hm;
  }
  return NULL;
}

// ISO_HEX strategy
template<>
void* __createHM<NUGA::ISO_HEX>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int zid)
{
  using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>;

  hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt));
  hm->zid = zid;
  return hm;
}

// DIR strategy
template<>
void* __createHM<NUGA::DIR>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int zid)
{
  if (typ != elt_t::HEXA)
  {
    PyErr_WarnEx(PyExc_Warning,
      "createHMesh: directionnal policy is only supported with Hexahedral mesh currently.", 1);
    return nullptr;
  }

  using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::DIR>;
  
  hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt));
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
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<int>>>& zone_to_zone_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  if (hmeshes.empty()) return 0;
  if (sensors.empty()) return 0;

  size_t nb_meshes = hmeshes.size();
  if (nb_meshes != sensors.size()) return 1;

  using adaptor_t = NUGA::adaptor_mpi<hmesh_t, sensor_t>;

  adaptor_t::run(hmeshes, sensors, zone_to_zone_to_list_owned, zonerank, COM, false/*do_agglo*/);

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
E_Int __adapt_wrapper
(E_Int elt_type, E_Int sensor_type,
std::vector<void*>&hookhmes, std::vector<void*>&hooksensors,
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<int>>>& zone_to_zone_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  E_Int err(0);
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

      //std::cout << "before __adapt " << std::endl;

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);

      //std::cout << "after __adapt " << std::endl;
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
  }
  else if (elt_type==elt_t::TETRA)
  {
    using ELT_type = K_MESH::Tetrahedron;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, STYPE>;

    std::vector<mesh_type*> hmeshes(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) hmeshes[i] = (mesh_type*)hookhmes[i];

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported
    
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
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
  }
  else if (elt_type==elt_t::PRISM3)
  {
    using ELT_type = K_MESH::Prism;
    using mesh_type = NUGA::hierarchical_mesh<ELT_type, STYPE>;

    std::vector<mesh_type*> hmeshes(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) hmeshes[i] = (mesh_type*)hookhmes[i];

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported
    
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_type>;

      std::vector<sensor_t*> sensors(nb_meshes);
      for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
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

      err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, varString, out);
    }
    // else if (sensor_type == 1) //xsensor
    // {
    // }
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
  }
  return err;
}

template <>
E_Int __adapt_wrapper<NUGA::ISO_HEX>
(E_Int elt_type/*dummy*/, E_Int sensor_type,
std::vector<void*>& hookhmes, std::vector<void*>& hooksensors,
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<int>>>& zone_to_zone_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  E_Int err(0);
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
E_Int __adapt_wrapper<NUGA::DIR>
(E_Int elt_type/*dummy*/, E_Int sensor_type,
std::vector<void*>& hookhmes, std::vector<void*>& hooksensors,
std::vector<int>& zonerank,
std::map<int, std::map<int, std::vector<int>>>& zone_to_zone_to_list_owned,
MPI_Comm COM,
const char* varString, PyObject *out)
{
  E_Int err(0);
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
/* Hierarchical Mesh Adaptation */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCells_mpi(PyObject* self, PyObject* args)
{

  //std::cout << "adaptCells : begin" << std::endl;
  PyObject *hook_hmeshes(nullptr), *hook_sensors(nullptr), *py_zone_to_zone_to_list_owned(nullptr);
  PyObject *py_zonerank(nullptr);
  MPI_Comm COM = MPI_COMM_WORLD;

  if (!PyArg_ParseTuple(args, "OOOO", &hook_hmeshes, &hook_sensors, &py_zone_to_zone_to_list_owned, &py_zonerank)) return NULL;
  //std::cout << "adaptCells : after parse tuple" << std::endl;

  // 1. GET MESHES AND SENSORS

  E_Int nb_meshes{1};
  E_Int nb_sensors{1};
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
  E_Int* elt_type{ nullptr }, *subdiv_type{ nullptr }, *hook_id{ nullptr };
  std::string* vString{ nullptr };
  //for unpacking sensors
  E_Int *hook_ss_id{ nullptr }, *sensor_type{ nullptr }, *smoothing_type{ nullptr }, *zid{nullptr};
  E_Int *subdiv_type_ss{ nullptr }, *elt_type_ss{ nullptr };

  //std::cout << "adaptCells : before loop" << std::endl;
  std::vector<E_Int> zids;
  for (E_Int m = 0; m < nb_meshes; ++m)
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

  // 3. GET zone_to_zone_to_list_owned

  std::map<int, std::map<int, std::vector<int>>> zone_to_zone_to_list_owned;

  if (PyDict_Check(py_zone_to_zone_to_list_owned))
  {
    E_Int nzid = PyDict_Size(py_zone_to_zone_to_list_owned);
    assert (nzid == nb_meshes);

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
        //std::cout << "est ce un numpy ??? " << isnumpy << std::endl;

        PyArrayObject* pyarr = reinterpret_cast<PyArrayObject*>(py_ptlist_owned);

        long ndims = PyArray_NDIM(pyarr);
        assert (ndims == 1); // vector
        long* dims = PyArray_SHAPE(pyarr);

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

  /*std::cout << "adaptCells : before __adapt_wrapper" << std::endl;
  std::cout << "sub type : " << *subdiv_type << std::endl;
  std::cout << "elt_type : " << *elt_type << std::endl;
  std::cout << "sensor_type : " << *sensor_type << std::endl;
  std::cout << "vString : " << vString << std::endl;
  std::cout << "hmeshes : " << hmeshes.size() << std::endl;
  std::cout << "sensors : " << sensors.size() << std::endl;*/

  // Adaptation
  // ==========
  PyObject *l(PyList_New(0));
  std::vector<E_Int> dummy;

  E_Int err(0);
  if (*subdiv_type == NUGA::ISO)
    err = __adapt_wrapper<ISO>(*elt_type, *sensor_type, hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, vString->c_str(), l);
  /*else if (*subdiv_type == NUGA::ISO_HEX)
    err = __adapt_wrapper<ISO_HEX>(*elt_type, *sensor_type, hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, vString->c_str(), l);
  else if (*subdiv_type == NUGA::DIR)
    err = __adapt_wrapper<DIR>(*elt_type, *sensor_type, hmeshes, sensors, zonerank, zone_to_zone_to_list_owned, COM, vString->c_str(), l);
*/
  //std::cout << "adaptCells : end" << std::endl;

  return (err) ? nullptr : l;
}
