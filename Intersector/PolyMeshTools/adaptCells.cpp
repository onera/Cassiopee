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
#include "Nuga/include/xsensor2.hxx"
#include "Nuga/include/nodal_sensor.hxx"
#include "Nuga/include/cell_sensor.hxx"
#include "Nuga/include/adaptor.hxx"
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

using ngon_type = ngon_t<K_FLD::IntArray>;
using elt_t = K_INTERSECTOR::eType;
using subdiv_t = NUGA::eSUBDIV_TYPE;

///
template <NUGA::eSUBDIV_TYPE STYPE>
void* __createHM(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, std::vector<std::vector<E_Int>>& bcptlists, E_Int zid, std::vector<std::pair<E_Int, std::vector<E_Int>>>& joinlists, void* pcom = nullptr);

// ISO strategy
template<>
void* __createHM<NUGA::ISO>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, std::vector<std::vector<E_Int>>& bcptlists, E_Int zid, std::vector<std::pair<E_Int, std::vector<E_Int>>>& joinlists, void* pcom)
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

    if (pcom == nullptr)
      return new hmesh_t(crd, cnt, 1/*idx start*/, bcptlists);
    else
    {
      using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
      using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

      com_t* com = (com_t*)pcom;
      //std::cout << "joinlists sz : " << joinlists.size() << std::endl;
      //std::cout << "com : " << com << std::endl;
      return new hmesh_t(zid, crd, cnt, bcptlists, joinlists, 1/*idx start*/, com);
    }
  }
  else if (typ == elt_t::TETRA)
  {
    using elt_type = K_MESH::Tetrahedron;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;

    if (pcom == nullptr)
      return new hmesh_t(crd, cnt, 1/*idx start*/, bcptlists);
    else
    {
      using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
      using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

      com_t* com = (com_t*)pcom;
      return new hmesh_t(zid, crd, cnt, bcptlists, joinlists, 1/*idx start*/, com);
    }
  }
  else if (typ == elt_t::PRISM3)
  {
    using elt_type = K_MESH::Prism;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;

    if (pcom == nullptr)
      return new hmesh_t(crd, cnt, 1/*idx start*/, bcptlists);
    else
    {
      using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
      using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

      com_t* com = (com_t*)pcom;
      return new hmesh_t(zid, crd, cnt, bcptlists, joinlists, 1/*idx start*/, com);
    }
  }
  else if (typ == elt_t::BASIC)
  {
    using elt_type = K_MESH::Basic;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;

    if (pcom == nullptr)
      return new hmesh_t(crd, cnt, 1/*idx start*/, bcptlists);
    else
    {
      using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
      using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

      com_t* com = (com_t*)pcom;
      return new hmesh_t(zid, crd, cnt, bcptlists, joinlists, 1/*idx start*/, com);
    }
  }
  return NULL;
}

// ISO_HEX strategy
template<>
void* __createHM<NUGA::ISO_HEX>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, std::vector<std::vector<E_Int>>& bcptlists, E_Int zid, std::vector<std::pair<E_Int, std::vector<E_Int>>>& joinlists, void* pcom)
{
  using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>;
  if (pcom == nullptr)
    return new hmesh_t(crd, cnt, 1/*idx start*/, bcptlists);
  else
  {
    /*using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
      using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

      com_t* com = (com_t*)pcom;
      return new hmesh_t(zid, crd, cnt, bcptlists, joinlists, 1, com);*/
    return nullptr;
  }
}

// DIR strategy
template<>
void* __createHM<NUGA::DIR>(E_Int typ, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, std::vector<std::vector<E_Int>>& bcptlists, E_Int zid, std::vector<std::pair<E_Int, std::vector<E_Int>>>& joinlists, void* pcom)
{
  if (typ != elt_t::HEXA)
  {
    PyErr_WarnEx(PyExc_Warning,
      "createHMesh: directionnal policy is only supported with Hexahedral mesh currently.", 1);
    return nullptr;
  }

  using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::DIR>;
  if (pcom == nullptr)
    return new hmesh_t(crd, cnt, 1/*idx start*/, bcptlists);
  else
  {
    using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
    using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

    com_t* com = (com_t*)pcom;
    return new hmesh_t(zid, crd, cnt, bcptlists, joinlists, 1/*idx start*/, com);
  }
}

//=============================================================================
/* get COM hook  */
//=============================================================================
void* unpackCOM(PyObject* hook, int *&hook_id, int *&subdiv_type, int *&elt_type, void **&packet)
{

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**)PyCObject_AsVoidPtr(hook);
#else
  packet = (void**)PyCapsule_GetPointer(hook, NULL);
#endif

  hook_id           = (int*)packet[0];        // type of hook
 
  if (*hook_id != COM_HOOK_ID)
  {
    PyErr_SetString(PyExc_TypeError,
      "unpackCOM: hook id checking failure.");
    return nullptr;
  }
  
  subdiv_type          = (int*)packet[2];        // subdivision type ISO, ISO_HEX, DIR...  
  elt_type             = (int*)packet[3];        // type of elements in hmesh
  
  return packet[1];
}

//============================================================================
/* Create a (list of) hmesh/hzone(s) and returns a (list of) hook(s) */
//============================================================================
PyObject* K_INTERSECTOR::createHMesh(PyObject* self, PyObject* args)
{
  PyObject* hook;
  void** packet = new void*[HMESH_PACK_SIZE];  // hook_ID, hmesh ptr, subdiv policy, elt type, varString

  int* hookid = new int;  packet[0] = hookid;
  //void* hmesh_ptr = nullptr;  packet[1] = hmesh_ptr;// templated hmesh type to build 
  int* subtype = new int; packet[2] = subtype;
  elt_t* etyp = new elt_t;    packet[3] = etyp;
  std::string* vString = new std::string; packet[4] = vString;
  int* zid = new int; packet[5] = zid;

  *hookid = HMESH_HOOK_ID;

  PyObject *arr;
  
  PyObject *pyBCptlitsts{nullptr}, *pyJzids{nullptr}, *pyJptlists{nullptr}, *hookCom{nullptr};
  if (!PYPARSETUPLEI(args, "OlOlOOO", "OiOiOOO", &arr, subtype, &pyBCptlitsts, zid, &pyJzids, &pyJptlists, &hookCom)) return nullptr;

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

  // BCs
  std::vector<std::vector<E_Int>> bcptlists;
  if (pyBCptlitsts != Py_None)
  {
    E_Int nb_ptl = PyList_Size(pyBCptlitsts);
    bcptlists.resize(nb_ptl);

    for (E_Int i=0; i < nb_ptl; ++i)
    {
      PyObject * pyBCptList = PyList_GetItem(pyBCptlitsts, i);
      E_Int *ptL, size, nfld;
      /*E_Int res2 = */K_NUMPY::getFromNumpyArray(pyBCptList, ptL, size, nfld, true/* shared*/);
      //std::cout << "res2/size/nfld : " << res2 << "/" << size << "/" << nfld << std::endl;

      std::vector<E_Int> vPtL(ptL, ptL+size);
      bcptlists[i] = vPtL;
    }
  }

  // joins
  void* com {nullptr};
  std::vector<std::pair<E_Int, std::vector<E_Int>>> joinlists;
  if (pyJzids != Py_None) 
  {
    E_Int nb_joins = PyList_Size(pyJzids);

    // JOINS

#ifdef DEBUG_2019
    E_Int nb_ptl = PyList_Size(pyJptlists);
    assert (nb_joins == nb_ptl);
#endif

    joinlists.resize(nb_joins);

    for (E_Int i=0; i < nb_joins; ++i)
    {
      PyObject* pyJzid = PyList_GetItem(pyJzids, i);
      E_Int jzid = PyLong_AsLong(pyJzid);
      
      PyObject*     pyJptList = PyList_GetItem(pyJptlists, i);
      E_Int *ptL, size, nfld;
      /*E_Int res = */K_NUMPY::getFromNumpyArray(pyJptList, ptL, size, nfld, true/* shared*/);
      //std::cout << "res2/size/nfld : " << res2 << "/" << size << "/" << nfld << std::endl;

      std::vector<E_Int> vPtL(ptL, ptL+size);
      //std::cout << "passed size for list : " << vPtL.size() << std::endl;
      joinlists[i] = std::make_pair(jzid, vPtL);
    }

    // COM
    assert (hookCom != Py_None);

    int* sub_type{ nullptr }, *elt_type{ nullptr }, *hook_id{ nullptr };
    void** packet{ nullptr };
    com = unpackCOM(hookCom, hook_id, sub_type, elt_type, packet);
  }

  if (*subtype == NUGA::ISO)
    packet[1] = __createHM<NUGA::ISO>(*etyp, crd, cnt, bcptlists, *zid, joinlists, com);
  else if (*subtype == NUGA::ISO_HEX)
    packet[1] = __createHM<NUGA::ISO_HEX>(*etyp, crd, cnt, bcptlists, *zid, joinlists, com);
  else if (*subtype == NUGA::DIR)
    packet[1] = __createHM<NUGA::DIR>(*etyp, crd, cnt, bcptlists, *zid, joinlists, com);

  if (packet[1] == nullptr) return Py_None;// the input mesh does not have basic elts

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif
  
  return hook;
}


//============================================================================
/* Creates a COM */
//============================================================================

template <NUGA::eSUBDIV_TYPE STYPE>
void* __createCOM(E_Int typ, E_Int NBZ);

// ISO strategy
template<>
void* __createCOM<NUGA::ISO>(E_Int typ, E_Int NBZ)
{
  if (typ == elt_t::UNKN)
  {
    PyErr_SetString(PyExc_ValueError,
      "createHMesh: input mesh to adapt must have basic elements and must be in NGON format.");
    return nullptr;
  }
  else if (typ == elt_t::HEXA)
  {
    using elt_type = K_MESH::Hexahedron;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;
    using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
    using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

    com_t* c = new com_t;
    c->agents.resize(NBZ, nullptr);
    return c;
  }
  else if (typ == elt_t::TETRA)
  {
    using elt_type = K_MESH::Tetrahedron;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;
    using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
    using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

    com_t* c = new com_t;
    c->agents.resize(NBZ, nullptr);
    return c;
  }
  else if (typ == elt_t::PRISM3)
  {
    using elt_type = K_MESH::Prism;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;
    using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
    using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

    com_t* c = new com_t;
    c->agents.resize(NBZ, nullptr);
    return c;
  }
  else if (typ == elt_t::BASIC)
  {
    using elt_type = K_MESH::Basic;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO>;
    using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
    using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

    com_t* c = new com_t;
    c->agents.resize(NBZ, nullptr);
    return c;
  }
  return NULL;
}

template<>
void* __createCOM<NUGA::ISO_HEX>(E_Int typ, E_Int NBZ)
{
    using elt_type = K_MESH::Polyhedron<0>;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::ISO_HEX>;
    using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
    using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

    com_t* c = new com_t;
    c->agents.resize(NBZ, nullptr);
    return c;
}

template<>
void* __createCOM<NUGA::DIR>(E_Int typ, E_Int NBZ)
{
    using elt_type = K_MESH::Hexahedron;
    using hmesh_t = NUGA::hierarchical_mesh<elt_type, NUGA::DIR>;
    using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
    using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;

    com_t* c = new com_t;
    c->agents.resize(NBZ, nullptr);
    return c;
}

///
PyObject* K_INTERSECTOR::createCom(PyObject* self, PyObject* args)
{
  PyObject* hook;
  void** packet = new void*[4];  // hook_ID, COM ptr   ////, subdiv policy, elt type

  E_Int* hookid = new E_Int;  packet[0] = hookid;
  E_Int* subtype = new E_Int; packet[2] = subtype;
  elt_t* etyp = new elt_t;    packet[3] = etyp;

  *hookid = COM_HOOK_ID;

  PyObject *arr;
  E_Int NBZ{1};
  if (!PYPARSETUPLEI(args, "Oll", "Oii", &arr, subtype, &NBZ)) return nullptr;

  // mesh
  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char *eltType, *varString;
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);

  std::unique_ptr<K_FLD::FloatArray> afmesh(f);  // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<K_FLD::IntArray> acmesh(cn); // to avoid to call explicit delete at several places in the code.
  if (err) return nullptr;
  //K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  *etyp = elt_t::UNKN;// Polyhedron
  if (*subtype == NUGA::ISO || *subtype == NUGA::DIR)
    *etyp = check_has_NGON_BASIC_ELEMENT(cnt);

  if (*subtype == NUGA::ISO)
    packet[1] = __createCOM<NUGA::ISO>(*etyp, NBZ);
  else if (*subtype == NUGA::ISO_HEX)
    packet[1] = __createCOM<NUGA::ISO_HEX>(*etyp, NBZ);
  else if (*subtype == NUGA::DIR)
    packet[1] = __createCOM<NUGA::DIR>(*etyp, NBZ);

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
void __conformizeHM(const void* hmesh_ptr, K_FLD::FloatArray& crdo, K_FLD::IntArray& cnto, 
                    std::vector<E_Int>& oids/*for BC*/, std::vector<E_Int>& jzids, 
                    std::vector<std::vector<E_Int>>& jptlists,
                    std::vector<std::vector<E_Int>>& bcptlists,
                    std::vector<std::vector<E_Float>>& fieldsC,
                    std::vector<std::vector<E_Float>>& fieldsN,
                    std::vector<std::vector<E_Float>>& fieldsF)
{
  using mesh_type = NUGA::hierarchical_mesh<ELT_t, STYPE>;

  mesh_type* hmesh = (mesh_type*)hmesh_ptr;

  // update the bc pointlists in hmesh
  hmesh->update_BCs();

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
  bcptlists = hmesh->BCptLists;
  if (!hmpgid_to_confpgid.empty())
  {
    for (size_t i=0; i < bcptlists.size(); ++i)
      for (size_t j=0; j < bcptlists[i].size(); ++j)
        bcptlists[i][j] = hmpgid_to_confpgid[bcptlists[i][j]-1]+1;
  }

  if (hmesh->join != nullptr){
    //std::cout << "updating join" << std::endl;
    hmesh->join->update();

    jptlists.resize(hmesh->join->jid_to_joinlist.size());

    size_t i{0};
    for (auto& j : hmesh->join->jid_to_joinlist)
    {
      E_Int jzid = j.first;

      jptlists[i] = j.second;

      if (!hmpgid_to_confpgid.empty())
      {
        for (size_t j=0; j < jptlists[i].size(); ++j)
          jptlists[i][j] = hmpgid_to_confpgid[jptlists[i][j]-1]+1;
      }

      jzids.push_back(jzid);
      ++i;
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
                    std::vector<E_Int>& oids/*for BC*/, std::vector<E_Int>& jzids, 
                    std::vector<std::vector<E_Int>>& jptlists,
                    std::vector<std::vector<E_Int>>& bcptlists,
                    std::vector<std::vector<E_Float>>& fieldsC,
                    std::vector<std::vector<E_Float>>& fieldsN,
                    std::vector<std::vector<E_Float>>& fieldsF)
{
  if (etype == elt_t::HEXA)
    __conformizeHM<K_MESH::Hexahedron, STYPE>(hmesh_ptr, crdo, cnto, oids, jzids, jptlists, bcptlists, fieldsC, fieldsN, fieldsF);
  else if (etype == (E_Int)elt_t::TETRA)
    __conformizeHM<K_MESH::Tetrahedron, STYPE>(hmesh_ptr, crdo, cnto, oids, jzids, jptlists, bcptlists, fieldsC, fieldsN, fieldsF);
  else if (etype == (E_Int)elt_t::PRISM3)
    __conformizeHM<K_MESH::Prism, STYPE>(hmesh_ptr, crdo, cnto, oids, jzids, jptlists, bcptlists, fieldsC, fieldsN, fieldsF);
  else if (etype == (E_Int)elt_t::BASIC)
    __conformizeHM<K_MESH::Basic, STYPE>(hmesh_ptr, crdo, cnto, oids, jzids, jptlists, bcptlists, fieldsC, fieldsN, fieldsF);
}

template <>
void __conformizeHM<NUGA::ISO_HEX>(E_Int etype/*dummy*/, const void* hmesh_ptr, K_FLD::FloatArray& crdo, K_FLD::IntArray& cnto,
                                   std::vector<E_Int>& oids/*for BC*/, 
                                   std::vector<E_Int>& jzids, std::vector<std::vector<E_Int>>& jptlists,
                                   std::vector<std::vector<E_Int>>& bcptlists,
                                   std::vector<std::vector<E_Float>>& fieldsC,
                                   std::vector<std::vector<E_Float>>& fieldsN,
                                   std::vector<std::vector<E_Float>>& fieldsF)
{
  __conformizeHM<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh_ptr, crdo, cnto, oids, jzids, jptlists, bcptlists, fieldsC, fieldsN, fieldsF);
}

template <>
void __conformizeHM<NUGA::DIR>(E_Int etype/*dummy*/, const void* hmesh_ptr, K_FLD::FloatArray& crdo, K_FLD::IntArray& cnto,
                                   std::vector<E_Int>& oids/*for BC*/, 
                                   std::vector<E_Int>& jzids, std::vector<std::vector<E_Int>>& jptlists,
                                   std::vector<std::vector<E_Int>>& bcptlists,
                                   std::vector<std::vector<E_Float>>& fieldsC,
                                   std::vector<std::vector<E_Float>>& fieldsN,
                                   std::vector<std::vector<E_Float>>& fieldsF)
{
  __conformizeHM<K_MESH::Hexahedron, NUGA::DIR>(hmesh_ptr, crdo, cnto, oids, jzids, jptlists, bcptlists, fieldsC, fieldsN, fieldsF);
}

//
PyObject* K_INTERSECTOR::conformizeHMesh(PyObject* self, PyObject* args)
{
  PyObject* hook, *pyfieldsC, *pyfieldsN, *pyfieldsF;
  if (!PyArg_ParseTuple(args, "OOOO", &hook, &pyfieldsC, &pyfieldsN, &pyfieldsF)) return nullptr;

  int* sub_type{ nullptr }, *elt_type{ nullptr }, *hook_id{ nullptr }, *zid(nullptr);
  std::string* vString{ nullptr };
  void** packet{ nullptr };
  void* hmesh = unpackHMesh(hook, hook_id, sub_type, elt_type, zid, vString, packet);


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

    /*std::cout << "res : " << res << std::endl;
    std::cout << "var : " << fvarStringsN << std::endl;
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
  std::vector<std::vector<E_Int>> jptlists, bcptlists;
  K_FLD::FloatArray crdo;

  if (*sub_type == NUGA::ISO)
    __conformizeHM<NUGA::ISO>(*elt_type, hmesh, crdo, cnto, dummy_oids, jzids, jptlists, bcptlists, fieldsC, fieldsN, fieldsF);
  else if (*sub_type == NUGA::ISO_HEX)
    __conformizeHM<NUGA::ISO_HEX>(*elt_type, hmesh, crdo, cnto, dummy_oids, jzids, jptlists, bcptlists, fieldsC, fieldsN, fieldsF);
  else if (*sub_type == NUGA::DIR)
    __conformizeHM<NUGA::DIR>(*elt_type, hmesh, crdo, cnto, dummy_oids, jzids, jptlists, bcptlists, fieldsC, fieldsN, fieldsF);

  //  0 : pushing out the mesh
  PyObject *tpl = K_ARRAY::buildArray(crdo, vString->c_str(), cnto, -1, "NGON", false);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // 1 : pushing an array to tell next ouptuts are join, bc or field stuff

  /*std::cout << "nb j ptlists : " << jptlists.size() << std::endl;
  std::cout << "nb bc ptlists: " << bcptlists.size() << std::endl;
  std::cout << "nb fields C:   " << fieldsC.size() << std::endl;
  std::cout << "nb fields F:   " << fieldsF.size() << std::endl;
  std::cout << "nb fields N:   " << fieldsN.size() << std::endl;*/

  E_Int ranges[] = {3,3,3,3,3,3};
  //ranges[0] = 3; // starting index in l for bcs
  if (jptlists.size()) ranges[1] = ranges[0] + jptlists.size(); //one-pas-the-end
  ranges[2] = ranges[1];
  if (bcptlists.size()) ranges[2] += bcptlists.size(); //one-pas-the-end
  ranges[3] = ranges[2];
  if (fieldsC.size()) ranges[3] += 1/*compacted in one array*/; //one-pas-the-end
  ranges[4] = ranges[3];
  if (fieldsN.size()) ranges[4] += 1;
  ranges[5] = ranges[4];
  if (fieldsF.size()) ranges[5] += 1;

  tpl = K_NUMPY::buildNumpyArray(&ranges[0], 6, 1, 0);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  assert (jzids.size() == jptlists.size());

  // 2 : joined zone ids
  tpl = Py_None;// append event if empty => to force to have l[3] as joins/bc pt list start
  if (!jzids.empty())
    tpl = K_NUMPY::buildNumpyArray(&jzids[0], jzids.size(), 1, 0);

  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  // then 

  // JOINS    : l[3]         -> l[ranges[1]-1]
  // BCs      : l[ranges[1]] -> l[ranges[2]-1]
  // FIELDS C : l[ranges[2]] -> l[ranges[3]-1]
  // FIELDS N : l[ranges[3]] -> l[ranges[4]-1]
  // FIELDS F : l[ranges[4]] -> l[ranges[5]-1]

  // JOINS
  for (size_t i=0; i < jzids.size(); ++i)
  {
    std::vector<E_Int>& new_ptl = jptlists[i];
    tpl = K_NUMPY::buildNumpyArray(&new_ptl[0], new_ptl.size(), 1, 0);

    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  // BCs
  for (size_t i=0; i < bcptlists.size(); ++i)
  {
    std::vector<E_Int>& new_ptl = bcptlists[i];
    tpl = K_NUMPY::buildNumpyArray(&new_ptl[0], new_ptl.size(), 1, 0);

    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  // FIELDS
  if (has_fieldsC)
  {
    K_FLD::FloatArray farr(nfields, fieldsC[0].size());
    for (size_t i=0; i < fieldsC.size(); ++i)
    {
      std::vector<double>& fld = fieldsC[i];
      for (size_t j = 0; j < fld.size(); ++j)farr(i, j) = fld[j];
    }

    PyObject* tpl = K_ARRAY::buildArray(farr, fvarStringsC, cnto, -1, feltType, false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  if (has_fieldsN)
  {
    K_FLD::FloatArray farr(nfields, fieldsN[0].size());
    for (size_t i=0; i < fieldsN.size(); ++i)
    {
      std::vector<double>& fld = fieldsN[i];
      for (size_t j = 0; j < fld.size(); ++j)farr(i, j) = fld[j];
    }

    PyObject* tpl = K_ARRAY::buildArray(farr, fvarStringsN, cnto, -1, feltType, false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  if (has_fieldsF) // fcadid only for now
  {
    PyObject* tpl = K_NUMPY::buildNumpyArray(&fieldsF[0][0], fieldsF[0].size(), 1, 0);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

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

template <>
void __deleteHM<NUGA::DIR>(E_Int etype/*dummy*/, const void* hmesh_ptr)
{
  __deleteHM<K_MESH::Hexahedron, NUGA::DIR>(hmesh_ptr);
}

PyObject* K_INTERSECTOR::deleteHMesh(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PyArg_ParseTuple(args, "O", &hook))
  {
      return NULL;
  }

  // recupere le hook
  int* sub_type{ nullptr }, *elt_type{ nullptr }, *hook_id{ nullptr }, *zid(nullptr);
  std::string* vString{ nullptr };
  void** packet{ nullptr };
  void* hmesh = unpackHMesh(hook, hook_id, sub_type, elt_type, zid, vString, packet);

  if (*sub_type == NUGA::ISO)
    __deleteHM<NUGA::ISO>(*elt_type, hmesh);
  else if (*sub_type == NUGA::ISO_HEX)
    __deleteHM<NUGA::ISO_HEX>(*elt_type, hmesh);
  else if (*sub_type == NUGA::DIR)
    __deleteHM<NUGA::DIR>(*elt_type, hmesh);

  delete hook_id;
  delete vString;
  delete sub_type;
  delete elt_type;
  delete zid;
  delete [] packet;
  Py_INCREF(Py_None);
  return Py_None;
}

//============================================================================
/* getDonnorPG */
//============================================================================
template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
E_Int getDonnorPG(NUGA::hierarchical_mesh<ELT_t, STYPE>& hmesh, std::vector<double>& fieldN, E_Int PG, std::vector<E_Int>& pth2a)
{
  E_Int PGparent = hmesh._PGtree.parent(PG);
  if (PGparent == IDX_NONE) return PG;
  
  E_Int nb_nodes = hmesh._ng.PGs.stride(PG);
  E_Int* p_nodes = hmesh._ng.PGs.get_facets_ptr(PG);
  
  E_Int check_PG = 0;
  for (E_Int n = 0; n < nb_nodes; ++n)
  {
    E_Int Ni = pth2a[p_nodes[n]-1]; //amesh
    //if (PGparent == 0 || PGparent == 6) printf(" Ni (hmesh): %3i \n", Ni);

    if (fieldN[Ni] > 1e99) ++check_PG;
  }
  
  
  if (check_PG == 0)
    return PG;
  else
    return getDonnorPG(hmesh, fieldN, PGparent, pth2a);
}

//============================================================================
/* MAJenfants */
//============================================================================
template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void MAJenfants(NUGA::hierarchical_mesh<ELT_t, STYPE>& hmesh, std::vector<double>& fieldN, E_Int PG, std::vector<E_Int>& pth2a)
{
    
  // Parent information:
  E_Int nb_nodes = hmesh._ng.PGs.stride(PG);
  E_Int* p_nodes = hmesh._ng.PGs.get_facets_ptr(PG);
  
  const E_Int* pchildren = hmesh._PGtree.children(PG);
  
  for (E_Int i = 0; i < nb_nodes; ++i)
  {
      E_Int Ni_parent = pth2a[p_nodes[i]-1]; //amesh (index of fieldN)
      
      //printf("Ni_parent: %3i,   fieldN[Ni_parent]: %3f \n", Ni_parent, fieldN[Ni_parent]);
      
      //PG child
      E_Int child = pchildren[i]; 
      E_Int nb_nodes_child = hmesh._ng.PGs.stride(child); 
      E_Int* p_nodes_child = hmesh._ng.PGs.get_facets_ptr(child);
      
      for (E_Int j = 0; j < nb_nodes_child; ++j)
      {
          E_Int Nj_child = pth2a[p_nodes_child[j]-1]; //amesh (index of fieldN)
          
          if (Ni_parent != Nj_child)
          {
            if (nb_nodes == 3) //TRI
            {
              if(i == 0 || (i == 1 && j == 2) ) fieldN[Nj_child] = fieldN[Ni_parent]/2;
              else fieldN[Nj_child] += fieldN[Ni_parent]/2;
            }
            else if (nb_nodes == 4) //QUAD
            {
              if(i == 0 || (i == 1 && j == 2) || ( i== 2 && j==3)) fieldN[Nj_child] = fieldN[Ni_parent]/2;
              else fieldN[Nj_child] += fieldN[Ni_parent]/2;
              if (i == 3 && j == 1)  fieldN[Nj_child] = fieldN[Nj_child]/2;
            }
          }
          //printf("Nj_child: %3i,   fieldN[Nj_child]: %3f \n", Nj_child, fieldN[Nj_child]);
      }
  }
  
  
  // Check if each child has other child -> Part recursive
  for (E_Int i = 0; i < hmesh._PGtree.nb_children(PG); ++i)
  {
    E_Int child = pchildren[i];
    E_Int nb_child = hmesh._PGtree.nb_children(child);
    if (nb_child > 0) MAJenfants(hmesh, fieldN, child, pth2a);
  }
  
}

//============================================================================
/* InterpolateHMeshNodalField */
//============================================================================
template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void __interpolateHMeshNodalField(const void* hmesh_ptrs, std::vector<double>& fieldN)
{
  using mesh_type = NUGA::hierarchical_mesh<ELT_t, STYPE>;
  mesh_type* hmesh = (mesh_type*)hmesh_ptrs;
  

  //wall_face_ids
  const auto & bcptlists = hmesh->BCptLists; 
  if (bcptlists.empty()) return;
  
  // hmesh2amesh
  std::vector<E_Int> pth2a;
  K_CONNECT::IdTool::reverse_indirection(hmesh->pthids0, pth2a);  

// PARTIE A
  std::set<E_Int> donnorFaces;  
  for(E_Int i = 0; i < bcptlists[0].size(); ++i)
  {
    E_Int PGi_hmesh = bcptlists[0][i] - 1;  

    E_Int PGparent = getDonnorPG(*hmesh, fieldN, PGi_hmesh, pth2a); 
    
    if (PGi_hmesh != PGparent)donnorFaces.insert(PGparent); 

  }
 
// PARTIE B
  
  for (auto it = donnorFaces.begin(); it != donnorFaces.end(); ++it) 
  {
    E_Int PGi_hmesh = *it;
    MAJenfants(*hmesh, fieldN, PGi_hmesh, pth2a);
  }  
}

template <NUGA::eSUBDIV_TYPE STYPE>
void __interpolateHMeshNodalField(E_Int etype, const void* hmesh_ptr, std::vector<double>& fieldN)
{
  if (etype == elt_t::HEXA)
    __interpolateHMeshNodalField<K_MESH::Hexahedron, STYPE>(hmesh_ptr, fieldN);
  else if (etype == (E_Int)elt_t::TETRA)
    __interpolateHMeshNodalField<K_MESH::Tetrahedron, STYPE>(hmesh_ptr, fieldN);
  else if (etype == (E_Int)elt_t::PRISM3)
    __interpolateHMeshNodalField<K_MESH::Prism, STYPE>(hmesh_ptr, fieldN);
  else if (etype == (E_Int)elt_t::BASIC)
    __interpolateHMeshNodalField<K_MESH::Basic, STYPE>(hmesh_ptr, fieldN);
}

template <>
void __interpolateHMeshNodalField<NUGA::ISO_HEX>(E_Int etype/*dummy*/, const void* hmesh_ptr, std::vector<double>& fieldN)
{
  __interpolateHMeshNodalField<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh_ptr, fieldN);
}

template <>
void __interpolateHMeshNodalField<NUGA::DIR>(E_Int etype/*dummy*/, const void* hmesh_ptr, std::vector<double>& fieldN)
{
  __interpolateHMeshNodalField<K_MESH::Hexahedron, NUGA::DIR>(hmesh_ptr, fieldN);
}

PyObject* K_INTERSECTOR::interpolateHMeshNodalField(PyObject* self, PyObject* args)
{
  PyObject *hook{nullptr}, *pyfieldN{nullptr};

  if (!PyArg_ParseTuple(args, "OO", &hook, &pyfieldN))
  {
      return NULL;
  }

  
  FloatArray* fi; E_Int ni, nj, nk;
  IntArray* c;
  char* varString;  char* eltType;
  E_Int ret = K_ARRAY::getFromArray(pyfieldN, varString, fi, ni, nj, nk, c, eltType);
  if (ret != 1 && ret != 2)  {
    PyErr_SetString(PyExc_TypeError, "getFromArray: invalid arrays input.");
    return NULL;
  }
  
  E_Float* pfieldN = fi->begin(0);
  E_Int npts = fi->getSize();
  
  

  std::vector<double> fieldN(npts);
  
  for(E_Int i=0; i<npts; ++i)
      fieldN[i] = pfieldN[i];

  // recupere le hook
  int* sub_type{ nullptr }, *elt_type{ nullptr }, *hook_id{ nullptr }, *zid{nullptr};
  std::string* vString{ nullptr };
  void** packet{ nullptr };
  void* hmesh = unpackHMesh(hook, hook_id, sub_type, elt_type, zid, vString, packet);
  
  if (*sub_type == NUGA::ISO)
    __interpolateHMeshNodalField<NUGA::ISO>(*elt_type, hmesh, fieldN);
  else if (*sub_type == NUGA::ISO_HEX)
    __interpolateHMeshNodalField<NUGA::ISO_HEX>(*elt_type, hmesh, fieldN);
  else if (*sub_type == NUGA::DIR)
    __interpolateHMeshNodalField<NUGA::DIR>(*elt_type, hmesh, fieldN);


  //Retourner le champ mis à jour
  K_FLD::FloatArray farr(1,fieldN.size());
  
  std::vector<double>& fld = fieldN;
  for (size_t j = 0; j < fieldN.size(); ++j)farr(0, j) = fld[j];

  PyObject* tpl = K_ARRAY::buildArray(farr, varString, *c, -1, eltType, false);
  
  delete fi;
  delete c;

  return tpl;
}

//============================================================================
/* Create a Geom sensor and returns a hook */
//============================================================================

template<typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void* __createSensor(void* hmesh, E_Int smoothing_type, E_Int itermax, E_Int metric_policy, E_Int sensor_type)
{
  if (sensor_type == 0)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t = NUGA::geom_sensor<hmesh_t>;
    return new sensor_t(*(hmesh_t*)hmesh, NUGA::eSmoother(smoothing_type), 1/*max_pts_per_cell*/, itermax);
  }
  else if (sensor_type == 2)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t = NUGA::nodal_sensor<hmesh_t>;
    return new sensor_t(*(hmesh_t*)hmesh);
  }
  else if (sensor_type == 3)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t = NUGA::cell_sensor<hmesh_t>;
    return new sensor_t(*(hmesh_t*)hmesh, NUGA::eSmoother(smoothing_type));
  }
  else if (sensor_type == 1 || sensor_type == 4)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t = NUGA::xsensor2<hmesh_t>;
    NUGA::eMetricPolicy policy = (NUGA::eMetricPolicy)metric_policy;
    return new sensor_t(*(hmesh_t*)hmesh, NUGA::eSmoother(smoothing_type), policy, itermax);
  }
  return nullptr;
}

///
template<>
void* __createSensor<K_MESH::Hexahedron, NUGA::DIR>(void* hmesh, E_Int smoothing_type, E_Int itermax, E_Int metric_policy, E_Int sensor_type)
{
  using ELT_t = K_MESH::Hexahedron;

  if (sensor_type == 0)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, NUGA::DIR>;
    using sensor_t = NUGA::geom_sensor<hmesh_t>;

    hmesh_t* hm = (hmesh_t*)(hmesh);

    return new sensor_t(*hm, NUGA::eSmoother(smoothing_type), 1/*max_pts_per_cell*/, itermax);
  }
  else if (sensor_type == 2)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, NUGA::DIR>;
    using sensor_t = NUGA::nodal_sensor<hmesh_t>;

    hmesh_t* hm = (hmesh_t*)(hmesh);

    return new sensor_t(*hm);
  }
  else if (sensor_type == 3)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, NUGA::DIR>;
    using sensor_t = NUGA::cell_sensor<hmesh_t>;

    hmesh_t* hm = (hmesh_t*)(hmesh);

    return new sensor_t(*hm, NUGA::eSmoother(smoothing_type));
  }
  else if (sensor_type == 4)
  {
    using hmesh_t = NUGA::hierarchical_mesh<ELT_t, NUGA::DIR>;
    using sensor_t = NUGA::xsensor2<hmesh_t>;

    hmesh_t* hm = (hmesh_t*)(hmesh);
    NUGA::eMetricPolicy policy = (NUGA::eMetricPolicy)metric_policy;

    return new sensor_t(*hm, NUGA::eSmoother(smoothing_type), policy, itermax);
  }
  return nullptr;
}

template<NUGA::eSUBDIV_TYPE STYPE>
void* __createSensor(E_Int etype, void* hmesh, E_Int sensor_type, E_Int smoothing_type, E_Int itermax, E_Int metric_policy);

template<>
void* __createSensor<NUGA::ISO>(E_Int elt_type, void* hmesh, E_Int sensor_type, E_Int smoothing_type, E_Int itermax, E_Int metric_policy)
{
  if (elt_type == elt_t::HEXA)
    return __createSensor<K_MESH::Hexahedron, NUGA::ISO>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
  else if (elt_type == elt_t::TETRA)
    return __createSensor<K_MESH::Tetrahedron, NUGA::ISO>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
  else if (elt_type == elt_t::PRISM3)
    return __createSensor<K_MESH::Prism, NUGA::ISO>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
  else if (elt_type == elt_t::BASIC)
    return __createSensor<K_MESH::Basic, NUGA::ISO>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
  else 
  {
    PyErr_SetString(PyExc_ValueError, "createSensor: wrong element type in hmesh for ISO strategy.");
    return nullptr;
  }
}

//
template<>
void* __createSensor<NUGA::ISO_HEX>(E_Int etype, void* hmesh, E_Int sensor_type, E_Int smoothing_type, E_Int itermax, E_Int metric_policy)
{
  return __createSensor<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
}


template<>
void* __createSensor<NUGA::DIR>(E_Int etype, void* hmesh, E_Int sensor_type, E_Int smoothing_type, E_Int itermax, E_Int metric_policy)
{
  return __createSensor<K_MESH::Hexahedron, NUGA::DIR>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
}

//============================================================================
/* Creates a Sensor */
//============================================================================
PyObject* K_INTERSECTOR::createSensor(PyObject* self, PyObject* args)
{
  PyObject *hook_sensor(nullptr);

  PyObject *hook_hmesh(nullptr);
  E_Int smoothing_type(0), sensor_type(0), itermax(0); 
  E_Int metric_policy(0); // O:ISO_MIN, 1:ISO_MEAN, 2:ISO_MAX, 3:ISO_MIN OR ISO_MAX 

  if (!PYPARSETUPLEI(args, "Ollll", "Oiiii", &hook_hmesh, &sensor_type, &smoothing_type, &itermax, &metric_policy)) return NULL;


  // // Unpack hmesh hook
  // // ==================
  int* subtype_hm{ nullptr }, *elttype_hm{ nullptr }, *hook_id{ nullptr }, *zid{nullptr};
  std::string* vString{ nullptr };
  void ** packet_h{ nullptr };
  void* hmesh = unpackHMesh(hook_hmesh, hook_id, subtype_hm, elttype_hm, zid, vString, packet_h);

  // // Create packet for sensor hook 
  // // ==============================
  void** packet_ss = new void*[6];  // hook ID, sensor type, hmesh ptr, smoothing type, elt type, subdiv type

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

  if (*sub_type == NUGA::ISO_HEX)
  {
    if (sensor_type == 2)
    {
      PyErr_SetString(PyExc_ValueError,
         "adaptCells: nodal sensor does not support ISO_HEX subdivision.");
      return nullptr;
    }
  }
  else if (*sub_type == NUGA::DIR)
  {
    if (sensor_type != 2 && sensor_type != 3)
    {
      PyErr_SetString(PyExc_ValueError,
       "adaptCells: DIR only works with cell/nodal sensor currently.");
      return nullptr;
    }
  }

  // HMESH PTR
  packet_ss[2] = nullptr;
  if (*subtype_hm == NUGA::ISO)
    packet_ss[2] = __createSensor<NUGA::ISO>(*elt_type, hmesh, sensor_type, smoothing_type, itermax, metric_policy);
  else if (*subtype_hm == NUGA::ISO_HEX)
    packet_ss[2] = __createSensor<NUGA::ISO_HEX>(*elt_type, hmesh, sensor_type, smoothing_type, itermax, metric_policy);
  else if (*subtype_hm == NUGA::DIR)
    packet_ss[2] = __createSensor<NUGA::DIR>(*elt_type, hmesh, sensor_type, smoothing_type, itermax, metric_policy);

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
  if (sensor_type == 0) //geom
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_type   = NUGA::geom_sensor<mesh_type>;
    sensor_type* sensor = (sensor_type*)sensor_ptrs;
    delete sensor;    
  }
  else if (sensor_type == 2) // nodal
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_type   = NUGA::nodal_sensor<mesh_type>;
    sensor_type* sensor = (sensor_type*)sensor_ptrs;
    delete sensor;    
  }
  else if (sensor_type == 3) // cell
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_type   = NUGA::cell_sensor<mesh_type>;
    sensor_type* sensor = (sensor_type*)sensor_ptrs;
    delete sensor;    
  }
  else if (sensor_type == 1 || sensor_type == 4) // xsensor2
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_type   = NUGA::xsensor2<mesh_type>;
    sensor_type* sensor = (sensor_type*)sensor_ptrs;
    delete sensor;    
  }
}

///
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

template <>
void __deleteSensor<NUGA::DIR>(E_Int etype/*dummy*/, E_Int sensor_type,  const void* sensor_ptr)
{
  __deleteSensor<K_MESH::Hexahedron, NUGA::DIR>(sensor_type, sensor_ptr);
}

PyObject* K_INTERSECTOR::deleteSensor(PyObject* self, PyObject* args)
{
  PyObject* hook_sensor;
  if (!PyArg_ParseTuple(args, "O", &hook_sensor))
  {
      return NULL;
  }

  // recupere le hook
  int *hook_ss_id{ nullptr }, *sensor_type{ nullptr }, *smoothing_type{ nullptr };
  int *subdiv_type{ nullptr }, *elt_type{ nullptr };
  void** packet{ nullptr };
  void* sensor = unpackSensor(hook_sensor, hook_ss_id, sensor_type, smoothing_type, subdiv_type, elt_type, packet);  

  if (*subdiv_type == NUGA::ISO)
    __deleteSensor<NUGA::ISO>(*elt_type, *sensor_type, sensor);
  else if (*subdiv_type == NUGA::ISO_HEX)
    __deleteSensor<NUGA::ISO_HEX>(*elt_type, *sensor_type, sensor);
  else if (*subdiv_type == NUGA::DIR)
    __deleteSensor<NUGA::DIR>(*elt_type, *sensor_type, sensor);

  delete hook_ss_id;
  delete sensor_type;
  delete smoothing_type;
  delete subdiv_type;
  delete elt_type;
  delete [] packet;
  Py_INCREF(Py_None);
  return Py_None;
}

// //============================================================================
// /* assign data to a sensor*/
// ============================================================================

template <NUGA::eSUBDIV_TYPE STYP> using incr_t = typename NUGA::adap_incr_type<STYP>::cell_incr_t; // E_Int (ISO) or int_tuple<3> (DIR)

template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void __assign_sensor_data
(int sensor_type, void* psensor, K_FLD::FloatArray& crdS, K_FLD::IntArray& cntS, std::vector<incr_t<STYPE>>& punctual_data)
{
  if (psensor == nullptr) return;

  if (sensor_type == 0) // geom
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t   = NUGA::geom_sensor<mesh_type>;
    sensor_t* sensor = (sensor_t*)psensor;

    sensor->assign_data(crdS);
  }
  else if (sensor_type == 2) // nodal
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t   = NUGA::nodal_sensor<mesh_type>;
    sensor_t* sensor = (sensor_t*)psensor;

    sensor->assign_data(punctual_data);
  }
  else if (sensor_type == 3) // cell
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t   = NUGA::cell_sensor<mesh_type>;
    sensor_t* sensor = (sensor_t*)psensor;

    sensor->assign_data(punctual_data);
  }
  else if (sensor_type == 1 || sensor_type == 4) // xsensor2
  {
    using mesh_type     = NUGA::hierarchical_mesh<ELT_t, STYPE>;
    using sensor_t   = NUGA::xsensor2<mesh_type>;
    sensor_t* sensor = (sensor_t*)psensor;
   
    ngon_type::eGEODIM geodim = ngon_type::get_ngon_geodim(cntS);
    //std::cout << "geodim : " << geodim << std::endl;    

    if (
      (geodim != ngon_type::eGEODIM::SURFACIC_CASSIOPEE) &&
      (geodim != ngon_type::eGEODIM::SURFACIC) &&
      (geodim != ngon_type::eGEODIM::VOLUMIC)
        )
      {
      //std::cout << "wrong type" << std::endl;
      return;
      }

    // so either SURFACIC, SURFACIC_CASSIOPEE or VOLUMIC

    if (geodim == ngon_type::eGEODIM::SURFACIC_CASSIOPEE)
    {
      //std::cout << "surfacic cassiopee" << std::endl;
      ngon_type ng(cntS);
      // convert to SURFACIC (NUGA)
      K_FLD::IntArray cnt1;
      ng.export_surfacic_view(cnt1);
      //std::cout << "exported" << std::endl;
      geodim = ngon_type::eGEODIM::SURFACIC;
      cntS=cnt1;
    }

    if (geodim == ngon_type::eGEODIM::VOLUMIC)
    {
      NUGA::ph_mesh_t mesh(crdS, cntS);
      NUGA::pg_smesh_t data;
      mesh.get_boundary<false>(data);
      sensor->assign_data(data);
    }
    else
    {
      NUGA::pg_smesh_t data(crdS, cntS);
      sensor->assign_data(data);
    }
  }
}

///
template <NUGA::eSUBDIV_TYPE STYPE>
void __assign_sensor_data
(int etype, int sensor_type, void* sensor_ptr, K_FLD::FloatArray& crdS, K_FLD::IntArray& cntS, std::vector<incr_t<STYPE>>& punctual_data)
{
  if (etype == elt_t::HEXA)
    __assign_sensor_data<K_MESH::Hexahedron, STYPE>(sensor_type, sensor_ptr, crdS, cntS, punctual_data);
  else if (etype == (E_Int)elt_t::TETRA)
    __assign_sensor_data<K_MESH::Tetrahedron, STYPE>(sensor_type, sensor_ptr, crdS, cntS, punctual_data);
  else if (etype == (E_Int)elt_t::PRISM3)
    __assign_sensor_data<K_MESH::Prism, STYPE>(sensor_type, sensor_ptr, crdS, cntS, punctual_data);
  else if (etype == (E_Int)elt_t::BASIC)
    __assign_sensor_data<K_MESH::Basic, STYPE>(sensor_type, sensor_ptr, crdS, cntS, punctual_data);
}

template <>
void __assign_sensor_data<NUGA::ISO_HEX>
(int etype, int sensor_type, void* sensor_ptr, K_FLD::FloatArray& crdS, K_FLD::IntArray& cntS, std::vector<E_Int>& punctual_data)
{
  __assign_sensor_data<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(sensor_type, sensor_ptr, crdS, cntS, punctual_data);
}

template <>
void __assign_sensor_data<NUGA::DIR>
(int etype, int sensor_type, void* sensor_ptr, K_FLD::FloatArray& crdS, K_FLD::IntArray& cntS, std::vector<int_tuple<3>>& punctual_data)
{
  __assign_sensor_data<K_MESH::Hexahedron, NUGA::DIR>(sensor_type, sensor_ptr, crdS, cntS, punctual_data);
}


PyObject* K_INTERSECTOR::assignData2Sensor(PyObject* self, PyObject* args)
{

  PyObject* hook_sensor{ nullptr }, *dataSensor(nullptr);
  if (!PyArg_ParseTuple(args, "OO", &hook_sensor, &dataSensor))
    return nullptr;
  if (dataSensor == Py_None) // nothing to do
    return Py_None;

  // recupere le hook
  int *hook_ss_id{ nullptr }, *sensor_type{ nullptr }, *smoothing_type{ nullptr };
  int *subdiv_type{ nullptr }, *elt_type{ nullptr };
  void** packet{ nullptr };
  void* sensor = unpackSensor(hook_sensor, hook_ss_id, sensor_type, smoothing_type, subdiv_type, elt_type, packet);

  std::vector<E_Int> sens_data;
  K_FLD::FloatArray fS;
  K_FLD::IntArray cnS;

  //geom/xsensor or nodal_sensor data ?

  if (PyList_Check(dataSensor)) // Array (mesh or coordinates)
  {
    E_Int ni, nj, nk;
    char* varString, *eltType;
    /*E_Int res = */K_ARRAY::getFromArray(dataSensor, varString, fS, ni, nj, nk, cnS, eltType);
    //std::cout << "res/eltType/fs sz/cns sz : " << res << "/" << eltType << "/" << fS.cols() << "/" << cnS.cols() << std::endl;
  }
  else // assuming numpy for nodal/cell
  {
     E_Int *nodv, size, nfld;
     /*E_Int res2 = */K_NUMPY::getFromNumpyArray(dataSensor, nodv, size, nfld, true/* shared*/);
     sens_data.resize(size);  
     for (E_Int i=0; i < size; ++i)
     {
       E_Int v = nodv[i];
       v = std::max((E_Int)0,v); //disable agglo currently
       sens_data[i] = v;
     }
  }

  if (*subdiv_type == NUGA::ISO)
    __assign_sensor_data<NUGA::ISO>(*elt_type, *sensor_type, sensor, fS, cnS, sens_data);
  else if (*subdiv_type == NUGA::ISO_HEX)
    __assign_sensor_data<NUGA::ISO_HEX>(*elt_type, *sensor_type, sensor, fS, cnS, sens_data);
  else if (*subdiv_type == NUGA::DIR)
  {
    //todo Imad : dealing with directional data set in the python
    std::vector<int_tuple<3>> sens_data_dir(sens_data.size(), int_tuple<3>(0));
    for (size_t u=0; u < sens_data.size(); ++u)
      sens_data_dir[u] = sens_data[u]; // pass a single value to the 3-tuple

    __assign_sensor_data<NUGA::DIR>(*elt_type, *sensor_type, sensor, fS, cnS, sens_data_dir);
  }

  Py_INCREF(Py_None);
  return Py_None;
}

//============================================================================
/* Deletes a COM */
//============================================================================
template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void __deleteCOM(const void* com_ptr)
{
  using hmesh_t = NUGA::hierarchical_mesh<ELT_t, STYPE>;
  using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
  using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;
  
  com_t* com = (com_t*)com_ptr;
  delete com;
}

template <>
void __deleteCOM<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(const void* com_ptr)
{
  using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>;
  using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
  using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;
  
  com_t* com = (com_t*)com_ptr;
  delete com;
}

template <>
void __deleteCOM<K_MESH::Hexahedron, NUGA::DIR>(const void* com_ptr)
{
  using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::DIR>;
  using exdata_t = typename NUGA::join_sensor<hmesh_t>::input_t;
  using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<hmesh_t, exdata_t>>;
  
  com_t* com = (com_t*)(com_ptr);
  delete com;
}

template <NUGA::eSUBDIV_TYPE STYPE>
void __deleteCOM(E_Int etype, const void* com_ptr)
{
  if (etype == elt_t::HEXA)
    __deleteCOM<K_MESH::Hexahedron, STYPE>(com_ptr);
  else if (etype == (E_Int)elt_t::TETRA)
    __deleteCOM<K_MESH::Tetrahedron, STYPE>(com_ptr);
  else if (etype == (E_Int)elt_t::PRISM3)
    __deleteCOM<K_MESH::Prism, STYPE>(com_ptr);
  else if (etype == (E_Int)elt_t::BASIC)
    __deleteCOM<K_MESH::Basic, STYPE>(com_ptr);
}

template <>
void __deleteCOM<NUGA::ISO_HEX>(E_Int etype, const void* com_ptr)
{
  __deleteCOM<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(com_ptr);
}

template <>
void __deleteCOM<NUGA::DIR>(E_Int etype, const void* com_ptr)
{
  __deleteCOM<K_MESH::Hexahedron, NUGA::DIR>(com_ptr);
}

PyObject* K_INTERSECTOR::deleteCom(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PyArg_ParseTuple(args, "O", &hook))
  {
      return NULL;
  }

  // recupere le hook
  int* sub_type{ nullptr }, *elt_type{ nullptr }, *hook_id{ nullptr };
  void** packet{ nullptr };
  void* com = unpackCOM(hook, hook_id, sub_type, elt_type, packet);
  
  if (*sub_type == NUGA::ISO)
    __deleteCOM<NUGA::ISO>(*elt_type, com);
  else if (*sub_type == NUGA::ISO_HEX)
    __deleteCOM<NUGA::ISO_HEX>(*elt_type, com);
  else if (*sub_type == NUGA::DIR)
    __deleteCOM<NUGA::DIR>(*elt_type, com);

  delete hook_id;
  delete sub_type;
  delete elt_type;
  delete [] packet;
  Py_INCREF(Py_None);
  return Py_None;
}

template <typename MESH_t, typename sensor_t>
E_Int __adapt
(std::vector<MESH_t*>& hmeshes, std::vector<sensor_t*>& sensors, const char* varString, PyObject *out)
{
  if (hmeshes.empty()) return 0;
  if (sensors.empty()) return 0;

  size_t nb_meshes = hmeshes.size();
  if (nb_meshes != sensors.size()) return 1;

  //std::cout << "nb_meshes ? " << nb_meshes << std::endl;

  if (nb_meshes == 1)
  {
    if (hmeshes[0] == nullptr) return 1;
    if (sensors[0] == nullptr) return 1;

    MESH_t* hmesh = hmeshes[0];

    hmesh->init();
    NUGA::adaptor<MESH_t, sensor_t>::run(*hmesh, *(sensors[0]));

    std::vector<E_Int> oids;
    K_FLD::IntArray cnto;

    hmesh->_ng.export_to_array(cnto);

    // pushing out the mesh
    PyObject *tpl = K_ARRAY::buildArray(hmesh->_crd, varString, cnto, -1, "NGON", false);
    PyList_Append(out, tpl);
    Py_DECREF(tpl);
  }
  else //mutliseq
  {
    using exdata_t = typename NUGA::join_sensor<MESH_t>::input_t;
    using com_t = typename NUGA::communicator<NUGA::jsensor_com_agent<MESH_t, exdata_t>>;

    //std::cout << "getting com" << std::endl;
    com_t* COM = hmeshes[0]->COM;
    //std::cout << "com : " << COM << std::endl;
    //std::cout << "MULTISEQ!!!" << std::endl;
    NUGA::adaptor<MESH_t, sensor_t>::run(hmeshes, sensors, false/*do_agglo*/, NUGA::ePara::SEQ, COM);
    //std::cout << "buildArray loop" << std::endl;
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
  }

  //std::cout << "__adapt : DONE" << std::endl;

  return 0;
}

///
template <typename ELT_t, subdiv_t STYPE>
E_Int __adapt_wrapper_lvl1
(E_Int sensor_type,
const char* varString, PyObject *out,
std::vector<void*>&hookhmes, std::vector<void*>&hooksensors)
{
  using mesh_type = NUGA::hierarchical_mesh<ELT_t, STYPE>;
  
  E_Int err = 0;
  size_t nb_meshes = hooksensors.size();

  std::vector<mesh_type*> hmeshes(nb_meshes);
  for (size_t i=0; i < nb_meshes; ++i) hmeshes[i] = (mesh_type*)hookhmes[i];

  if (sensor_type == 0) // geom sensor
  {
    using sensor_t = NUGA::geom_sensor<mesh_type>;

    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

    err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, varString, out);
  }
  else if (sensor_type == 2) //nodal sensor
  {
    using sensor_t = NUGA::nodal_sensor<mesh_type>;

    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

    err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, varString, out);
  }
  else if (sensor_type == 3) //cell sensor
  {
    using sensor_t = NUGA::cell_sensor<mesh_type>;

    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

    //std::cout << "before __adapt " << std::endl;

    err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, varString, out);

    //std::cout << "after __adapt " << std::endl;
  }
  else if (sensor_type == 1 || sensor_type == 4) // xsensor2
  {
    using sensor_t = NUGA::xsensor2<mesh_type>;

    std::vector<sensor_t*> sensors(nb_meshes);
    for (size_t i=0; i < nb_meshes; ++i) sensors[i] = (sensor_t*)hooksensors[i];

    err = __adapt<mesh_type, sensor_t>(hmeshes, sensors, varString, out);
  }

  return err;
}

///
template <subdiv_t STYPE>
E_Int __adapt_wrapper_lvl0
(E_Int elt_type, E_Int sensor_type,
const char* varString, PyObject *out,
std::vector<void*>&hookhmes, std::vector<void*>&hooksensors)
{
  E_Int err(0);

  if (elt_type==elt_t::HEXA)
    err = __adapt_wrapper_lvl1<K_MESH::Hexahedron, STYPE>(sensor_type, varString, out, hookhmes, hooksensors);
  else if (elt_type==elt_t::TETRA)
    err = __adapt_wrapper_lvl1<K_MESH::Tetrahedron, STYPE>(sensor_type, varString, out, hookhmes, hooksensors);
  else if (elt_type==elt_t::PRISM3)
    err = __adapt_wrapper_lvl1<K_MESH::Prism, STYPE>(sensor_type, varString, out, hookhmes, hooksensors);
  else if (elt_type==elt_t::BASIC)
    err = __adapt_wrapper_lvl1<K_MESH::Basic, STYPE>(sensor_type, varString, out, hookhmes, hooksensors);

  return err;
}

template <>
E_Int __adapt_wrapper_lvl0<NUGA::ISO_HEX>
(E_Int elt_type, E_Int sensor_type,
const char* varString, PyObject *out,
std::vector<void*>&hookhmes, std::vector<void*>&hooksensors)
{
  E_Int err(0);
  assert (elt_type == elt_t::UNKN);

  err = __adapt_wrapper_lvl1<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(sensor_type, varString, out, hookhmes, hooksensors);

  return err;
}

template <>
E_Int __adapt_wrapper_lvl0<NUGA::DIR>
(E_Int elt_type, E_Int sensor_type,
const char* varString, PyObject *out,
std::vector<void*>&hookhmes, std::vector<void*>&hooksensors)
{
  E_Int err(0);
  assert (elt_type == elt_t::HEXA);

  err = __adapt_wrapper_lvl1<K_MESH::Hexahedron, NUGA::DIR>(sensor_type, varString, out, hookhmes, hooksensors);

  return err;
}

//=============================================================================
/* Hierarchical Mesh Adaptation */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCells(PyObject* self, PyObject* args)
{
  //std::cout << "adaptCells : begin" << std::endl;
  PyObject *hook_hmeshes(nullptr), *hook_sensors(nullptr);

  if (!PyArg_ParseTuple(args, "OO", &hook_hmeshes, &hook_sensors)) return NULL;

  //std::cout << "adaptCells : after parse tuple" << std::endl;

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
  int* elt_type{ nullptr }, *subdiv_type{ nullptr }, *hook_id{ nullptr }, *zid{nullptr};
  std::string* vString{ nullptr };
  //for unpacking sensors
  int *hook_ss_id{ nullptr }, *sensor_type{ nullptr }, *smoothing_type{ nullptr };
  int *subdiv_type_ss{ nullptr }, *elt_type_ss{ nullptr };

  //std::cout << "adaptCells : before loop" << std::endl;

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
    err = __adapt_wrapper_lvl0<ISO>(*elt_type, *sensor_type, vString->c_str(), l, hmeshes, sensors);
  else if (*subdiv_type == NUGA::ISO_HEX)
    err = __adapt_wrapper_lvl0<ISO_HEX>(*elt_type, *sensor_type, vString->c_str(), l, hmeshes, sensors);
  else if (*subdiv_type == NUGA::DIR)
    err = __adapt_wrapper_lvl0<DIR>(*elt_type, *sensor_type, vString->c_str(), l, hmeshes, sensors);

  //std::cout << "adaptCells : end" << std::endl;

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
  
  std::vector<std::vector<E_Int>> bcs;

  mesh_type hmesh(crd, cnt, 1/*idx start*/, bcs);
  hmesh.init();
  sensor_type sensor(hmesh, NUGA::eSmoother(smoothing_type), 1/*max pts per cell*/, itermax);
  /*E_Int err = */sensor.assign_data(crdS);
  
  NUGA::adaptor<mesh_type, sensor_type>::run(hmesh, sensor);

  //std::cout << "output leaves..." << std::endl
  ngon_type ngo;
  std::vector<E_Int> pghids1, phhids1, hmpgid_to_confpgid;
  hmesh.conformize(ngo, hmpgid_to_confpgid, pghids1, phhids1);

  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);

  //std::cout << "output ..." << std::endl;
  
  PyObject* tpl = K_ARRAY::buildArray(hmesh._crd, varString, cnto, -1, "NGON", false);

  delete f; delete cn;
  return tpl;
}


//=======================  Intersector/PolyMeshTools/splitCells.cpp ====================
