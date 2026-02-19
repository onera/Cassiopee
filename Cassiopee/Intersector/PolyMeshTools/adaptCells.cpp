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

using ngon_type = ngon_t<K_FLD::IntArray>;
using elt_t = K_INTERSECTOR::eType;
using subdiv_t = NUGA::eSUBDIV_TYPE;

//=============================================================================
/* Initialize the mesh (roerient & shift_geom) */
//=============================================================================
PyObject* K_INTERSECTOR::initForAdaptCells(PyObject* self, PyObject* args)
{
  PyObject *arr, *py_dict_transfo_to_list;

  if (!PYPARSETUPLE_(args, OO_, &arr, &py_dict_transfo_to_list)) return nullptr;

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
    delete f; delete cn;
    return nullptr;
  }


  // We reorient the PG of our NGON
  ngi.flag_externals(1);

  DELAUNAY::Triangulator dt;
  bool has_been_reversed;
  err = ngon_type::reorient_skins(dt, crd, ngi, has_been_reversed); //orientate normal outwards
  if (err) { delete f; delete cn; return nullptr; }

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
  delete f; delete cn;
  return m;
}

///
template <typename ELT_t, NUGA::eSUBDIV_TYPE STYPE>
void* __createHM(K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int zid)
{
  using hmesh_t = NUGA::hierarchical_mesh<ELT_t, STYPE>;
  hmesh_t* hm = new hmesh_t(crd, ngon_type(cnt));
  hm->zid = zid;
  return hm;
}

//============================================================================
/* Create a (list of) hmesh/hzone(s) and returns a (list of) hook(s) */
//============================================================================
PyObject* K_INTERSECTOR::createHMesh(PyObject* self, PyObject* args)
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

  if (!PYPARSETUPLE_(args, O_ II_, &arr, subtype, zid)) return nullptr;

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
  if (*subtype == NUGA::ISO || *subtype == NUGA::DIR_PROTO || *subtype == NUGA::DIR)
    *etyp = check_has_NGON_BASIC_ELEMENT(cnt);

  if (*subtype == NUGA::ISO)
  {
    if (*etyp  == elt_t::UNKN)
    {
      PyErr_WarnEx(PyExc_Warning,
        "createHMesh: input mesh to adapt must have basic elements and must be in NGON format : no adaptation for this mesh.", 1);
      return nullptr;
    }
    else if (*etyp == elt_t::HEXA)
      packet[1] = __createHM<K_MESH::Hexahedron, NUGA::ISO>(crd, cnt, *zid);
    else if (*etyp == elt_t::TETRA)
      packet[1] = __createHM<K_MESH::Tetrahedron, NUGA::ISO>(crd, cnt, *zid);
    else if (*etyp == elt_t::PRISM3)
      packet[1] = __createHM<K_MESH::Prism, NUGA::ISO>(crd, cnt, *zid);
    else if (*etyp == elt_t::BASIC)
      packet[1] = __createHM<K_MESH::Basic, NUGA::ISO>(crd, cnt, *zid);
  }
  else if (*subtype == NUGA::ISO_HEX)
  {
    if (*etyp != elt_t::UNKN)
    {
      PyErr_WarnEx(PyExc_Warning,
        "createHMesh: ISO_HEX policy is only supported Polyhedral mesh.", 1);
      return nullptr;
    }
    packet[1] = __createHM<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(crd, cnt, *zid);
  }
  else if (*subtype == NUGA::DIR_PROTO)
  {
    if (*etyp != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "createHMesh: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    packet[1] = __createHM<K_MESH::Hexahedron, NUGA::DIR_PROTO>(crd, cnt, *zid);
  }
  else if (*subtype == NUGA::DIR)
  {
    if (*etyp != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "createHMesh: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    packet[1] = __createHM<K_MESH::Hexahedron, NUGA::DIR>(crd, cnt, *zid);
  }

  if (packet[1] == nullptr) { Py_INCREF(Py_None); return Py_None; } // the input mesh does not have basic elts

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif
  
  return hook;
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

PyObject* K_INTERSECTOR::deleteHMesh(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PYPARSETUPLE_(args, O_, &hook))
  {
    return NULL;
  }

  // recupere le hook
  int* sub_type{ nullptr }, *elt_type{ nullptr }, *hook_id{ nullptr }, *zid(nullptr);
  std::string* vString{ nullptr };
  void** packet{ nullptr };
  void* hmesh_ptr = unpackHMesh(hook, hook_id, sub_type, elt_type, zid, vString, packet);

  if (*sub_type == NUGA::ISO)
  {
    if (*elt_type == elt_t::HEXA)
      __deleteHM<K_MESH::Hexahedron, NUGA::ISO>(hmesh_ptr);
    else if (*elt_type == (E_Int)elt_t::TETRA)
      __deleteHM<K_MESH::Tetrahedron, NUGA::ISO>(hmesh_ptr);
    else if (*elt_type == (E_Int)elt_t::PRISM3)
      __deleteHM<K_MESH::Prism, NUGA::ISO>(hmesh_ptr);
    else if (*elt_type == (E_Int)elt_t::BASIC)
      __deleteHM<K_MESH::Basic, NUGA::ISO>(hmesh_ptr);
  }
  else if (*sub_type == NUGA::ISO_HEX)
  {
    if (*elt_type != elt_t::UNKN)
    {
      PyErr_WarnEx(PyExc_Warning,
        "deleteHMesh: ISO_HEX policy is only supported Polyhedral mesh.", 1);
      return nullptr;
    }
    __deleteHM<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh_ptr);
  }
  else if (*sub_type == NUGA::DIR_PROTO)
  {
    if (*elt_type != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "deleteHMesh: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    __deleteHM<K_MESH::Hexahedron, NUGA::DIR_PROTO>(hmesh_ptr);
  }
  else if (*sub_type == NUGA::DIR)
  {
    if (*elt_type != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "deleteHMesh: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    __deleteHM<K_MESH::Hexahedron, NUGA::DIR>(hmesh_ptr);
  }

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
void __interpolateHMeshNodalField(const void* hmesh_ptrs, std::vector<double>& fieldN, std::vector<std::vector<E_Int>>& bcptlists)
{
  using mesh_type = NUGA::hierarchical_mesh<ELT_t, STYPE>;
  mesh_type* hmesh = (mesh_type*)hmesh_ptrs;
  
  if (bcptlists.empty()) return;
  
  // hmesh2amesh
  std::vector<E_Int> pth2a;
  K_CONNECT::IdTool::reverse_indirection(hmesh->pthids0, pth2a);  

// PARTIE A
  std::set<E_Int> donnorFaces;  
  for(size_t i = 0; i < bcptlists[0].size(); ++i)
  {
    E_Int PGi_hmesh = bcptlists[0][i] - 1;  
    PGi_hmesh = hmesh->pghids0[PGi_hmesh];


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

///
PyObject* K_INTERSECTOR::interpolateHMeshNodalField(PyObject* self, PyObject* args)
{
  PyObject *hook{nullptr}, *pyfieldN{nullptr}, *py_bcptlists{nullptr};

  if (!PYPARSETUPLE_(args, OOO_, &hook, &pyfieldN, &py_bcptlists))
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
      /*E_Int res2 = */K_NUMPY::getFromPointList(pyBCptList, ptL, size, nfld);
      //std::cout << "res2/size/nfld : " << res2 << "/" << size << "/" << nfld << std::endl;

      std::vector<E_Int> vPtL(ptL, ptL+size);
      bcptlists[i] = vPtL;
    }
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
  void* hmesh_ptr = unpackHMesh(hook, hook_id, sub_type, elt_type, zid, vString, packet);
  
  if (*sub_type == NUGA::ISO)
  {
    if (*elt_type == elt_t::HEXA)
      __interpolateHMeshNodalField<K_MESH::Hexahedron, NUGA::ISO>(hmesh_ptr, fieldN, bcptlists);
    else if (*elt_type == (E_Int)elt_t::TETRA)
      __interpolateHMeshNodalField<K_MESH::Tetrahedron, NUGA::ISO>(hmesh_ptr, fieldN, bcptlists);
    else if (*elt_type == (E_Int)elt_t::PRISM3)
      __interpolateHMeshNodalField<K_MESH::Prism, NUGA::ISO>(hmesh_ptr, fieldN, bcptlists);
    else if (*elt_type == (E_Int)elt_t::BASIC)
      __interpolateHMeshNodalField<K_MESH::Basic, NUGA::ISO>(hmesh_ptr, fieldN, bcptlists);
    }
  else if (*sub_type == NUGA::ISO_HEX)
  {
    if (*elt_type != elt_t::UNKN)
    {
      PyErr_WarnEx(PyExc_Warning,
        "interpolateHMeshNodalField: ISO_HEX policy is only supported Polyhedral mesh.", 1);
      return nullptr;
    }
    __interpolateHMeshNodalField<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh_ptr, fieldN, bcptlists);
  }
  else if (*sub_type == NUGA::DIR_PROTO)
  {
    if (*elt_type != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "interpolateHMeshNodalField: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    __interpolateHMeshNodalField<K_MESH::Hexahedron, NUGA::DIR_PROTO>(hmesh_ptr, fieldN, bcptlists);
  }
  else if (*sub_type == NUGA::DIR)
  {
    if (*elt_type != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "interpolateHMeshNodalField: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    __interpolateHMeshNodalField<K_MESH::Hexahedron, NUGA::DIR>(hmesh_ptr, fieldN, bcptlists);
  }


  //Retourner les champs mis à jour
  K_FLD::FloatArray farr(1,fieldN.size());
  
  std::vector<double>& fld = fieldN;
  for (size_t j = 0; j < fieldN.size(); ++j)farr(0, j) = fld[j];

  PyObject* tpl = K_ARRAY::buildArray(farr, varString, *c, -1, eltType, false);
  
  delete fi;
  delete c;

  return tpl;
}

//============================================================================
/* Creates a Sensor */
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
  if (sensor_type == 5)
  {
    using hmesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::DIR>;
    using sensor_t = NUGA::metric_sensor<hmesh_t>;

    hmesh_t* hm = (hmesh_t*)(hmesh);

    return new sensor_t(*hm);
  }
  return nullptr;
}

///
PyObject* K_INTERSECTOR::createSensor(PyObject* self, PyObject* args)
{
  PyObject *hook_sensor(nullptr);

  PyObject *hook_hmesh(nullptr);
  E_Int smoothing_type(0), sensor_type(0), itermax(0); 
  E_Int metric_policy(0); // O:ISO_MIN, 1:ISO_MEAN, 2:ISO_MAX, 3:ISO_MIN OR ISO_MAX 

  if (!PYPARSETUPLE_(args, O_ IIII_, &hook_hmesh, &sensor_type, &smoothing_type, &itermax, &metric_policy)) return NULL;


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
  else if (*sub_type == NUGA::DIR_PROTO)
  {
    /*if (sensor_type != 2 && sensor_type != 3)
    {
      PyErr_SetString(PyExc_ValueError,
       "adaptCells: DIR_PROTO only works with cell/nodal sensor currently.");
      return nullptr;
    }*/
  }
  else if (*sub_type == NUGA::DIR)
  {
    if (sensor_type != 5)
    {
      PyErr_SetString(PyExc_ValueError,
       "adaptCells: DIR only works with metric sensor.");
      return nullptr;
    }
  }

  // HMESH PTR
  packet_ss[2] = nullptr;
  if (*subtype_hm == NUGA::ISO)
  {
    if (*elt_type == elt_t::HEXA)
      packet_ss[2] = __createSensor<K_MESH::Hexahedron, NUGA::ISO>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
    else if (*elt_type == elt_t::TETRA)
      packet_ss[2] = __createSensor<K_MESH::Tetrahedron, NUGA::ISO>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
    else if (*elt_type == elt_t::PRISM3)
      packet_ss[2] = __createSensor<K_MESH::Prism, NUGA::ISO>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
    else if (*elt_type == elt_t::BASIC)
      packet_ss[2] = __createSensor<K_MESH::Basic, NUGA::ISO>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
    else 
    {
      PyErr_SetString(PyExc_ValueError, "createSensor: wrong element type in hmesh for ISO strategy.");
      return nullptr;
    }
  }
  else if (*subtype_hm == NUGA::ISO_HEX)
  {
    if (*elt_type != elt_t::UNKN)
    {
      PyErr_WarnEx(PyExc_Warning,
        "createSensor: ISO_HEX policy is only supported Polyhedral mesh.", 1);
      return nullptr;
    }
    packet_ss[2] = __createSensor<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
  }
  else if (*subtype_hm == NUGA::DIR_PROTO)
  {
    if (*elt_type != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "createSensor: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    packet_ss[2] = __createSensor<K_MESH::Hexahedron, NUGA::DIR_PROTO>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
  }
  else if (*subtype_hm == NUGA::DIR)
  {
    if (*elt_type != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "createSensor: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    packet_ss[2] = __createSensor<K_MESH::Hexahedron, NUGA::DIR>(hmesh, smoothing_type, itermax, metric_policy, sensor_type);
  }

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
PyObject* K_INTERSECTOR::deleteSensor(PyObject* self, PyObject* args)
{
  PyObject* hook_sensor;
  if (!PYPARSETUPLE_(args, O_, &hook_sensor))
  {
    return NULL;
  }

  // recupere le hook
  int *hook_ss_id{ nullptr }, *sensor_type{ nullptr }, *smoothing_type{ nullptr };
  int *subdiv_type{ nullptr }, *elt_type{ nullptr };
  void** packet{ nullptr };
  void* sensor = unpackSensor(hook_sensor, hook_ss_id, sensor_type, smoothing_type, subdiv_type, elt_type, packet);  

  if (*subdiv_type == NUGA::ISO)
  {
    if (*elt_type == elt_t::HEXA)
      __deleteSensor<K_MESH::Hexahedron, NUGA::ISO>(*sensor_type, sensor);
    else if (*elt_type == (E_Int)elt_t::TETRA)
      __deleteSensor<K_MESH::Tetrahedron, NUGA::ISO>(*sensor_type, sensor);
    else if (*elt_type == (E_Int)elt_t::PRISM3)
      __deleteSensor<K_MESH::Prism, NUGA::ISO>(*sensor_type, sensor);
    else if (*elt_type == (E_Int)elt_t::BASIC)
      __deleteSensor<K_MESH::Basic, NUGA::ISO>(*sensor_type, sensor);
  }
  else if (*subdiv_type == NUGA::ISO_HEX)
  {
    if (*elt_type != elt_t::UNKN)
    {
      PyErr_WarnEx(PyExc_Warning,
        "deleteSensor: ISO_HEX policy is only supported Polyhedral mesh.", 1);
      return nullptr;
    }
    __deleteSensor<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(*sensor_type, sensor);
  }
  else if (*subdiv_type == NUGA::DIR_PROTO)
  {
    if (*elt_type != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "deleteSensor: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    __deleteSensor<K_MESH::Hexahedron, NUGA::DIR_PROTO>(*sensor_type, sensor);
  }
  else if (*subdiv_type == NUGA::DIR)
  {
    if (*elt_type != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "deleteSensor: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    __deleteSensor<K_MESH::Hexahedron, NUGA::DIR>(*sensor_type, sensor);
  }

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

template <NUGA::eSUBDIV_TYPE STYP> using incr_t = typename NUGA::adap_incr_type<STYP>::cell_incr_t; // E_Int (ISO) or int_tuple<3> (DIR_PROTO)

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
      mesh.get_boundary(data);
      sensor->assign_data(data);
    }
    else
    {
      NUGA::pg_smesh_t data(crdS, cntS);
      sensor->assign_data(data);
    }
  }
  else if (sensor_type == 5)
  {
    assert(crdS.rows() == 6);

    using mesh_type  = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::DIR>;
    using sensor_t   = NUGA::metric_sensor<mesh_type>;
    sensor_t* sensor = (sensor_t*)psensor;

    sensor->assign_data(crdS);
  }
}


PyObject* K_INTERSECTOR::assignData2Sensor(PyObject* self, PyObject* args)
{

  PyObject* hook_sensor{ nullptr }, *dataSensor(nullptr);
  if (!PYPARSETUPLE_(args, OO_, &hook_sensor, &dataSensor))
    return nullptr;
  if (dataSensor == Py_None) // nothing to do
    { Py_INCREF(Py_None); return Py_None; }

  // recupere le hook
  int *hook_ss_id{ nullptr }, *sensor_type{ nullptr }, *smoothing_type{ nullptr };
  int *subdiv_type{ nullptr }, *elt_type{ nullptr };
  void** packet{ nullptr };
  void* sensor = unpackSensor(hook_sensor, hook_ss_id, sensor_type, smoothing_type, subdiv_type, elt_type, packet);

  std::vector<E_Int> sens_data;
  K_FLD::FloatArray fS;
  K_FLD::IntArray cnS;

  //geom/xsensor or nodal_sensor data ?

  if (PyList_Check(dataSensor)) // Array (mesh or coordinates) OR metric_field as 6 numpies
  {
    if (PyList_Size(dataSensor) == 6) // metric_field as 6 numpies
    {
      E_Float *mxx, *mxy, *mxz, *myy, *myz, *mzz;
      E_Int size = -1, nfld = -1;

      PyObject *py_mxx = PyList_GetItem(dataSensor, 0);
      K_NUMPY::getFromNumpyArray(py_mxx, mxx, size, nfld);
      
      PyObject *py_mxy = PyList_GetItem(dataSensor, 1);
      K_NUMPY::getFromNumpyArray(py_mxy, mxy, size, nfld);
      
      PyObject *py_mxz = PyList_GetItem(dataSensor, 2);
      K_NUMPY::getFromNumpyArray(py_mxz, mxz, size, nfld);
      
      PyObject *py_myy = PyList_GetItem(dataSensor, 3);
      K_NUMPY::getFromNumpyArray(py_myy, myy, size, nfld);
      
      PyObject *py_myz = PyList_GetItem(dataSensor, 4);
      K_NUMPY::getFromNumpyArray(py_myz, myz, size, nfld);
      
      PyObject *py_mzz = PyList_GetItem(dataSensor, 5);
      K_NUMPY::getFromNumpyArray(py_mzz, mzz, size, nfld);

      // remplir fS
      fS.resize(6, nfld);
      
      for (E_Int i = 0; i < nfld; i++) {
        fS(0,i) = mxx[i];
        fS(1,i) = mxy[i];
        fS(2,i) = mxz[i];
        fS(3,i) = myy[i];
        fS(4,i) = myz[i];
        fS(5,i) = mzz[i]; 
      }
    }
    else { // Array (mesh or coordinates)
     E_Int ni, nj, nk;
    char* varString, *eltType;
    /*E_Int res = */K_ARRAY::getFromArray(dataSensor, varString, fS, ni, nj, nk, cnS, eltType);
    //std::cout << "res/eltType/fs sz/cns sz : " << res << "/" << eltType << "/" << fS.cols() << "/" << cnS.cols() << std::endl;
    }
  }
  else // assuming numpy for nodal/cell
  {
     E_Int *nodv, size, nfld;
     /*E_Int res2 = */K_NUMPY::getFromNumpyArray(dataSensor, nodv, size, nfld);
     sens_data.resize(size);  
     for (E_Int i=0; i < size; ++i)
     {
       E_Int v = nodv[i];
       v = std::max((E_Int)0,v); //disable agglo currently
       sens_data[i] = v;
     }
  }

  if (*subdiv_type == NUGA::ISO)
  {
    if (*elt_type == elt_t::HEXA)
    __assign_sensor_data<K_MESH::Hexahedron, NUGA::ISO>(*sensor_type, sensor, fS, cnS, sens_data);
  else if (*elt_type == (E_Int)elt_t::TETRA)
    __assign_sensor_data<K_MESH::Tetrahedron, NUGA::ISO>(*sensor_type, sensor, fS, cnS, sens_data);
  else if (*elt_type == (E_Int)elt_t::PRISM3)
    __assign_sensor_data<K_MESH::Prism, NUGA::ISO>(*sensor_type, sensor, fS, cnS, sens_data);
  else if (*elt_type == (E_Int)elt_t::BASIC)
    __assign_sensor_data<K_MESH::Basic, NUGA::ISO>(*sensor_type, sensor, fS, cnS, sens_data);
  }
  else if (*subdiv_type == NUGA::ISO_HEX)
    __assign_sensor_data<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(*sensor_type, sensor, fS, cnS, sens_data);
  else if (*subdiv_type == NUGA::DIR_PROTO)
  {
    std::vector<int_tuple<3>> sens_data_dir(sens_data.size(), int_tuple<3>(0));
    for (size_t u=0; u < sens_data.size(); ++u)
      sens_data_dir[u] = sens_data[u]; // pass a single value to the 3-tuple

    __assign_sensor_data<K_MESH::Hexahedron, NUGA::DIR_PROTO>(*sensor_type, sensor, fS, cnS, sens_data_dir);
  }
  else if (*subdiv_type == NUGA::DIR)
  {
    // Imad
    std::vector<int_tuple<3>> dummy;
    __assign_sensor_data<K_MESH::Hexahedron, NUGA::DIR>(*sensor_type, sensor, fS, cnS, dummy);
  }

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* generic adatptation function : SEQ or MULTI-SEQ*/
//=============================================================================
template <typename hmesh_t, typename sensor_t>
E_Int __adapt_lvl0
(std::vector<void*>& hookhmeshes, std::vector<void*>& hooksensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_owned,
const char* varString, PyObject *out)
{
  if (hookhmeshes.empty()) return 0;
  if (hooksensors.empty()) return 0;

  size_t nb_meshes = hookhmeshes.size();
  if (nb_meshes != hooksensors.size()) return 1;

  using para_algo_t = omp_algo<hmesh_t, E_Int>; // SEQ or multi-SEQ

  using adaptor_t = NUGA::adaptor_para<para_algo_t, hmesh_t, sensor_t>;

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

  a.run(hmeshes, zids, zone_to_rid_to_list_owned, rid_to_zones);

  for (size_t i=0; i < nb_meshes; ++i)
  {
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
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_owned,
const char* varString, PyObject *out)
{
  int err(0);

  if (elt_type==elt_t::HEXA)
  {
    using mesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, STYPE>;
   
    if (sensor_type == 0) // geom sensor
    {
      using sensor_t = NUGA::geom_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 1 || sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
  }
  else if (elt_type==elt_t::TETRA)
  {
    using mesh_t = NUGA::hierarchical_mesh<K_MESH::Tetrahedron, STYPE>;

    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
  }
  else if (elt_type==elt_t::PRISM3)
  {
    using mesh_t = NUGA::hierarchical_mesh<K_MESH::Prism, STYPE>;
    
    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
  }
  else if (elt_type==elt_t::BASIC)
  {
    using mesh_t = NUGA::hierarchical_mesh<K_MESH::Basic, STYPE>;

    if (sensor_type == 1) sensor_type = 0; //currently xsensor not supported

    if (sensor_type == 0)
    {
      using sensor_t = NUGA::geom_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 2) //nodal sensor
    {
      using sensor_t = NUGA::nodal_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 3) //cell sensor
    {
      using sensor_t = NUGA::cell_sensor<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
    else if (sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }
  }
  return err;
}

template <>
int __adapt_wrapper<NUGA::ISO_HEX>
(int elt_type/*dummy*/, int sensor_type,
std::vector<void*>& hmeshes, std::vector<void*>& sensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_owned,
const char* varString, PyObject *out)
{
  E_Int err(0);
  
  using mesh_t = NUGA::hierarchical_mesh<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>;
  
  if (sensor_type == 0) // geom sensor
  {
    using sensor_t = NUGA::geom_sensor<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
  }
  else if (sensor_type == 2) //nodal sensor
  {
    using sensor_t = NUGA::nodal_sensor<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
  }
  else if (sensor_type == 3) //cell sensor
  {
    using sensor_t = NUGA::cell_sensor<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
  }
  else if (sensor_type == 1 || sensor_type == 4) // xsensor2
    {
      using sensor_t = NUGA::xsensor2<mesh_t>;
      err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
    }

  return err;
}

template <>
int __adapt_wrapper<NUGA::DIR_PROTO>
(int elt_type/*dummy*/, int sensor_type,
std::vector<void*>& hmeshes, std::vector<void*>& sensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_owned,
const char* varString, PyObject *out)
{
  E_Int err(0);

  using mesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::DIR_PROTO>;

  if (sensor_type == 0) // geom sensor
  {
    using sensor_t = NUGA::geom_sensor<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
  }
  else if (sensor_type == 2) //nodal sensor
  {
    using sensor_t = NUGA::nodal_sensor<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
  }
  else if (sensor_type == 3) //cell sensor
  {
    using sensor_t = NUGA::cell_sensor<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
  }
  else if (sensor_type == 1 || sensor_type == 4) // xsensor2
  {
    using sensor_t = NUGA::xsensor2<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
  }

  return err;
}

template <>
int __adapt_wrapper<NUGA::DIR>
(int elt_type/*dummy*/, int sensor_type,
std::vector<void*>& hmeshes, std::vector<void*>& sensors,
std::map<int, std::pair<int,int>>& rid_to_zones,
std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_owned,
const char* varString, PyObject *out)
{
  E_Int err(0);

  using mesh_t = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::DIR>;

  if (sensor_type == 5) // metric sensor
  {
    using sensor_t = NUGA::metric_sensor<mesh_t>;
    err = __adapt_lvl0<mesh_t, sensor_t>(hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, varString, out);
  }

  return err;
}

//=============================================================================
/* Hierarchical Mesh Adaptation */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCells(PyObject* self, PyObject* args)
{
  PyObject *hook_hmeshes(nullptr), *hook_sensors(nullptr), *py_zone_to_rid_to_list_owned(nullptr);
  PyObject *py_rid_to_zones(nullptr);
  if (!PYPARSETUPLE_(args, OOOO_, &hook_hmeshes, &hook_sensors, &py_zone_to_rid_to_list_owned, &py_rid_to_zones)) return NULL;
  
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


  // 3. GET zone_to_rid_to_list_owned
  std::map<int, std::map<int, std::vector<E_Int>>> zone_to_rid_to_list_owned;
  convert_dico_to_map___int_int_vecint(py_zone_to_rid_to_list_owned, zone_to_rid_to_list_owned);
  //assert (zone_to_rid_to_list_owned == nb_meshes);

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
  // assert (zone_to_rid_to_list_owned == nb_meshes);

  // Adaptation
  // ==========
  PyObject *l(PyList_New(0));
  std::vector<E_Int> dummy;

  E_Int err(0);
  if (*subdiv_type == NUGA::ISO)
    err = __adapt_wrapper<ISO>(*elt_type, *sensor_type, hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, vString->c_str(), l);
  else if (*subdiv_type == NUGA::ISO_HEX)
    err = __adapt_wrapper<ISO_HEX>(*elt_type, *sensor_type, hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, vString->c_str(), l);
  else if (*subdiv_type == NUGA::DIR_PROTO)
    err = __adapt_wrapper<DIR_PROTO>(*elt_type, *sensor_type, hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, vString->c_str(), l);
  else if (*subdiv_type == NUGA::DIR)
    err = __adapt_wrapper<DIR>(*elt_type, *sensor_type, hmeshes, sensors, rid_to_zones, zone_to_rid_to_list_owned, vString->c_str(), l);

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

  if (!PYPARSETUPLE_(args, O_ R_ II_, &arrS, &bratio, &smoothing_type, &itermax)) return NULL;

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

  //std::cout << "adapt box..." << std::endl;

  using mesh_type = NUGA::hierarchical_mesh<K_MESH::Hexahedron, NUGA::ISO>;
  using sensor_type = NUGA::geom_sensor/*geom_static_sensor*/<mesh_type>;

  // We reorient the PG of our NGON
  ngi.flag_externals(1);

  {
    DELAUNAY::Triangulator dt;
    bool has_been_reversed;
    int err = ngon_type::reorient_skins(dt, crd, ngi, has_been_reversed);
    if (err)
    {
      delete f; delete cn;
      return nullptr;
    }
  }
  
  mesh_type hmesh(crd, ngi);

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

void convert_ptl_to_hmesh_ids(const std::vector<E_Int>& indir, std::vector<E_Int>& ptl)
{
  for (size_t i=0; i < ptl.size(); ++i)
    ptl[i]=indir[ptl[i]-1]+1;
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
                    std::vector<std::vector<E_Float>>& fieldsF,
					const bool conformize)
{
  using mesh_type = NUGA::hierarchical_mesh<ELT_t, STYPE>;

  mesh_type* hmesh = (mesh_type*)hmesh_ptr;

  if (hmesh->pghids0.empty()) //the 3 are empty/full at the same time
  {
    K_CONNECT::IdTool::init_inc(hmesh->pthids0, hmesh->_nb_pts0);
    K_CONNECT::IdTool::init_inc(hmesh->pghids0, hmesh->_nb_pgs0);
    K_CONNECT::IdTool::init_inc(hmesh->phhids0, hmesh->_nb_phs0);
  }

  // update the bc pointlists in hmesh
  for (size_t b=0; b < bcptlists.size(); ++b)
  {
    convert_ptl_to_hmesh_ids(hmesh->pghids0, bcptlists[b]);
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

      convert_ptl_to_hmesh_ids(hmesh->pghids0, ptlL);
      hmesh->update_pointlist(ptlL, false);

      convert_ptl_to_hmesh_ids(hmesh->pghids0, ptlR);
      hmesh->update_pointlist(ptlR, true);

      r_to_ptl.second = ptlL;
      r_to_ptl.second.insert(r_to_ptl.second.end(), ALL(ptlR));
    }
    else
    {
      convert_ptl_to_hmesh_ids(hmesh->pghids0, r_to_ptl.second);
      hmesh->update_pointlist(r_to_ptl.second, reverse);
    }
  }

  ngon_type ngo;
  std::vector<E_Int> pghids1, phhids1, hmpgid_to_confpgid;
  if (conformize)
  	hmesh->conformize(ngo, hmpgid_to_confpgid, pghids1, phhids1);
  else
  	hmesh->no_conformize(ngo, hmpgid_to_confpgid, pghids1, phhids1);


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

///
PyObject* K_INTERSECTOR::conformizeHMesh(PyObject* self, PyObject* args)
{ 
  //hooks[i], bcptlists, jzone_to_ptlist, fieldsC, fieldsN, fieldsF
  PyObject* hook, *py_bcptlists, *py_rid_to_ptlist, *py_rid_to_zones;
  PyObject* pyfieldsC, *pyfieldsN, *pyfieldsF;
  E_Int conformize;
  
  if (!PYPARSETUPLE_(args, OOOO_ OOO_ I_, &hook, &py_bcptlists, &py_rid_to_ptlist, &py_rid_to_zones, 
    &pyfieldsC, &pyfieldsN, &pyfieldsF, &conformize)) return NULL;

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
      /*E_Int res2 = */K_NUMPY::getFromPointList(pyBCptList, ptL, size, nfld);
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

      npy_intp* dims = PyArray_SHAPE(pyarr);

      E_Int ptl_sz = dims[0];
      
      //long* dataPtr = static_cast<long*>(PyArray_DATA(pyarr));
      E_Int* dataPtr = (E_Int*)PyArray_DATA(pyarr);

      std::vector<E_Int> ptl(ptl_sz);
      for (size_t u=0; u < (size_t)ptl_sz; ++u) ptl[u] = dataPtr[u];

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
      K_NUMPY::getFromNumpyArray(fieldFi, Fid);
      pFid[j] = Fid->begin(1);
      fieldsF[j].resize(Fid->getSize());
    }

    for (E_Int j = 0; j < nfieldsF; ++j)
      for (size_t i = 0; i < fieldsF[j].size(); ++i)
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
  {
    //std::cout << "conformize: " << conformize << std::endl;
    if (*elt_type == elt_t::HEXA)
      __conformizeHM<K_MESH::Hexahedron, NUGA::ISO>(hmesh, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF, conformize);
    else if (*elt_type == (E_Int)elt_t::TETRA)
      __conformizeHM<K_MESH::Tetrahedron, NUGA::ISO>(hmesh, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF, conformize);
    else if (*elt_type == (E_Int)elt_t::PRISM3)
      __conformizeHM<K_MESH::Prism, NUGA::ISO>(hmesh, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF, conformize);
    else if (*elt_type == (E_Int)elt_t::BASIC)
      __conformizeHM<K_MESH::Basic, NUGA::ISO>(hmesh, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF, conformize);
  }
  else if (*sub_type == NUGA::ISO_HEX)
  {
    if (*elt_type != elt_t::UNKN)
    {
      PyErr_WarnEx(PyExc_Warning,
        "conformizeHMesh: ISO_HEX policy is only supported Polyhedral mesh.", 1);
      return nullptr;
    }
    __conformizeHM<K_MESH::Polyhedron<0>, NUGA::ISO_HEX>(hmesh, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF, conformize);
  }
  else if (*sub_type == NUGA::DIR_PROTO)
  {
    if (*elt_type != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "conformizeHMesh: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    __conformizeHM<K_MESH::Hexahedron, NUGA::DIR_PROTO>(hmesh, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF, conformize);
  }
  else if (*sub_type == NUGA::DIR)
  {
    if (*elt_type != elt_t::HEXA)
    {
      PyErr_WarnEx(PyExc_Warning,
        "conformizeHMesh: directionnal policy is only supported with Hexahedral mesh currently.", 1);
      return nullptr;
    }
    __conformizeHM<K_MESH::Hexahedron, NUGA::DIR>(hmesh, crdo, cnto, rid_to_ptlist, rid_to_zones, bcptlists, fieldsC, fieldsN, fieldsF, conformize);
  }


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
      farr.resize(fieldsC.size(), fieldsC[0].size());
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
      farr.resize(fieldsN.size(), fieldsN[0].size());
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

static void transpose_pointlists
(
  const std::map<int, std::pair<int,int>>& rid_to_zones,
  const std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_owned,
  std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_opp
)
{
  zone_to_rid_to_list_opp.clear();

  // OMP : just transpose
  if (!zone_to_rid_to_list_owned.empty())
  {
    for (const auto& it : zone_to_rid_to_list_owned)
    {
      int zid = it.first;
      const auto& rid_to_list = it.second;
      for (const auto & z2L : rid_to_list)
      {
        int rid = z2L.first;
        auto it = rid_to_zones.find(rid);

        int jzid = (it->second.first == zid) ? it->second.second : it->second.first;
        const auto & ptL = z2L.second;

        zone_to_rid_to_list_opp[jzid][rid] = ptL;
      }
    }
  }

  assert(zone_to_rid_to_list_owned.size() == zone_to_rid_to_list_opp.size());
}

//=============================================================================
/* Transpose the owned PointLists : to update PointListDonor */
//=============================================================================
PyObject* K_INTERSECTOR::transposePointLists(PyObject* self, PyObject* args)
{
  // zonerank, Cmpi.rank, Cmpi.size, zone_to_rid_to_list_owned
  PyObject *py_rid_to_zones(nullptr), *py_zone_to_rid_to_list_owned(nullptr);

  if (!PYPARSETUPLE_(args, OO_, &py_rid_to_zones, &py_zone_to_rid_to_list_owned)) return nullptr;

  // 2. GET POINTLISTS MAP 
  std::map<int, std::map<int, std::vector<E_Int>>> zone_to_rid_to_list_owned;
  convert_dico_to_map___int_int_vecint(py_zone_to_rid_to_list_owned, zone_to_rid_to_list_owned);
  //assert (zone_to_rid_to_list_owned.size() == nb_meshes);

  // 3. GET RID_TO_ZONES MAP 
  //todo VD : py_rid_to_zones => rid_to_zones
  std::map<int, std::pair<int,int>> rid_to_zones;
  convert_dico_to_map__int_pairint(py_rid_to_zones, rid_to_zones);

  // 3. EXCHANGE : JUST TRANSPOSE (because OMP context)
  std::map<int, std::map<int, std::vector<E_Int>>> zone_to_rid_to_list_opp;
  transpose_pointlists(rid_to_zones, zone_to_rid_to_list_owned, zone_to_rid_to_list_opp);

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


//=======================  Intersector/PolyMeshTools/adaptCells.cpp ====================
