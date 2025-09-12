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

# include "intersector.h"
# include <string>
# include <sstream> 

#include "Nuga/include/mesh_t.hxx"

//#include <iostream>
#include <memory>
#include "dico_to_stl.h"
#include "Nuga/include/close_cells.hxx"

#include <map>

# include "mpi.h"
# include "mpi4py/mpi4py.h"

using namespace std;
using namespace K_FLD;
using namespace NUGA;

//=============================================================================
//
//=============================================================================
PyObject* K_INTERSECTOR::closeCells_mpi(PyObject* self, PyObject* args)
{
  PyObject *py_arrs(nullptr);
  PyObject *py_zids(nullptr), *py_zid_to_rid_to_list(nullptr), *py_zonerank(nullptr);
  PyObject *py_rid_to_zones(nullptr), *mpi4pyCom(nullptr);
  if (!PYPARSETUPLE_(args, OOOO_ OO_, &py_arrs, &py_zids, &py_zid_to_rid_to_list, &py_zonerank, &py_rid_to_zones, &mpi4pyCom)) return NULL;

  void* pt_comm = (void*)&(((PyMPICommObject*)mpi4pyCom)->ob_mpi);
  MPI_Comm COM = *((MPI_Comm*) pt_comm);

  // 1. GET MESHES 

  E_Int nb_meshes{0};
  if (PyList_Check(py_arrs))
    nb_meshes = PyList_Size(py_arrs);

  //std::cout << "nb_meshes : " << nb_meshes << std::endl;

  if (nb_meshes == 0) return nullptr;

  std::vector<K_FLD::FloatArray*> f(nb_meshes);
  std::vector<K_FLD::IntArray*> cn(nb_meshes);
  char* varString, *eltType;

  std::vector<NUGA::ph_mesh_t> ph_meshes(nb_meshes);
  std::vector<NUGA::ph_mesh_t*> ptr_ph_meshes(nb_meshes);

  using ngon_type = ngon_t<K_FLD::IntArray>;
  
  for (E_Int m = 0; m < nb_meshes; ++m)
  {
    PyObject* arr = PyList_GetItem(py_arrs, m);
    // Check array
    E_Int err = check_is_NGON(arr, f[m], cn[m], varString, eltType);
    if (err) return nullptr;
    
    // conversion to the generic mesh interface
    ph_meshes[m].crd = *(f[m]);
    ph_meshes[m].cnt = std::move(ngon_type(*(cn[m])));

    ptr_ph_meshes[m] = &ph_meshes[m];
  }

  // 2. GET ZIDS 
  std::vector<int> zids(nb_meshes);
  assert (nb_meshes == PyList_Size(py_zids));
  for (E_Int m = 0; m < nb_meshes; ++m)
  {
    PyObject* pyz = PyList_GetItem(py_zids, m);
    E_Int zid = (int) PyInt_AsLong(pyz);
    zids[m]=zid;
  }

  // 3. GET RID TO LIST MAP 
  std::map<int, std::map<int, std::vector<E_Int>>> zid_to_rid_to_list;
  convert_dico_to_map___int_int_vecint(py_zid_to_rid_to_list, zid_to_rid_to_list);
  //assert (zid_to_rid_to_list.size() == nb_meshes);

  // 4. GET ZONE RANK
  std::vector<int> zonerank;
  if (PyDict_Check(py_zonerank))
  {
    E_Int nranks = PyDict_Size(py_zonerank);
    for (E_Int r = 0; r < nranks; ++r)
    {
      PyObject* key = PyInt_FromLong ((long) r);
      PyObject* py_rank = PyDict_GetItem(py_zonerank,key);
      assert (py_rank);
      E_Int rank = (int) PyInt_AsLong(py_rank);
      //std::cout << "key/item : " << r << "/" << my_val <<std::endl;//<< *item << std::endl;
      zonerank.push_back(rank);
    }
  }

  // 5. GET RID_TO_ZONES MAP 
  //todo VD : py_rid_to_zones => rid_to_zones
  std::map<int, std::pair<int,int>> rid_to_zones;
  convert_dico_to_map__int_pairint(py_rid_to_zones, rid_to_zones);

  // 6. CLOSE
  using para_algo_t = NUGA::hybrid_para_algo<NUGA::ph_mesh_t, E_Float>; // SEQ or multi-SEQ
  using closecell_t = NUGA::close_cells< para_algo_t, NUGA::ph_mesh_t>;

  closecell_t cc;
  cc.run(ptr_ph_meshes, zids, zid_to_rid_to_list, rid_to_zones, zonerank, COM);

  // pushing out the result : the set of closed meshes
  PyObject *l(PyList_New(0));
  for (size_t i=0; i < (size_t)nb_meshes; ++i)
  {
    K_FLD::IntArray cnto;
    ptr_ph_meshes[i]->cnt.export_to_array(cnto);

    // pushing out the mesh
    PyObject *tpl = K_ARRAY::buildArray(ptr_ph_meshes[i]->crd, varString, cnto, -1, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  for (int m = 0; m < nb_meshes; ++m)
  {
    delete f[m]; delete cn[m];
  }

  return l;
}



//=======================  Intersector/PolyMeshTools/utils_mpi.cpp ====================
