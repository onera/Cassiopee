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

# include "converter.h"
# include <string>
# include <sstream> 

//#define DEBUG_GHOST
//#define DEBUG_ZONE_T

# include "Fld/ngon_t.hxx"
# include "Nuga/include/zone_t.hxx"

//#include <iostream>
#include <memory>


using namespace K_FLD;

using ngon_type = ngon_t<K_FLD::IntArray>;
using crd_t = K_FLD::FloatArray;
using zone_type = NUGA::zone_t<crd_t, ngon_type>;

E_Int check_is_NGON(PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType)
{
  E_Int ni, nj, nk;
  
  E_Int res = K_ARRAY::getFromArray(arr, varString, f1, ni, nj, nk,
                                    cn1, eltType);
     
  bool err = (res !=2);
  err |= (strcmp(eltType, "NGON") != 0);
  if (err)
  {
    //std::cout << "input error : err => " << err << std::endl;
    //std::cout << "input error : eltType => " << eltType << std::endl;
    PyErr_SetString(PyExc_TypeError, "input error : invalid array, must be a unstructured NGON array.");//fixme triangulateExteriorFaces : PASS A STRING AS INPUT
    delete f1; delete cn1;
    return 1;
  }

  // Check coordinates.
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if ((posx == -1) || (posy == -1) || (posz == -1))
  {
    PyErr_SetString(PyExc_TypeError, "input error : can't find coordinates in array.");//fixme  conformUnstr
    delete f1; delete cn1;
    return 1;
  }
  
  return 0;
}

void add_n_topo_layers(std::vector<zone_type>& zones, E_Int i0, E_Int NLAYERS)
{
  E_Int nb_zones = zones.size();
  zone_type& Zi0 = zones[i0];

  // Coloring all zones polygons with a unique color keeping the priority INNER < JOIN < BC
  E_Int color(0);
  
  Zi0.init_pg_types(color);
    
  for (E_Int i=0; i < nb_zones; ++i)
  {
    if (i == i0) continue;
    zones[i].init_pg_types(color);
  }
  
  // Set PG types for current zone
  Zi0.set_BC_types(color); // to put current zone BC at the tail
      
  // Reset PHs types to 0 for all zones (incl. current) 
  for (E_Int j=0; j < nb_zones; ++j) zones[j].reset_ph_types(UNSET_COL);
  Zi0.reset_ph_types(PH_INNER_COL);
      
  // Flag first layer : traverse the joins
  for (auto& itrz : Zi0._joins)
  {
    E_Int jid = itrz.first;
    zone_type* rz = itrz.second.first;
    Vector_t<E_Int>& f2e = rz->get_f2e();
     
    const Vector_t<E_Int>& racids = rz->get_join_pg_list(jid);  
    E_Int siz = racids.size(); 
           
    for (E_Int i = 0; i < siz; ++i)
    {
      E_Int PGi = racids[i] - 1;
      E_Int leftPH = f2e[PGi] - 1;
        
      rz->set_ph_type(leftPH, 1);
    }
  }

  // Build layers from second to last
  Vector_t<E_Int> zones4next;
  std::set<zone_type*> zones4this;
  if (NLAYERS > 1)
  {
    zones4next.resize(nb_zones, false);
    
    for (auto& itrz : Zi0._joins)
      zones4this.insert(itrz.second.first) ;
  }

  //
  for (E_Int i = 1; i < NLAYERS; ++i)
  {
    zones4next.clear();
    zones4next.resize(nb_zones, false);
    
    for (auto z = zones4this.begin(); z != zones4this.end(); ++z)        
      (*z)->flag_connected_to_type(i+1, i, zones4next);
    
    if (i == NLAYERS - 1) break;
    
    zones4this.clear();
    for (E_Int j=0; j < nb_zones; ++j) if (zones4next[j]) zones4this.insert(&zones[j]);
   }

   // Agglomeration to Zi
  
  for (E_Int i = 0; i < nb_zones; ++i)
  {
    if (i == i0) continue;
    
    E_Int nbe = zones[i].reduce_to_positive_types();
    if (nbe == 0) continue;

#ifdef DEBUG_ZONE_T
    //zones[i].draw_boundaries();
    // K_FLD::IntArray cnto;
    // zones[i]._ng.export_to_array(cnto);
    // std::ostringstream o;
    // o << "p_" << i << ".plt";
    // MIO::write(o.str().c_str(), zones[i]._crd, cnto, "NGON");
    // Zi0.draw_boundaries();
#endif

    Zi0.append(zones[i], true/*keep joins*/);

#ifdef DEBUG_ZONE_T
    //Zi0.draw_boundaries();
#endif
  }

  Zi0.set_pg_colors();
  Zi0.sort_by_type();

}

//=============================================================================
/* Add ghost cells to a polyhedral multibloc mesh*/
//=============================================================================
PyObject* K_CONVERTER::addGhostCellsNG(PyObject* self, PyObject* args)
{
  PyObject *arrs, *if2es, *z_j_ptlist, *z_j_ptlistD, *z_j_ptlist_sizes, *z_j_donnorIDs, *z_bc_ptlist_sizes, *z_bc_ptlist;
  E_Int NLAYERS(1);
    
  if (!PyArg_ParseTuple(args, "OOOOOOOOi", &arrs, &if2es, &z_j_ptlist_sizes, &z_j_ptlist, &z_j_ptlistD, &z_j_donnorIDs, &z_bc_ptlist_sizes, &z_bc_ptlist, &NLAYERS)) return NULL;

  E_Int nb_zones = PyList_Size(arrs);
  //std::cout << "K_CONVERTER::addGhostCellsNG : " << nb_zones << std::endl;

  /*std::cout << "nb of arrs : " << PyList_Size(arrs) << std::endl;
  std::cout << "nb of if2es : " << PyList_Size(if2es) << std::endl;
  std::cout << "nb of z_j_ptlist_sizes : " << PyList_Size(z_j_ptlist_sizes) << std::endl;
  std::cout << "nb of z_j_ptlist : " << PyList_Size(z_j_ptlist) << std::endl;
  std::cout << "nb of z_j_ptlistD : " << PyList_Size(z_j_ptlistD) << std::endl;
  std::cout << "nb of z_j_donnorIDs : " << PyList_Size(z_j_donnorIDs) << std::endl;
  std::cout << "nb of z_bc_ptlist_sizes : " << PyList_Size(z_bc_ptlist_sizes) << std::endl;
  std::cout << "nb of z_bc_ptlist : " << PyList_Size(z_bc_ptlist) << std::endl;
  std::cout << "nb of layers : " << NLAYERS << std::endl;*/

  std::vector<zone_type*> zones(nb_zones, nullptr);
  std::vector<FloatArray*> crds(nb_zones, nullptr);
  std::vector<IntArray*> cnts(nb_zones, nullptr);
  std::vector<ngon_type*> ngs(nb_zones, nullptr);
  std::vector<E_Int*> f2es(nb_zones, nullptr);
  std::vector<FldArrayI> F2Es(nb_zones);

  std::vector<E_Int*> j_pointLists, j_pointListsD;
  char* varString, *eltType;
  
  E_Int err(0);

  //std::cout << "gathering inputs... " << std::endl;

  //E_Int nb_f2E = PyList_Size(if2es);
    
  E_Int ok = !err;

  for (E_Int i = 0; (i < nb_zones) && ok; ++i)
  {
  	//std::cout << "gathering inputs : check_is_NGON for " << i << std::endl;

  	PyObject* o = PyList_GetItem(arrs, i);

    err = check_is_NGON(o, crds[i], cnts[i], varString, eltType);
    if (err) break;


    //std::cout << "gathering inputs : creating ngon_type " << i << std::endl;

    ngs[i] = new ngon_type(*cnts[i]);

    // Get or build F2E

    PyObject* oo = PyList_GetItem(if2es, i);

    E_Int c, r;
    bool buildit = true;
      
    if (oo != Py_None){
      ok = K_NUMPY::getFromNumpyArray(oo, f2es[i], c, r, true/*shared*/);
      buildit = !ok;
    }

    //std::cout << f2e << std::endl;
    if (buildit) //error or missing : build it 
    {   
      //std::cout << "building F2E..." << std::endl;
      ngs[i]->build_noF2E(F2Es[i]);
      f2es[i] = F2Es[i].begin();
    }    

    zones[i] = new zone_type(i, *crds[i], *ngs[i], f2es[i], 0/*F2E_NONE*/);
  }

  E_Int jcount(0), bccount(0);

  // JOINS & BCS

  std::set<std::pair<E_Int, E_Int> > jdone;
  std::vector<std::pair<E_Int, E_Int> > j_to_pair; // donnor association

  for (E_Int i = 0; (i < nb_zones) && ok; ++i)
  {
  	  	
    PyObject* pyo_j_ptLs = PyList_GetItem(z_j_ptlist, i);
    PyObject* pyo_j_ptLs_sz = PyList_GetItem(z_j_ptlist_sizes, i);
    
    PyObject* pyo_j_ptLs_D = PyList_GetItem(z_j_ptlistD, i);
    PyObject* pyo_j_donIds = PyList_GetItem(z_j_donnorIDs, i);

    PyObject* pyo_bc_ptLs = PyList_GetItem(z_bc_ptlist, i);
    //PyObject* pyo_bc_ptLs_sz = PyList_GetItem(z_bc_ptlist_sizes, i);

    E_Int r, c, *donIds;

    ok = K_NUMPY::getFromNumpyArray(pyo_j_donIds, donIds, c, r, true/*shared*/);
    if (!ok)
    {
    	std::cout << "ERROR : could not get donnor ids" << std::endl; 
    	break;
    }
     
    E_Int *ptL_sz;
    ok = K_NUMPY::getFromNumpyArray(pyo_j_ptLs_sz, ptL_sz, c, r, true/*shared*/);
    if (!ok) 
    {
    	std::cout << "ERROR : could not get point list sizes" << std::endl; 
    	break;
    }

    // JOINS

    E_Int nb_joins = PyList_Size(pyo_j_ptLs);

    //std::cout << "nb of joins in c side  (before): " << nb_joins << std::endl;

    j_pointLists.reserve(j_pointLists.size() + nb_joins);
    j_pointListsD.reserve(j_pointListsD.size() + nb_joins);


    for (E_Int j=0; j < nb_joins; ++j)
    {

      std::pair<E_Int, E_Int> p1(i, donIds[j]), op(donIds[j], i); 

      if (jdone.find(op) != jdone.end()) continue;
      jdone.insert(p1);

      E_Int *ptL(nullptr), *ptL_D(nullptr);

      {
        PyObject* pyo_j_ptL = PyList_GetItem(pyo_j_ptLs, j);	
        E_Int c,r;
        ok = K_NUMPY::getFromNumpyArray(pyo_j_ptL, ptL, c, r, true/*shared*/);
        
        if (!ok) 
        {
          std::cout << "ERROR : could not get current point list" << std::endl; 
          break;
        }
        
        j_pointLists.push_back(ptL);
      }

      {
      	
        PyObject* pyo_ptL_D = PyList_GetItem(pyo_j_ptLs_D, j);	
        
        E_Int c,r;
        ok = K_NUMPY::getFromNumpyArray(pyo_ptL_D, ptL_D, c, r, true/*shared*/);
        
        if (!ok)
        {
          std::cout << "ERROR : could not get current point list Donnor" << std::endl; 
          break;
        }
        
        j_pointListsD.push_back(ptL_D);
      }
      
      //std::cout << "nb racs : " << ptL_sz[j] << std::endl;
      zone_type::join_zones(*zones[i], ptL, *zones[donIds[j]], ptL_D, ptL_sz[j], jcount++);

      j_to_pair.push_back(std::make_pair(i, donIds[j]));

    }

    // BCS

    E_Int nb_bcs = PyList_Size(pyo_bc_ptLs);
    for (E_Int b=0; b < nb_bcs; ++b)
    {
      PyObject* pyo_bc_ptL = PyList_GetItem(pyo_bc_ptLs, b);  

      E_Int nb_ids,r;
      E_Int *ptL(nullptr);    
      ok = K_NUMPY::getFromNumpyArray(pyo_bc_ptL, ptL, nb_ids, r, true/*shared*/);
      if (!ok) 
      {
        std::cout << "ERROR : could not get current point list" << std::endl; 
        break;
      }

      zones[i]->add_boundary(bccount++, ptL, nb_ids);
    }

  }

  err = !ok;

  PyObject *node_list(PyList_New(0));

  if (!err)
  {

    // MAJ des F2E pour orienter correctement les elements frontière
    for (E_Int i=0; i < nb_zones; ++i) zones[i]->reorient_f2e_boundaries();
    
    // MAJ des F2E pour orienter correctement les elements frontière
    for (E_Int i=0; i < nb_zones; ++i) zones[i]->set_f2e_across_joins();

    
    //std::cout << "Adding layers.." << std::endl;
     
    // Adding layers
    Vector_t<zone_type> tmp_zones;

    for (E_Int i = 0; i < nb_zones; ++i)
    {

      tmp_zones.clear();

      for (E_Int j = 0; j < nb_zones; ++j) tmp_zones.push_back(*zones[j]);

      for (E_Int j = 0; j < nb_zones; ++j)
      {
      	for (E_Int k = 0; k < nb_zones; ++k)
      		tmp_zones[j].change_joins_zone(zones[k], &tmp_zones[k]);
      }
      
      add_n_topo_layers(tmp_zones, i, NLAYERS);
      zone_type& Zghost = tmp_zones[i];

      if (!err)
      {
        PyObject *onode(PyList_New(0));

        K_FLD::IntArray cnto;
        Zghost._ng.export_to_array(cnto);
        PyObject *mesh = K_ARRAY::buildArray(Zghost._crd, varString, cnto, 8, "NGON", false);
        PyList_Append(onode, mesh);
        Py_DECREF(mesh);

        PyObject *joins(PyList_New(0));

        for (auto j : Zghost._joins)
        {
          // (oid, zid, ids)

          //std::cout << "adding join " << j->first << std::endl;

          PyObject * jids = K_NUMPY::buildNumpyArray(&j.second.second[0], j.second.second.size(), 1, 0);
          E_Int iid = j.first;
          auto it1 = Zghost._rac_inId2outId.find(iid);
          if (it1 == Zghost._rac_inId2outId.end()) continue;
          E_Int j_oid = Zghost._rac_inId2outId[iid];

          std::pair<E_Int, E_Int>& p = j_to_pair[j_oid];
          E_Int zd = (p.first == i) ? p.second : p.first; // donnor zone

          //std::cout << "adding join : donnor zone : " << zd << std::endl;

          PyObject* oid_ids = PyTuple_New(3);
          
          PyObject* jid = PyInt_FromLong(j_oid);
          PyObject* zid = PyInt_FromLong(zd);

          PyTuple_SetItem(oid_ids, 0, jid);
          PyTuple_SetItem(oid_ids, 1, zid);
          PyTuple_SetItem(oid_ids, 2, jids);

          //std::cout << "adding join b" << std::endl;

          PyList_Append(joins, oid_ids);
          //Py_DECREF(oid_ids);
          //Py_DECREF(jids);

          //std::cout << "adding join c" << std::endl;
        }

        PyList_Append(onode, joins);
        Py_DECREF(joins);

        PyObject *bcs(PyList_New(0));

        for (auto j : Zghost._bcs)
        {
          PyObject * obc = K_NUMPY::buildNumpyArray(&j.second[0], j.second.size(), 1, 0);

          E_Int bc_oid = Zghost._bc_inId2outId[j.first];

          PyObject* oid_ids = PyTuple_New(2);
          PyObject* id = PyInt_FromLong(bc_oid);
          PyTuple_SetItem(oid_ids, 0, id);
          PyTuple_SetItem(oid_ids, 1, obc);
          PyList_Append(bcs, oid_ids);
          Py_DECREF(oid_ids);
        }

        PyList_Append(onode, bcs);
        PyList_Append(node_list, onode);

      }
    }
  }

  for (E_Int i=0; i < nb_zones; ++i)
  {
    delete ngs[i];
    delete cnts[i];
    delete crds[i];
    delete zones[i];
  }
   
  return node_list;
}

//=======================  Intersector/PolyMeshTools/utils.cpp ====================
