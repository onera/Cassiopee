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

//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_CLOSE_CELLS_HXX
#define NUGA_CLOSE_CELLS_HXX

#include "Nuga/include/hybrid_para_algo.hxx"
#include "Nuga/include/ngon_t.hxx"

namespace NUGA
{


  template <typename para_algo_t, typename mesh_t>
  class close_cells : public para_algo_t
  {
  public:

    using parent_t = para_algo_t;
    using plan_t = K_FLD::FloatArray;
    using zid_to_rid_to_ptlist_t = typename parent_t::zid_to_rid_to_ptlist_t;
    using id_to_PG_to_plan_t = typename parent_t::id_to_PG_to_plan_t;
    using rid_to_zones_t = typename parent_t::rid_to_zones_t;

    inline int get_data_stride() override { return 3; }

    ///
    bool prepare_data_to_send
    (
      const mesh_t & mesh,
      const std::map<int, std::vector<E_Int>>& rid_to_list,
      id_to_PG_to_plan_t & rid_to_PG_to_plan
    ) override
    {
      bool has_packs{ false };
      //std::cout << "pack : loop on joins : " << join << std::endl;
      for (auto& it : rid_to_list)
      {
        E_Int rid = it.first;
        auto& ptlist = it.second;

        for (size_t i = 0; i < ptlist.size(); ++i)
        {
          plan_t p; //only implemented for IntArray
          E_Int PGi = ptlist[i] - 1;
          //std::cout << "i/PGi/sz : " << i << "/" << PGi << "/" << ptlist.size() << std::endl;
          
          extract_loop_plan(mesh.crd, mesh.cnt.PGs, PGi, true/*reverse*/, p);
          
          if (p.getSize())
          {
            rid_to_PG_to_plan[rid][i/*PGi*/] = p;
            has_packs = true;
          }
        }
      }
      return has_packs;
    }

    ///
    void autonomous_run(const std::vector<mesh_t*>& meshes, int i) override
    {

      {
        //0. store (in _type attrib.) nb of nodes before close_phs for each PG
        // this is used in prepare_data_to_send to select only those that have been refined
        E_Int npgs = meshes[i]->cnt.PGs.size();
        meshes[i]->cnt.PGs._type.clear();
        meshes[i]->cnt.PGs._type.resize(npgs);
        for (E_Int k = 0; k < npgs; ++k)
        {
          int nnodes = meshes[i]->cnt.PGs.stride(k);
          meshes[i]->cnt.PGs._type[k] = nnodes;
        }
      }

      //1. merge nodes : coincident vertices (appearing when refining faces sharing an edge by Polygon::imprint) need to be merge before close_phs call
      auto& nodal_metric2 = meshes[i]->get_nodal_metric2(eMetricType::ISO_MIN, true/*because coincident points exist here*/); //fixme hpc : currently nodal_metric2 recomputed each time. should be extended instead 

#ifdef DEBUG_CLOSECELLS
      E_Float minm = *std::min_element(ALL(nodal_metric2));
      assert (minm > 0.); // mandatory : to ensure coincident vertices will be merged
#endif
      
      /*E_Int nmerges = */meshes[i]->cnt.join_phs(meshes[i]->crd, nodal_metric2, EPSILON);//small tol to deal only with coincident vertices generated after Polygon::imprint call
      
      //2. close
      //std::cout << "close_cell : autonomous_run : close_phs" << std::endl;
      ngon_type::close_phs(meshes[i]->cnt, meshes[i]->crd);
    }

    ///
    bool run_with_data
    (
      const std::vector<mesh_t*>& meshes,
      const std::vector<int>& zids,
      const id_to_PG_to_plan_t & zid_to_PG_to_plan
    ) override
    {
      bool has_changes = false;

      E_Int NBZ{ E_Int(meshes.size()) };

//#pragma omp parallel for reduction ( || : has_changes) if(PARA == COARSE_OMP)
      //no reduction mode #pragma omp parallel for if(PARA == COARSE_OMP)          
      for (E_Int i = 0; i < NBZ; ++i)
      {
        int zid = zids[i];
        mesh_t& mesh = *meshes[i];
        auto& crd = mesh.crd;

        auto it_PG_to_plan = zid_to_PG_to_plan.find(zid);
        if (it_PG_to_plan == zid_to_PG_to_plan.end()) continue;

        auto & PG_to_plan = it_PG_to_plan->second;

        E_Int npgs = mesh.cnt.PGs.size();
        E_Int npts0 = mesh.crd.cols();

        ngon_unit new_pgs;
        std::vector<E_Int> molecPG;
        bool modified = false;

        for (E_Int j = 0; j < npgs; ++j)
        {
          E_Int nnodes = mesh.cnt.PGs.stride(j);
          const E_Int* pnodes = mesh.cnt.PGs.get_facets_ptr(j);

          auto it = PG_to_plan.find(j);
          if (it == PG_to_plan.end())
          {
            new_pgs.add(nnodes, pnodes);
            continue;
          }

          auto & p1 = it->second;
          
          molecPG.clear();
          K_MESH::Polygon::imprint(crd, pnodes, nnodes, p1, 1.e-6/*RTOL*/, molecPG);

          if (molecPG.empty())
          {
            new_pgs.add(nnodes, pnodes);
            continue;
          }

          new_pgs.add(molecPG.size(), &molecPG[0]);
        }

        E_Int npts1 = crd.cols();

        modified = (npts1 > npts0);

        if (!modified) continue;

        mesh.cnt.PGs = new_pgs;
        mesh.cnt.PGs.updateFacets();

        autonomous_run(meshes, i); // close_phs. fixme hpc : should pass PH_to_process to close_phs..

        has_changes = true;
      }

      return has_changes;
    }
    
    ///
    void extract_loop_plan(const K_FLD::FloatArray& crd, const ngon_unit& PGs, E_Int PGi, bool reverse, plan_t& plan) const
    {
      plan.clear(); // plan est un IntArray (DynArray<E_Int>) : une 'matrice' nnodes x 1

      E_Int nnodes = PGs.stride(PGi);

      if (PGs._type[PGi] == nnodes) return; // i.e. if this PG has not been refined, skip

      const E_Int* pnodes = PGs.get_facets_ptr(PGi);

      if (!reverse)
      {
        for (E_Int n = 0; n < nnodes; ++n)
        {
          const E_Float* p = crd.col(pnodes[n] - 1);
          plan.pushBack(p, p + 3);
        }
      }
      else
      {
        for (E_Int n = 0; n < nnodes; ++n)
        {
          const E_Float* p = crd.col(pnodes[nnodes - n - 1] - 1);
          plan.pushBack(p, p + 3);
        }
      }
    }


  };

}

#endif
