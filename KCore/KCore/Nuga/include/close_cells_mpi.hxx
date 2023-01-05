/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_CLOSE_CELLS_MPI_HXX
#define NUGA_CLOSE_CELLS_MPI_HXX

#include "Nuga/include/hybrid_para_algo.hxx"

namespace NUGA
{


  template <typename mesh_t>
  class close_cells_mpi : public hybrid_para_algo<mesh_t, atomdata<1>>
  {
  public:

    using data_t = typename atomdata<1>::type;// == K_FLD::FloatArray;

    inline int get_data_stride() override { return 3; }

    ///
    bool prepare_data_to_send
    (
      const mesh_t & mesh,
      const std::map<E_Int, std::vector<E_Int>>& rid_to_list,
      std::map<int, std::map<int, data_t>> & rid_to_PG_to_plan
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
          data_t p; //only implemented for IntArray
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

    void autonomous_run(const std::vector<mesh_t*>& meshes, int i) override
    {
      ngon_type::close_phs(meshes[i]->cnt, meshes[i]->crd);
    }

    bool run_with_data
    (
      const std::vector<mesh_t*>& meshes,
      const std::map<int, std::map<int, data_t>> & zid_to_data
    ) override
    {
      bool has_changes = false;

      E_Int NBZ{ E_Int(meshes.size()) };

//#pragma omp parallel for reduction ( || : has_changes) if(PARA == COARSE_OMP)
      //no reduction mode #pragma omp parallel for if(PARA == COARSE_OMP)          
      for (E_Int i = 0; i < NBZ; ++i)
      {
        mesh_t& mesh = *meshes[i];
        auto& crd = mesh.crd;

        auto it_PG_to_plan = zid_to_data.find(i);
        if (it_PG_to_plan == zid_to_data.end()) continue;

        auto & PG_to_plan = it_PG_to_plan->second;

        E_Int npgs = mesh.cnt.PGs.size();

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
          K_MESH::Polygon::imprint(crd, pnodes, nnodes, p1, molecPG);

          if (molecPG.empty())
          {
            new_pgs.add(nnodes, pnodes);
            continue;
          }

          new_pgs.add(molecPG.size(), &molecPG[0]);
          modified = true;
        }

        if (!modified) continue;

        mesh.cnt.PGs = new_pgs;
        mesh.cnt.PGs.updateFacets();
        has_changes = true;

        ngon_type::close_phs(mesh.cnt, mesh.crd); //fixme : fill PHto_process
      }

      return has_changes;
    }
    
    ///
    void extract_loop_plan(const K_FLD::FloatArray& crd, const ngon_unit& PGs, E_Int PGi, bool reverse, data_t& plan) const
    {
      plan.clear(); // plan est un IntArray (DynArray<E_Int>) : une 'matrice' nnodes x 1

      E_Int nnodes = PGs.stride(PGi);

      if (nnodes < 4) return; //fixme

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