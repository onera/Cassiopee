/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_JOIN_HXX
#define NUGA_JOIN_HXX

#include <vector>

namespace NUGA
{

  ///
  template <typename mesh_t>
  struct join_t
  {
    E_Int                               id;
    mesh_t&                             mesh;
    std::map<E_Int, std::vector<E_Int>> jid_to_joinlist;
    E_Int                               idx_start;

    join_t(E_Int i, mesh_t& m, E_Int idx_strt) :id(i), mesh(m), idx_start(idx_strt) {}

    void link(E_Int jzid, std::vector<E_Int>& pointlist);

    void update(); // refresh jid_to_joinlist
  };

  ///
  template <typename mesh_t>
  void join_t<mesh_t>::link(E_Int jid, std::vector<E_Int>& pointlist)
  {
    jid_to_joinlist[jid] = pointlist;

    // shift_geom on joins (and therefore reapply reorder_pgs on attached PHs)
    for (size_t i = 0; i < pointlist.size(); ++i)
    {
      E_Int PGi = pointlist[i] - idx_start;
      E_Int* nodes = mesh._ng.PGs.get_facets_ptr(PGi);
      E_Int nnodes = mesh._ng.PGs.stride(PGi);
      K_MESH::Polygon::shift_geom(mesh._crd, nodes, nnodes, idx_start);

      E_Int PHL = mesh._F2E(0, PGi);
      E_Int PHR = mesh._F2E(1, PGi); //should be IDX_NONE for oriented joins

      if (PHL != IDX_NONE)
        K_MESH::Basic::reorder_pgs(mesh._ng, mesh._F2E, PHL);
      if (PHR != IDX_NONE)
        K_MESH::Basic::reorder_pgs(mesh._ng, mesh._F2E, PHR);
    }
  }

  ///
  template <typename mesh_t>
  void join_t<mesh_t>::update()
  {
    //std::cout << "join_t<mesh_t>::update() : begin" << std::endl;
    std::vector<E_Int> ids;

    // update pointlists
    for (auto& j : jid_to_joinlist)
    {
      E_Int jzid = j.first;
      auto& ptlist = j.second;

      std::vector<E_Int> new_ptlist;
      for (size_t i = 0; i < ptlist.size(); ++i)
      {
        E_Int PGi = ptlist[i] - idx_start;
        //std::cout << "PGi : " << PGi << std::endl;        
        if (mesh._PGtree.is_enabled(PGi))
          new_ptlist.push_back(PGi);
        else // look in the genealogy where are the enabled
        {
          ids.clear();
          bool reverse = (jzid < id); // the lower id dictate the sorting
          mesh.extract_enabled_pgs_descendance(PGi, reverse, ids);

          if (!ids.empty()) //refinement
            new_ptlist.insert(new_ptlist.end(), ALL(ids));
          else //agglo : get the enabled ancestor
          {
            E_Int pid{IDX_NONE};
            mesh._PGtree.get_enabled_parent(PGi, pid);
            assert (pid != IDX_NONE);
            new_ptlist.push_back(pid);
          }
        }
      }

      ptlist = new_ptlist;
      // remove duplicates due to agglo
      K_CONNECT::IdTool::compress_unic(ptlist);
      // for the outside world
      K_CONNECT::IdTool::shift(ptlist, idx_start);

      //std::cout << "join_t<mesh_t>::update() : end" << std::endl;

      /*for (size_t i=0; i< ptlist.size(); ++i)
        std::cout << ptlist[i] << "/";
      std::cout << std::endl;*/
    }
  }
}
#endif