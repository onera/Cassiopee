/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef NUGA_XCELLN_HXX
#define NUGA_XCELLN_HXX

#include "Nuga/include/classifyer.hxx"
#include "Nuga/include/clipper.hxx"

#ifdef DEBUG_XCELLN
#include "Nuga/include/medit.hxx"
#endif

namespace NUGA
{
  ///
  template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
  class xcellnv : public classifyer<XCELLN_VAL, zmesh_t, bound_mesh_t>
  {
  public:
    using parent_t  = classifyer<XCELLN_VAL, zmesh_t, bound_mesh_t>;
    using wdata_t   = typename parent_t::wdata_t;
    using outdata_t = typename parent_t::outdata_t;
    using aelt_t    = typename zmesh_t::aelt_t;

    xcellnv(double RTOL) : parent_t(RTOL) {}

    outdata_t __process_X_cells(zmesh_t const & z_mesh, std::vector< bound_mesh_t*> const & mask_bits, wdata_t & wdata)
    {
      E_Int ncells = z_mesh.ncells();
      assert(ncells == wdata.size());

      z_mesh.get_nodal_tolerance();

      outdata_t xcelln(ncells, OUT);

      std::vector<int> cands;
      std::vector<aelt_t> bits, tmpbits; // one clip can produce several bits

      for (E_Int i = 0; i < ncells; ++i)
      {
        color_t const & idata = wdata[i];

        double v = double(idata);
        xcelln[i] = v;
        if (v != X) continue;

        bits.clear();
        bits.push_back(z_mesh.aelement(i)); // autonomous element directly stolen by bits (rvalue ref) 
        
        double v0 = bits[0].metrics(); // initial surface&normal (surf)/volume 

        NUGA::CLIP::poly_clip<zmesh_t, bound_mesh_t>(bits, v0, idata.masks, mask_bits, parent_t::_RTOL);

        // accumulated volume
        double vcur = 0.;
        for (size_t b = 0; b < bits.size(); ++b)
          vcur += bits[b].extent();

        xcelln[i] = vcur / v0;
      }

      return std::move(xcelln);
    }
    
  };

  ///
  template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
  class xcellno : public classifyer<XCELLN_OUT, zmesh_t, bound_mesh_t>
  {
  public:
    using parent_t  = classifyer<XCELLN_OUT, zmesh_t, bound_mesh_t>;
    using wdata_t   = typename parent_t::wdata_t;
    using outdata_t = typename parent_t::outdata_t; // zmesh_t
    using aelt_t    = typename zmesh_t::aelt_t;

    xcellno(double RTOL) : parent_t(RTOL) {}

    void finalize(outdata_t& outdata)
    {
      //fixme : hack for fully inside non prior
      // inferior but uncolored => assume it means fully in so IN
      if (outdata.full_out && parent_t::_is_inferior)
        outdata.mesh.clear();
    }

    outdata_t __process_X_cells(zmesh_t const & z_mesh, std::vector< bound_mesh_t*> const & mask_bits, wdata_t & wdata)
    {
      E_Int ncells = z_mesh.ncells();
      assert(ncells == wdata.size());

      z_mesh.get_nodal_tolerance();

      outdata_t xmesh;

      std::vector<bool> keep(ncells, false);
      for (E_Int i = 0; i < ncells; ++i)
      {
        color_t const & idata = wdata[i];
        double v = double(idata);
        keep[i] = (v == OUT);
      }

      xmesh.mesh = z_mesh;
      K_CONNECT::IdTool::init_inc(xmesh.mesh.flag, xmesh.mesh.ncells()); // for history
      xmesh.mesh.compress(keep);

      xmesh.full_out = (ncells == xmesh.mesh.ncells());

      std::vector<int> cands;
      std::vector<aelt_t> bits, tmpbits; // one clip can produce several bits

      for (E_Int i = 0; i < ncells; ++i)
      {
        color_t const & idata = wdata[i];

        double v = double(idata);
        if (v != X) continue;

        bits.clear();
        bits.push_back(z_mesh.aelement(i)); // autonomous element directly stolen by bits (rvalue ref) 

        double v0 = bits[0].metrics(); // initial surface&normal (surf)/volume 

        NUGA::CLIP::poly_clip<zmesh_t, bound_mesh_t>(bits, v0, idata.masks, mask_bits, parent_t::_RTOL);

        for (size_t b = 0; b < bits.size(); ++b)
          xmesh.mesh.add(bits[b]);
        xmesh.mesh.flag.resize(xmesh.mesh.ncells(), i); //use flag for history
      }

      return std::move(xmesh);
    }

  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  using IntVec = std::vector<E_Int>;
  using prior_t = std::map< E_Int, IntVec>;

  using ii_pair_t = std::pair<E_Int, E_Int>;
  using no_ii_pair_t = K_MESH::NO_Edge;

  struct aIsLessThanb : std::binary_function<int, int, bool>
  {

    aIsLessThanb(const prior_t p) :_p(p) {}

    bool operator()(int a, int b) const {

      auto it = _p.find(a);
      if (it == _p.end()) return false;

      for (size_t i = 0; i < it->second.size(); ++i)
        if (it->second[i] == b) return true;

      return false;
    }

    prior_t _p;
  };

  inline void comp_priorities(const std::vector<ii_pair_t> & priority, prior_t & decrease_prior_per_comp, IntVec& rank_wnps)
  {
    // WARNING DECREASING PRIORITY UPON EXIT

    // 0. get the max comp id
    E_Int max_comp_id = -1;
    for (auto it = priority.begin(); it != priority.end(); ++it)
    {
      max_comp_id = std::max(max_comp_id, it->first);
      max_comp_id = std::max(max_comp_id, it->second);
    }

    E_Int nb_comps = max_comp_id + 1;

    decrease_prior_per_comp.clear();
    rank_wnps.resize(nb_comps, 0);

    // // 1. CLEANING INPUT (discard redundancy & bilateral)
    std::set<no_ii_pair_t> upairs;
    for (auto it = priority.begin(); it != priority.end(); ++it) // make pairs unique and unilateral
      upairs.insert(no_ii_pair_t(it->first, it->second));
    // // in case of bilateral input keep the first
    std::vector<ii_pair_t> clean_priority;
    for (auto it = priority.begin(); it != priority.end(); ++it)
    {
      no_ii_pair_t p(it->first, it->second);
      if (upairs.find(p) == upairs.end()) continue; // p is redundant here and first occ has been remove in upairs so discard p
      upairs.erase(p);
      clean_priority.push_back(*it);
    }

    // // 2. 
    // // store the priors per comp
    for (size_t i = 0; i < clean_priority.size(); ++i)
    {
      E_Int Lcompid = clean_priority[i].first;
      E_Int Rcompid = clean_priority[i].second;

      //std::cout << " pair : " << Lcompid << "/" << Rcompid << std::endl;

      decrease_prior_per_comp[Lcompid].push_back(Rcompid);
      ++rank_wnps[Lcompid];

      decrease_prior_per_comp[Rcompid].push_back(Lcompid); //opposite pair : for considering WNPs
    }

    // IntVec all_comps;
    // for (auto it = decrease_prior_per_comp.begin(); it != decrease_prior_per_comp.end(); ++it)
    //   all_comps.push_back(it->first);

    aIsLessThanb pred(decrease_prior_per_comp);
    rank_wnps.resize(1 + decrease_prior_per_comp.rbegin()->first, 0); //sized as the max comp num

    E_Int i = 0;
    for (auto it = decrease_prior_per_comp.begin(); it != decrease_prior_per_comp.end(); ++it, ++i)
    {
      E_Int compid = it->first;
      std::sort(ALL(it->second), pred);
      std::reverse(ALL(it->second)); // WARNING DECREASING PRIORITY DONE HERE
                                     //std::cout << "WNP rank for comp " << it->first << " is " << rank_wnps[compid] << std::endl;
    }
  }

  using zmesh_t = NUGA::pg_smesh_t;
  using bmesh_t = NUGA::edge_mesh_t;
  //using zmesh_t = NUGA::ph_mesh_t;
  //using bmesh_t = NUGA::pg_smesh_t;

  ///
  template <typename classifyer_t> inline 
  void MOVLP_xcelln_1zone(K_FLD::FloatArray& z_crd, K_FLD::IntArray& z_cnt,
    const IntVec& z_priorities, E_Int rank_wnp,
    const std::vector<K_FLD::FloatArray*> &mask_crds, const std::vector<K_FLD::IntArray*>& mask_cnts,
    std::vector< std::vector<E_Int>> &mask_wall_ids,
    typename classifyer_t::outdata_t& z_xcelln, E_Float RTOL)
  {
    zmesh_t  z_mesh(z_crd, z_cnt);  // polygonal surface mesh

    classifyer_t classs(RTOL);

    std::vector<bmesh_t*> mask_meshes;
    classs.prepare(z_mesh, mask_crds, mask_cnts, mask_wall_ids, z_priorities, rank_wnp, mask_meshes);
    classs.compute(z_mesh, mask_meshes, z_xcelln);
    classs.finalize(z_xcelln); //replace IN/OUT with 0./1.

    for (size_t u = 0; u < mask_meshes.size(); ++u)
      delete mask_meshes[u];
  }

  ///
  template <typename classifyer_t> inline
  void MOVLP_xcelln_zones(const std::vector<K_FLD::FloatArray*> &crds, const std::vector<K_FLD::IntArray*>& cnts,
    const std::vector<E_Int>& comp_id, std::vector<std::pair<E_Int, E_Int>> & priority,
    const std::vector<K_FLD::FloatArray*> &mask_crds, const std::vector<K_FLD::IntArray*>& mask_cnts,
    std::vector< std::vector<E_Int>> &mask_wall_ids,
    std::vector<typename classifyer_t::outdata_t > & xcelln, E_Float RTOL)
  {
    //priority map
    prior_t sorted_comps_per_comp;
    IntVec rank_wnps;
    comp_priorities(priority, sorted_comps_per_comp, rank_wnps);

    E_Int nb_zones = crds.size();

    xcelln.resize(nb_zones); // xcelln can beclipped mesh

#ifndef NETBEANSZ
//#pragma omp parallel for
#endif
    for (E_Int z = 0; z < nb_zones; ++z)
    {
      //std::cout << "processing zone : " << z << " from comp " << comp_id[z] << std::endl;

      auto it = sorted_comps_per_comp.find(comp_id[z]);
      if (it == sorted_comps_per_comp.end()) continue;

      E_Int z_rank_wnp = rank_wnps[it->first];
      IntVec& z_priorities = it->second;
      
      MOVLP_xcelln_1zone<classifyer_t>(*crds[z], *cnts[z], z_priorities, z_rank_wnp, mask_crds, mask_cnts, mask_wall_ids, xcelln[z], RTOL);
    }
  }

}

#endif
