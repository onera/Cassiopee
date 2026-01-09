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
//Authors : Sam Landier (sam.landier@onera.fr)

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

    xcellnv(E_Float RTOL) : parent_t(RTOL) {}

    outdata_t __process_X_cells(zmesh_t const & z_mesh, std::vector< bound_mesh_t*> const & mask_bits, wdata_t & wdata)
    {
      size_t ncells = z_mesh.ncells();
      assert(ncells == wdata.size());

      z_mesh.get_nodal_metric2();

      outdata_t xcelln(ncells, OUT);
      if (mask_bits.empty()) return xcelln;// completely visible

      std::vector<E_Int> cands;
      std::vector<aelt_t> bits, tmpbits; // one clip can produce several bits

      for (size_t i = 0; i < ncells; ++i)
      {
        color_t const & idata = wdata[i];

        E_Float v = E_Float(idata);
        xcelln[i] = v;
        if (v != X) continue;

        bits.clear();
        bits.push_back(z_mesh.aelement(i)); // autonomous element directly stolen by bits (rvalue ref) 
        
        E_Float v0 = bits[0].metrics(); // initial surface&normal (surf)/volume 

        NUGA::CLIP::poly_clip<zmesh_t, bound_mesh_t>(bits, v0, idata.masks, mask_bits, parent_t::_RTOL);

        // accumulated volume
        E_Float vcur = 0.;
        for (size_t b = 0; b < bits.size(); ++b)
          vcur += bits[b].extent();

        xcelln[i] = vcur / v0;

#ifdef DEBUG_XCELLN
        if (xcelln[i] >= 1.1)
        {
          std::cout << "valeur bizarre pour : " << i << std::endl;
          std::cout << "v0 : " << v0 << std::endl;
          std::cout << "vcur : " << vcur << std::endl;
          std::ostringstream o; o << "bizarre_" << i;
          medith::write(o.str().c_str(), z_mesh.crd, z_mesh.cnt, i);
        }
#endif
        assert(xcelln[i] < 1.1);
        xcelln[i] = std::min(1., xcelln[i]);
      }

      return xcelln;
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

    xcellno(E_Float RTOL) : parent_t(RTOL) {}

    void finalize(const zmesh_t& m, outdata_t& outdata)
    {
      eClassify loc = OUT;
      K_SEARCH::BBox3D zbx;
      outdata.mesh.bbox(zbx); // use truncated mesh instead of m
      if (outdata.full_out) {
        //std::cout << "xcellno testing" << std::endl;
        loc = NUGA::CLASSIFY::classify(zbx, parent_t::_comp_boxes, parent_t::_z_priorities);
      }
      if (loc == IN)
        outdata.mesh.clear();
    }

    outdata_t __process_X_cells(zmesh_t const & z_mesh, std::vector< bound_mesh_t*> const & mask_bits, wdata_t & wdata)
    {
      size_t ncells = z_mesh.ncells();
      assert(ncells == wdata.size());

      outdata_t xmesh;
      bool docomp(false);
      std::vector<bool> keep(ncells, false);
      for (size_t i = 0; i < ncells; ++i)
      {
        color_t const & idata = wdata[i];
        E_Float v = E_Float(idata);

        keep[i] = (v == OUT);
        docomp |= (v != OUT);
      }

      xmesh.mesh = z_mesh;
      K_CONNECT::IdTool::init_inc(xmesh.mesh.flag, xmesh.mesh.ncells()); // for history
      if (docomp) xmesh.mesh.compress(keep);//contains now only non-X elements

      xmesh.full_out = (ncells == (size_t)xmesh.mesh.ncells());
      if (xmesh.full_out) return xmesh;

      z_mesh.get_nodal_metric2();

      std::vector<E_Int> cands;
      std::vector<aelt_t> bits, tmpbits; // one clip can produce several bits

      for (size_t i = 0; i < ncells; ++i)
      {
        color_t const & idata = wdata[i];

        E_Float v = E_Float(idata);
        if (v != X) continue;

        bits.clear();
        bits.push_back(z_mesh.aelement(i)); // autonomous element directly stolen by bits (rvalue ref) 

        E_Float v0 = bits[0].metrics(); // initial surface&normal (surf)/volume 

        NUGA::CLIP::poly_clip<zmesh_t, bound_mesh_t>(bits, v0, idata.masks, mask_bits, parent_t::_RTOL);

        for (size_t b = 0; b < bits.size(); ++b)
        {
#ifdef DEBUG_XCELLN
          std::ostringstream o;
          o << "cut_" << i;
         // medith::write(o.str().c_str(), bits[b]);
#endif

          xmesh.mesh.add(bits[b], true/*append coords*/);
        }
        xmesh.mesh.flag.resize(xmesh.mesh.ncells(), i); //use flag for history
      }

      xmesh.mesh.cnt.updateFacets();
      std::vector<E_Int> nids;
      zmesh_t::trait::compact_to_used_nodes(xmesh.mesh.cnt, xmesh.mesh.crd, nids);
      return xmesh;
    }

  };

 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  using IntVec = std::vector<E_Int>;
  using prior_t = std::map< E_Int, IntVec>;

  using ii_pair_t = std::pair<E_Int, E_Int>;
  using no_ii_pair_t = K_MESH::NO_Edge;

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
    rank_wnps.resize(max_comp_id + 1, 0);

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
    decrease_prior_per_comp.clear();
  
    for (size_t i = 0; i < clean_priority.size(); ++i)
    {
      E_Int Lcompid = clean_priority[i].first;
      E_Int Rcompid = clean_priority[i].second;

      //std::cout << " pair : " << Lcompid << "/" << Rcompid << std::endl;

      decrease_prior_per_comp[Lcompid].push_back(Rcompid);
      ++rank_wnps[Lcompid];
    }

    E_Int i = 0;
    for (auto it = decrease_prior_per_comp.begin(); it != decrease_prior_per_comp.end(); ++it, ++i)
    {
      //E_Int compid = it->first;
      std::sort(ALL(it->second), 
        [&decrease_prior_per_comp] (E_Int a, E_Int b)
        {
        auto it = decrease_prior_per_comp.find(a);
        if (it == decrease_prior_per_comp.end()) return false;

        for (size_t i = 0; i < it->second.size(); ++i)
          if (it->second[i] == b) return true;

        return false;
        });
      std::reverse(ALL(it->second)); // WARNING DECREASING PRIORITY DONE HERE
      //std::cout << "WNP rank for comp " << compid << " is " << rank_wnps[compid] << std::endl;
    }

    // append with opposite pairs for considering WNPS
    for (size_t j = 0; j < clean_priority.size(); ++j)
    {
      E_Int Lcompid = clean_priority[j].first;
      E_Int Rcompid = clean_priority[j].second;
      decrease_prior_per_comp[Rcompid].push_back(Lcompid);
    }
  }

  ///
  template <typename classifyer_t> inline
  E_Int MOVLP_xcelln_1zone(typename classifyer_t::zmesh_t& z_mesh,
    const IntVec& z_priorities, E_Int rank_wnp,
    const std::vector<K_FLD::FloatArray> &mask_crds, const std::vector<K_FLD::IntArray>& mask_cnts,
    std::vector< std::vector<E_Int>> &mask_wall_ids,
    const std::vector<K_SEARCH::BBox3D>& comp_boxes,
    typename classifyer_t::outdata_t& z_xcelln, E_Float RTOL)
  {
    //using zmesh_t = typename classifyer_t::zmesh_t;
    using bmesh_t = typename classifyer_t::bmesh_t;

    E_Int err(0);

    classifyer_t classs(RTOL);

#ifdef DEBUG_XCELLN
    std::cout << "MOVLP_xcelln_1zone : PREPARE " << std::endl;
#endif

    std::vector<bmesh_t*> mask_meshes;
    bmesh_t WP, WNP;
    classs.prepare(z_mesh, mask_crds, mask_cnts, mask_wall_ids, comp_boxes, z_priorities, rank_wnp, mask_meshes, WP, WNP);

#ifdef DEBUG_XCELLN
    std::cout << "MOVLP_xcelln_1zone : COMPUTE " << std::endl;
#endif

    err = classs.compute(z_mesh, mask_meshes, WP, WNP, z_xcelln);

#ifdef DEBUG_XCELLN
    std::cout << "MOVLP_xcelln_1zone : FINALIZE " << std::endl;
#endif

    if(!err) classs.finalize(z_mesh, z_xcelln); //replace IN/OUT with 0./1.

    for (size_t u = 0; u < mask_meshes.size(); ++u)
      delete mask_meshes[u];

    return err;
  }

  ///
  template <typename classifyer_t> inline
  void MOVLP_xcelln_zones(const std::vector<K_FLD::FloatArray> &crds, const std::vector<K_FLD::IntArray>& cnts,
                          const std::vector< std::vector<E_Int>> &zone_wall_ids,
                          const std::vector<E_Int>& comp_id, std::vector<std::pair<E_Int, E_Int>> & priority,
                          const std::vector<K_FLD::FloatArray> &mask_crds, const std::vector<K_FLD::IntArray>& mask_cnts,
                          std::vector< std::vector<E_Int>> &mask_wall_ids,
                          std::vector<typename classifyer_t::outdata_t > & xcelln, E_Float RTOL)
  {
    //priority map
    prior_t sorted_comps_per_comp;
    IntVec rank_wnps;
    comp_priorities(priority, sorted_comps_per_comp, rank_wnps);

    E_Int nb_zones = crds.size();
    //std::cout << "total nb of zones : " << nb_zones << std::endl;
    //std::cout << "nb of masks : " << mask_crds.size() << std::endl;

    assert(zone_wall_ids.size() == (size_t)nb_zones);

    xcelln.resize(nb_zones); // xcelln can be clipped mesh

    using zmesh_t = typename classifyer_t::zmesh_t;

    // compute component boxes
    E_Int nb_comps = *std::max_element(ALL(comp_id)) + 1;
    std::vector<K_SEARCH::BBox3D> comp_boxes(nb_comps);

    if (typeid(zmesh_t) == typeid(ph_mesh_t))
    {
      for (E_Int c = 0; c < nb_comps; ++c)
        comp_boxes[c].compute(mask_crds[c]);//skin is sufficient to compute bbox in 3D
    }
    else
    {
      for (E_Int z = 0; z < nb_zones; ++z)
      {
        K_SEARCH::BBox3D zbx;
        zbx.compute(crds[z]); //a close surface gives a null mask so compute box on patch which is light to compute
        comp_boxes[comp_id[z]].merge(zbx);
      }
    }


#ifndef NETBEANSZ
//#pragma omp parallel for
#endif
    for (E_Int z = 0; z < nb_zones; ++z)
    {

#ifdef DEBUG_XCELLN
      //if (z != 29) continue;
      std::cout << "processing zone : " << z << " from comp " << comp_id[z] << std::endl;
#endif

      auto it = sorted_comps_per_comp.find(comp_id[z]);
      if (it == sorted_comps_per_comp.end()) continue;

      E_Int z_rank_wnp = rank_wnps[it->first];
      IntVec& z_priorities = it->second;

#ifdef DEBUG_XCELLN
      //NUGA::chrono c;
      //c.start();

      //std::cout << "priorities : ";
      //for (size_t u = 0; u < it->second.size(); ++u)std::cout << it->second[u] << "/";
      //std::cout << std::endl;
      //std::cout << "rank WNP : " << z_rank_wnp << std::endl;
#endif

      using zmesh_t = typename classifyer_t::zmesh_t;

      zmesh_t zm(crds[z], cnts[z], 1/*ASSUME ORIENTED*/);

      //std::cout << "nb of walls : " << zone_wall_ids[z].size() << std::endl;

      zm.set_boundary_type(BCWALL, zone_wall_ids[z]);// only implemented for volumes due to double wall management
            
      E_Int err = MOVLP_xcelln_1zone<classifyer_t>(zm, z_priorities, z_rank_wnp, mask_crds, mask_cnts, mask_wall_ids, comp_boxes, xcelln[z], RTOL);
      if (err) break;

      //std::cout << "processed in : " << c.elapsed() << " s." << std::endl;
    }
  }

}

#endif
