

/*
 
 
 
              NUGA 
 
 
 
 */
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_CLIPPER_HXX
#define NUGA_CLIPPER_HXX

#include <vector>
#include "Nuga/Boolean/BAR_Conformizer.h"
#include "Nuga/Delaunay/T3Mesher.h"

#ifdef DEBUG_CLIPPER
#include "Nuga/include/medit.hxx"
#endif

namespace NUGA
{
  namespace CLIP
  {
    using crd_t = K_FLD::FloatArray;
    using cnt_t = K_FLD::IntArray;

    inline E_Int isolated_clip(const crd_t& crd1, const crd_t& crd2, const cnt_t& cutter, E_Float RTOL, std::vector<crd_t>& res, bool& true_clip)
    {
      true_clip = false;
      res.clear();

      E_Int nb_nodes = crd1.cols();
      cnt_t subj(2, nb_nodes);
      for (E_Int i = 0; i < crd1.cols(); ++i)
      {
        subj(0, i) = i;
        subj(1, i) = (i + 1) % nb_nodes;
      }

      int err(0);
      int nb_edges1(nb_nodes);

      //gather subj & cutter 
      K_FLD::FloatArray crd(crd1);
      E_Float W[3];
      //fixme //////////////////////////////////////
      {
        NUGA::aPolygon aPGsubj(std::move(crd));//fixme
        aPGsubj.normal<3>(W);                     //fixme
        crd = std::move(aPGsubj.m_crd);//give it back
      }
      //////////////////////////////////////////////
      E_Int nb_pts1 = crd.cols();
      crd.pushBack(crd2);
      K_FLD::IntArray cnt(subj), cnt2(cutter);

      for (size_t i = 0; i < cnt2.cols(); ++i) std::swap(cnt2(0, i), cnt2(1, i));//boolean diff mode

      cnt2.shift(nb_pts1);
      cnt.pushBack(cnt2);

      // got 2D (in sub ref frame)
      K_FLD::FloatArray P(3, 3), iP(3, 3);
      
      FittingBox::computeAFrame(W, P);
      iP = P;
      K_FLD::FloatArray::inverse3(iP);
      FittingBox::transform(crd, iP);
      E_Float z = crd(2, 0);
      crd.resize(2, crd.cols());

      // compute an overall abstol
      E_Float min_d, max_d, ABSTOL;
      K_CONNECT::MeshTool::computeMinMaxEdgeSqrLength<2>(crd, cnt, min_d, max_d);
      ABSTOL = ::sqrt(min_d) * RTOL;

      // conformize this cloud
      std::vector<E_Int> ancE2;
      NUGA::BAR_Conformizer<2> conformizer;
      conformizer._brute_force = true;

#ifdef DEBUG_CLIPPER
      conformizer._silent_errors = false;
#else
      conformizer._silent_errors = true;
#endif

      E_Int nb_edges0 = cnt.cols();
      err = conformizer.run(crd, cnt, ancE2, nullptr/*&priority*/, ABSTOL, nb_edges1, 1 /*one iter only*/);
      if (err)
        return err;

#ifdef DEBUG_CLIPPER
      {
        K_FLD::FloatArray tcrd(crd);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("in", tcrd, cnt, "BAR");
      }
#endif

      if (cnt.cols() == nb_edges0) return 0;// no intersections => fully visible or hidden

      true_clip = true;

      DELAUNAY::MeshData data;
      data.pos = &crd;
      data.connectB = &cnt;

      DELAUNAY::MesherMode mode;
      mode.mesh_mode = mode.TRIANGULATION_MODE;
      mode.remove_holes = true;
      mode.silent_errors = conformizer._silent_errors;
      DELAUNAY::T3Mesher<E_Float> mesher(mode);

      mesher.run(data);

      if (data.connectM.cols() == 0) return 0; //empty here means IN

#ifdef DEBUG_CLIPPER
      {
        K_FLD::FloatArray tcrd(*data.pos);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("out", tcrd, data.connectM, "TRI");
      }
#endif

      // go back 3D
      data.pos->resize(3, data.pos->cols(), z);
      FittingBox::transform(*data.pos, P); // back to original ref frame  

#ifdef DEBUG_CLIPPER
      //medith::write("boolean", *data.pos, data.connectM, "TRI");
#endif

      int minc = *std::min_element(ALL(data.colors));
      int maxc = *std::max_element(ALL(data.colors));

      if (minc != maxc) // non-connex
      {
        std::map<int, K_FLD::IntArray> col_to_cntB;
        E_Int nbe = data.connectM.cols();
        for (size_t i = 0; i < nbe; ++i)
        {
          int coli = data.colors[i];
          K_FLD::IntArray::const_iterator pS = data.connectM.col(i);
          for (size_t n = 0; n < 3; ++n)
          {
            E_Int& Kv = data.neighbors(n, i);
            if (Kv == E_IDX_NONE) //border
            {
              int E[] = { *(pS + (n + 1) % 3),*(pS + (n + 2) % 3) };
              col_to_cntB[coli].pushBack(E, E + 2);
            }
          }
        }

        std::vector<int> PGi;
        res.resize(col_to_cntB.size());
        int i = 0;
        std::vector<E_Int> nids;
        for (auto it = col_to_cntB.begin(); it != col_to_cntB.end(); ++it, ++i)
        {
          //
          nids.clear();
          K_FLD::FloatArray crd(*data.pos); //should not hurt as meshes are small
          K_CONNECT::MeshTool::compact_to_mesh(crd, it->second, nids);

          PGi.clear();
          // sort the nodes
          BARSplitter::getSortedNodes(it->second, PGi);

          E_Int str = crd.rows();
          res[i].reserve(str, PGi.size());
          for (size_t j = 0; j < PGi.size(); ++j)
            res[i].pushBack(crd.col(PGi[j]), crd.col(PGi[j]) + str);
        }
      }
      else //connex
      {
        K_FLD::IntArray cB;
        K_CONNECT::MeshTool::getBoundaryT3Mesh(data.connectM, data.neighbors, cB);

        std::vector<E_Int> nids;
        K_FLD::FloatArray& crd = *data.pos;
        K_CONNECT::MeshTool::compact_to_mesh(crd, cB, nids);

        std::vector<E_Int> PGi;
        // sort the nodes
        BARSplitter::getSortedNodes(cB, PGi);

        res.resize(1);
        E_Int str = crd.rows();
        res[0].reserve(str, PGi.size());
 
        for (size_t j = 0; j < PGi.size(); ++j)
          res[0].pushBack(crd.col(PGi[j]), crd.col(PGi[j]) + str);
      }

      return err;
    }

    template <typename aELT1, typename aELT2>
    bool compute(aELT1 const & subj, aELT2 const & cutter, E_Float RTOL, std::vector<aELT1>& bits)
    {
      bool true_clip(false);
      std::vector<crd_t> res;
      /*E_Int err = */isolated_clip(subj.m_crd, cutter.crd, cutter.cnt, RTOL, res, true_clip);

      bits.clear();
      bits.reserve(res.size());
      //for (size_t i = 0; i < res.size(); ++i)
        //bits.push_back(std::move(res[i]));
      std::move(ALL(res), std::back_inserter(bits));
 
      return true_clip;
    }

    template<typename aelt_t>
    void __replace_append_and_next_iter(std::vector<aelt_t>& bits, E_Int& b, std::vector<aelt_t>& toadd)
    {
      size_t sz = toadd.size();

      if (sz != 0) // not empty
      {
        auto it = toadd.begin();
        bits[b] = std::move(*it);                               // replace the current with the first bit to add
        if (sz == 1) return;
        bits.reserve(bits.size() + sz - 1);
        std::move(++it, toadd.end(), std::back_inserter(bits)); // push back remainers
        toadd.clear();
      }
      else
      {
        if (b < bits.size() - 1) // put the last (if exist) in the current, pop back the last
        {
          bits[b] = std::move(bits.back());
          --b; // to treat it at next iter
        }

        bits.pop_back(); // erase if single, remove the replaced one otherwise
      }
    }

    ///
    template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
    void poly_clip
    (std::vector<typename zmesh_t::aelt_t>& bits, double v0, 
     std::vector<E_Int> const & mask_ids, std::vector< bound_mesh_t*> const & mask_bits, double RTOL)
    {
      using loc_t = typename zmesh_t::loc_t;

      double vcur(v0), veps(v0*E_EPSILON);
      K_SEARCH::BBox3D bx0;
      bx0.compute(bits[0].m_crd);
      bx0.enlarge(RTOL);

      std::vector<E_Int> cands;
      std::vector<typename zmesh_t::aelt_t> tmpbits;

      size_t nmask = mask_ids.size();
      for (size_t m = 0; (m < nmask) && (vcur > veps); ++m)
      {
        E_Int mi = mask_ids[m];
        const bound_mesh_t& mask_bit = *mask_bits[mi];
        const loc_t& mask_loc = *(mask_bit.get_localizer());

#ifdef DEBUG_XCELLN
        medith::write("mask", mask_bit.crd, mask_bit.cnt);
#endif
        E_Int idx_start = mask_bit.index_start;
        //
        E_Int nbits = (E_Int)bits.size(); // bits can be appended when split give multiple bits
        for (E_Int b = 0; b < nbits; ++b)
        {
          auto& ae1 = bits[b];

          cands.clear();
          mask_loc.get_candidates(ae1, ae1.m_crd, cands, idx_start, RTOL);

          bool just_io = cands.empty();
          if (just_io) // try with inital box to catch the bounds and do i/o test
            mask_loc.get_candidates(bx0, cands, idx_start, RTOL);

          // autonomous cutter front
          bound_mesh_t acut_front(mask_bit, cands);

#ifdef DEBUG_XCELLN
          medith::write("subj", ae1.m_crd, &ae1.m_nodes[0], ae1.m_nodes.size(), 0);// or medith::write<ngon_type>("subj", z_mesh.crd, z_mesh.cnt, i);
          medith::write("cutter_front", mask_bit.crd, mask_bit.cnt, cands, idx_start);
#endif
          //CLIPPING
          tmpbits.clear();
          if (!just_io)
            just_io = !NUGA::CLIP::compute(ae1, acut_front, 1.e-4/*parent_t::_RTOL*/, tmpbits); //robust_clip returns true if true clip

          // IO : current bit does not intersect front.
          if (just_io)
            if (NUGA::CLASSIFY::classify(ae1, acut_front) == OUT) continue;

          // COMPRESS STRATEGY (with MOVE SEMANTICS) for bits => b can be reduced of 1 to treat the replaced at next iter 
          // if tmpbits is empty (IN) => compress (erase bits if single, put the last instead of current otherwise)
          // else replace the first bit, append the other bits  . 
          __replace_append_and_next_iter(bits, b, tmpbits);
          if (bits.size() < nbits) nbits = bits.size(); //update iff the last have replaced the current 
        }
      }
    }
  } // CLIP
} // NUGA

#endif /* NUGA_INTERSECT_HXX */

