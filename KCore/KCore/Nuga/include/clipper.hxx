/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_CLIPPER_HXX
#define NUGA_CLIPPER_HXX

#include <vector>
#include "Nuga/include/BAR_Conformizer.h"
#include "Nuga/include/T3Mesher.h"

#include "Nuga/include/ph_clipper.hxx"
#include "Nuga/include/classifyer.hxx"


#ifdef DEBUG_CLIPPER
#include "Nuga/include/medit.hxx"
#endif

namespace NUGA
{
  namespace CLIP
  {
    using crd_t = K_FLD::FloatArray;

    template <typename aELT1, typename aELT2> //sub & cut are non const as they can be reoriented/reversed
    inline E_Int isolated_clip(aELT1& sub, aELT2& cut, NUGA::INTERSECT::eOPER oper, E_Float RTOL, std::vector<aELT1>& res, bool& true_clip);

    template <>
    inline E_Int isolated_clip<aPolygon, edge_mesh_t>(aPolygon& sub, edge_mesh_t& cut, NUGA::INTERSECT::eOPER oper, E_Float RTOL, std::vector<aPolygon>& bits, bool& true_clip)
    {
      using cnt_t = K_FLD::IntArray;
      const crd_t& crd1 = sub.m_crd;
      const crd_t& crd2 = cut.crd;
      const cnt_t& cutter = cut.cnt;

      true_clip = false;
      bits.clear();

      std::vector<crd_t> crd_res; //intermediate format as crds

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
        NUGA::aPolygon aPGsubj(std::move(crd));
        aPGsubj.normal<3>(W);                   
        crd = std::move(aPGsubj.m_crd);//give it back
      }
      //////////////////////////////////////////////
      E_Int nb_pts1 = crd.cols();
      crd.pushBack(crd2);
      K_FLD::IntArray cnt(subj), cnt2(cutter);

      if (oper == NUGA::INTERSECT::DIFFERENCE)
        connect_trait<LINEIC, true>::reverse_orient(cnt2);

      cnt2.shift(nb_pts1);
      cnt.pushBack(cnt2);

      // got 2D (in sub ref frame)
      K_FLD::FloatArray P(3, 3), iP(3, 3);
      
      NUGA::computeAFrame(W, P);
      iP = P;
      K_FLD::FloatArray::inverse3(iP);
      NUGA::transform(crd, iP);

      std::vector<double> zs;
      crd.extract_field(2, zs); //keep 3rd coord appart (altitude)
      //compute zmean
      double zmean(0.);
      for (E_Int k = 0; k < nb_pts1; ++k) zmean += zs[k];
      zmean /= nb_pts1;
      
      // compute an overall 3D abstol
      // using meL : NOT WORKING because after first cut, sub is overdefined so underestimate Lref, hence some good cut edges are discarded
      /*E_Float min_d, max_d, ABSTOL;
      NUGA::MeshTool::computeMinMaxEdgeSqrLength<3>(crd, cnt, min_d, max_d);
      double Lref = ::sqrt(min_d);*/
      // using sub bbox
      /* NOT WORKING neither when it's a planar case where the plane is axi-aligned => one corrdinate is the same => Lref gets null
      K_SEARCH::BBox3D bx;
      sub.bbox(sub.m_crd, bx);
      double Lref = std::min(bx.maxB[0] - bx.minB[0], std::min(bx.maxB[1] - bx.minB[1], bx.maxB[2] - bx.minB[2]));*/
      // using sub bbox in 2D frame
      K_SEARCH::BBox2D bx(crd, nb_pts1);
      double Lref = std::min(bx.maxB[0] - bx.minB[0], bx.maxB[1] - bx.minB[1]);
 
      double ABSTOL = Lref * RTOL;
      ABSTOL = std::max(ABSTOL, ZERO_M);

      // discard false overlaps among front (now we are in 2D, those with big altitudes)
      {
        std::vector<bool> keep(cnt.cols(), true);
        bool do_compact(false);
        for (E_Int i = nb_edges1; i < cnt.cols(); ++i)
        {
          double z1 = zs[cnt(0, i)] - zmean;
          double z2 = zs[cnt(1, i)] - zmean;

          if (z1*z2 < 0.) continue; // means crossing 

          double mz = std::min(::fabs(z1), ::fabs(z2));
          keep[i] = (mz < Lref); //at least one inside interf zone
          do_compact |= !keep[i];
        }
        if (do_compact)
        {
          std::vector<E_Int> nids;
          E_Int nb_edgesi = cnt.cols();
          K_FLD::IntArray::compact(cnt, keep, nids);
          // IMPORTANT : discard also those bits from cutter to synchronize it, in case of a I/O classify test afterwards
          assert (cnt2.cols() == (nb_edgesi - nb_edges1));
          std::vector<bool> kp(cnt2.cols(), true);// make keep info  relative to cutter
          for (E_Int i = nb_edges1; i < nb_edgesi; ++i) kp[i-nb_edges1]=keep[i];
          cut.compress(kp);
        }
      }
      // now apply zmean to front points such remaining ones at the end will be roughly on subj supporting surface
      for (size_t k = nb_pts1; k < zs.size(); ++k) zs[k] = zmean;

      crd.resize(2, crd.cols());//now pure 2D

      // conformize this cloud
      std::vector<E_Int> ancE2;
      NUGA::BAR_Conformizer<2> conformizer(true/* keep track of nodes history*/);
      conformizer._brute_force = true;

#ifdef DEBUG_CLIPPER
      conformizer._silent_errors = false;
#else
      conformizer._silent_errors = true;
#endif

#ifdef DEBUG_CLIPPER
      {
        K_FLD::FloatArray tcrd(crd);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("before_conf", tcrd, cnt, "BAR");
      }
#endif

      E_Int nb_edges0 = cnt.cols();
      err = conformizer.run(crd, cnt, ancE2, nullptr/*&priority*/, ABSTOL, nb_edges1, 1 /*one iter only*/);
      if (err)
        return err;

#ifdef DEBUG_CLIPPER
      {
        K_FLD::FloatArray tcrd(crd);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("after_conf", tcrd, cnt, "BAR");
      }
#endif

      if (cnt.cols() == nb_edges0) return 0;// no intersections => fully visible or hidden

      // use the T3mesher to classify

      DELAUNAY::MeshData data;
      data.pos = &crd;
      data.connectB = &cnt;

      DELAUNAY::MesherMode mode;
      mode.mesh_mode = mode.TRIANGULATION_MODE;
      mode.remove_holes = true;
      mode.silent_errors = true; // conformizer._silent_errors;
      DELAUNAY::T3Mesher<E_Float> mesher(mode);
      
      mesher.seed_random(1);
      err = mesher.run(data);
      if (err) return err;

      if (data.connectM.cols() == 0) return 0; //empty here means IN

      // 1 region after triangulation => so something might have happened on its boundary but not inside so consider as untouched
      if (data.mono_connex) return 0;

      true_clip = true;

#ifdef DEBUG_CLIPPER
      {
        K_FLD::FloatArray tcrd(*data.pos);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("out", tcrd, data.connectM, "TRI");
      }
#endif

      // transfer altitude
      std::vector<double> zsnew(data.pos->cols(), NUGA::FLOAT_MAX);
      const auto& pnids = conformizer.get_node_history();
      for (size_t i = 0; i < pnids.size(); ++i)
        if (pnids[i] != IDX_NONE) zsnew[pnids[i]] = zs[i];
      zs = std::move(zsnew);

      // interpolate missings using x history
      auto xedge = conformizer.get_x_history();
      for (size_t i = 0; i < zs.size(); ++i)
      {
        if (zs[i] != NUGA::FLOAT_MAX) continue;
        int xe = xedge[i];
        
        assert(xe != IDX_NONE); // i is an x point and has its histo
        
        if (xe >= subj.cols()) continue;

        int N0 = subj(0, xe);
        double * P0 = crd.col(N0);
        int N1 = subj(1, xe);
        double * P1 = crd.col(N1);
        double * Px = crd.col(i);

        double L2 = NUGA::sqrDistance(P0, P1, 2);
        double lam2 = NUGA::sqrDistance(P0, Px, 2);
        double l = ::sqrt(lam2 / L2);

        zs[i] = (1. - l) * zs[N0] + l * zs[N1];
      }

      
      for (size_t k = 0; k < zs.size(); ++k)
      {
        if (zs[k] == NUGA::FLOAT_MAX) zs[k] = zmean;
      }

      // go back 3D
      data.pos->resize(3, data.pos->cols());
      data.pos->set_field(2, zs); //
      NUGA::transform(*data.pos, P); // back to original ref frame  

#ifdef DEBUG_CLIPPER
      //medith::write("boolean", *data.pos, data.connectM, "TRI");
#endif

      int minc = *std::min_element(ALL(data.colors));
      int maxc = *std::max_element(ALL(data.colors));

      if (minc != maxc) // non-connex or glued on some edges due to a big tol
      {
        std::map<int, K_FLD::IntArray> col_to_cntB;
        E_Int nbe = data.connectM.cols();
        for (E_Int i = 0; i < nbe; ++i)
        {
          int coli = data.colors[i];
          K_FLD::IntArray::const_iterator pS = data.connectM.col(i);
          for (size_t n = 0; n < 3; ++n)
          {
            E_Int& Kv = data.neighbors(n, i);
            if (Kv == IDX_NONE || data.colors[Kv] != coli) //color border
            {
              int E[] = { *(pS + (n + 1) % 3),*(pS + (n + 2) % 3) };
              col_to_cntB[coli].pushBack(E, E + 2);
            }
          }
        }

        std::vector<int> PGi;
        crd_res.resize(col_to_cntB.size());
        int i = 0;
        std::vector<E_Int> nids;
        for (auto it = col_to_cntB.begin(); it != col_to_cntB.end(); ++it, ++i)
        {
          //
          nids.clear();
          K_FLD::FloatArray crd1(*data.pos); //should not hurt as meshes are small
          NUGA::MeshTool::compact_to_mesh(crd1, it->second, nids);

          PGi.clear();
          // sort the nodes
          BARSplitter::getSortedNodes(it->second, PGi);

          E_Int str = crd1.rows();
          crd_res[i].reserve(str, PGi.size());
          for (size_t j = 0; j < PGi.size(); ++j)
            crd_res[i].pushBack(crd1.col(PGi[j]), crd1.col(PGi[j]) + str);
        }
      }
      else //connex
      {
        K_FLD::IntArray cB;
        NUGA::MeshTool::getBoundaryT3Mesh(data.connectM, data.neighbors, cB);

        std::vector<E_Int> nids;
        K_FLD::FloatArray& crd1 = *data.pos;
        NUGA::MeshTool::compact_to_mesh(crd1, cB, nids);

        std::vector<E_Int> PGi;
        // sort the nodes
        BARSplitter::getSortedNodes(cB, PGi);

        crd_res.resize(1);
        E_Int str = crd1.rows();
        crd_res[0].reserve(str, PGi.size());
 
        for (size_t j = 0; j < PGi.size(); ++j)
          crd_res[0].pushBack(crd1.col(PGi[j]), crd1.col(PGi[j]) + str);
      }

      bits.reserve(crd_res.size());
      std::move(ALL(crd_res), std::back_inserter(bits));

      return err;
    }

    template <>
    inline E_Int isolated_clip<aPolygon, aPolygon>(aPolygon& sub, aPolygon& cut, NUGA::INTERSECT::eOPER oper, E_Float RTOL, std::vector<aPolygon>& bits, bool& true_clip)
    {
      edge_mesh_t ecut(cut);
      return isolated_clip<aPolygon, edge_mesh_t>(sub, ecut, oper, RTOL, bits, true_clip);
    }

    template <>
    inline E_Int isolated_clip<aPolyhedron<0>, pg_smesh_t>(aPolyhedron<0>& sub, pg_smesh_t& cut, NUGA::INTERSECT::eOPER oper, E_Float RTOL, std::vector<aPolyhedron<0>>& bits, bool& true_clip)
    {
      bool dbg(true);
      using acrd_t = K_FLD::ArrayAccessor<K_FLD::FloatArray>;

      bool inward1 = (cut.oriented == -1);
      sub.reorient(false);

      E_Float ps_min(1. - EPSILON);
      E_Int contact;
      
      haPolyhedron<0> result;
      aPolyhedron<0> cutter(cut.cnt, cut.crd);

      acrd_t acrd1(sub.m_crd), acrd2(cutter.m_crd);

      E_Int ret = NUGA::INTERSECT::isolated_clip(acrd1, sub, inward1, acrd2, cutter, ps_min, RTOL, result, contact, true_clip, dbg);
      if (!result.empty())bits.push_back(result);
      return ret;
    }

    template <>
    inline E_Int isolated_clip<aPolyhedron<0>, aPolyhedron<0>>(aPolyhedron<0>& sub, aPolyhedron<0>& cut, NUGA::INTERSECT::eOPER oper, E_Float RTOL, std::vector<aPolyhedron<0>>& bits, bool& true_clip)
    {
      if (oper == NUGA::INTERSECT::INTERSECTION)
        cut.reorient(false);

      pg_smesh_t ecut(cut);
      return isolated_clip<aPolyhedron<0>, pg_smesh_t>(sub, ecut, oper, RTOL, bits, true_clip);
    }


    template <typename aELT1, typename aELT2>
    bool compute(aELT1 & subj, aELT2 & cutter, NUGA::INTERSECT::eOPER oper, std::vector<aELT1>& bits)
    {
      bool true_clip(false);
      double RTOL = 1.e-9;
      E_Int err(1);
      for (size_t i = 0; i < 8 && err; ++i)
      {
        bits.clear();
        true_clip = false;
        err = isolated_clip(subj, cutter, oper, RTOL, bits, true_clip);
        RTOL *= 10;
      }

      if (err || !true_clip)
      {
#ifdef DEBUG_CLIPPER
        if (err) std::cout << "error clipping => just IO test" << std::endl;
        else if (!true_clip) std::cout << "singular config (contact..) => just IO test" << std::endl;
#endif
        bits.clear();
        return false; // to classify with a I/O test (even in case of error => roughly good approx)
      }
 
      return true;// true_clip;
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
        if (b < (E_Int)bits.size() - 1) // put the last (if exist) in the current, pop back the last
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

      double vcur(v0), veps(v0*EPSILON);
      K_SEARCH::BBox3D bx0;
      bx0.compute(bits[0].m_crd);
      bx0.enlarge(RTOL);

      std::vector<E_Int> cands;
      std::vector<typename zmesh_t::aelt_t> tmpbits;

      size_t nmask = mask_ids.size();
      for (size_t m = 0; (m < nmask) && (vcur > veps); ++m)
      {
        E_Int mi = mask_ids[m];
        if (mask_bits[mi] == nullptr) continue; // fixme : this is walls mask that does not exist anymore : clearing it from mask_ids before loop ?
        const bound_mesh_t& mask_bit = *mask_bits[mi];
        const loc_t& mask_loc = *(mask_bit.get_localizer());

#ifdef DEBUG_CLIPPER
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

          if (cands.empty()) continue; //fixme : why this happens ?
          //assert(!cands.empty());

          // autonomous cutter front
          bound_mesh_t acut_front(mask_bit, cands, idx_start);

#ifdef DEBUG_CLIPPER
          medith::write("subj", ae1);// or medith::write<ngon_type>("subj", z_mesh.crd, z_mesh.cnt, i);
          medith::write("cutter_front", acut_front.crd, acut_front.cnt);
#endif
          //CLIPPING
          tmpbits.clear();
          if (!just_io)
            just_io = !NUGA::CLIP::compute(ae1, acut_front, NUGA::INTERSECT::DIFFERENCE, tmpbits); //robust_clip returns true if true clip

          // IO : current bit does not intersect front (or only along its boundaries)
          if (just_io)
          {
            assert(tmpbits.empty());
            NUGA::eClassify wher = NUGA::CLASSIFY::classify(ae1, acut_front, true/*unambiguous*/);

#ifdef DEBUG_CLIPPER
            medith::write("cutter_front", acut_front.crd, acut_front.cnt);
            std::cout << "where OUT ? " << (wher == OUT) << std::endl;
#endif
            if (wher == OUT)
              continue; //will be considered with next front
          }

          // COMPRESS STRATEGY (with MOVE SEMANTICS) for bits => b can be reduced of 1 to treat the replaced at next iter 
          // if tmpbits is empty (IN) => compress (erase bits if single, put the last instead of current otherwise)
          // else replace the first bit, append the other bits  . 
          __replace_append_and_next_iter(bits, b, tmpbits);
          if ((E_Int)bits.size() < nbits) nbits = bits.size(); //update iff the last have replaced the current 
        }

        // update vcur
        vcur = 0.;
        for (size_t u = 0; u < bits.size(); ++u) vcur += bits[u].metrics();

      }
    }
  } // CLIP
} // NUGA

#endif /* NUGA_INTERSECT_HXX */

