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
//Authors : Sam Landier (sam.landier@onera.fr)

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

    ///////////////////// private functions ////////////////////////////////////////////////////////////////////////////////////////////////

    inline void __filtrate_opposite_edges(K_FLD::IntArray& cnt, std::vector<bool>& keep)
    {
      std::map<K_MESH::Edge, std::vector<int>> w_E2ids;
      bool do_compress = false;
      //
      for (E_Int e = 0; e < cnt.cols(); ++e)
      {
        auto ee = K_MESH::Edge(cnt(0, e), cnt(1, e));
        auto oee = K_MESH::Edge(cnt(1, e), cnt(0, e));

        auto it1 = w_E2ids.find(ee);
        auto it = w_E2ids.find(oee);
        if (it != w_E2ids.end())
          w_E2ids[oee].push_back(e);
        if (it == w_E2ids.end() && it1 == w_E2ids.end())
          w_E2ids[ee].push_back(e);

      }

      keep.resize(cnt.cols(), true);

      for (auto ii : w_E2ids)
      {
        if (ii.second.size() == 1) continue;
        do_compress = true;
        for (size_t e = 0; e < ii.second.size(); ++e)
          keep[ii.second[e]] = false;
      }

      if (do_compress)
      {
        K_CONNECT::keep<bool> pred(keep);
        K_CONNECT::IdTool::compress(cnt, pred);
      }
    }

    template <typename aELT>
    inline void __filtrate_outside_edges(aELT& sub, const K_FLD::FloatArray& crd2D, K_FLD::IntArray& cnt, int id_start, std::vector<std::pair<E_Int,E_Int>>& xedge, std::vector<bool>& keep)
    {
      // WARNING : sub cnt ids must refer to crd2D

      keep.resize(cnt.cols(), true);

      std::vector<bool> w_keep(crd2D.cols(), true);
      DELAUNAY::Triangulator dt;
      bool is_in = true;
      //
      for (E_Int k = id_start; k < crd2D.cols(); ++k)
      {
        if (size_t(k) < xedge.size() && xedge[k].first != IDX_NONE)
        {
          w_keep[k] = true;
          continue; // X node => lying on (impacted) sub : no test required (which is good it as it might fail)
        }

        const double * P = crd2D.col(k);
        sub.template fast_is_in_pred<DELAUNAY::Triangulator, 2>(dt, crd2D, P, is_in, 1.e-9);
        w_keep[k] = is_in;
      }

      bool do_compress = false;
      keep.resize(cnt.cols(), true);
      //
      for (E_Int k = 0; k < cnt.cols(); ++k)
      {
        auto pK = cnt.col(k);
        keep[k] = (w_keep[*pK] && w_keep[*(pK + 1)]);
        do_compress |= !keep[k];
      }

      if (do_compress)
      {
        K_CONNECT::keep<bool> pred(keep);
        K_CONNECT::IdTool::compress(cnt, pred);
      }
    }

    template <typename aELT>
    inline void __filtrate_outside_edges2(aELT& cut, const K_FLD::FloatArray& crd2D, K_FLD::IntArray& cnt, int id_end, std::vector<std::pair<int, int>>& xedge, std::vector<bool>& keep)
    {
      // WARNING : sub cnt ids must refer to crd2D

      keep.resize(cnt.cols(), true);

      std::vector<bool> w_keep(crd2D.cols(), true);
      DELAUNAY::Triangulator dt;
      bool is_in = true;
      //
      for (size_t k = 0; k < id_end; ++k)
      {
        if (k < xedge.size() && xedge[k].first != IDX_NONE)
        {
          w_keep[k] = true;
          continue; // X node => lying on (impacted) sub : no test required (which is good it as it might fail)
        }

        const double * P = crd2D.col(k);
        cut.template fast_is_in_pred<DELAUNAY::Triangulator, 2>(dt, crd2D, P, is_in, 1.e-9);
        w_keep[k] = is_in;
      }

      bool do_compress = false;
      keep.resize(cnt.cols(), true);
      //
      for (size_t k = 0; k < cnt.cols(); ++k)
      {
        auto pK = cnt.col(k);
        keep[k] = (w_keep[*pK] && w_keep[*(pK + 1)]);
        do_compress |= !keep[k];
      }

      if (do_compress)
      {
        K_CONNECT::keep<bool> pred(keep);
        K_CONNECT::IdTool::compress(cnt, pred);
      }
    }

    template <short DIM>
    void __normalize(K_FLD::FloatArray& crd)
    {
      //try again with a normalized contour
      K_SEARCH::BoundingBox<DIM> box;
      box.compute(crd);

      double dX[DIM];
      for (size_t k = 0; k < DIM; ++k)
        dX[k] = box.maxB[k] - box.minB[k];

      for (int u = 0; u < crd.cols(); ++u)
      {
        for (size_t k = 0; k < DIM; ++k)
          if (dX[k] != 0.) crd(k, u) = (crd(k, u) - box.minB[k]) / dX[k];
      }
    }

    inline void go_2D(bool proj_on_first, const E_Float* W, K_FLD::FloatArray& crd, std::vector<double>& zs, double& zmean, K_FLD::FloatArray& P, const K_FLD::IntArray& cnt, E_Int nb_pts1, E_Int nb_edges1)
    {      
      K_FLD::FloatArray iP(3, 3);

      NUGA::computeAFrame(W, P);
      iP = P;
      K_FLD::FloatArray::inverse3(iP);
      NUGA::transform(crd, iP);

      crd.extract_field(2, zs); //keep 3rd coord appart (altitude)
      
      //compute zmean
      zmean = 0.;
      if (proj_on_first) {
        for (E_Int k = 0; k < nb_pts1; ++k) zmean += zs[k];
        zmean /= nb_pts1;
      }
      else {
        for (E_Int k = nb_pts1; k < crd.cols(); ++k) zmean += zs[k];
        zmean /= (crd.cols() - nb_pts1);
      }

      if (proj_on_first)  // now apply zmean to front points such remaining ones at the end will be roughly on subj supporting surface
        for (size_t k = nb_pts1; k < zs.size(); ++k) zs[k] = zmean;
      else
        for (E_Int k = 0; k < nb_pts1; ++k) zs[k] = zmean;

      crd.resize(2, crd.cols());//now pure 2D
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template <typename aELT1, typename aELT2> //sub & cut are non const as they can be reoriented/reversed
    inline E_Int isolated_clip(aELT1& sub, aELT2& cut, NUGA::INTERSECT::eOPER oper, E_Float ARTOL, std::vector<aELT1>& res, bool& true_clip);
    
    ///
    template <>
    inline E_Int isolated_clip<aPolygon, edge_mesh_t>(aPolygon& sub, edge_mesh_t& cut, NUGA::INTERSECT::eOPER oper, E_Float ARTOL, std::vector<aPolygon>& bits, bool& true_clip)
    {
      using cnt_t = K_FLD::IntArray;
      const crd_t& crd1 = sub.m_crd;
      const crd_t& crd2 = cut.crd;
      const cnt_t& cutter = cut.cnt;

      true_clip = false;
      bits.clear();

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

#ifdef DEBUG_CLIPPER
      {
        medith::write("before_conf3D", crd, cnt, "BAR");
      }
#endif

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
      
      // DISCARD FALSE OVERLAPS among fronts (now we are in 2D, those with big altitudes)
      // using meL : NOT WORKING because after first cut, sub is overdefined so underestimate Lref, hence some good cut edges are discarded
      // using bbox in 3D frame : NOT WORKING neither when it's a planar case where the plane is axi-aligned => one corrdinate is the same => Lref gets null
      // using bbox in 2D frame
      K_SEARCH::BBox2D bx(crd, nb_pts1);
      std::vector<E_Int> new_edge_ids; //in case of compacting, need to propagate original ids in xedge
      double Lref = std::min(bx.maxB[0] - bx.minB[0], bx.maxB[1] - bx.minB[1]);

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
          E_Int nb_edgesi = cnt.cols();
          K_FLD::IntArray::compact(cnt, keep, new_edge_ids);
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
        __normalize<2>(tcrd);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("before_conf2D", tcrd, cnt, "BAR");
      }
#endif

      double ABSTOL2D = ARTOL;
      if (ARTOL < 0.) //relative
      {
        //2D ABS TOLERANCE
        double Lref2 = NUGA::FLOAT_MAX;
        for (int k = 0; k < cnt.cols(); ++k)
        {
          double d2 = NUGA::sqrDistance(crd.col(cnt(0, k)), crd.col(cnt(1, k)), 3);
          Lref2 = std::min(d2, Lref2);
        }

        ABSTOL2D = -::sqrt(Lref2) * ARTOL;
        ABSTOL2D = std::max(ABSTOL2D, ZERO_M);
      }

      //
      E_Int nb_edges0 = cnt.cols();
      err = conformizer.run(crd, cnt, ancE2, nullptr/*&priority*/, ABSTOL2D, nb_edges1, 1 /*one iter only*/);
      if (err) return err;

#ifdef DEBUG_CLIPPER
      {
        K_FLD::FloatArray tcrd(crd);
        __normalize<2>(tcrd);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("after_conf2D", tcrd, cnt, "BAR");
      }
#endif

      auto xedge = conformizer.get_x_history();

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

#ifdef DEBUG_CLIPPER
      {
        K_FLD::FloatArray tcrd(*data.pos);
        __normalize<2>(tcrd);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("out", tcrd, data.connectM, "TRI");
      }
#endif

      true_clip = true;

      // transfer altitude
      std::vector<double> zsnew(data.pos->cols(), NUGA::FLOAT_MAX);
      const auto& pnids = conformizer.get_node_history();
      for (size_t i = 0; i < pnids.size(); ++i)
        if (pnids[i] != IDX_NONE) zsnew[pnids[i]] = zs[i];
      zs = std::move(zsnew);

      // interpolate missings using x history
      
      for (size_t i = 0; i < zs.size(); ++i)
      {
        if (zs[i] != NUGA::FLOAT_MAX) continue;            // non valuated node
        if (i < pnids.size() && pnids[i] != E_Int(i) && pnids[i] != IDX_NONE) continue;   // merged node

        auto xe = xedge[i];
        
        assert(xe.first != IDX_NONE); // i is an x point and has its histo
        
        if (xe.first >= subj.cols()) continue;

        int N0 = subj(0, xe.first);
        double * P0 = crd.col(N0);
        int N1 = subj(1, xe.first);
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

      std::vector<crd_t> crd_res;              //intermediate format as crds
      int minc = *std::min_element(ALL(data.colors));
      int maxc = *std::max_element(ALL(data.colors));

      std::vector<int> PGi;
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
              E_Int E[] = { *(pS + (n + 1) % 3),*(pS + (n + 2) % 3) };
              col_to_cntB[coli].pushBack(E, E + 2);
            }
          }
        }

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
        crd_res.resize(1);
        K_FLD::IntArray cB;
        NUGA::MeshTool::getBoundaryT3Mesh(data.connectM, data.neighbors, cB);

        std::vector<E_Int> nids;
        K_FLD::FloatArray& crd1 = *data.pos;
        NUGA::MeshTool::compact_to_mesh(crd1, cB, nids);

        PGi.clear();
        // sort the nodes
        BARSplitter::getSortedNodes(cB, PGi);

        E_Int str = crd1.rows();
        crd_res[0].reserve(str, PGi.size());
 
        for (size_t j = 0; j < PGi.size(); ++j)
          crd_res[0].pushBack(crd1.col(PGi[j]), crd1.col(PGi[j]) + str);
      }

      int nbits = crd_res.size();
      bits.reserve(nbits);
      std::move(ALL(crd_res), std::back_inserter(bits));

      return err;
    }

    ///
    template <>
    inline E_Int isolated_clip<aPolygon, aPolygon>(aPolygon& sub, aPolygon& cut, NUGA::INTERSECT::eOPER oper, E_Float ARTOL, std::vector<aPolygon>& bits, bool& true_clip)
    {
      using cnt_t         = K_FLD::IntArray;
      const crd_t& crd1   = sub.m_crd;
      const crd_t& crd2   = cut.m_crd;
      edge_mesh_t e_cut(cut); 
      const cnt_t& cutter = e_cut.cnt;

      // always project on sub unless cut normal is computed before getting here
      bool proj_on_first = (cut.m_normal[0] == NUGA::FLOAT_MAX);

      true_clip = false;
      bits.clear();

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
        
      E_Int nb_pts1 = crd.cols();
      crd.pushBack(crd2);
      K_FLD::IntArray cnt(subj), cnt2(cutter);

      if (oper == NUGA::INTERSECT::DIFFERENCE)
        connect_trait<LINEIC, true>::reverse_orient(cnt2);

      cnt2.shift(nb_pts1);
      cnt.pushBack(cnt2);

#ifdef DEBUG_CLIPPER
      {
        medith::write("before_conf3D", crd, cnt, "BAR");
      }
#endif

      // project one on the other
      const E_Float* W = proj_on_first ? sub.get_normal() : cut.get_normal();

      std::vector<double> zs;
      double zmean;
      K_FLD::FloatArray P(3, 3);
      go_2D(proj_on_first, W, crd, zs, zmean, P, cnt, nb_pts1, nb_edges1);

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
        __normalize<2>(tcrd);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("before_conf2D", tcrd, cnt, "BAR");
      }
#endif

      double ABSTOL2D = ARTOL;
      if (ARTOL < 0.) //relative
      {
        //2D ABS TOLERANCE
        double Lref2 = NUGA::FLOAT_MAX;
        for (int k = 0; k < cnt.cols(); ++k)
        {
          double d2 = NUGA::sqrDistance(crd.col(cnt(0, k)), crd.col(cnt(1, k)), 3);
          Lref2 = std::min(d2, Lref2);
        }

        ABSTOL2D = -::sqrt(Lref2) * ARTOL;
        ABSTOL2D = std::max(ABSTOL2D, ZERO_M);
      }

      //
      err = conformizer.run(crd, cnt, ancE2, nullptr/*&priority*/, ABSTOL2D, -nb_nodes, 1 /*one iter only*/);
      if (err)
        return err;

#ifdef DEBUG_CLIPPER
      {
        K_FLD::FloatArray tcrd(crd);
        __normalize<2>(tcrd);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("after_conf2D", tcrd, cnt, "BAR");
      }
#endif

      auto xedge = conformizer.get_x_history();

      // Remove identical edges
      {
        std::vector<E_Int> dupIds;
        int nb_removed = MeshTool::removeDuplicated(cnt, dupIds, true/*stict orient*/);
        if (nb_removed)
        {

          std::vector<E_Int> new_ancE2;
          new_ancE2.reserve(ancE2.size());
          for (E_Int k = 0; k < E_Int(dupIds.size()); ++k)
          {
            if (dupIds[k] == k)
              new_ancE2.push_back(ancE2[k]);
          }
          ancE2 = new_ancE2;
        }
      }

      // FILTERS
      std::vector<bool> keep;
      
      // INTERECTION & DIFFERENCE
      {
        // Filter #1 : remove oppposite edges 
        __filtrate_opposite_edges(cnt, keep);

        {
          K_CONNECT::keep<bool> pred1(keep);
          K_CONNECT::IdTool::compress(ancE2, pred1);
        }

        if (cnt.cols() == 0)
        {
          // Means that sub & cut are exactly matching, with opposite orientation
          assert(oper != NUGA::INTERSECT::INTERSECTION); // bacause sub & cut ahve same orientation, so we cannot have ALL edges with opposite orientation
          assert(oper == NUGA::INTERSECT::DIFFERENCE); // this function is currently only tested for INTERSECTION/DIFFERNCE cases only.
          // Empty answer (solved topologically)
          sub.clear();
          return 0;
        }

#ifdef DEBUG_CLIPPER
        {
          K_FLD::FloatArray tcrd(crd);
          __normalize<2>(tcrd);
          tcrd.resize(3, tcrd.cols(), 0.);
          medith::write("after_conf2D_filter1", tcrd, cnt, "BAR");
        }
#endif
      }

      // INTERECTION & DIFFERENCE
      {
        // Filter #2 : burn edges oustide subj
        __filtrate_outside_edges(sub, crd, cnt, nb_pts1, xedge, keep);

        {
          K_CONNECT::keep<bool> pred1(keep);
          K_CONNECT::IdTool::compress(ancE2, pred1);
        }

#ifdef DEBUG_CLIPPER
        {
          K_FLD::FloatArray tcrd(crd);
          __normalize<2>(tcrd);
          tcrd.resize(3, tcrd.cols(), 0.);
          medith::write("after_conf2D_filter2", tcrd, cnt, "BAR");
        }
#endif
      }
      
      /*if (oper == NUGA::INTERSECT::INTERSECTION)
      {
        // Filter #3 : burn edges oustide cut
        keep.clear();
        K_CONNECT::IdTool::shift(cut.m_nodes, nb_pts1);
        __filtrate_outside_edges2(cut, crd, cnt, nb_pts1, xedge, keep);

        {
          K_CONNECT::keep<bool> pred1(keep);
          K_CONNECT::IdTool::compress(ancE2, pred1);
        }

#ifdef DEBUG_CLIPPER
        {
          K_FLD::FloatArray tcrd(crd);
          __normalize<2>(tcrd);
          tcrd.resize(3, tcrd.cols(), 0.);
          medith::write("after_conf2D_filter3", tcrd, cnt, "BAR");
        }
#endif
        if (cnt.cols() == 0)
        {
          // Means that sub & cut are not overlapping
          true_clip = true;
          return 0;
        }
      }*/

      //
      {
        // Filter #4 : burn edges oustide sub
        std::map<E_Int, E_Int> node_to_count;
        NUGA::MeshTool::build_node_arity(cnt, node_to_count);
        NUGA::MeshTool::burn_free_branches(cnt, node_to_count, keep);

        {
          K_CONNECT::keep<bool> pred1(keep);
          K_CONNECT::IdTool::compress(ancE2, pred1);
        }

#ifdef DEBUG_CLIPPER
        {
          K_FLD::FloatArray tcrd(crd);
          __normalize<2>(tcrd);
          tcrd.resize(3, tcrd.cols(), 0.);
          medith::write("after_conf2D_filter4", tcrd, cnt, "BAR");
        }
#endif
      }

      //assert(cnt.cols() != 0);
      if (cnt.cols() == 0)
      {
        // Empty answer
        sub.clear();
        return 0;
      }

      int maxAnc = *std::max_element(ALL(ancE2));
      if (maxAnc < nb_nodes)
      {
        return 0; // "no" intersections => fully visible or hidden ;
      }

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
      if (err != 0 && err != 77/*open contour*/)
      {
        return 1;
      }

      if (data.connectM.cols() == 0) //checkme
      {
        if (oper == NUGA::INTERSECT::DIFFERENCE)
          sub.clear();
        //checme : is it always meanin OUT for INTERSECTION ?
        return 0; 
      }

      // check if the result is just impacted subj
      {
        int maxAnc = -1;
        std::map<K_MESH::NO_Edge, int> edge_to_anc;
        for (int e = 0; e < cnt.cols(); ++e)
          edge_to_anc[K_MESH::NO_Edge(cnt(0, e), cnt(1, e))] = ancE2[e];

        K_FLD::IntArray cB;
        NUGA::MeshTool::getBoundaryT3Mesh(data.connectM, data.neighbors, cB);

        for (int e = 0; e < cB.cols(); ++e)
        {
          auto pK = cB.col(e);
          K_MESH::NO_Edge E(*pK, *(pK + 1));
          auto it = edge_to_anc.find(E);
          assert(it != edge_to_anc.end());
          //if (it == edge_to_anc.end()) continue; //checkme : can happen ?
          maxAnc = std::max(maxAnc, it->second);
        }
        assert(maxAnc != -1);
        if (maxAnc < nb_nodes) return 0; // impacted subj => IO test
      }

#ifdef DEBUG_CLIPPER
      {
        K_FLD::FloatArray tcrd(*data.pos);
        __normalize<2>(tcrd);
        tcrd.resize(3, tcrd.cols(), 0.);
        medith::write("out", tcrd, data.connectM, "TRI");
      }
#endif

      true_clip = true;

      // transfer altitude
      std::vector<double> zsnew(data.pos->cols(), NUGA::FLOAT_MAX);
      const auto& pnids = conformizer.get_node_history();
      for (size_t i = 0; i < pnids.size(); ++i)
        if (pnids[i] != IDX_NONE) zsnew[pnids[i]] = zs[i];
      zs = std::move(zsnew);

      // interpolate missings using x history
      
      for (size_t i = 0; i < zs.size(); ++i)
      {
        if (zs[i] != NUGA::FLOAT_MAX) continue;            // non valuated node
        if (i < pnids.size() && pnids[i] != E_Int(i) && pnids[i] != IDX_NONE) continue;   // merged node

        auto xe = xedge[i];
        
        assert(xe.first != IDX_NONE); // i is an x point and has its histo
        
        if (xe.first >= subj.cols()) continue;

        int N0 = subj(0, xe.first);
        double * P0 = crd.col(N0);
        int N1 = subj(1, xe.first);
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

      // build history
      std::vector<E_Int> poids;
      K_CONNECT::IdTool::reverse_indirection(pnids, poids);
      poids.resize(data.pos->cols(), IDX_NONE);

      // rule : 
      // sub id =>
      // cut id =>
      // X node => SZUDOR code of both xe (intersecting edges)
      std::vector<long> l_m_poids(poids.size(), IDX_NONE);
      for (size_t k = 0; k < poids.size(); ++k)
      {
        if (poids[k] < nb_nodes)       // sub
          l_m_poids[k] = sub.m_poids[poids[k]];
        else if (poids[k] != IDX_NONE) // cutter
          l_m_poids[k] = cut.m_poids[poids[k]- nb_nodes];
        else if (k < xedge.size())
        {
          int e0 = xedge[k].first;
          int e1 = xedge[k].second;

          if (e0 == IDX_NONE || e1 == IDX_NONE) continue; // needs both ids 

          e1 -= nb_nodes;

          l_m_poids[k] = -(NUGA::szudzik_pairing(e0, e1) + 1);
        }
      }

      std::vector<crd_t> crd_res;              //intermediate format as crds
      std::vector<std::vector<long>> poids_res; // to migrate nodes history

      int minc = *std::min_element(ALL(data.colors));
      int maxc = *std::max_element(ALL(data.colors));

      std::vector<int> PGi;
      std::vector<long> new_poids;

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
              E_Int E[] = { *(pS + (n + 1) % 3),*(pS + (n + 2) % 3) };
              col_to_cntB[coli].pushBack(E, E + 2);
            }
          }
        }

        crd_res.resize(col_to_cntB.size());
        poids_res.resize(col_to_cntB.size());
        
        int i = 0;
        std::vector<E_Int> nids;
        for (auto it = col_to_cntB.begin(); it != col_to_cntB.end(); ++it, ++i)
        {
          //
          nids.clear();
          K_FLD::FloatArray crd1(*data.pos); //should not hurt as meshes are small
          NUGA::MeshTool::compact_to_mesh(crd1, it->second, nids);

          new_poids.clear();
          new_poids.resize(crd1.cols(), IDX_NONE);
          for (size_t k = 0; k < l_m_poids.size(); ++k)
          {
            if (nids[k] == IDX_NONE) continue;
            new_poids[nids[k]] = l_m_poids[k];
          }

          PGi.clear();
          // sort the nodes
          BARSplitter::getSortedNodes(it->second, PGi);

          E_Int str = crd1.rows();
          crd_res[i].reserve(str, PGi.size());
          poids_res[i].resize(crd1.cols(), IDX_NONE);
          for (size_t j = 0; j < PGi.size(); ++j)
          {
            crd_res[i].pushBack(crd1.col(PGi[j]), crd1.col(PGi[j]) + str);
            poids_res[i][j] = new_poids[PGi[j]];
          }
        }
      }
      else //connex
      {
        crd_res.resize(1);
        poids_res.resize(1);

        K_FLD::IntArray cB;
        NUGA::MeshTool::getBoundaryT3Mesh(data.connectM, data.neighbors, cB);

        std::vector<E_Int> nids;
        K_FLD::FloatArray& crd1 = *data.pos;
        NUGA::MeshTool::compact_to_mesh(crd1, cB, nids);

        new_poids.clear();
        new_poids.resize(crd1.cols(), IDX_NONE);
        for (size_t k = 0; k < l_m_poids.size(); ++k)
        {
          if (nids[k] == IDX_NONE) continue;
          new_poids[nids[k]] = l_m_poids[k];
        }

        PGi.clear();
        // sort the nodes
        BARSplitter::getSortedNodes(cB, PGi);

        E_Int str = crd1.rows();

        crd_res[0].reserve(str, PGi.size());
        poids_res[0].resize(crd1.cols(), IDX_NONE);

        for (size_t j = 0; j < PGi.size(); ++j)
        {
          crd_res[0].pushBack(crd1.col(PGi[j]), crd1.col(PGi[j]) + str);
          poids_res[0][j] = new_poids[PGi[j]];
        }
      }

      size_t nbits = crd_res.size();
      bits.reserve(nbits);
      std::move(ALL(crd_res), std::back_inserter(bits));
      // pass the history
      for (size_t k = 0; k < nbits; ++k)
      {
        assert(poids_res[k].size() == (size_t)bits[k].m_crd.cols());
        bits[k].m_poids = poids_res[k];
      }

      return err;
    }

    
    ///
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
      double RTOL = -1.e-9;
      E_Int err(1);
      for (size_t i = 0; i < 8 && err; ++i)
      {
        bits.clear();
        true_clip = false;
        err = isolated_clip(subj, cutter, oper, RTOL, bits, true_clip);
        RTOL *= 10;
      }

      // singularities/errors management

      if (err || !true_clip)
      {
#ifdef DEBUG_CLIPPER
        if (err) std::cout << "error clipping => just IO test" << std::endl;
        else if (!true_clip) std::cout << "singular config (contact..) => just IO test" << std::endl;
#endif
        bits.clear();
        return false; // to classify with a I/O test (even in case of error => roughly good approx)
      }
      else if (true_clip && oper == NUGA::INTERSECT::INTERSECTION) //fixme
      {
        // singularity I : 2 faces side by side with one overlapping edge => should be empty intersection
        if (bits.size() == 2)// singularity I is equivalent to that test IFF both faces are convex : cannot have 2 pieces. Real test should use surface : sum (surf(bits)) < min(surf(subj), surf(cutter))
        {
          bits.clear(); 
          return false;
        }
      }
 
      return true;// true_clip;
    }

    template<typename aelt_t>
    void __replace_append_and_next_iter(std::vector<aelt_t>& bits, int& b, std::vector<aelt_t>& toadd)
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
        if (b < (int)bits.size() - 1) // put the last (if exist) in the current, pop back the last
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
        int idx_start = mask_bit.index_start;
        //
        int nbits = (int)bits.size(); // bits can be appended when split give multiple bits
        for (int b = 0; b < nbits; ++b)
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
          //medith::write("subj", ae1);// or medith::write<ngon_type>("subj", z_mesh.crd, z_mesh.cnt, i);
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
          if ((int)bits.size() < nbits) nbits = bits.size(); //update iff the last have replaced the current 
        }

        // update vcur
        vcur = 0.;
        for (size_t u = 0; u < bits.size(); ++u) vcur += bits[u].metrics();

      }
    }
  } // CLIP
} // NUGA

#endif /* NUGA_INTERSECT_HXX */

