

/*
 
 
 
              NUGA 
 
 
 
 */
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_CLIPPER_HXX
#define NUGA_CLIPPER_HXX

#include <vector>
#include "Nuga/Boolean/BAR_Conformizer.h"
#include "Nuga/Delaunay/T3Mesher.h"
namespace NUGA
{
  namespace CLIP
  {
    using crd_t = K_FLD::FloatArray;
    using cnt_t = K_FLD::IntArray;

    E_Int isolated_clip(const crd_t& crd1, const crd_t& crd2, const cnt_t& cutter, E_Float RTOL, std::vector<crd_t>& res, bool& true_clip)
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
      E_Int nb_pts1 = crd.cols();
      crd.pushBack(crd2);
      K_FLD::IntArray cnt(subj), cnt2(cutter);

      for (size_t i = 0; i < cnt2.cols(); ++i) std::swap(cnt2(0, i), cnt2(1, i));

      cnt2.shift(nb_pts1);
      cnt.pushBack(cnt2);

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
        medith::write("D:\\slandier\\DATA\\dev\\cases\\in.mesh", tcrd, cnt, "BAR");
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

#ifdef DEBUG_CLIPPER
      {
        data.pos->resize(3, data.pos->cols(), 0.);
        medith::write("D:\\slandier\\DATA\\dev\\cases\\toto.mesh", *data.pos, data.connectM, "TRI");
      }
#endif

      int minc = *std::min_element(ALL(data.colors));
      int maxc = *std::max_element(ALL(data.colors));

      if (minc != maxc) // non-connex
      {
        std::map<int, K_FLD::IntArray> col_to_cntB;
        for (size_t i = 0; i < data.connectM.cols(); ++i)
        {
          int coli = data.colors[i];
          K_FLD::IntArray::const_iterator pS = data.connectM.col(i);
          for (size_t n = 0; n < 3; ++n)
          {
            if (data.neighbors(n, i) == E_IDX_NONE)
            {
              int E[] = { *(pS + n),*(pS + (n + 1) % 3) };
              col_to_cntB[coli].pushBack(E, E + 2);
            }
          }
        }

        std::vector<int> PGi;
        res.resize(col_to_cntB.size());
        int i = 0;
        for (auto it = col_to_cntB.begin(); it != col_to_cntB.end(); ++it, ++i)
        {
          //
          K_FLD::FloatArray crd(*data.pos);
          std::vector<E_Int> nids;
          K_CONNECT::MeshTool::compact_to_mesh(crd, it->second, nids);

          PGi.clear();
          // sort the nodes
          BARSplitter::getSortedNodes(it->second, PGi);

          for (size_t j = 0; j < PGi.size(); ++j)
            res[i].pushBack(crd.col(PGi[j]), crd.col(PGi[j]) + crd.rows());
        }
      }
      else //connex
      {

        K_FLD::IntArray cB;
        K_CONNECT::MeshTool::getBoundaryT3Mesh(data.connectM, data.neighbors, cB);

        K_FLD::FloatArray crd(*data.pos);
        std::vector<E_Int> nids;
        K_CONNECT::MeshTool::compact_to_mesh(crd, cB, nids);

        std::vector<E_Int> PGi;
        // sort the nodes
        BARSplitter::getSortedNodes(cB, PGi);

        res.resize(1);
        for (size_t j = 0; j < PGi.size(); ++j)
          res[0].pushBack(crd.col(PGi[j]), crd.col(PGi[j]) + crd.rows());
      }

      return err;
    }

    template <typename aELT1, typename aELT2>
    bool compute(aELT1 const & subj, aELT2 const & cutter, E_Float RTOL, std::vector<aELT1>& bits)
    {
      bool true_clip(false);
      std::vector<crd_t> res;
      isolated_clip(subj.m_crd, cutter.crd, cutter.cnt, RTOL, res, true_clip);

      bits.clear();
      bits.reserve(res.size());
      //for (size_t i = 0; i < res.size(); ++i)
        //bits.push_back(std::move(res[i]));
      std::move(std::begin(res), std::end(res), std::back_inserter(bits));
 
      return true_clip;
    }
  } // CLIP
} // NUGA

#endif /* NUGA_INTERSECT_HXX */

