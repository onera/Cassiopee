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

#ifndef NUGA_DISPLACEMENT_HXX
#define NUGA_DISPLACEMENT_HXX

#include "Nuga/include/vertex.h"
#include "Nuga/include/mesh_t.hxx"
#include "Nuga/include/collider.hxx"

//#define DISPLACEMENT_DBG

#ifdef DISPLACEMENT_DBG
#include "Nuga/include/medit.hxx"
#endif

namespace NUGA
{

  template <typename mesh_t1, typename mesh_t2>
  static std::vector<direction> immerse_nodes(mesh_t1& m1, const mesh_t2& m2, double ARTOL)
  {
    eMetricType mtype = NUGA::ISO_MIN;

    if (m2.ARG1 == SURFACIC)
      assert(m2.oriented == 1);

    // get vertices from mesh1
    std::vector<NUGA::vertex> vertices;

    auto sqrmetric1 = m1.get_nodal_metric2(mtype);

    vertices.reserve(m1.crd.cols());
    for (E_Int i = 0; i < m1.crd.cols(); ++i)
      vertices.push_back(NUGA::vertex(m1.crd.col(i), sqrmetric1[i]));

    // identify singular points
    std::vector<std::vector<E_Int>> pt_to_faces;
    NUGA::COLLIDE::get_colliding(vertices, m2, ARTOL, mtype, pt_to_faces);

    //std::cout << "compute displacements directions" << std::endl;
    std::vector<NUGA::direction> dirs;
    dirs.resize(vertices.size());
    // default : using faces normal
    for (size_t i = 0; i < pt_to_faces.size(); ++i)
    {
      dirs[i].vec[0] = dirs[i].vec[1] = dirs[i].vec[2] = 0.;

      if (pt_to_faces[i].empty()) continue; // regular point

      E_Int nbpgs = pt_to_faces[i].size();

      if (nbpgs < 3) // 1 or 2 : PG normal or average of the 2
      {
        for (size_t p = 0; p < pt_to_faces[i].size(); ++p)
        {
          const E_Int& fid = pt_to_faces[i][p];
          auto PG = m2.element(fid);
          E_Float nn[3];
          PG.template ndS<3>(m2.crd, nn);
          NUGA::normalize<3>(nn);

  #ifdef DEBUG_IMMERSION
          E_Float l2 = ::sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
          assert(::fabs(l2 - 1.) < EPSILON);// DEGEN
  #endif

          // immersion => inward => -=
          dirs[i].vec[0] -= nn[0];
          dirs[i].vec[1] -= nn[1];
          dirs[i].vec[2] -= nn[2];
        }
      }
      else // this vertex is near a m2's node => compute nodal normal : weighted with angle between each pair of ray in the node shell (sum is 2PI)
      {
        //std::cout << "compute node normal" << std::endl;
        double TOL2 = ARTOL*ARTOL;
        if (ARTOL < 0.) //relative
          TOL2 *= vertices[i].val2;

        // find near node : 
        E_Int PGi0 = pt_to_faces[i][0];
        auto PG0 = m2.element(PGi0);
        double dmin2= FLOAT_MAX;
        E_Int Ni{IDX_NONE};
        for (size_t u=0; u < (size_t)PG0.nb_nodes(); ++u) 
        {
          double d2 = NUGA::sqrDistance(m2.crd.col(PG0.node(u)), vertices[i].vec, 3);
          if (d2 < dmin2)
          {
            dmin2 = d2;
            Ni = PG0.node(u);
          }
        }
        assert (dmin2 < TOL2);
        assert (Ni != IDX_NONE);

        //std::cout << "Near node : " << Ni << std::endl;

        E_Int count{0};

        //std::cout << "compute the rays angles" << std::endl;

        for (size_t p = 0; p < pt_to_faces[i].size(); ++p)
        {
          const E_Int& fid = pt_to_faces[i][p];
          auto PG = m2.element(fid);

          //std::cout << p << "-th face : " << fid << std::endl;

          E_Int nnodes = PG.nb_nodes();

          E_Int pos = K_CONNECT::IdTool::get_pos(PG.begin(), nnodes, Ni+1);
          assert (pos != -1);
          if (pos == -1) continue;
          ++count;

          //std::cout << "pos ? : " << pos << std::endl;

          E_Int Nip1 = PG.node((pos+1)%nnodes);
          E_Int Nim1 = PG.node((pos+nnodes-1)%nnodes);

          //std::cout << "before angular" << std::endl;

          E_Float nn[3];
          NUGA::angular_weighted_normal(m2.crd.col(Nim1), m2.crd.col(Ni), m2.crd.col(Nip1), nn);

          //std::cout << "after angular" << std::endl;

          // immersion => inward => -=
          dirs[i].vec[0] -= nn[0];
          dirs[i].vec[1] -= nn[1];
          dirs[i].vec[2] -= nn[2];

        }

        if ((size_t)count != pt_to_faces[i].size())
          std::cout << "WRONG LOGIC" << std::endl;
      }

      // normalize
      NUGA::normalize<3>(dirs[i].vec);
    }

    // compute displacements norms

    m2.get_nodal_metric2(mtype);

    for (size_t i = 0; i < vertices.size(); ++i)
    {
      dirs[i].val2 = 0.;

      if (pt_to_faces[i].empty()) continue; // regular point

      // compute max distance among sticking PGs to satisfy all of them
      for (size_t p = 0; p < pt_to_faces[i].size(); ++p)
      {
        const E_Int& fid = pt_to_faces[i][p];
        auto PG = m2.element(fid);

        E_Int P0 = PG.node(0);         // first face node
        double nPG[3];                 // normal to PG
        PG.template normal<K_FLD::FloatArray, 3>(m2.crd, nPG);

        const double* Pi = vertices[i].vec;

        double I[3]; // intersection of pt in computed dir on PG plane.
        double h = NUGA::project(m2.crd.col(P0), nPG, Pi, dirs[i].vec, I); // == PiP0.nPG
        
        if (h < 0.) continue; //already immersed

        double w = NUGA::dot<3>(nPG, dirs[i].vec);
        assert(w < 0.);
        w = 1. / w;

        double dmax = ::sqrt(std::min(vertices[i].val2, PG.Lref2(m2.nodal_metric2)));

        double TOLi = ARTOL;
        if (ARTOL < 0.) // relative
          TOLi *= -dmax;

        TOLi = std::min(0.95*dmax, TOLi);
        
        if (::fabs(h) >= TOLi) 
        {
          continue;//too far
        }

        double depth = 0.1*(dmax - TOLi);
        double lambdaMove = -w * (h + depth);

        // following assert because dir must be well oriented : inward the surface if above, outward otherwise
        // so attractive if above, repulsive if below (to put it far enough to exit interference zone)
        
        dirs[i].val2 = std::max(lambdaMove*lambdaMove, dirs[i].val2);
      }
    }

    // apply displacements to m1
#ifdef DEBUG_IMMERSION
    auto crd10 = m1.crd;
#endif
    E_Int nb_disp{0};
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      if (pt_to_faces[i].empty()) continue; // regular point

      if (::sqrt(dirs[i].val2) > EPSILON) ++nb_disp;

      NUGA::sum<3>(::sqrt(dirs[i].val2), dirs[i].vec, 1., m1.crd.col(i), m1.crd.col(i));
    }

    std::cout << "NB OF DISPLACED POINTS : " << nb_disp << std::endl;
#ifdef DEBUG_IMMERSION
    

    ngon_unit ngumoved;
    for (size_t i = 0; i < m1.cnt.PGs.size(); ++i)
    {
      const E_Int* nodes = m1.cnt.PGs.get_facets_ptr(i);
      E_Int nnodes = m1.cnt.PGs.stride(i);
      E_Int count{ 0 };
      for (size_t j = 0; j < nnodes; ++j)
      {
        E_Int Ni = nodes[j] - m1.index_start;

        if (!pt_to_faces[Ni].empty())
          ++count;
      }

      if (count != 0)
        ngumoved.add(nnodes, nodes);
    }

    ngumoved.updateFacets();
    medith::write("pgmoved", m1.crd, ngumoved);
    medith::write("pgtomove", crd10, ngumoved);
#endif

    // returns moved vertices and opposite displacements (for geometrically get back to original pos in boolean)
    std::vector<direction> odirs;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      if (pt_to_faces[i].empty()) continue; // regular point

      dirs[i].flag = i; //store node id
      // reverse for applying reverse displacement
      dirs[i].vec[0] = - dirs[i].vec[0];
      dirs[i].vec[1] = - dirs[i].vec[1];
      dirs[i].vec[2] = - dirs[i].vec[2];

      odirs.push_back(dirs[i]);
    }

    return odirs;
  }

  template <typename mesh_t1, typename mesh_t2>
  static std::vector<direction> compute_displacement_for_singular_nodes(mesh_t1& m1, const mesh_t2& m2, double ARTOL, eMetricType mtype, bool inward)
  {
    // inward = True  => immersion
    // inward = False => repulsion

    if (m2.ARG1 == SURFACIC)
      assert(m2.oriented == 1);

    // get vertices from mesh1
    std::vector<NUGA::vertex> vertices;

    auto sqrmetric1 = m1.get_nodal_metric2(/*mtype*/);

    vertices.reserve(m1.crd.cols());
    for (E_Int i = 0; i < m1.crd.cols(); ++i)
    {
      vertices.push_back(NUGA::vertex(m1.crd.col(i), sqrmetric1[i]));
    }

    // identify singular points
    std::vector<std::vector<E_Int>> pt_to_faces;
    NUGA::COLLIDE::get_colliding(vertices, m2, ARTOL, mtype, pt_to_faces);

    //std::cout << "compute displacements directions" << std::endl;
    std::vector<NUGA::direction> dirs;
    dirs.resize(vertices.size());
    for (size_t i = 0; i < dirs.size(); ++i)
    {
      dirs[i].val2 = dirs[i].vec[0] = dirs[i].vec[1] = dirs[i].vec[2] = 0.;
    }
    // default : using faces normal
    for (size_t i = 0; i < pt_to_faces.size(); ++i)
    {
      //dirs[i].vec[0] = dirs[i].vec[1] = dirs[i].vec[2] = 0.;

      if (pt_to_faces[i].empty()) continue; // regular point

      E_Int nbpgs = pt_to_faces[i].size();

      if (nbpgs < 3) // 1 or 2 : PG normal or average of the 2
      {
        for (size_t p = 0; p < pt_to_faces[i].size(); ++p)
        {
          const E_Int& fid = pt_to_faces[i][p];
          auto PG = m2.element(fid);
          E_Float nn[3];
          PG.template ndS<3>(m2.crd, nn);
          NUGA::normalize<3>(nn);

#ifdef DEBUG_IMMERSION
          E_Float l2 = ::sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
          assert(::fabs(l2 - 1.) < EPSILON);// DEGEN
#endif
          dirs[i].vec[0] += nn[0];
          dirs[i].vec[1] += nn[1];
          dirs[i].vec[2] += nn[2];
        }
      }
      else // this vertex is near a m2's node => compute nodal normal : weighted with angle between each pair of ray in the node shell (sum is 2PI)
      {
        //std::cout << "compute node normal" << std::endl;
        double TOL2 = ARTOL * ARTOL;
        if (ARTOL < 0.) //relative
          TOL2 *= vertices[i].val2;

        // find near node : 
        E_Int PGi0 = pt_to_faces[i][0];
        auto PG0 = m2.element(PGi0);
        double dmin2 = FLOAT_MAX;
        E_Int Ni{ IDX_NONE };
        for (E_Int u = 0; u < PG0.nb_nodes(); ++u)
        {
          double d2 = NUGA::sqrDistance(m2.crd.col(PG0.node(u)), vertices[i].vec, 3);
          if (d2 < dmin2)
          {
            dmin2 = d2;
            Ni = PG0.node(u);
          }
        }
        assert(dmin2 < TOL2);
        assert(Ni != IDX_NONE);

        //std::cout << "Near node : " << Ni << std::endl;

        size_t count{ 0 };

        //std::cout << "compute the rays angles" << std::endl;

        for (size_t p = 0; p < pt_to_faces[i].size(); ++p)
        {
          const E_Int& fid = pt_to_faces[i][p];
          auto PG = m2.element(fid);

          //std::cout << p << "-th face : " << fid << std::endl;
          
          E_Int nnodes = PG.nb_nodes();

          E_Int pos = K_CONNECT::IdTool::get_pos(PG.begin(), nnodes, Ni + 1);
          assert(pos != -1);
          if (pos == -1) continue;
          ++count;

          //std::cout << "pos ? : " << pos << std::endl;

          E_Int Nip1 = PG.node((pos + 1) % nnodes);
          E_Int Nim1 = PG.node((pos + nnodes - 1) % nnodes);

          //std::cout << "before angular" << std::endl;

          E_Float nn[3];
          NUGA::angular_weighted_normal(m2.crd.col(Nim1), m2.crd.col(Ni), m2.crd.col(Nip1), nn);

          dirs[i].vec[0] += nn[0];
          dirs[i].vec[1] += nn[1];
          dirs[i].vec[2] += nn[2];
        }

        if (count != pt_to_faces[i].size())
          std::cout << "WRONG LOGIC" << std::endl;
      }

      // normalize
      NUGA::normalize<3>(dirs[i].vec);
    }

    // compute displacements norms

    m2.get_nodal_metric2(mtype);
    DELAUNAY::Triangulator dt;

    for (size_t i = 0; i < vertices.size(); ++i)
    {
      dirs[i].val2 = 0.;

      if (pt_to_faces[i].empty()) continue; // regular point

      const double* Pi = vertices[i].vec;
 
      // compute max distance among sticking PGs to satisfy all of them
      for (size_t p = 0; p < pt_to_faces[i].size(); ++p)
      {
        const E_Int& fid = pt_to_faces[i][p];
        auto PG = m2.element(fid);

        double TOLi = ARTOL;
        if (ARTOL < 0.) // relative
        {
          double PGLref2 = PG.Lref2(m2.nodal_metric2);
          TOLi *= -::sqrt(std::min(vertices[i].val2, PGLref2));
        }
        double hmove = (1. + EPSILON) * TOLi;

        PG.triangulate(dt, m2.crd);

        double lambdaImin = FLOAT_MAX;
        double wmin{ 0 };
        // find the closest triangle
        for (E_Int t = 0; t < PG.nb_tris(); ++t)
        {
          E_Int T[3];
          PG.triangle(t, T);

          const E_Float * P0 = m2.crd.col(T[0]);
          const E_Float * P1 = m2.crd.col(T[1]);
          const E_Float * P2 = m2.crd.col(T[2]);

          E_Float W[3];
          K_MESH::Triangle::normal(P0, P1, P2, W);

          double w = NUGA::dot<3>(W, dirs[i].vec);
          if (::fabs(w) < ZERO_M) continue; // should not happen : dir must not be ortho to nPG

          double I[3]; // intersection of pt in computed dir on PG plane.
          double lambdaI = NUGA::project(P0, W, Pi, dirs[i].vec, I); // == Pi + lambdaI*dir = I

          if (::fabs(lambdaI) < ::fabs(lambdaImin))
          {
            lambdaImin = lambdaI;
            wmin = w;
          }
        }

        assert(wmin != 0.);
        assert(lambdaImin != FLOAT_MAX);

        wmin = 1. / wmin;

        bool reverse = ((inward && (wmin > 0.)) || (!inward && (wmin < 0.)));

        if (reverse)
        {
          dirs[i].vec[0] = -dirs[i].vec[0];
          dirs[i].vec[1] = -dirs[i].vec[1];
          dirs[i].vec[2] = -dirs[i].vec[2];
          lambdaImin = -lambdaImin;
        }

        double lambdaMove = (::fabs(wmin) * hmove) + lambdaImin;

        // following assert because dir must be well oriented : inward the surface if above, outward otherwise
        // so attractive if above, repulsive if below (to put it far enough to exit interference zone)
        assert(lambdaMove > 0.);
        //lambdaMove = w * (h + hmove);
        dirs[i].val2 = std::max(lambdaMove*lambdaMove, dirs[i].val2);
      }
    }

    // returns moved vertices and displacements
    std::vector<direction> odirs;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      if (pt_to_faces[i].empty()) continue; // regular point

      dirs[i].flag = i; //store node id
      dirs[i].flag2 = pt_to_faces[i][0];
      odirs.push_back(dirs[i]);
    }

    return odirs;
  }

  
  template <typename mesh_t>
  static bool move_matching_nodes(const K_FLD::FloatArray& crd, const std::vector<direction>& dirs, mesh_t& m)
  {
    using acrd_t = K_FLD::ArrayAccessor<K_FLD::FloatArray>;
    acrd_t acrd(crd);
    K_SEARCH::KdTree<> tree(acrd, EPSILON);
    bool found{ false };
    double d2;
    std::vector<bool> targeted(crd.cols(), false);
    for (E_Int i = 0; i < m.crd.cols(); ++i)
    {
      int N = tree.getClosest(m.crd.col(i), ZERO_M*ZERO_M, d2);

      if (N == IDX_NONE) continue;
      assert(targeted[N] == false);

      targeted[N] = true;
      NUGA::sum<3>(::sqrt(dirs[N].val2), dirs[N].vec, 1., m.crd.col(i), m.crd.col(i));
      found = true;
    }
    return found;
  }


  enum eBCType { BCNONE = 0, BCWALL = 1 };

  ///
  template <typename bound_mesh_t> inline
    void move_double_walls
    (bound_mesh_t* bit, const bound_mesh_t& zbound, double ARTOL, eMetricType mtype, double AMAX, const std::vector<E_Int>& zmask_wall_ids, std::vector<bool>& is_dw)
  {
    is_dw.clear();

    // detect singular nodes on mask walls and compute their displacement (outward of zbound)
    bound_mesh_t only_mask_walls(*bit);
    std::vector<E_Int> ptoids;

    int nbcells = bit->ncells();

    bound_mesh_t zwalls(zbound);
    //E_Int nb_occ1 = std::count(ALL(zbound.cnt._type), 1);
    //std::cout << "nb of walls in zbound : " << nb_occ1 << std::endl;
    std::vector<bool> keep(zwalls.ncells(), false);
    for (size_t u = 0; u < keep.size(); ++u)
      keep[u] = (zwalls.cnt._type[u] == BCWALL);

    zwalls.compress(keep);

    if (zwalls.ncells() == 0) return;

    {
      std::vector<bool> keep(nbcells, false);
      for (size_t u = 0; u < zmask_wall_ids.size(); ++u)
      {
        E_Int wid = zmask_wall_ids[u];
        keep[wid] = true;
      }
      std::vector<E_Int> ptnids;
      only_mask_walls.compress(keep, ptnids);
      K_CONNECT::IdTool::reverse_indirection(ptnids, ptoids);

#ifdef DISPLACEMENT_DBG
      {
        std::ostringstream o;
        o << "mask_walls_";// << i;
        medith::write<>(o.str().c_str(), only_mask_walls.crd, only_mask_walls.cnt);
        medith::write<>("zwalls", zwalls.crd, zwalls.cnt);
      }
#endif
    }
    std::vector<direction> nodes_to_move = NUGA::compute_displacement_for_singular_nodes(only_mask_walls, zwalls, ARTOL, mtype, false /*so outward*/);

    // convert this compact info
    std::vector<direction> all_nodes_disp(bit->crd.cols());
    std::vector<E_Int> is_singular_node(bit->crd.cols(), IDX_NONE);
    for (size_t i = 0; i < nodes_to_move.size(); ++i)
    {
      E_Int Ni = nodes_to_move[i].flag;
      E_Int oNi = ptoids[Ni];
      all_nodes_disp[oNi].val2 = nodes_to_move[i].val2;
      all_nodes_disp[oNi].vec[0] = nodes_to_move[i].vec[0];
      all_nodes_disp[oNi].vec[1] = nodes_to_move[i].vec[1];
      all_nodes_disp[oNi].vec[2] = nodes_to_move[i].vec[2];

      is_singular_node[oNi] = nodes_to_move[i].flag2; // for each node to move, gives face id in zbound to evaluate colinearity
    }

    // flag double-walls (DW) on mask to emerge them and then discard them 
    // (those with singular nodes and overlapping zboundaries)
    is_dw.clear();
    is_dw.resize(nbcells, false);

    for (E_Int i = 0; i < bit->cnt.size(); ++i)
    {
      const E_Int* nodes = bit->cnt.get_facets_ptr(i);
      E_Int nnodes = bit->cnt.stride(i);

      auto PGi = bit->element(i);
      E_Float nni[3];
      PGi.template ndS<3>(bit->crd, nni);
      NUGA::normalize<3>(nni);

      NUGA::vecval average_n;
      average_n.vec[0] = average_n.vec[1] = average_n.vec[2] = 0.;
      bool singular{ false };
      for (E_Int j = 0; j < nnodes; ++j)
      {
        E_Int Ni = nodes[j] - bit->index_start;

        if (is_singular_node[Ni] != IDX_NONE)
        {
          singular = true;

          E_Int fid = is_singular_node[Ni];
          auto PG = zwalls.element(fid);
          E_Float nn[3];
          PG.template ndS<3>(zwalls.crd, nn);
          NUGA::normalize<3>(nn);

          average_n.vec[0] += nn[0];
          average_n.vec[1] += nn[1];
          average_n.vec[2] += nn[2];
        }
      }

      if (singular)
      {
        NUGA::normalize<3>(average_n.vec);
        /*if (i == 14183)
        {
        K_FLD::FloatArray crdtmp(bit->crd);
        E_Float nP[3];
        nP[0] = crdtmp(0, nodes[0] - 1) + average_n.vec[0];
        nP[1] = crdtmp(1, nodes[0] - 1) + average_n.vec[1];
        nP[2] = crdtmp(2, nodes[0] - 1) + average_n.vec[2];
        crdtmp.pushBack(nP, nP + 3);
        E_Int E[] = { nodes[0] - 1, crdtmp.cols() - 1 };
        K_FLD::IntArray cntmp;
        cntmp.pushBack(E, E + 2);
        std::ostringstream o;
        o << "norm2";
        medith::write(o.str().c_str(), crdtmp, cntmp, "BAR");
        }*/

        //double val = ::fabs(NUGA::dot<3>(nni, average_n.vec));
        double angle = NUGA::normals_angle(nni, average_n.vec);
        if (angle < AMAX)
          is_dw[i] = true;
      }
    }

#ifdef DISPLACEMENT_DBG
    {
      bound_mesh_t tmp(*bit);
      tmp.compress(is_dw);

      std::ostringstream o;
      o << "double_walls_";// << i;
      medith::write<>(o.str().c_str(), tmp.crd, tmp.cnt);
    }
#endif

    // compute displacement for missing nodes on double walls by extrapolation
    for (E_Int i = 0; i < bit->cnt.size(); ++i)
    {
      if (!is_dw[i]) continue;

      const E_Int* nodes = bit->cnt.get_facets_ptr(i);
      E_Int nnodes = bit->cnt.stride(i);
      E_Int nb_missing{ 0 };

      NUGA::vecval v;
      v.val2 = 0.;
      v.vec[0] = v.vec[1] = v.vec[2] = 0.;

      for (E_Int j = 0; j < nnodes; ++j)
      {
        E_Int Ni = nodes[j] - 1;
        if (all_nodes_disp[Ni].val2 == FLOAT_MAX)
          ++nb_missing;
        else
        {
          v.val2 += ::sqrt(all_nodes_disp[Ni].val2);
          v.vec[0] += all_nodes_disp[Ni].vec[0];
          v.vec[1] += all_nodes_disp[Ni].vec[1];
          v.vec[2] += all_nodes_disp[Ni].vec[2];
        }
      }

      if (nb_missing)
      {
        assert(nb_missing != nnodes);    // at least one singular node
        v.val2 /= (nnodes - nb_missing); // averaging the contributions
        v.val2 *= v.val2;
        NUGA::normalize<3>(v.vec);

        for (E_Int j = 0; j < nnodes; ++j)
        {
          E_Int Ni = nodes[j] - 1;

          if (all_nodes_disp[Ni].val2 == FLOAT_MAX)
          {
            all_nodes_disp[Ni].val2 = v.val2;
            all_nodes_disp[Ni].vec[0] = v.vec[0];
            all_nodes_disp[Ni].vec[1] = v.vec[1];
            all_nodes_disp[Ni].vec[2] = v.vec[2];
          }
        }
      }
    }

    // reset degenerating displacements : fixme : use neighborhood instead of looping on the whole mesh
    {
      bool carry_on = true;
      E_Int iter{0}, itermax{10};
      while (carry_on)
      {
#ifdef DISPLACEMENT_DBG
        std::cout << "reseting iter : " << iter << std::endl;
#endif
        DELAUNAY::Triangulator dt;
        K_FLD::FloatArray crd = bit->crd;
        E_Int bad{0};

        // move singular nodes
        for (size_t i = 0; i < all_nodes_disp.size(); ++i)
        {
          if (all_nodes_disp[i].val2 == FLOAT_MAX) continue;

          E_Float * Pi = crd.col(i);
          NUGA::sum<3>(1., Pi, ::sqrt(all_nodes_disp[i].val2), all_nodes_disp[i].vec, Pi);
        }

        std::vector<bool> validpt(crd.cols(), true);
        for (E_Int i=0; i < bit->ncells(); ++i)
        {
          auto elt = bit->element(i);
          //std::cout << "triangulating" << std::endl;
          E_Int err = elt.triangulate(dt, crd);
          if (err)
          {
#ifdef DISPLACEMENT_DBG
            std::cout << "error for face : " << i << " : need to unmove" << std::endl;
#endif
            ++bad;
            for (E_Int n = 0; n < elt.nb_nodes(); ++n)
            {
              E_Int Ni = elt.node(n);
              assert (Ni > -1 && (size_t)Ni < all_nodes_disp.size());
              //std::cout << "disabling node : " << Ni << std::endl;
              all_nodes_disp[Ni].val2 = FLOAT_MAX;
            }
          }
        }

        carry_on = (bad > 0) && (iter++ < itermax);
      }
    }

    // move singular nodes
    for (size_t i = 0; i < all_nodes_disp.size(); ++i)
    {
      if (all_nodes_disp[i].val2 == FLOAT_MAX) continue;

      E_Float * Pi = bit->crd.col(i);
      NUGA::sum<3>(1., Pi, ::sqrt(all_nodes_disp[i].val2), all_nodes_disp[i].vec, Pi);
    }

#ifdef DISPLACEMENT_DBG
    {
      bound_mesh_t tmp(*bit);
      tmp.compress(is_dw);

      std::ostringstream o;
      o << "moved_double_walls_";// << i;
      medith::write<>(o.str().c_str(), tmp.crd, tmp.cnt);
    }
#endif
  }

  /// 
  template <> inline
    void move_double_walls<edge_mesh_t>
    (edge_mesh_t* bit, const edge_mesh_t& zbound, double ARTOL, eMetricType mtype, double AMAX, const std::vector<E_Int>& zmask_wall_ids, std::vector<bool>& is_dw)
  {
    // not implemented for surfaces
  }

}

#endif

