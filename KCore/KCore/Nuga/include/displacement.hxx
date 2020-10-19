/*
 
 
 
              NUGA 
 
 
Author : Sam Landier (sam.landier@onera.fr) 
 */

#ifndef NUGA_DISPLACEMENT_HXX
#define NUGA_DISPLACEMENT_HXX

#include "Nuga/include/vertex.h"
#include "Nuga/include/mesh_t.hxx"
#include "Nuga/include/collider.hxx"

namespace NUGA
{

  template <typename mesh_t1, typename mesh_t2>
  static std::vector<direction> immerse_nodes(mesh_t1& m1, const mesh_t2& m2, double ARTOL)
  {
    if (m2.ARG1 == SURFACIC)
      assert(m2.oriented == 1);

    K_FLD::FloatArray weirds;

    // get vertices from mesh1
    std::vector<NUGA::vertex> vertices;

    auto sqrmetric1 = m1.get_nodal_metric2();

    vertices.reserve(m1.crd.cols());
    for (size_t i = 0; i < m1.crd.cols(); ++i)
      vertices.push_back(NUGA::vertex(m1.crd.col(i), sqrmetric1[i]));

    // identify singular points
    std::vector<std::vector<E_Int>> pt_to_faces;
    NUGA::COLLIDE::get_colliding(vertices, m2, ARTOL, pt_to_faces);

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

        //weirds.push_back()

        // find near node : 
        E_Int PGi0 = pt_to_faces[i][0];
        auto PG0 = m2.element(PGi0);
        double dmin2= FLOAT_MAX;
        E_Int Ni{IDX_NONE};
        for (size_t u=0; u < PG0.nb_nodes(); ++u) 
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

        if (count != pt_to_faces[i].size())
          std::cout << "WRONG LOGIC" << std::endl;
      }

      // normalize
      NUGA::normalize<3>(dirs[i].vec);
    }

    // compute displacements norms

    m2.get_nodal_metric2();

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
        h = ::fabs(h); // if h < 0, already immersed

        double w = NUGA::dot<3>(nPG, dirs[i].vec);
        assert(w < 0.);
        w = 1. / w;

        double TOLi = ARTOL;
        if (ARTOL < 0.) // relative
        {
          double PGLref2 = PG.Lref2(m2.nodal_metric2);
          TOLi *= -::sqrt(std::min(vertices[i].val2, PGLref2));
        }

        double hmove = (1. + EPSILON) * TOLi;
        double lambdaMove = w * (h - hmove);

        // following assert because dir must be well oriented : inward the surface if above, outward otherwise
        // so attractive if above, repulsive if below (to put it far enough to exit interfernce zone
        assert(lambdaMove > 0.);
        //lambdaMove = w * (h + hmove);
        dirs[i].val2 = std::max(lambdaMove*lambdaMove, dirs[i].val2);
      }
    }

    // apply displacements to m1
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      if (pt_to_faces[i].empty()) continue; // regular point

      NUGA::sum<3>(::sqrt(dirs[i].val2), dirs[i].vec, 1., m1.crd.col(i), m1.crd.col(i));
    }

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
      double d2;
      int N = tree.getClosest(m.crd.col(i), ZERO_M*ZERO_M, d2);

      if (N == IDX_NONE) continue;
      assert(targeted[N] == false);

      targeted[N] = true;
      NUGA::sum<3>(::sqrt(dirs[N].val2), dirs[N].vec, 1., m.crd.col(i), m.crd.col(i));
      found = true;
    }
    return found;
  }

}

#endif

