/*    
    Copyright 2013-2019 Onera.

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

#include "MeshTool.h"
#include "MeshElement/Triangle.h"
#include <math.h>
#include "Connect/EltAlgo.h"
#include "Connect/GeomAlgo.h"
#include "Connect/BARSplitter.h"
#include "Search/BbTree.h"
#ifdef DEBUG_MESHTOOL
#include <sstream>
#include "IO/io.h"
#endif
//#include <hash_set>

using namespace K_CONST;

K_CONNECT::MeshTool::MeshTool(const tree_type& tree, E_Float tolerance):
_tree(tree), _tolerance(tolerance)
{
}

K_CONNECT::MeshTool::~MeshTool(void)
{
}

/*****************  private ***********************************/

E_Int
K_CONNECT::MeshTool::__getContainingElement
(const E_Float* point, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
 const K_FLD::IntArray& neighbors, const K_CONT_DEF::bool_vector_type& mask, size_type Kseed) const
{
  E_Bool      brute_force(Kseed == E_IDX_NONE), found(false);
  size_type   K(Kseed), nb_tris(connect.cols()), b(0), visited(-1);

  if (brute_force)
  {
    for (size_type i = 0; (i < nb_tris) && (b >= 0); ++i)
    {
      K = i;
      if (!mask[K]) //skip invalidated elements
        continue; 
      b = __searchDirection(K, point, pos, connect, /*random selection*/false);
    }
  }
  else
  {
    while (!found && ++visited < nb_tris)
    {
      b = __searchDirection(K, point, pos, connect, /*random selection*/true);
      found = (b == -1);
      if (b >= 0)
        K = neighbors(b, K);
    }
  }

  return (b == -1) ? K : E_IDX_NONE;
}

E_Int
K_CONNECT::MeshTool::__searchDirection
(size_type K, const E_Float* point, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Bool random) const
{
  size_type                         b(0), Ni, Nj;
  K_FLD::IntArray::const_iterator    pK = connect.col(K);

  if (random) // Random value between 0,1 and 2.
    b = std::rand() % K_MESH::Triangle::NB_NODES;//fixme : not a good random...

  // The right direction is the edge without any visbility on P.
  for (size_type i = 0; i < K_MESH::Triangle::NB_NODES; ++i)
  {
    Ni = *(pK + (b+1) % K_MESH::Triangle::NB_NODES);
    Nj = *(pK + (b+2) % K_MESH::Triangle::NB_NODES);

    if (K_MESH::Triangle::surface<2>(pos.col(Ni), pos.col(Nj), point) < 0.)
      return b;

    b = (b+1) % K_MESH::Triangle::NB_NODES;
  }

  return (-1); // P is inside K : the 3 surfaces are positive i.e the edges are all visible from P.
}

void
K_CONNECT::MeshTool::getBoundary(const K_FLD::IntArray& connect, K_FLD::IntArray& connectB)
{
  std::map<K_MESH::NO_Edge, E_Int> edge_to_count;
  std::map<K_MESH::NO_Edge, E_Int>::iterator it;
  E_Int COLS(connect.cols()), n, Ni, Nj;
  K_MESH::NO_Edge b;
  E_Int E[2];
  E_Int NBNODES = connect.rows();
  K_FLD::IntArray::const_iterator pS;

  connectB.clear();

  for (E_Int Si = 0; Si < COLS; ++Si)
  {
    pS = connect.col(Si);

    for (n = 0; n < NBNODES; ++n)
    {
      Ni = *(pS+n);
      Nj = *(pS+(n+1)%NBNODES);
      b.setNodes(Ni, Nj);
      it = edge_to_count.find(b);

      if (it == edge_to_count.end())
        edge_to_count.insert(std::make_pair(b, 1));
      else
        ++it->second;
    }
  }

  for (E_Int Si = 0; Si < COLS; ++Si)
  {
    pS = connect.col(Si);

    for (n = 0; n < NBNODES; ++n)
    {
      Ni = *(pS+n);
      Nj = *(pS+(n+1)%NBNODES);
      b.setNodes(Ni, Nj);
      
      if (edge_to_count[b] == 1)
      {
        E[0] = Ni; E[1] = Nj;
        connectB.pushBack(E, E+2);
      }
    }
  }
}

///
void
K_CONNECT::MeshTool::getBoundaryT3Mesh
(const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors, K_CONT_DEF::int_pair_vector_type& boundaries)
{
  K_FLD::IntArray::const_iterator pVi;
  boundaries.clear();
  
  for (E_Int i = 0; i < neighbors.cols(); ++i)
  {
    pVi = neighbors.col(i);
    for (E_Int n = 0; n < 3; ++n)
    {
      if (*(pVi+n) == E_IDX_NONE)
        boundaries.push_back(std::make_pair(i, n));
    }
  }
}

///
void
K_CONNECT::MeshTool::getBoundaryT3Mesh
(const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors, K_FLD::IntArray& cB)
{
  K_FLD::IntArray::const_iterator pVi, pKi;
  cB.clear();

  for (E_Int i = 0; i < neighbors.cols(); ++i)
  {
    pVi = neighbors.col(i);
    pKi = connect.col(i);
    for (E_Int n = 0; n < 3; ++n)
    {
      if (*(pVi + n) == E_IDX_NONE)
      {
        E_Int E[] = { *(pKi + (n + 1) % 3), *(pKi + (n + 2) % 3) };
        cB.pushBack(E, E + 2);
      }
    }
  }
}
///
void
K_CONNECT::MeshTool::getNonManifoldEdges
(const K_FLD::IntArray& connect, K_FLD::IntArray& connectB)
{
  std::map<K_MESH::NO_Edge, E_Int> edge_to_count;
  std::map<K_MESH::NO_Edge, E_Int>::iterator it;
  E_Int COLS(connect.cols()), n, Ni, Nj;
  K_MESH::NO_Edge b;
  E_Int E[2];
  E_Int NBNODES = connect.rows();
  K_FLD::IntArray::const_iterator pS;

  connectB.clear();

  for (E_Int Si = 0; Si < COLS; ++Si)
  {
    pS = connect.col(Si);

    for (n = 0; n < NBNODES; ++n)
    {
      Ni = *(pS+n);
      Nj = *(pS+(n+1)%NBNODES);
      b.setNodes(Ni, Nj);

      it = edge_to_count.find(b);
      if (it == edge_to_count.end())
        edge_to_count.insert(std::make_pair(b, 1));
      else
        ++it->second;
    }
  }

  for (E_Int Si = 0; Si < COLS; ++Si)
  {
    pS = connect.col(Si);

    for (n = 0; n < NBNODES; ++n)
    {
      Ni = *(pS+n);
      Nj = *(pS+(n+1)%NBNODES);
      b.setNodes(Ni, Nj);
      
      if (edge_to_count[b] > 2)
      {
        E[0] = Ni; E[1] = Nj;
        connectB.pushBack(E, E+2);
      }
    }
  }
}

///
void
K_CONNECT::MeshTool::getAttachedElements
(const K_FLD::IntArray& connectB, const K_FLD::IntArray& connectS, 
 K_FLD::IntArray& connectElts)
{
  K_CONT_DEF::int_vector_type     nodes;
  K_CONT_DEF::int_set_type        bnodes;
  E_Int                           nb_elts(connectS.cols()), nb_nodes(connectS.rows());
  K_FLD::IntArray::const_iterator pS;

  connectElts.clear();

  connectB.uniqueVals(nodes);
  bnodes.insert(nodes.begin(), nodes.end());

  for (E_Int i = 0; i < nb_elts; ++i)
  {
    pS = connectS.col(i);
    for (E_Int n = 0; n < nb_nodes; ++n)
      if (bnodes.find(*(pS+n)) != bnodes.end())
      {
        connectElts.pushBack(pS, pS+nb_nodes);
        break;
      }
  }
}

void
K_CONNECT::MeshTool::computeNodeNormals
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, K_FLD::FloatArray& normals)
{
  K_CONT_DEF::int_vector_type     nodes, count;
  E_Float                         zero(0), normal[3], *p;
  K_FLD::IntArray                 connectElts;
  K_FLD::IntArray::const_iterator pS;
  E_Int                           nb_elts;

  if (connect.cols() == 0)
    return;

  normals.clear();
  normals.resize(3, pos.cols(), &zero);
  count.resize(pos.cols(), 0);

  nb_elts = connect.cols();
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    pS = connect.col(i);
    element_type::normal(pos.col(*pS), pos.col(*(pS+1)), pos.col(*(pS+2)), normal);
    for (E_Int n = 0; n < element_type::NB_NODES; ++n)
    {
      for (E_Int k = 0; k < 3; ++k)
        normals(k, *(pS+n)) += normal[k];

      ++count[*(pS+n)];
    }
  }
  for (E_Int n = 0; n < pos.cols(); ++n)
  {
    if (count[n] > 1)
    {
      for (E_Int k = 0; k < 3; ++k)
        normals(k,n) /= count[n];
    }

    if (count[n] != 0)
    {
      p = normals.col(n);
      K_FUNC::normalize<3>(p);
    }
  }
}

E_Int
K_CONNECT::MeshTool::computeNodeNormalsFromPGNormals
(const K_FLD::FloatArray& PG_normals, const ngon_unit& pgs, K_FLD::FloatArray& node_normals)
{
  K_CONT_DEF::int_vector_type     count;
  E_Float                         *p;
  E_Int nb_pgs = pgs.size();
  
  if (nb_pgs == 0)
    return 0;
  
  E_Int idmaxp1 = pgs.get_facets_max_id();

  node_normals.clear();
  node_normals.resize(3, idmaxp1, 0.);
  count.resize(idmaxp1, 0);

  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    const E_Int* nodes = pgs.get_facets_ptr(i);
    E_Int nb_nodes = pgs.stride(i);
    
    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      E_Int Ni = *(nodes + n) - 1;
      p = node_normals.col(Ni);
      for (E_Int k = 0; k < 3; ++k)
        p[k] += PG_normals(k,i);
       ++count[Ni];
    }
  }
  
  E_Int c = 0;
  for (E_Int n = 0; n < idmaxp1; ++n)
  {
    p = node_normals.col(n);
    
    if (count[n] > 1)
    {
      for (E_Int k = 0; k < 3; ++k)
        p[k] /= count[n];
      K_FUNC::normalize<3>(p);
      ++c;
    }
  }
  
  return idmaxp1;
}

E_Int
K_CONNECT::MeshTool::computeNodeNormals
(const K_FLD::FloatArray& crd, const ngon_unit& pgs, K_FLD::FloatArray& normals, E_Int smooth_iters)
{
  // First Compute A Valid PG normal field : if a normal failed, look around
  E_Int nb_pgs = pgs.size();
  K_FLD::FloatArray PGnormals(3, nb_pgs, 0.);
  E_Float normal[3];
  std::vector<bool> is_degen(nb_pgs, false);
  E_Int nb_degen=0;
  
  typedef K_FLD::ArrayAccessor<K_FLD::FloatArray > acrd_t;
  acrd_t acrd(crd);
  
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    const E_Int* nodes = pgs.get_facets_ptr(i);
    E_Int nb_nodes = pgs.stride(i);
    K_MESH::Polygon::normal<acrd_t, 3>(acrd, nodes, nb_nodes, 1, normal);
    E_Float l2 = ::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    
    is_degen[i] = (::fabs(l2 - 1.) >= E_EPSILON); // DEGEN
    
    if (!is_degen[i]) 
    {
      PGnormals(0,i) = normal[0];
      PGnormals(1,i) = normal[1];
      PGnormals(2,i) = normal[2];
    }
    else
      ++nb_degen;
  }
  
  ngon_unit neighbors;
  if (nb_degen)
    K_MESH::Polygon::build_pg_neighborhood(pgs, neighbors);

  while (nb_degen > 0)
  {
    K_CONT_DEF::int_vector_type     count(nb_pgs);
    E_Int nb_degen0=nb_degen;

    // Smooth with nbeighbors
    for (E_Int i=0; i < nb_pgs; ++i)
    {
      if (!is_degen[i]) continue;
      
      E_Int nb_neighs = neighbors.stride(i);
      E_Int* neighs = neighbors.get_facets_ptr(i);
      E_Int nb_contrib=0;
      
      for (E_Int n=0; n < nb_neighs; ++n)
      {
        E_Int PGn = *(neighs + n);
        if (PGn == E_IDX_NONE) continue;
        if (is_degen[PGn]) continue;
        
        PGnormals(0,i) += PGnormals(0,PGn);
        PGnormals(1,i) += PGnormals(1,PGn);
        PGnormals(2,i) += PGnormals(2,PGn);
        
        ++nb_contrib;
      }
      
      if (nb_contrib)
      {
        --nb_degen;
        is_degen[i]=false;
        for (E_Int k = 0; k < 3; ++k)
          PGnormals(k,i) /= nb_contrib;
        K_FUNC::normalize<3>(PGnormals.col(i));
      }
    }
    
    if (nb_degen0 == nb_degen) break; // error ? some island are there.
  };

  // average it at nodes
  E_Int idmax = computeNodeNormalsFromPGNormals(PGnormals, pgs, normals);
  
  if (smooth_iters > 0)//smooth it a bit
    smoothNodeNormals(pgs, normals, smooth_iters);

  return idmax;
}

E_Int K_CONNECT::MeshTool::smoothNodeNormals(const ngon_unit& pgs, K_FLD::FloatArray& normals, E_Int smooth_iters)
{
  if (smooth_iters <= 0) return 0;
  
  //std::cout << "smoothing..." << std::endl;
  
  E_Int nb_pgs = pgs.size();
  
  while (smooth_iters--)
  {
    for (E_Int i = 0; i < nb_pgs; ++i)
  {
    const E_Int* nodes = pgs.get_facets_ptr(i);
    E_Int nb_nodes = pgs.stride(i);
    
    E_Float incr[3], FACTOR = 0.5;
    
    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      E_Int Nim1 = (n == 0) ? *(nodes + nb_nodes - 1) - 1 : *(nodes + n - 1) - 1;
      E_Int Ni = *(nodes + n) - 1;
      E_Int Nip1 = *(nodes + (n+1)%nb_nodes) - 1;
      
      K_FUNC::sum<3>(0.5, normals.col(Nim1), 0.5, normals.col(Nip1), incr);
      
      
      K_FUNC::normalize<3>(incr);
      
      E_Float l2 = ::sqrt(incr[0] * incr[0] + incr[1] * incr[1] + incr[2] * incr[2]);
      if(::fabs(l2 - 1.) >= E_EPSILON) // DEGEN
        continue;
      
      K_FUNC::sum<3>(1. - FACTOR, incr, FACTOR, normals.col(Ni), normals.col(Ni));
      K_FUNC::normalize<3>(normals.col(Ni));
    }
  }
  }
  return 1;  
}

void
K_CONNECT::MeshTool::compact_to_mesh
(K_FLD::FloatArray& pos, K_FLD::IntArray& connect, std::vector<E_Int>& new_IDs, E_Int* N0)
{
  // Compact
  E_Int ROWS(connect.rows());
  std::vector<bool> flag;
  if (N0 && (*N0 > 0) && (*N0 < pos.cols()))
    flag.resize(*N0, true);
  flag.resize(pos.cols(), false);
  for (E_Int i = 0; i < connect.cols(); ++i)
  {
    for (E_Int k = 0; k < ROWS; ++k)
      flag[connect(k,i)] = true;
  }

  K_FLD::FloatArray::compact(pos, flag, new_IDs);  // Remove any node not in the mesh and above N0.
  K_FLD::IntArray::changeIndices(connect, new_IDs);// Change mesh indices accordingly.
}

void
K_CONNECT::MeshTool::compact_to_mesh
(const K_FLD::FloatArray& pos0, const K_FLD::IntArray& connect0,
 K_FLD::FloatArray& pos1, K_FLD::IntArray& connect1,
 std::vector<E_Int>& oids)
{
  pos1.clear();
  connect1.clear();
  oids.clear();
  
  size_t id=0, ROWS(connect0.rows()), COLS(connect0.cols());
  const E_Int* pS0;
  E_Int* pS1;
  std::map<E_Int, E_Int> nids;
  std::map<E_Int, E_Int>::const_iterator it;
  
  connect1.resize(ROWS, COLS);
  
  for (size_t i=0; i < COLS; ++i)
  {
    pS0 = connect0.col(i);
    pS1 = connect1.col(i);
    for (size_t n=0; n < ROWS; ++n)
    {
      const E_Int& Ni = *(pS0+n);
      it = nids.find(Ni);
      if (it == nids.end()){
        nids[Ni]=id++;
        oids.push_back(Ni);
      }
        
      *(pS1+n)=nids[Ni];
    }
  }
  
  ROWS=pos0.rows();
  pos1.reserve(ROWS, id);
  for (size_t i=0; i < id; ++i)
    pos1.pushBack(pos0.col(oids[i]), pos0.col(oids[i])+ROWS);
}

///
void
K_CONNECT::MeshTool::boundingBox
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectE2,
 E_Float* minB, E_Float* maxB)
{
  E_Int COLS(connectE2.cols());
  K_FLD::IntArray::const_iterator pS;
  E_Int Ni;
  const E_Float* Pi;

  for (E_Int k = 0; k < 3; ++k)
  {
    minB[k] = K_CONST::E_MAX_FLOAT;
    maxB[k] = -K_CONST::E_MAX_FLOAT;
  }

  // Bounding Box.
  for (E_Int i = 0; i < COLS; ++i)
  {
    pS = connectE2.col(i);
    for (E_Int n = 0; n < 2; ++n)
    {
      Ni = *(pS+n);
      Pi = pos.col(Ni);
      for (E_Int k = 0; k < 3; ++k)
      {
        minB[k] = (minB[k] > Pi[k]) ? Pi[k] : minB[k];
        maxB[k] = (maxB[k] < Pi[k]) ? Pi[k] : maxB[k];
      }
    }
  }
}

///
void
K_CONNECT::MeshTool::boundingBox
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
 E_Int* minN, E_Int* maxN)
{
  E_Int COLS(connect.cols()), ROWS(connect.rows());
  K_FLD::IntArray::const_iterator pS;
  E_Int Ni;
  const E_Float* Pi;
  E_Float minB[3], maxB[3];

  for (E_Int k = 0; k < 3; ++k)
  {
    minB[k] = K_CONST::E_MAX_FLOAT;
    maxB[k] = -K_CONST::E_MAX_FLOAT;
    minN[k] = E_IDX_NONE;
    maxN[k] = E_IDX_NONE;
  }

  // Bounding Box.
  for (E_Int i = 0; i < COLS; ++i)
  {
    pS = connect.col(i);
    for (E_Int n = 0; n < ROWS; ++n)
    {
      Ni = *(pS+n);
      Pi = pos.col(Ni);
      for (E_Int k = 0; k < 3; ++k)
      {
        if (minB[k] > Pi[k])
        {
          minB[k] = Pi[k];
          minN[k] = Ni;
        }
        if (maxB[k] < Pi[k])
        {
          maxB[k] = Pi[k];
          maxN[k] = Ni;
        }
      }
    }
  }
}

///
void
K_CONNECT::MeshTool::flipT3(K_FLD::IntArray& connectT3)
{
  K_FLD::IntArray::iterator pS;
  E_Int COLS(connectT3.cols());
  for (E_Int c = 0; c < COLS; ++c)
  {
    pS = connectT3.col(c);
    std::swap(*pS, *(pS+1));
  }
}

void
K_CONNECT::MeshTool::metricAnisoE22D
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectE2, double kn, K_FLD::FloatArray& M)
{
  M.clear();
  M.resize(3, pos.cols());

  K_FLD::FloatArray L, Me(2, 2), Mi(2, 2), tmp(2,2);
  computeEdgesSqrLengths<2>(pos, connectE2, L);
  
  E_Int Ni, Nj, NBEDGES(connectE2.cols());

  // Initialize as iso.
  for (E_Int i = 0; i < NBEDGES; ++i)
  {
    Ni = connectE2(0,i);
    Nj = connectE2(1,i);
    M(0,Ni) = M(2, Ni) = M(0,Nj) = M(2, Nj) = 1. / L(0, i);
    M(1, Ni) = M(1, Nj) =0.;
  }

  
  E_Float h, m11, m22, Vt[2], Vn[2];
  for (E_Int i = 0; i < NBEDGES; ++i)
  {
    Ni = connectE2(0,i);
    Nj = connectE2(1,i);

    K_FUNC::diff<2>(pos.col(Nj), pos.col(Ni), Vt);
    h = K_FUNC::normalize<2>(Vt);
    Vn[0] = -Vt[1];
    Vn[1] = Vt[0];

    m11 = 1./(h*h);
    if (kn < 0.)
      m22 = m11/(kn*kn);
    else if (kn > 0.)
      m22 = 1./(kn*kn);
    else
      m22 = m11;
    
    //edge metric Me
    Me(0,0) = m11*Vt[0]*Vt[0] + m22*Vn[0]*Vn[0];
    Me(0,1) = Me(1,0) = m11*Vt[0]*Vt[1] + m22*Vn[0]*Vn[1];
    Me(1,1) = m11*Vt[1]*Vt[1] + m22*Vn[1]*Vn[1];

    for (E_Int n = 0; n < 2; ++n)
    {
      Ni = connectE2(n,i);
      // Mi = Me inter Mi
      Mi(0,0) = M(0,Ni);
      Mi(1,0) = Mi(0,1) = M(1,Ni);
      Mi(1,1) = M(2, Ni);
      K_LINEAR::DelaunayMath::intersect(Mi, Me, tmp);
      M(0,Ni) = Me/*tmp*/(0,0);
      M(1, Ni) = Me/*tmp*/(0,1);
      M(2, Ni) = Me/*tmp*/(1,1);
    }
  }
}

///
bool
K_CONNECT::MeshTool::detectDuplicated
(const K_FLD::IntArray& connect, std::vector<E_Int>& dupIds, bool strict_orient)
{
  if (strict_orient)
    return detectDuplicated<K_MESH::O_Triangle>(connect, dupIds);
  else
    return detectDuplicated<K_MESH::NO_Triangle>(connect, dupIds);
}

///
E_Int
K_CONNECT::MeshTool::removeDuplicated
(K_FLD::IntArray& connect, std::vector<E_Int>& dupIds, bool strict_orient)
{
  if (connect.rows() == 3)
  {
    if (strict_orient)
      return removeDuplicated<K_MESH::O_Triangle>(connect, dupIds);
    else
      return removeDuplicated<K_MESH::NO_Triangle>(connect, dupIds);
  }
  else if (connect.rows() == 2)
  {
    if (strict_orient)
      return removeDuplicated<K_MESH::Edge>(connect, dupIds);
    else
      return removeDuplicated<K_MESH::NO_Edge>(connect, dupIds);
  }
  return 0;
}

///
void
K_CONNECT::MeshTool::computeMinMaxIndices
(const K_FLD::IntArray& connect, E_Int& min_i, E_Int& max_i)
{
  E_Int                           COLS(connect.cols()), NB_NODES(connect.rows());
  K_FLD::IntArray::const_iterator pS;

  min_i = E_IDX_NONE;
  max_i = - 1;

  // Fast returns
  if (COLS == 0)
    return;

  for (E_Int c = 0; c < COLS; ++c)
  {
    pS = connect.col(c);

    for (E_Int n = 0; n < NB_NODES; ++n)
    {
      const E_Int& Ni = *(pS+n); 
      min_i = (Ni < min_i) ? Ni : min_i;
      max_i = (Ni > max_i) ? Ni : max_i;
    }
  }
}

E_Int K_CONNECT::MeshTool::get_polygonal_boundary
(const K_FLD::IntArray& cT3, const std::vector<E_Int>& ids, std::deque<E_Int>& PGi, 
 std::set<K_MESH::Edge>& w_oe_set, std::map<E_Int, E_Int>& w_n_map)
{
  // Assume consistent orientation and connex set

  w_oe_set.clear();
  w_n_map.clear();
  PGi.clear();

  std::set<K_MESH::Edge>::iterator it;
  K_MESH::Edge E, revE;
  std::map<E_Int, E_Int>::iterator itN;
  K_FLD::IntArray::const_iterator pKi;

  for (size_t i = 0; i < ids.size(); ++i)
  {
    E_Int ti = ids[i];
    pKi = cT3.col(ti);
    for (E_Int n = 0; n < 3; ++n)
    {
      E_Int Ni = *(pKi + n);
      E_Int Nj = *(pKi + (n + 1) % 3);

      revE.setNodes(Nj, Ni);
      it = w_oe_set.find(revE);
      if (it != w_oe_set.end())//already in
        w_oe_set.erase(*it);
      else
        w_oe_set.insert(K_MESH::Edge(Ni, Nj));
    }
  }

  if (w_oe_set.empty())
    return 1;

  for (it = w_oe_set.begin(); it != w_oe_set.end(); ++it)
  {
    E_Int Ni = it->node(0);
    E_Int Nj = it->node(1);

#ifdef DEBUG_IT
    itN = w_n_map.find(Ni);
    assert(itN == w_n_map.end());
#endif

    w_n_map[Ni] = Nj;//Nj is next to Ni
  }

  if (w_n_map.empty())
    return 1;

  E_Int Nbegin = w_n_map.begin()->first;
  E_Int count = w_n_map.size();
  E_Int Ncur = Nbegin;
  E_Int Nnext;

  do
  {
    Nnext = w_n_map[Ncur];
    PGi.push_back(Ncur);
    Ncur = Nnext;
  } 
  while (Nnext != Nbegin && count-- > 0);

  return ((PGi.size() != (ids.size()+2)) || (Nnext != Nbegin)); // otherwise multiple loops or error
}

E_Int  K_CONNECT::MeshTool::aggregate_convex
(const K_FLD::FloatArray& crd, const K_FLD::IntArray& connectT3, const E_Float* normal, ngon_unit& agg_pgs, E_Float convexity_tol)
{
  typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
  acrd_t acrd(crd);

  // WARNING : ASSUME THAT THE INPUT T3 set IS CONNEX
  K_FLD::IntArray neighbors;
  std::vector<K_FLD::IntArray*> pool, trash;

  // Init
  K_FLD::IntArray* connectM = new K_FLD::IntArray(connectT3);
  

#ifdef DEBUG_MESHTOOL
  //enabling_write("cM.mesh", crd, *connectM, "TRI");
#endif

  pool.push_back(connectM);
  trash.push_back(connectM);//to delete at the end

  std::vector<E_Int> PGi;
  K_FLD::IntArray connectB;

  E_Int err = 0;
  E_Int iter(0), iterMax(connectT3.cols() * 100);
  while (!pool.empty() && !err)
  {
    if (iter++ > iterMax)
    {
      err=1;
      break;
    }
    K_FLD::IntArray* ci = pool.back(); pool.pop_back();
    K_FLD::IntArray& connectMi = *ci;

    if (connectMi.cols() == 0)
      continue;
    
    if (connectMi.cols() == 1)
    {
      E_Int T[] = {connectMi(0,0)+1, connectMi(1,0) + 1, connectMi(2,0) + 1};
      
      agg_pgs.add(3, &T[0]);
      continue;
    }

#ifdef DEBUG_MESHTOOL
    std::ostringstream o;
    o << "cM.mesh";
    MIO::write(o.str().c_str(), crd, connectMi, "TRI");
#endif

    //Build the neighbourhood matrix.
    K_CONNECT::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(connectMi, neighbors);

#ifdef DEBUG_MESHTOOL  
    std::cout << neighbors << std::endl;
#endif

    // Get the contour.
    K_CONNECT::MeshTool::getBoundaryT3Mesh(connectMi, neighbors, connectB);

    // sort the nodes
    BARSplitter::getSortedNodes(connectB, PGi);
    
#ifdef DEBUG_MESHTOOL  
    std::cout << connectB << std::endl;
    for (size_t i=0; i < PGi.size(); ++i)
      std::cout << PGi[i] << std::endl;
#endif

    // get worst concavity
    E_Int iworst, ibest;
    E_Int nb_nodes = PGi.size();
    bool convex = K_MESH::Polygon::is_convex(crd, &PGi[0], nb_nodes, 0, normal, convexity_tol, iworst, ibest/*unused here*/);

    if (convex)
    {
      K_CONNECT::IdTool::shift(PGi, 1);
      agg_pgs.add(nb_nodes, &PGi[0]);
    }
    else // do the split starting at Nworst
    {
      // Find element and local id of the concavity
      E_Int Nim1 = PGi[(iworst > 0) ? iworst-1 : nb_nodes - 1];
      E_Int Ni   = PGi[iworst];
      
      E_Int K0, n0(E_IDX_NONE);
      //find element info for iworst : K0, n0 : WARNING : we assume that all connectMi T3 are consistently oriented between them and with PGi
      for (E_Int k = 0; k < connectMi.cols(); ++k)
      {
        const E_Int* pK = connectMi.col(k);
        for (E_Int n = 0; n < 3; ++n)
          if (*(pK + n) == Nim1 && *(pK + (n + 1) % 3) == Ni)
          {
            K0 = k; n0 = (n+2)%3; break;
          }
      }

      assert (n0 != E_IDX_NONE);
      if (n0 != E_IDX_NONE)//fixme : how the contrary can happen ?
      {
        K_FLD::IntArray *c1(new K_FLD::IntArray), *c2(new K_FLD::IntArray);

        err = cut_mesh_convex_line(crd, connectMi, normal, connectB, K0, n0, neighbors, *c1, *c2);
        if (!err)
        {
          pool.push_back(c1); pool.push_back(c2);
        }
        trash.push_back(c1);
        trash.push_back(c2);
      }
    }
  }

  //cleaning
  for (size_t i = 0; i < trash.size(); ++i)
    if (trash[i])delete trash[i];

  return err;
}

///
E_Int K_CONNECT::MeshTool::cut_mesh_convex_line
(const K_FLD::FloatArray& crd, const K_FLD::IntArray& connectT3, const E_Float* normal, const K_FLD::IntArray& connectB,
E_Int K0, E_Int n0, K_FLD::IntArray& neighbors,
K_FLD::IntArray& connectT31, K_FLD::IntArray& connectT32)
{
  assert(K0 != E_IDX_NONE);

  std::set<E_Int> bnodes;
  bnodes.clear();
  for (E_Int i = 0; i < connectB.cols(); ++i) bnodes.insert(connectB(0, i));

  E_Float Eip1Aip1[3], Ek[3], Dir[3], ps;
  E_Int Ai = 0, Aip1, Ki(K0), Kip1, Kim1(K0), ni((n0 + 1) % 3), nim1 = 0, Ei(connectT3((n0 + 1) % 3, K0)), Eip1(connectT3((n0 + 2) % 3, K0));
  K_FLD::IntArray::const_iterator pKi;

  //
  while (1)
  {
    K_FUNC::diff<3>(crd.col(Eip1), crd.col(Ei), &Dir[0]);
    K_FUNC::normalize<3>(Dir);

    while (1)
    {
      pKi = connectT3.col(Ki);
      Aip1 = *(pKi + (ni + 2) % 3);
      Kip1 = neighbors(ni, Ki);

      if (Kip1 == E_IDX_NONE) break;

      K_FUNC::diff<3>(crd.col(Aip1), crd.col(Eip1), &Eip1Aip1[0]);
      K_FUNC::normalize<3>(Eip1Aip1);
      K_FUNC::crossProduct<3>(Dir, Eip1Aip1, Ek);
      ps = K_FUNC::dot<3>(Ek, normal);

      if (ps < -E_EPSILON) break;

      // current best
      Kim1 = Ki;
      Ai = Aip1;
      nim1 = ni;

      //nex step
      Ki = Kip1;
      pKi = connectT3.col(Ki);
      ni = K_MESH::Triangle::getLocalNodeId(pKi, Aip1);
    };

    // cut edge : (Eip1, Ai)

    //do the cut for this edge by updating neighbors.
    neighbors(nim1, Kim1) = E_IDX_NONE;
    neighbors((ni + 2) % 3, Ki) = E_IDX_NONE;
    /*if (Kip1 != E_IDX_NONE)
    {
    pKi = connectT3.col(Kip1);
    ni = K_MESH::Triangle::getLocalNodeId(pKi, Ai);
    neighbors((ni+2)%3, Kip1) = E_IDX_NONE;
    }*/

    if (bnodes.find(Ai) != bnodes.end()) break;

    Ei = Eip1;
    Eip1 = Ai;
    Ki = Kim1;
    ni = (nim1 + 1) % 3; //new Ei
  };

  Vector_t<E_Int> colors;
  K_CONNECT::EltAlgo<K_MESH::Triangle>::coloring(neighbors, colors);//hpc

  E_Int colmax = *std::max_element(colors.begin(), colors.end());
  if (colmax != 1)
  {
#ifdef DEBUG_BOOLEAN
    //MIO::write("part.plt", _coord, connectT3, "TRI");
#endif
    return 1;
  }

  for (size_t i = 0; i < colors.size(); ++i)
  {
    if (colors[i] == 0)
      connectT31.pushBack(connectT3.col(i), connectT3.col(i) + 3);
    else
    {
      assert(colors[i] == 1); //we must have 2 parts only.
      connectT32.pushBack(connectT3.col(i), connectT3.col(i) + 3);
    }
  }

  return 0;
}

E_Int K_CONNECT::MeshTool::starify_from_node
(const K_FLD::FloatArray& crd, const E_Int* nodes, E_Int nb_nodes, E_Int index_start, E_Int istar, K_FLD::IntArray& connectT3, K_FLD::IntArray& neighbors, const E_Float* normal, E_Float convexity_tol)
{
  E_Float tol2(1.e-4);
  E_Int Nw = nodes[istar] - index_start; // worst concavity
  //get its ancestors
  std::vector<E_Int> local_id(connectT3.cols(), E_IDX_NONE); //tell if it's an ancestor and Nw local node id in it
  K_FLD::IntArray::iterator pS, pSn;
  
  for (E_Int S=0; S < connectT3.cols(); ++S)
  {
    pS = connectT3.col(S);

    if (*pS == Nw)
      local_id[S]=0;
    else if (*(pS+1) == Nw)
      local_id[S]=1;
    else if (*(pS+2) == Nw)
      local_id[S]=2;
  }
  
  K_SEARCH::BoundingBox<3> box;
  
  for (E_Int S=0; S < connectT3.cols(); ++S)
  {
    const E_Int& b = local_id[S];
    if (b == E_IDX_NONE) // this element doesn't have Nw as a node
      continue;
    
    E_Int Sn = neighbors(b, S); 
    
    if (Sn == E_IDX_NONE)
      continue;

    pS = connectT3.col(S);
    E_Int Ni = *(pS + b);
    E_Int Nj = *(pS + (b+1) % 3);
    E_Int Nl = *(pS + (b+2) % 3);
    
    pSn = connectT3.col(Sn);
    E_Int bn = element_type::getOppLocalNodeId(S, b, connectT3, neighbors);

    E_Int Nk = *(pSn + bn);
    
    // Lref for Adim
    box.compute(crd, nodes, nb_nodes, 1);
    E_Float Lref2 = std::max(box.maxB[0]-box.minB[0], box.maxB[1]-box.minB[1]);
    Lref2 = std::max(box.maxB[2]-box.minB[2], Lref2);
    Lref2 *= Lref2;
    
    // Convexity test : skip if the pair of element (S, Sn) is not convex (i.e. edge swapping is not doable).
    // Test done by visibility criterion.
    E_Float s1 = element_type::surface(crd.col(Ni), crd.col(Nj), crd.col(Nk), 3);
    E_Float s2 = element_type::surface(crd.col(Ni), crd.col(Nk), crd.col(Nl), 3);
    bool isConvex  = (s1 > tol2*Lref2) && (s2 > tol2*Lref2);
    if (!isConvex)
      continue;

    // Neighbors to modify
    E_Int S1 = neighbors((b+1) % 3, S);
    E_Int S2 = neighbors((bn+1) % 3, Sn);
    E_Int b1 = element_type::getOppLocalNodeId(S, (b+1) % 3, connectT3, neighbors);
    E_Int b2 = element_type::getOppLocalNodeId(Sn, (bn+1) % 3, connectT3, neighbors);

    // Update elements S and Sn (connect)
    *(pS  + (b+2)  % 3) = Nk;
    *(pSn + (bn+2) % 3) = Ni;

    // Update the neighboring (neighbors)
    neighbors((b+1)  % 3, S)  = Sn;
    neighbors((bn+1) % 3, Sn) = S;
    if ((S1 != E_IDX_NONE) && (b1 != E_IDX_NONE))
      neighbors(b1, S1)              = Sn;
    neighbors(bn, Sn)                = S1;
    if ((S2 != E_IDX_NONE) && (b2 != E_IDX_NONE))
      neighbors(b2, S2)              = S;
    neighbors(b, S)                  = S2;
  }
  
  return 0;
}

E_Int 
K_CONNECT::MeshTool::get_edges_lying_on_plane(const K_FLD::FloatArray& crd, E_Int index_start, const std::set<K_MESH::NO_Edge>& edges, E_Int Np, const E_Float* normal, E_Float abstol, std::set<K_MESH::NO_Edge>& lyingEs)
{
  typedef std::set<K_MESH::NO_Edge> eset_t;
  eset_t::const_iterator itE;
  E_Float Eij[3], u[3], abstol2(abstol*abstol);

  lyingEs.clear();

  for (itE = edges.begin(); itE != edges.end(); ++itE)
  {
    E_Int Ni = itE->node(0) - index_start;
    E_Int Nj = itE->node(1) - index_start;
    
    K_FUNC::diff<3>(crd.col(Nj), crd.col(Ni), Eij);
    K_FUNC::normalize<3>(Eij);
    
    bool quasi_aligned = (::fabs(K_FUNC::dot<3>(Eij, normal)) < 0.25);
    if (!quasi_aligned)
      continue;
    //E_Float L2 = K_FUNC::sqrNorm<3>(Eij);
    //if (L2 < E_EPSILON) L2 = 1.;

    // is Ni close enoough to the plane ?
    K_FUNC::diff<3>(crd.col(Ni), crd.col(Np), u);

    E_Float h2 = K_FUNC::dot<3>(u, normal);
    h2 *= h2;
    //E_Float seuil = L2*tol_rel*tol_rel;
    //if (h2 > L2*tol_rel*tol_rel)
    //continue;
    if (h2 > abstol2)
      continue;

    // is Nj close enoough to the plane ?
    K_FUNC::diff<3>(crd.col(Nj), crd.col(Np), u);

    h2 = K_FUNC::dot<3>(u, normal);
    h2 *= h2;

    //if ((h2) > L2*tol_rel*tol_rel)
      //continue;
    if (h2 > abstol2)
      continue;

    lyingEs.insert(*itE);
  }
  
  // ctach any missing edge : thos for which bot ends are part of edges considered as lying
  std::set <E_Int> lyingnodes;
  for (itE = lyingEs.begin(); itE != lyingEs.end(); ++itE)
  {
    lyingnodes.insert(itE->node(0));
    lyingnodes.insert(itE->node(1));
  }
  
  for (itE = edges.begin(); itE != edges.end(); ++itE)
  {
    E_Int Ni = itE->node(0);
    
    if (lyingnodes.find(Ni) == lyingnodes.end())
      continue;
    
    E_Int Nj = itE->node(1);
    
    if (lyingnodes.find(Nj) == lyingnodes.end())
      continue;
    
    lyingEs.insert(*itE); 
  }
  
  return 0;
}

void K_CONNECT::MeshTool::build_node_arity
(const std::set<K_MESH::NO_Edge>& edges, std::map<E_Int, E_Int>& node_to_count)
{
  typedef std::set<K_MESH::NO_Edge> eset_t;
  
  node_to_count.clear();
  for (eset_t::iterator it = edges.begin(); it != edges.end(); ++it)
  {
    if (node_to_count.find(it->node(0)) == node_to_count.end())
      node_to_count[it->node(0)] = 1;
    else
      ++node_to_count[it->node(0)];
    if (node_to_count.find(it->node(1)) == node_to_count.end())
      node_to_count[it->node(1)] = 1;
    else
      ++node_to_count[it->node(1)];
  }
}
///
void K_CONNECT::MeshTool::burn_free_branches(std::set<K_MESH::NO_Edge>& edges, std::map<E_Int, E_Int>& node_to_count)
{
  bool an_edge_has_been_removed = true;
  
  typedef std::set<K_MESH::NO_Edge> eset_t;
  eset_t tmp = edges;
  
  while (an_edge_has_been_removed && !edges.empty())
  {
    an_edge_has_been_removed = false;
    for (eset_t::iterator it = tmp.begin(); it != tmp.end(); ++it)
    {
      if (node_to_count[it->node(0)] == 1)
      {
        edges.erase(*it); an_edge_has_been_removed = true;
        --node_to_count[it->node(0)];
        --node_to_count[it->node(1)];
      }
      else if (node_to_count[it->node(1)] == 1)
      {
        edges.erase(*it); an_edge_has_been_removed = true;
        --node_to_count[it->node(0)];
        --node_to_count[it->node(1)];
      }
    }
  };
}

///
E_Float K_CONNECT::MeshTool::get_max_deviation
(const K_FLD::FloatArray& crd, const K_FLD::IntArray& cT3, const K_FLD::IntArray& neighT3)
{
  E_Float amax=0.;
  E_Float ni[3], nj[3];

  K_FLD::IntArray::const_iterator pKi;

  for (E_Int i = 0; i < cT3.cols(); ++i)
  {
    pKi = cT3.col(i);
    K_MESH::Triangle::normal(crd, pKi, ni);

    for (E_Int n = 0; n < 3; ++n)
    {
      E_Int e0 = *(pKi + n);
      E_Int e1 = *(pKi + (n + 1) % 3);
      const E_Float* E0 = crd.col(e0);
      const E_Float* E1 = crd.col(e1);

      E_Int j = neighT3(n, i);
      if (j == E_IDX_NONE)
        continue;

      K_MESH::Triangle::normal(crd, cT3.col(j), nj);

      E_Float alpha = K_CONNECT::GeomAlgo<K_MESH::Triangle>::angle_measure(ni, nj, E0, E1);
      amax = std::max(amax, ::fabs(K_CONST::E_PI - alpha)); // max deviation from flatness
    }
  }

  return amax;
}

///
void K_CONNECT::MeshTool::computeIncidentEdgesSqrLengths
(const K_FLD::FloatArray& crd, const ngon_unit& pgs, K_FLD::FloatArray& L)
{
  E_Int                           nb_pgs(pgs.size());
  
  L.clear();
  
  if (nb_pgs == 0) return;
  
  E_Int idmaxp1 = pgs.get_facets_max_id();
  
  L.resize(1, idmaxp1, K_CONST::E_MAX_FLOAT);  //mins
  L.resize(2, idmaxp1, -K_CONST::E_MAX_FLOAT); // maxs
  
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    const E_Int* nodes = pgs.get_facets_ptr(i);
    E_Int nb_nodes = pgs.stride(i);
    
    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      E_Int Ni = *(nodes + n) - 1;
      E_Int Nj = *(nodes + (n+1) % nb_nodes) - 1;
      E_Float l = K_FUNC::sqrDistance(crd.col(Ni), crd.col(Nj), 3);
      
      L(0, Ni) = (l < L(0, Ni)) ? l : L(0, Ni);
      L(1, Ni) = (l > L(1, Ni)) ? l : L(1, Ni);
    }
  }
}
