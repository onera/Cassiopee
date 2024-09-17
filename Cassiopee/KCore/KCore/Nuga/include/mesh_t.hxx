/*    
    Copyright 2013-2024 Onera.

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

#include "Nuga/include/DynArray.h"
#include "Nuga/include/ngon_t.hxx"
#include "Nuga/include/localizer.hxx"
# include "Nuga/include/BbTree.h"
#include "Nuga/include/polygon.hxx"
#include "Nuga/include/polyhedron.hxx"
#include "Nuga/include/metric.hxx"

#ifndef NUGA_MESH_T_HXX
#define NUGA_MESH_T_HXX

using ngon_type = ngon_t<K_FLD::IntArray>;

namespace NUGA
{

enum eGEODIM { LINEIC = 1, SURFACIC = 2, VOLUMIC = 3};

// : default impl : fixed stride (Basic element in any DIM)
template <eGEODIM GEODIM, bool fixed_stride>
struct connect_trait;


// LINEIC : BAR
template <>
struct connect_trait<LINEIC, true>
{
  using cnt_t = K_FLD::IntArray;
  using elt_t = K_MESH::Edge; using aelt_t = K_MESH::aEdge;
  using neighbor_t = cnt_t;
  using construct_elt_t = NUGA::aPolygon; // a polygon is a polyline (closed) => define a constructor

  static const E_Int index_start=0;

  static const bool BOUND_STRIDE = true;
  
  using             bound_elt_t  = int; //dummy :  need to define it to compile, but unused for BARs so put int as a dummy type;
  using             bound_trait  = int; //dummy : idem

  static E_Int ncells(const cnt_t& c) {return c.cols();}
  static void unique_indices(const cnt_t&c, std::vector<E_Int>& uinds) { c.uniqueVals(uinds);}

  static void shift(cnt_t& c, E_Int v) { c.shift(v);}
  static void append(cnt_t& c, const cnt_t& to_append){c.pushBack(to_append);}
  
  static void compress(cnt_t&c, const std::vector<bool>& keep)
  {
    std::vector<E_Int> nids;
    K_CONNECT::keep<bool> pred_keep(keep);
    K_CONNECT::IdTool::compress(c, pred_keep, nids);
  }
  static void compress(cnt_t&c, const std::vector<E_Int>& keepids, int idx_start)
  {
    K_CONNECT::IdTool::compress(c, keepids, idx_start);
  }

  static cnt_t compress_(cnt_t const& c, const std::vector<E_Int>& keepids,int idx_start)
  {
    return K_CONNECT::IdTool::compress_(c, keepids, idx_start);
  }

  static void compact(cnt_t&c, const std::vector<bool>& keep){ std::vector<E_Int> nids; K_FLD::IntArray::compact(c, keep, nids);}
  static void compact_to_used_nodes(cnt_t& c, K_FLD::FloatArray& crd, std::vector<E_Int>& nids)
  {
    NUGA::MeshTool::compact_to_mesh(crd, c, nids);
  }

  static K_FLD::FloatArray compact_to_used_nodes(cnt_t&c, K_FLD::FloatArray const & crdi, std::vector<E_Int>& nids)
  {
    std::vector<E_Int> oids;
    K_FLD::FloatArray lcrd;
    K_FLD::IntArray lcnt;
    NUGA::MeshTool::compact_to_mesh(crdi, c, lcrd, lcnt, oids); //use this version to minimize mem print
    c = std::move(lcnt);

    // build nids, sized as crdi
    K_CONNECT::IdTool::reverse_indirection(oids, nids);
    nids.resize(crdi.cols(), IDX_NONE);

    return lcrd;
  }

  static void compute_nodal_metric2(const K_FLD::FloatArray& crd, const cnt_t& cnt, std::vector<E_Float>& nodal_metric2, eMetricType metric_type)
  {
    // should be factorized in a behaviour class
    K_FLD::FloatArray L;
    NUGA::MeshTool::computeIncidentEdgesSqrLengths(crd, cnt, L);

    if (metric_type == eMetricType::ISO_MIN)
      L.extract_field(0, nodal_metric2); // 0: extract the min
    else if (metric_type == eMetricType::ISO_MAX)
      L.extract_field(1, nodal_metric2); // 0: extract the min
    else if (metric_type == eMetricType::ISO_MEAN)
    {
      nodal_metric2.resize(crd.cols(), 0.);
      for (E_Int i = 0; i < cnt.cols(); ++i)
      {
        nodal_metric2[i] = 0.5* (L(0, i) + L(1, i));
      }
    }
  }
  
  static void reverse_orient(cnt_t&c)
  {
    for (E_Int i=0; i<c.cols(); ++i)std::swap(c(0,i), c(1,i));
  }

  // polygon => closed polyline
  static void contruct_from_elt(const NUGA::aPolygon& e, K_FLD::FloatArray& crd, cnt_t& cnt)
  {
    crd = e.m_crd;
    E_Int nnodes = crd.cols();
    cnt.resize(2, nnodes);
    for (E_Int i = 0; i < nnodes; ++i)
    {
      cnt(0, i) = e.m_nodes[i];
      cnt(1, i) = e.m_nodes[(i + 1) % nnodes];
    }
  }
  
};

// SURF
template <>
struct connect_trait<SURFACIC, false>
{
  using cnt_t = ngon_unit;
  using elt_t = K_MESH::Polygon;  using aelt_t = NUGA::aPolygon;
  using neighbor_t = cnt_t;
  using construct_elt_t = NUGA::aPolyhedron<UNKNOWN>;  // a polyhedron is a surface (closed) => define a constructor

  static const eGEODIM BOUND_GEODIM = LINEIC;
  static const bool    BOUND_STRIDE = true;
 
  using                bound_elt_t = K_MESH::Edge;
  using                bound_trait = int; //dummy : need to define it to compile, but unused for BARs so put int as a dummy type;

  static const E_Int index_start=1;

  static E_Int ncells(const cnt_t& c) {return c.size();}
  static void shift(cnt_t& c, E_Int v) { c.shift(v);}
  static void unique_indices(const cnt_t& c, std::vector<E_Int>& uinds) { c.unique_indices(uinds);}

  static void append(cnt_t& c, const cnt_t& to_append){c.append(to_append);}

  template <typename ELT_t> static void add(K_FLD::FloatArray& crd, cnt_t& c, ELT_t const& e, bool capitalize_coords);
  
  static void build_neighbors(const cnt_t& c, neighbor_t& neighbors)
  {
    K_MESH::Polygon::build_pg_neighborhood(c, neighbors);
  }

  static void get_boundary(const cnt_t& c, K_FLD::IntArray& bc)
  {
    NUGA::MeshTool::getBoundary(c, bc);
    bc.shift(-1); //NGON store 1-based, IntArray is 0-based
  }
  static void get_boundary(const cnt_t& c, K_FLD::IntArray& bc, std::vector<E_Int>& ancestors)
  {
    NUGA::MeshTool::getBoundary(c, bc, &ancestors);
    bc.shift(-1); //NGON store 1-based, IntArray is 0-based
  }
  
  static void get_boundaries_by_color(const cnt_t& c, const neighbor_t& neighborz, std::map<E_Int, K_FLD::IntArray>& col_to_bound)
  {
    // WARNING : neighborz must contain some colors : one for the inner faces, and some other for the boundary
    for (E_Int i=0; i < c.size(); ++i)
    {
      E_Int nneigh = neighborz.stride(i);
      const E_Int* pneigh = neighborz.get_facets_ptr(i); //0-based
      
      E_Int nnodes = c.stride(i);
      const E_Int* pnodes = c.get_facets_ptr(i); //1-based
      
      assert (nnodes == nneigh);
      
      for (E_Int j=0; j < nneigh; ++j)
      {
        E_Int const & col = pneigh[j];
        E_Int E[] = {pnodes[j]-1, pnodes[(j+1)%nnodes]-1};
        
        col_to_bound[col].pushBack(E, E+2);
      }
    }
  }

  template <typename T> static void set_boundary_type(cnt_t const& c, T typ, const std::vector<E_Int>& ids) { /*todo ?*/ }
  
  static cnt_t compress_(cnt_t const& c, const std::vector<E_Int>& keepids, int idx_start)
  {
    std::vector<E_Int> oids;
    ngon_unit pgs;
    c.extract(keepids, pgs, oids, idx_start);
    return pgs;
  }

  static void compress(cnt_t&c, const std::vector<bool>& ikeep)
  {
    K_CONNECT::keep<bool> pred(ikeep);
    Vector_t<E_Int> pgnids, pgoids;
    ngon_unit pgs;
    c.extract_by_predicate(pred, pgs, pgoids, pgnids);
    c = pgs;
  }
  
  static void compact_to_used_nodes(cnt_t& c, K_FLD::FloatArray& crd, std::vector<E_Int>& nids)
  {
    ngon_type::compact_to_used_nodes(c, crd, nids);
  }

  static K_FLD::FloatArray compact_to_used_nodes(cnt_t&c, K_FLD::FloatArray const & crdi, std::vector<E_Int>& nids)
  {
    return ngon_type::compress_to_used_nodes(c, crdi, nids);
  }

  static void compute_nodal_metric2(const K_FLD::FloatArray& crd, const cnt_t& cnt, std::vector<E_Float>& nodal_metric2, eMetricType metric_type)
  {
    K_FLD::FloatArray L;
    NUGA::MeshTool::computeIncidentEdgesSqrLengths(crd, cnt, L);

    E_Int col{ 0 };
    if (metric_type == eMetricType::ISO_MIN) col = 0;
    else if (metric_type == eMetricType::ISO_MAX) col = 1;

    L.extract_field(col, nodal_metric2); // 0: extract the min
  }

  static void reverse_orient(cnt_t&c)
  {
    for (E_Int PGi = 0; PGi < c.size(); ++PGi)
    {
      E_Int s = c.stride(PGi);
      E_Int* p = c.get_facets_ptr(PGi);
      std::reverse(p, p + s);
    }
  }

  // polyhedron => closed surface
  static void contruct_from_elt(const NUGA::aPolyhedron<0>& e, K_FLD::FloatArray& crd, cnt_t& cnt)
  {
    crd = e.m_crd;
    cnt = e.m_pgs;
  }

  static void build_global_edge_ids(const cnt_t& cnt, cnt_t& glob_edge_ids)
  {
    std::map<K_MESH::NO_Edge, E_Int> e2id;
    std::vector<int> molecule;
    K_MESH::NO_Edge noe;
    E_Int count(0);

    E_Int npgs = cnt.size();
    for (E_Int i = 0; i < npgs; ++i)
    {
      const E_Int* nodes = cnt.get_facets_ptr(i);
      int nnodes = cnt.stride(i);

      molecule.clear();
      molecule.push_back(nnodes);

      for (int n = 0; n < nnodes; ++n)
      {
        noe.setNodes(nodes[n], nodes[(n+1)%nnodes]);
        auto ite = e2id.find(noe);

        if (ite == e2id.end())
        {
          molecule.push_back(count);
          e2id[noe] = count++;
        }
        else
          molecule.push_back(ite->second);
      }

      glob_edge_ids.add(molecule);

    }

    glob_edge_ids.updateFacets();
  }
};

// Specialisation aelt_t
template <> inline
void
connect_trait<SURFACIC, false>::add<NUGA::aPolygon>(K_FLD::FloatArray& crd, ngon_unit& c, aelt_t const& e, bool append_vertices)
{

  if (append_vertices || e.m_poids.empty()) // fixme : add a condition : or inconsistent with crd
  {
    // append the vertices

    E_Int shift = crd.cols() + index_start - e.shift();
    crd.pushBack(e.m_crd);
    c.add(e.nb_nodes(), e.begin(), shift);
  }
  else
  {
    // capitalize the nodes  to ease conformity
    // here we use the original node ids
    // and update crd

    assert(crd.rows() == 3);
    assert(e.m_poids.size() == e.m_nodes.size());

    auto nods = e.m_nodes; //lazy way of construct
    E_Int shift = index_start - e.shift();
    //
    for (size_t n = 0; n < nods.size(); ++n)
    {
      ASSERT_IN_VECRANGE(e.m_poids, nods[n]);

      E_Int oid = e.m_poids[nods[n]];

      ASSERT_IN_DYNARANGE(crd, oid);

      crd(0, oid) = e.m_crd(0, nods[n]);
      crd(1, oid) = e.m_crd(1, nods[n]);
      crd(2, oid) = e.m_crd(2, nods[n]);
      nods[n] = oid;

      assert(nods[n] != IDX_NONE);
    }

    c.add(nods.size(), &nods[0], shift);
  }
}

// VOL 
template <>
struct connect_trait<VOLUMIC, false>
{
  using cnt_t = ngon_type;
  using elt_t = K_MESH::Polyhedron<UNKNOWN>; using aelt_t = NUGA::aPolyhedron<UNKNOWN>;
  using neighbor_t = ngon_unit;

  using construct_elt_t = E_Int;//DUMMY : non-sense in volumic

  static const eGEODIM BOUND_GEODIM = SURFACIC;
  static const bool    BOUND_STRIDE = false;
 
  using                bound_elt_t = K_MESH::Polygon;
  using                bound_trait = connect_trait<SURFACIC, false>;
  
  static const E_Int index_start=1;

  static E_Int ncells(const cnt_t& c) {return c.PHs.size();}
  static void shift(cnt_t& c, E_Int v) { c.PGs.shift(v);}
  static void unique_indices(const cnt_t& c, std::vector<E_Int>& uinds) { c.PGs.unique_indices(uinds);}

  template <typename ELT_t> static void add(K_FLD::FloatArray& crd, cnt_t& c, ELT_t const& e, bool capitalize_coords/*not used yet*/)
  {
    E_Int ptshift = crd.cols();
    E_Int npgs0 = c.PGs.size();
    //E_Int nphs0 = c.PHs.size();

    crd.pushBack(e.m_crd);
    c.PGs.append(e.m_pgs);
    c.PGs.shift(ptshift, npgs0);
    c.PHs.add(e.m_faces.size(), &e.m_faces[0], npgs0);
  }

  ///
  static void build_neighbors(const cnt_t& c, neighbor_t& neighbors)
  {
    c.build_ph_neighborhood(neighbors);
  }

  ///
  static void get_boundary(cnt_t& c, ngon_unit& bc)
  {
    //for transferring type
    //auto ctype_cpy = c.PHs._type;
    auto ftype_cpy = c.PGs._type; //save it before flag_externals

    c.flag_externals(INITIAL_SKIN);
    Vector_t<E_Int> oids;
    c.PGs.extract_of_type(INITIAL_SKIN, bc, oids);

    if (ftype_cpy.empty()) return;

    bc._type.resize(oids.size());
    for (size_t i = 0; i < oids.size(); ++i)
      bc._type[i] = ftype_cpy[oids[i]];
  }

  ///
  static void get_boundary(cnt_t& c, ngon_unit& bc, std::vector<E_Int>& oids, std::vector<E_Int>& ancestors)
  {
    oids.clear();

    c.flag_externals(INITIAL_SKIN);
    c.PGs.extract_of_type(INITIAL_SKIN, bc, oids);

    ancestors.clear();
    ancestors.resize(c.PGs.size(), IDX_NONE);
    for (E_Int i = 0; i < c.PHs.size(); ++i)
    {
      E_Int s = c.PHs.stride(i);
      const E_Int* f = c.PHs.get_facets_ptr(i);

      for (E_Int j = 0; j < s; ++j)
        ancestors[f[j] - 1] = i;
    }

    std::vector<E_Int> new_anc(bc.size(), IDX_NONE);
    for (size_t i = 0; i < new_anc.size(); ++i) new_anc[i] = ancestors[oids[i]];
    ancestors = new_anc;
  }

  ///
  static void get_boundary(cnt_t& c, ngon_unit& bc, std::vector<E_Int>& ancestors)
  {
    std::vector<E_Int> oids;
    get_boundary(c, bc, oids, ancestors);
  }

  ///
  static E_Int get_nth_neighborhood (cnt_t const& c, neighbor_t const& neighbors, E_Int N, const std::vector<E_Int>& iList, std::vector<E_Int>& oList)
  {
    oList.clear();

    N = std::max(E_Int(0), N);
    if (N == 0) return 0;

    std::set<E_Int> front, iset, oset;

    iset.insert(ALL(iList));

    do 
    {
      for (auto PHi : iset)
      {
        E_Int nneighs = neighbors.stride(PHi);
        const E_Int* pneighs = neighbors.get_facets_ptr(PHi);

        for (E_Int n = 0; n < nneighs; ++n)
        {
          if (pneighs[n] == IDX_NONE) continue;
          front.insert(pneighs[n]);
        }
      }

      for (auto e : iset)
        front.erase(e);

      oset.insert(ALL(front));
      iset = front;
      front.clear();
    }
    while (--N > 0);

    for (size_t i = 0; i < iList.size(); ++i)
      oset.erase(iList[i]);

    oList.insert(oList.end(), ALL(oset));

    return 0;

  }

  //fixme : slow
  /*static void get_nth_neighborhood
  (cnt_t const& c, neighbor_t const& neighbors, E_Int i, E_Int N, std::set<E_Int>& neighs)
  {
    neighs.clear();

    // assumption : N >= 1
    if (N < 1) return;

    E_Int nneighs = neighbors.stride(i);
    const E_Int* pneighs = neighbors.get_facets_ptr(i);

    std::set<E_Int> processed, pool, pool2;
    
    pool.insert(i);

    while (N-- > 0)
    {
      pool2 = pool;
      for (auto e : pool)
      {
        processed.insert(e);
        pool2.erase(e);
        neighs.insert(e);

        E_Int nneighs = neighbors.stride(e);
        const E_Int* pneighs = neighbors.get_facets_ptr(e);

        for (E_Int n = 0; n < nneighs; ++n)
        {
          if (pneighs[n] == IDX_NONE) continue;

          if (processed.find(pneighs[n]) == processed.end())
          {
            pool2.insert(pneighs[n]);
            neighs.insert(pneighs[n]);
          }
        }
      }
      pool = pool2;
    }
    neighs.erase(i);
  }*/

  template <typename T> static void set_boundary_type(cnt_t const& c, T typ, const std::vector<E_Int>& ids)
  {
    c.PGs._type.clear();
    c.PGs._type.resize(c.PGs.size(), 0/*NUGA::ANY*/);
    for (size_t u = 0; u < ids.size(); ++u)
      c.PGs._type[ids[u]] = E_Int(typ);

    //std::cout << "wall ids list size : " << ids.size() << std::endl;
  }

  static cnt_t compress_(cnt_t const& c, const std::vector<E_Int>& keepids, int idx_start)
  {
    std::vector<E_Int> oids;
    cnt_t ng;
    c.PHs.extract(keepids, ng.PHs, oids, idx_start);

    ng.PGs = c.PGs;

    Vector_t<E_Int> pgnids, phnids;
    ng.remove_unreferenced_pgs(pgnids, phnids);

    return ng;
  }

  static void compress(cnt_t&c, const std::vector<bool>& ikeep)
  {
    K_CONNECT::keep<bool> pred(ikeep);
    Vector_t<E_Int> pgnids, pgoids;
    ngon_unit phs;
    c.PHs.extract_by_predicate(pred, phs, pgoids, pgnids);
    c.PHs = phs;

    Vector_t<E_Int> phnids;
    c.remove_unreferenced_pgs(pgnids, phnids);
  }

  static void compress(ngon_unit& neigh, const std::vector<bool>& ikeep)
  {
    connect_trait<SURFACIC, false>::compress(neigh, ikeep);
  }

  static void compact_to_used_nodes(cnt_t& c, K_FLD::FloatArray& crd, std::vector<E_Int>& nids)
  {
    nids.clear();
    ngon_type::compact_to_used_nodes(c.PGs, crd);
  }

  static K_FLD::FloatArray compact_to_used_nodes(cnt_t&c, K_FLD::FloatArray const & crdi, std::vector<E_Int>& nids)
  {
    return ngon_type::compress_to_used_nodes(c.PGs, crdi, nids);
  }

  static void compute_nodal_metric2(const K_FLD::FloatArray& crd, const cnt_t& c, std::vector<E_Float>& nodal_metric2, eMetricType metric_type)
  {
    // should be factorized in a behaviour class
    K_FLD::FloatArray L;
    NUGA::MeshTool::computeIncidentEdgesSqrLengths(crd, c.PGs, L);

    E_Int col{ 0 };
    if (metric_type == eMetricType::ISO_MIN) col = 0;
    else if (metric_type == eMetricType::ISO_MAX) col = 1;

    L.extract_field(col, nodal_metric2); // 0: extract the min
  }
};


///
template <eGEODIM GEODIM, bool FIXSTRIDE>
struct mesh_t
{
  static const  eGEODIM ARG1 = GEODIM;
  //static const  bool ARG2 = FIXSTRIDE;

  using trait = connect_trait<GEODIM, FIXSTRIDE>;

  static const E_Int BOUND_STRIDE = trait::BOUND_STRIDE;

  using bound_mesh_t = mesh_t<eGEODIM(GEODIM - 1), BOUND_STRIDE>; // recursive definition! => require dummy terminal struct : mesh_t<0,true>
  using bound_elt_t  = typename bound_mesh_t::elt_t;
  using bound_trait  = typename bound_mesh_t::trait;  
 
  static const E_Int index_start  = trait::index_start;

  using cnt_t           = typename trait::cnt_t;
  
  using elt_t           = typename trait::elt_t;
  using aelt_t          = typename trait::aelt_t;
  using construct_elt_t = typename trait::construct_elt_t;
  
  using loc_t = localizer<K_SEARCH::BbTree3D>;
  using neighbor_t = typename trait::neighbor_t;
  
  K_FLD::FloatArray   crd;
  cnt_t               cnt;
  mutable std::vector<E_Int>   e_type, flag;
  mutable std::vector<E_Float> nodal_metric2; // square dist
  mutable loc_t*               localiz;
  mutable neighbor_t*          neighbors;
  int                          oriented;
  mutable eMetricType          metric_type;

  // CONSTRUCTORS / DESTRUCTOR //

  mesh_t():localiz(nullptr), neighbors(nullptr), oriented(0), metric_type(eMetricType::ISO_MIN){}

  mesh_t(const K_FLD::FloatArray &crd, const K_FLD::IntArray& cnt, int orient=0):crd(crd),cnt(cnt), localiz(nullptr), neighbors(nullptr), oriented(orient), metric_type(eMetricType::ISO_MIN) {}
  
  mesh_t& operator=(const mesh_t&m)
  {
    crd = m.crd; cnt = m.cnt;
    oriented = m.oriented;
    
    // these attribute are rferring to a previous state
    if (localiz != nullptr) delete localiz;
    if (neighbors != nullptr) delete neighbors;
    localiz   = nullptr;
    neighbors = nullptr;

    if (m.neighbors != nullptr)
    {
      neighbors  = new neighbor_t;
      *neighbors = *m.neighbors;
    }
    
    e_type = m.e_type;
    flag = m.flag;
    nodal_metric2 = m.nodal_metric2;
    metric_type = m.metric_type;

    return *this;
  }

  ///
  mesh_t& operator=(const mesh_t &&m)
  {
    crd = std::move(m.crd); //call to DynArray::operator=(DynArray&&)
    cnt = std::move(m.cnt);
    oriented = m.oriented;

    // these attribute are rferring to a previous state
    if (localiz != nullptr) delete localiz;
    if (neighbors != nullptr) delete neighbors;
    // steal new ones
    localiz = m.localiz;
    neighbors = m.neighbors;
    
    e_type = std::move(m.e_type);
    flag = std::move(m.flag);
    nodal_metric2 = std::move(m.nodal_metric2);
    metric_type = m.metric_type;

    m.localiz = nullptr;
    m.neighbors = nullptr;

    return *this;
  }
  
  mesh_t(const mesh_t& m) :localiz(nullptr), neighbors(nullptr) { *this = m; }
  mesh_t(const mesh_t&& m):localiz(nullptr), neighbors(nullptr) { *this = m; }

  mesh_t(const mesh_t& m, std::vector<E_Int>& ids, int idx_start) :cnt(trait::compress_(m.cnt, ids, idx_start)), localiz(nullptr), neighbors(nullptr)
  {
    // fixme : hpc issue for more than lineic
    std::vector<E_Int> nids;
    crd = trait::compact_to_used_nodes(cnt, m.crd, nids);
    oriented = m.oriented;

    nodal_metric2 = K_CONNECT::IdTool::compact_(m.nodal_metric2, nids); //sync the metric
    metric_type = m.metric_type;

    //sync the neighborhood
    if (m.neighbors != nullptr)
    {
      neighbors = new neighbor_t;
      *neighbors = trait::compress_(*m.neighbors, ids, idx_start);
    }
  }

  // to create a bound mesh from a "parent" mesh
  template <bool USTRIDE>
  mesh_t(const mesh_t<eGEODIM(GEODIM+1), USTRIDE>& parent_mesh): 
    crd(crd), localiz(nullptr), neighbors(nullptr), oriented(parent_mesh.oriented)
  {
    parent_mesh.template get_boundary<FIXSTRIDE>(*this);
  }

  mesh_t(const construct_elt_t& e):localiz(nullptr), neighbors(nullptr), oriented(0), metric_type(eMetricType::ISO_MIN)
  {
    trait::contruct_from_elt(e, crd, cnt);// polygon/polyhedron => polyline/surface 
  }

  ~mesh_t()
  {
    if (localiz != nullptr) { delete localiz; localiz = nullptr; }
    if (neighbors != nullptr) { delete neighbors; neighbors = nullptr; }
  }

  // end constructors / destructor //

  // ACCESSORS //

  elt_t element(E_Int i) const { return elt_t(cnt,i);}
  aelt_t aelement(E_Int i) const 
  {
    elt_t e(cnt, i);

    // Lref2
    // if nodal_metric2 exist and is 'valid' (same size as crd) => use it
    // otherwise compute Lref2 based only on min edge length of the element (so might by over estimated)
    E_Float Lr2 = (E_Int(nodal_metric2.size()) == crd.cols()) ? e.Lref2(nodal_metric2) : e.Lref2(crd);
    
    return aelt_t(e, crd, Lr2);
  }
    
  void get_boundary(E_Int i, int j, bound_elt_t& b) const
  {
    elt_t e(cnt, i);
    e.getBoundary(j, b);
  }
  
  E_Int ncells() const {return trait::ncells(cnt);}

  // end accessors

  //
  void get_boundary(bound_mesh_t& bound_mesh) 
  {
    trait::get_boundary(cnt, bound_mesh.cnt);
    bound_mesh.crd = crd;
    bound_mesh.oriented = oriented;
    //compact crd to boundary only
    std::vector<E_Int> nids;
    bound_trait::compact_to_used_nodes(bound_mesh.cnt, bound_mesh.crd, nids);
    //tranfer metric if any
    if (!nodal_metric2.empty())
    {
      bound_mesh.nodal_metric2.resize(bound_mesh.crd.cols(), NUGA::FLOAT_MAX);
      for (size_t i = 0; i < nids.size(); ++i)
      {
        if (nids[i] != IDX_NONE)
          bound_mesh.nodal_metric2[nids[i]] = nodal_metric2[i];
      }
    }
  }

  ///
  void get_boundary(bound_mesh_t& bound_mesh, std::vector<E_Int>& ancestors) 
  {
    trait::get_boundary(cnt, bound_mesh.cnt, ancestors);
    bound_mesh.crd = crd;
    bound_mesh.oriented = oriented;
    //compact crd to boundary only
    std::vector<E_Int> nids;
    bound_trait::compact_to_used_nodes(bound_mesh.cnt, bound_mesh.crd, nids);
  }

  ///
  void get_boundary(bound_mesh_t& bound_mesh, std::vector<E_Int>& oids, std::vector<E_Int>& ancestors)
  {
    trait::get_boundary(cnt, bound_mesh.cnt, oids, ancestors);
    bound_mesh.crd = crd;
    bound_mesh.oriented = oriented;
    //compact crd to boundary only
    std::vector<E_Int> nids;
    bound_trait::compact_to_used_nodes(bound_mesh.cnt, bound_mesh.crd, nids);
  }
  
  ///
  void get_boundaries_by_color(std::map<E_Int, bound_mesh_t>& col_to_boundmesh) const
  {
    using bcnt_t = typename bound_mesh_t::cnt_t;
    
    const neighbor_t* neighborz = get_neighbors(); // should be computed and set with colors before to have more than IDX_NONE
    
    std::map<E_Int, bcnt_t> col_to_boundcnt;
    trait::get_boundaries_by_color(cnt, *neighborz, col_to_boundcnt);
    
    for (auto it : col_to_boundcnt)
    {
      bcnt_t& bcn = it.second;
      bound_mesh_t& bmesh = col_to_boundmesh[it.first] = bound_mesh_t(crd, bcn);
      
      std::vector<E_Int> nids;
      bound_trait::compact_to_used_nodes(bmesh.cnt, bmesh.crd, nids);
    }
    
  }
  
  /// 
  void reverse_orient()
  {
    trait::reverse_orient(cnt);
    if (oriented != 0) oriented = -oriented;
  }
  
  ///
  void append(const mesh_t& m)
  {
    cnt_t mcnt = m.cnt;//fixme copy just because of the non-constness of the shift : shift should also have an option to start from an inex to make it at the end
    trait::shift(mcnt, crd.cols());
    
    crd.pushBack(m.crd);
    trait::append(cnt, mcnt);

    e_type.insert(e_type.end(), ALL(m.e_type));
  }

  ///
  void add(aelt_t const & e, bool capitalize_coords)
  {
    trait::add(crd, cnt, e, capitalize_coords);
  }

  ///
  void clear()
  {
    crd.clear();
    cnt.clear();
  }

  ///
  template <typename T>
  void set_type(T typ, const std::vector<E_Int>& ids)
  {
    //std::cout << "set_type : 1" << std::endl;
    E_Int nbcells = ncells();
    e_type.resize(nbcells, IDX_NONE);
    //std::cout << "set_type : nb cell : " << nbcells << std::endl;
    for (size_t i=0; i < ids.size(); ++i)
    {
      //std::cout << "id/ncell : " << ids[i] << "/" << nbcells << std::endl;
      if (ids[i] >= nbcells) continue; 
      e_type[ids[i]] = (E_Int)typ;
    }
    //std::cout << "set_type : 3" << std::endl;
  }

  template <typename T>
  void set_boundary_type(T typ, const std::vector<E_Int>& ids)
  {
    trait::set_boundary_type(cnt, typ, ids);
  }
  
  void set_flag(E_Int i, int val) const // fixme : constness here is bad design
  {
    if (i >= flag.size()) flag.resize(ncells(), IDX_NONE);
    flag[i] = val;
  }

  void build_global_edge_ids(ngon_unit& glob_edge_ids) const
  {
    trait::build_global_edge_ids(cnt, glob_edge_ids);
  }

  void bbox(K_SEARCH::BBox3D& box) const 
  {
    std::vector<E_Int> uinds;
    trait::unique_indices(cnt, uinds);
    box.compute(crd, uinds, index_start);
  }
  
  void bbox(E_Int i, K_SEARCH::BBox3D& box) const 
  {
    assert(i < ncells());

    elt_t e = element(i);
    e.bbox(crd, box);
  }

  double Lref2() const
  {
    if (! nodal_metric2.empty())
      return *std::min_element(ALL(nodal_metric2));

    // Lref2 est calcule uniquement avec les aretes de la cellule => Possiblement surestime 
    // car des cellules voisines pourraient avoir des aretes plus petites et donc definir une valeur nodale plus petite
    double minLref2 = NUGA::FLOAT_MAX;
    for (E_Int i=0; i < ncells(); ++i)
    {
      elt_t e =element(i);
      minLref2 = std::min(e.Lref2(crd), minLref2);
      //std::cout << "minLref2/i/ncel" << minLref2 << "/" << i << "/" << ncells() << std::endl;
    }
    return minLref2;
  }
  
  double Lref2(E_Int i) const
  { 
    elt_t e =element(i);
    
    if (! nodal_metric2.empty())
      return e.Lref2(nodal_metric2); // taking into acount surrounding metric field (based on neighbors)
    // otherwise only based on this element
    return e.Lref2(crd);               // isolated val
  }

  double Lref2(const std::vector<E_Int>& cands, E_Int idx_start) const
  {
    double val = NUGA::FLOAT_MAX;
    for (size_t i = 0; i < cands.size(); ++i)
    {
      val = std::min(val, Lref2(cands[i] - idx_start));
    }
    return val;
  }
  
  // WARNING : if T=bool, use predicate compress, if T=int, use keep as a list of indices to keep
  template <typename T>
  void compress(const std::vector<T>& keep)
  {
    std::vector<E_Int> ptnids;
    compress(keep, ptnids);
  }

  template <typename T>
  void compress(const std::vector<T>& keep, std::vector<E_Int>& ptnids)
  {
    ptnids.clear();

    trait::compress(cnt, keep);
    trait::compact_to_used_nodes(cnt, crd, ptnids);

    if (ptnids.empty())
      nodal_metric2.clear();
    else
      K_CONNECT::IdTool::compact(nodal_metric2, ptnids); //sync the metric 

    if (neighbors != nullptr)
      trait::compress(*neighbors, keep); //sync the neighborhood

    if (e_type.size() == keep.size())
      K_CONNECT::IdTool::compact(e_type, keep);
    else e_type.clear();

    if (flag.size() == keep.size())
      K_CONNECT::IdTool::compact(flag, keep);
    else flag.clear();
  }

  void build_localizer() const 
  {
    if (ncells() == 0) return;

    if (localiz != nullptr) delete localiz;

    using tree_t = K_SEARCH::BbTree3D;

    tree_t* t = new tree_t(crd, cnt);
    localiz = new loc_t(t, EPSILON/*fixme*/);   //not a leak : mask_loc owes t
  }

  const loc_t* get_localizer(void) const
  {
    if (localiz == nullptr) build_localizer();
    return localiz; //can be null if cnt is empty
  }

  void build_neighbors() const
  {
    if (ncells() == 0) return;
    if (neighbors != nullptr) delete neighbors;

    neighbors = new neighbor_t;
    trait::build_neighbors(cnt, *neighbors);
  }

  const neighbor_t* get_neighbors() const
  {
    if (neighbors != nullptr) { delete neighbors; neighbors = nullptr; }
    if (neighbors == nullptr) build_neighbors();
    return neighbors;
  }
  neighbor_t* get_neighbors()
  {
    if (neighbors == nullptr) build_neighbors();
    return neighbors;
  }

  void get_nth_neighborhood(E_Int N, std::vector<E_Int>& iList, std::vector<E_Int>& oList) const
  {
    get_neighbors();//lazy build
    trait::get_nth_neighborhood(cnt, *neighbors, N, iList, oList);
  }

  /* fixme : slow impl : to review
  void get_nth_neighborhood(E_Int i, E_Int N, std::set<E_Int>& neighs) const
  {
    get_neighbors();//lazy build
    trait::get_nth_neighborhood(cnt, *neighbors, i, N, neighs);
  }*/

  const std::vector<E_Float>& get_nodal_metric2(eMetricType mtype = eMetricType::ISO_MIN, bool based_on_cnt_only=false) const // based_on_cnt_only is only relevant woth ISO_MIN (cf. build_nodal_metric2)
  {
    if ( (metric_type != mtype) || nodal_metric2.empty() || ((E_Int)nodal_metric2.size() != crd.cols()))
       build_nodal_metric2(mtype, based_on_cnt_only); // already computed and correctly sized
    return nodal_metric2;
  }

  void build_nodal_metric2(eMetricType mtype, bool based_on_cnt_only=false) const
  {
    metric_type = mtype;
    trait::compute_nodal_metric2(crd, cnt, nodal_metric2, metric_type);

    if (mtype == ISO_MIN && !based_on_cnt_only)
    {
      // check if close point other than cnt, specially with folded surfaces. WARNING : better to use on watertight mesh i.e without doublons
      using acrd_t = K_FLD::ArrayAccessor<K_FLD::FloatArray>;
      acrd_t acrd(crd);
      K_SEARCH::KdTree<> tree(acrd);

      double d2;
      for (size_t i = 0; i < nodal_metric2.size(); ++i)
      {
        double r2 = (1. - EPSILON) * nodal_metric2[i]; //reduce it to discard nodes connected to i.
        E_Int N = tree.getClosest(i, r2, d2);
        if (N != IDX_NONE)nodal_metric2[i] = d2;
      }
    }
  }

};

template <>
struct mesh_t<NUGA::eGEODIM(0), true>
{
  using elt_t = int; //dummy
  using trait = int; //dummy
};

///
template <typename MT>
using boundary_t = mesh_t<MT::trait::BOUND_GEODIM, MT::trait::BOUND_STRIDE>;

using ph_mesh_t = mesh_t<VOLUMIC, false>;
using pg_smesh_t  = mesh_t<SURFACIC,false>;
//using pure_basic_smesh_t = mesh_t<25, true>; 
using edge_mesh_t = mesh_t<LINEIC,true>; // == boundary_t< mesh_t<25, FIXSTRIDE>> qqsoit FIXSTRIDE

}

#endif
