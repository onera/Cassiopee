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

#ifndef NUGA_SPLITTER_H
#define NUGA_SPLITTER_H

#include "Nuga/include/DynArray.h"
#include "Nuga/include/ngon_t.hxx"
#include <vector>
#include <set>
#include <map>
#include <deque>
#include "Nuga/include/defs.h"


namespace NUGA
{
  // For applying transformations
  struct transfo_t
  {
    enum ePHset { CONCAVES = 0, NONSTARS };
    enum eSplitPolicy {CONVEXIFY = 0, STARIFY_CONCAVES, STARIFY_FROM_CHAIN_ENDS};

    transfo_t() :convexity_tol(CONVEXITY_TOL){}

    E_Float concavity_threshold; //between faces
    E_Float convexity_threshold; //between faces
    E_Float convexity_tol;       //between edges

    bool improve_qual;

    std::set<E_Int> nodes;       //optional, used when starifying at specifed nodes

  };

  typedef E_Int(*transfo_func)(const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& tri_pgs, Vector_t<E_Int>& oids, const transfo_t& params, const Vector_t<bool>* to_process);



  class Splitter
  {

#ifdef DEBUG_PHSPLITTER
    static const E_Int faultyPH = 1;
#endif

  public:
    typedef ngon_t<K_FLD::IntArray>                 ngon_type;
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    typedef std::set<K_MESH::NO_Edge>               eset_t;
    typedef std::set<E_Int>                         iset_t;
    typedef std::vector<E_Int>                      ivec_t;
    typedef std::map<E_Int, iset_t >                id_to_ids_t;
    typedef std::map<E_Int, E_Int>                  id_to_id_t;
    typedef std::vector<std::deque<E_Int> >         idlists_t;
    typedef std::map<K_MESH::NO_Edge, E_Float>      edge_angle_t;
  
  public: //PG SPLITTING

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// POLYGON SPLIT WRAPPERS /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    inline static E_Int split_pgs
      (transfo_func func, const K_FLD::FloatArray& crd, const ngon_type& ngi, E_Float convexity_tol, ngon_type& ngo, const Vector_t<E_Int>* PHlist = 0);
    ///
    template <typename TriangulatorType>
    inline static E_Int triangulate_external_pgs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Int in_or_out, const transfo_t& qual_param, ngon_type& ngo); //in_or_out : 0(internal only), 1(extrnal only), 2(both)
    ///
    template <typename TriangulatorType>
    inline static E_Int triangulate_specified_pgs(const K_FLD::FloatArray& crd, ngon_type& ngi, const E_Int* PGlist, E_Int sz, const transfo_t& qual_param, ngon_type& ngo); 
    ///
    template <typename TriangulatorType>
    inline static E_Int starify_pgs_at_chain_nodes
      (const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orient, E_Float concave_threshold, E_Float convex_threshold, E_Float convexity_tol, ngon_type& ngo, const Vector_t<E_Int>* PHlist = 0);
    ///
    template <typename TriangulatorType>
    inline static E_Int prepareCellsSplit(const K_FLD::FloatArray& crd, ngon_type& ngi,
      transfo_t::ePHset phset, transfo_t::eSplitPolicy policy, E_Float PH_conc_threshold, E_Float PH_cvx_threshold,
      E_Float PG_cvx_threshold, ngon_type& ngo);
 
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // SPLIT ACTIVE FUNCTIONS (fitting expected signature for split_pgs) ////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Triangulates specified PGs
    template<typename TriangulatorType>
    inline static E_Int triangulate_pgs
      (const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& tri_pgs, Vector_t<E_Int>& oids, const transfo_t& qual_param, const Vector_t<bool>* to_process = 0);
    ///
    template <typename TriangulatorType>
    inline static E_Int convexify_pgs
      (const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& convex_pgs, Vector_t<E_Int>& colors, const transfo_t& angular_params, const Vector_t<bool>* process = 0);
    ///
    template <typename TriangulatorType>
    inline static E_Int starify_pgs_on_reflex
      (const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& pgso, Vector_t<E_Int>& oids,
      const transfo_t& params, const Vector_t<bool>* process = 0);
    ///
    template <typename TriangulatorType>
    inline static E_Int starify_pgs_at_specified_nodes
      (const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& pgso, Vector_t<E_Int>& oids,
      const transfo_t& params, const Vector_t<bool>* process = 0);

  public: // PH SPLITTING
    
     //
    template<typename TriangulatorType>
    inline static E_Int split_non_star_phs
    (K_FLD::FloatArray& crd, ngon_type& ngi, E_Float concave_threshold, E_Float convex_threshold, E_Float RTOL, E_Float GRmin, E_Float Fluxmax, ngon_type& ngo)
    {
      // Create orientation info
      ngon_unit orient;
      E_Int err = ngon_type::build_orientation_ngu<TriangulatorType>(crd, ngi, orient);
      if (err) return 0;
  
      std::vector<ngon_type::ePathoPH> PHtypes;
      ngon_type::detect_pathological_PHs<TriangulatorType>(crd, ngi, orient, PHtypes);
      
      std::vector<E_Int> PHlist;
      for (size_t i=0; i < PHtypes.size(); ++i)
        if (PHtypes[i] == ngon_type::CONCAVITY_TO_SPLIT) PHlist.push_back(i);
      
      return split_phs<TriangulatorType>(crd, ngi, orient, concave_threshold, convex_threshold, RTOL, GRmin, Fluxmax, ngo, PHlist);
    }
    
    //
    template<typename TriangulatorType>
    inline static E_Int split_phs
    (K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orient, E_Float concave_threshold, E_Float convex_threshold, E_Float RTOL, E_Float GRmin, E_Float Fluxmax, ngon_type& ngo, const Vector_t<E_Int>& PHlist)
    {
      ngo.clear();//important
      
      if (PHlist.empty())
      {
        std::cout << "Nothing to split" << std::endl;
        ngo=ngi;
        return 0;
      }
      
      ngo.PGs = ngi.PGs;
      //E_Int nb_pgs = ngo.PGs.size();

#ifdef DEBUG_SPLITTER 
      E_Int faultyPH = -1;
#endif
      
      E_Int nb_err(0),nb_fixed(0), nb_false(0), nb_ok(0);
      E_Int nb_phs = ngi.PHs.size();
      
      std::vector<bool> split_it(nb_phs, false);
      for (size_t i=0; i < PHlist.size() ; ++i) split_it[PHlist[i]]=true;
      
      for (E_Int i=0; i < nb_phs ; ++i)
      {
        E_Int PHi = i;
        
        if (!split_it[i])
        {
          ++nb_ok;
          ngo.addPH(ngi, PHi);
          continue;
        }

#ifdef DEBUG_SPLITTER
        std::ostringstream o;
        o << "PH_" << PHi << ".plt";
        NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH(o.str().c_str(), crd, ngi, PHi);
#endif

        ngon_type ngsplit;

        E_Int err = NUGA::Splitter::single_concave_split<DELAUNAY::Triangulator>
                        (crd, 1/*index_start*/, ngi.PGs, ngi.PHs.get_facets_ptr(PHi), ngi.PHs.stride(PHi), orient.get_facets_ptr(PHi), 
                         concave_threshold, convex_threshold, RTOL, GRmin, Fluxmax, ngsplit);

        if (err == -1)
        {
          ngo.clear();
          std::cout << "Polyhedron " << PHi << " is open ! input mesh is corrupted " << std::endl;
          return err;
        }
        else if (err == 1)
        {
          ++nb_err;
          //std::cout << "cannot split " << PHi << std::endl;
          ngo.addPH(ngi, PHi);

#ifdef DEBUG_SPLITTER
          std::ostringstream o;
#ifndef WIN32
          o << "undoable_PH_" << PHi << ".plt";
#else
          o << "undoable_PH_" << PHi << ".tp";
#endif
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH(o.str().c_str(), crd, ngi, PHi);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PHT3(crd, ngi, PHi);
#endif
        }
        else if (err == 0 && ngsplit.PHs.size() == 2)
        {
          ++nb_fixed;
          ngo.append(ngsplit);

#ifdef DEBUG_SPLITTER 

          if (PHi == faultyPH)
          {
            for (size_t b = 0; b< ngsplit.PHs.size(); ++b)
            {
              /*std::set<K_MESH::NO_Edge> free_edges;
              bool ok = K_MESH::Polyhedron<0>::is_closed(ngsplit.PGs, ngsplit.PHs.get_facets_ptr(b), ngsplit.PHs.stride(b), free_edges);
              assert(ok);*/
              std::ostringstream o;
              o << "bit_" << PHi << "_" << b << ".tp";
              NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH(o.str().c_str(), crd, ngsplit, b);
              NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PHT3(crd, ngsplit, b);
            }
          }
#endif
        }
        else if (err == 0)
        {
          ++nb_false;
          ngo.addPH(ngi, PHi);
        }
      }  

      ngo.clean_connectivity(ngo, crd, -1/*ngon_dim*/, 0./*tolerance*/, false/*remove_dup_phs*/, false/*do_omp*/);      
      
      E_Int nb_to_split = nb_phs - nb_ok;
      std::cout << "nb of phs requiring a split : " << nb_to_split << " over " << nb_phs << std::endl;
      
      std::cout << "nb of successfull splits : " << nb_fixed << " over " << nb_to_split << " elements." << std::endl;
      std::cout << "nb of non-necessary splits : " << nb_false << " over " << nb_phs << " elements." << std::endl;
      std::cout << "nb of non-handled (unsplittable) elements : " << nb_err << " over " << nb_to_split << " elements." << std::endl;
      
      return 0;
    }

    //
    template<typename TriangulatorType>
    inline static E_Int single_concave_split(const K_FLD::FloatArray& crd, E_Int index_start,
                                             const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient,
                                             E_Float alpha_conc, E_Float alpha_cvx, E_Float rtol, E_Float GRmin, E_Float Fluxmax, ngon_type& twoPHs,
                                             const E_Float** normals = 0);

  
    ///
    inline static E_Int __split_pgs
      (const K_FLD::FloatArray& crd, const ngon_type& ngi, ngon_type& ngo, ngon_unit& split_graph, transfo_func transFunc, const transfo_t& params, Vector_t<bool>* to_process = 0);
    
    ///
    inline static E_Int __split_pgs
      (const K_FLD::FloatArray& crd, ngon_type& ngio, const ngon_unit& splitpgs, const Vector_t<E_Int>& oids);

   private:
    ///
    inline static void __apply_pgs_splits(ngon_unit& PHs, const ngon_unit& split_graph);
    ///
    template <typename TriangulatorType>
    inline static E_Int __get_chain_nodes
      (const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orient, transfo_t& params, const Vector_t<E_Int>* PHlist = 0);

    inline static E_Int __append_chain_nodes(const K_FLD::FloatArray& crd, E_Int index_start,
                                             const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient,
                                             E_Float alpha_conc, E_Float alpha_cvx, iset_t& chain_nodes, const E_Float** normals = 0);

    // valid = with no node of arity more than 2 (with current strategy)
    inline static void __get_valid_concave_chains(const eset_t& reflex_edges, idlists_t& chains, ivec_t& chains_type);
    //
    inline static void __append_all_topo_paths(const K_FLD::FloatArray& crd, const std::deque<E_Int>& chain, E_Int index_start, const id_to_ids_t& nodes_graph, const iset_t& bad_nodes, idlists_t& paths);
    //
    inline static void __append_all_geom_paths(const K_FLD::FloatArray& crd, const std::deque<E_Int>& chain, E_Int index_start, const id_to_ids_t& nodes_graph, const iset_t& bad_nodes, 
                                               const eset_t all_edges, E_Float abstol, idlists_t& paths);
    //
    inline static void __get_nodes_graph(const eset_t& edges, id_to_ids_t& nodes_graph);
    //
    inline static void __get_bad_nodes(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const std::deque<E_Int>& chain, const ivec_t& simple_colors, iset_t& bad_nodes);
    //
    template <typename TriangulatorType>
    inline static bool __valid_split(const K_FLD::FloatArray& crd, const std::vector<E_Int>& chain, E_Int index_start,
                              const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, E_Float v, E_Float GRmin, E_Float Fluxmax,
                              const ivec_t& simple_colors, const eset_t& reflex_edges, E_Float angle_max, E_Float chain_angle, E_Float & angle, E_Float minA, E_Float maxA, ngon_type& twoPHs);

    inline static E_Float __get_abs_tol(const K_FLD::FloatArray& crd, E_Float rtol, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs);
    
    Splitter(){}
    ~Splitter(){}
  };


  ///////////////////////////////////////////////////////////////////////////////////

  /// Triangulates specified PGs
  template<typename TriangulatorType>
  E_Int NUGA::Splitter::triangulate_pgs
    (const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& tri_pgs, Vector_t<E_Int>& oids, const transfo_t& qual_param, const Vector_t<bool>* to_process)
  {
    tri_pgs.clear();
    oids.clear();

    // triangulate specified PGs
    TriangulatorType dt;
    K_FLD::IntArray connectT3;

    bool shuffle = true;//qual_param.improve_qual ? true : false;
    bool imp_qual = qual_param.improve_qual ? true : false;

    E_Int err = ngon_type::triangulate_pgs<TriangulatorType>(PGs, crd, connectT3, oids, shuffle, imp_qual, to_process);
    // convert triangulation to ngon_unit
    ngon_unit::convert_fixed_stride_to_ngon_unit(connectT3, 1/*shift : convert to 1-based*/, tri_pgs);
    //update history
    if (PGs._ancEs.cols())
    {
      for (size_t i = 0; i < oids.size(); ++i)
        tri_pgs._ancEs.pushBack(PGs._ancEs.col(oids[i]), PGs._ancEs.col(oids[i]) + 2);
    }
    if (!PGs._type.empty())
    {
      for (size_t i = 0; i < oids.size(); ++i)
        tri_pgs._type.push_back(PGs._type[oids[i]]);
    }
    return err;
  }
  
  /// Convexify specified PGs
  template <typename TriangulatorType>
  E_Int NUGA::Splitter::convexify_pgs
    (const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& convex_pgs, Vector_t<E_Int>& oids, const transfo_t& params, const Vector_t<bool>* process)
  {
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);

    K_FLD::IntArray connectT3;
    ngon_unit cvx_pgs;
    E_Float normal[3];
    TriangulatorType dt;

    oids.clear();

    PGs.updateFacets();

    E_Int sz = PGs.size(), err(0);
    for (E_Int i = 0; (i<sz); ++i)
    {
      if (process && (*process)[i] == false)
        continue;

      const E_Int* nodes = PGs.get_facets_ptr(i);
      E_Int nb_nodes = PGs.stride(i);
      
      if (nb_nodes == 3) continue;

      connectT3.clear();
      err = K_MESH::Polygon::triangulate(dt, crd, nodes, nb_nodes, 1/*index start*/, connectT3, false, true);

      if (err) //pass this PG as it is
      {
        convex_pgs.add(nb_nodes, nodes);
        oids.push_back(i);
        continue;
      }

      K_MESH::Polygon::normal<acrd_t, 3>(acrd, nodes, nb_nodes, 1, normal);

      cvx_pgs.clear();
      err = NUGA::MeshTool::aggregate_convex(crd, connectT3, normal, cvx_pgs, params.convexity_tol);

      if (!err)
      {
        convex_pgs.append(cvx_pgs);
        oids.resize(convex_pgs.size(), i);
      }
      else
      {
        //std::cout << "WARNING : failed on aggregate_convex for polygon : " << i << std::endl;
        ngon_unit& ngu = cvx_pgs;
        ngon_unit::convert_fixed_stride_to_ngon_unit(connectT3, 1, ngu);
        convex_pgs.append(ngu);
        oids.resize(convex_pgs.size(), i);
      }
    }

    convex_pgs.spread_attributes(PGs, oids);

    return err;
  }

  /// Starify specified PGs at their worst reflex node
  template <typename TriangulatorType>
  E_Int NUGA::Splitter::starify_pgs_on_reflex
    (const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& pgso, Vector_t<E_Int>& oids,
    const transfo_t& params, const Vector_t<bool>* process)
  {
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);

    K_FLD::IntArray connectT3, neighbors;
    ngon_unit ngstar;
    E_Float normal[3];
    TriangulatorType dt;

    oids.clear();

    PGs.updateFacets();

    E_Int sz = PGs.size(), err(0);

    for (E_Int i = 0; (i<sz) && !err; ++i)
    {
      if (process && (*process)[i] == false)
        continue;

      const E_Int* nodes = PGs.get_facets_ptr(i);
      E_Int nb_nodes = PGs.stride(i);

      if (nb_nodes == 3)
        continue;

      err = K_MESH::Polygon::triangulate(dt, crd, nodes, nb_nodes, 1/*index start*/, connectT3, neighbors); //connectT3, neighbors are overwritten
      if (err)
        continue;

#ifdef DEBUG_NGON_T
      medith::write("before_star.mesh", crd, connectT3, "TRI");
#endif

      K_MESH::Polygon::normal<acrd_t, 3>(acrd, nodes, nb_nodes, 1, normal);

      E_Int iworst, ibest;
      bool convex = K_MESH::Polygon::is_convex(crd, nodes, nb_nodes, 1/*index_start*/, normal, params.convexity_tol, iworst, ibest);

      if (convex) //nothing to do
        continue;

      NUGA::MeshTool::starify_from_node(crd, nodes, nb_nodes, 1/*index start*/, iworst, connectT3, neighbors, normal, params.convexity_tol);

#ifdef DEBUG_NGON_T
      medith::write("after_star.mesh", crd, connectT3, "TRI");
#endif
      //in case of error, connecT3 is not changed
      ngstar.clear();
      ngon_unit::convert_fixed_stride_to_ngon_unit(connectT3, 1, ngstar);
      pgso.append(ngstar);
      oids.resize(pgso.size(), i);
    }

    pgso.spread_attributes(PGs, oids);

    return err;
  }

  /// Starify specified PGs at specified nodes in params
  template <typename TriangulatorType>
  E_Int NUGA::Splitter::starify_pgs_at_specified_nodes
    (const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& pgso, Vector_t<E_Int>& oids,
    const transfo_t& params, const Vector_t<bool>* process)
  {
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);

    TriangulatorType dt;

    ngon_unit ngstar;
    E_Float normal[3];

    K_FLD::IntArray connectT3, neighbors;

    E_Int sz(PGs.size()), err(0);

    for (E_Int i = 0; (i < sz); ++i)
    {
      if (process && (*process)[i] == false)
        continue;

      const E_Int* nodes = PGs.get_facets_ptr(i);
      E_Int nb_nodes = PGs.stride(i);

      if (nb_nodes == 3)
        continue;

      bool starify = false;
      for (E_Int n = 0; (n < nb_nodes) && !starify; ++n)
      {
        E_Int Ni = *(nodes + n);
        starify = params.nodes.find(Ni) != params.nodes.end();
      }
      if (!starify)
        continue;

      err = K_MESH::Polygon::triangulate(dt, crd, nodes, nb_nodes, 1/*index start*/, connectT3, neighbors); //connectT3, neighbors are overwritten
      if (err)
        continue;

      K_MESH::Polygon::normal<acrd_t, 3>(acrd, nodes, nb_nodes, 1, normal);

      for (E_Int n = 0; (n < nb_nodes); ++n)
      {
        E_Int Ni = *(nodes + n);
        if (params.nodes.find(Ni) == params.nodes.end())
          continue;

        NUGA::MeshTool::starify_from_node(crd, nodes, nb_nodes, 1/*index start*/, n, connectT3, neighbors, normal, params.convexity_tol);
      }

      //in case of error, connecT3 is not changed
      ngstar.clear();
      ngon_unit::convert_fixed_stride_to_ngon_unit(connectT3, 1, ngstar);
      pgso.append(ngstar);
      oids.resize(pgso.size(), i);
    }

    pgso.spread_attributes(PGs, oids);

#ifdef DEBUG_NGON_T
    assert(pgso.attributes_are_consistent());
#endif

    return 0;
  }
}

///
template <typename TriangulatorType>
E_Int NUGA::Splitter::triangulate_specified_pgs(const K_FLD::FloatArray& crd, ngon_type& ngi, const E_Int* PGlist, E_Int sz, const transfo_t& qual_param, ngon_type& ngo)
{
  E_Int nb_pgs = ngi.PGs.size();
  Vector_t<bool> toprocess(nb_pgs, (sz > 0) ? false : true);
  for (E_Int i = 0; i < sz; ++i)
    toprocess[PGlist[i]] = true;
  
  // triangulate any specified PG
  ngon_unit split_graph;
  return __split_pgs(crd, ngi, ngo, split_graph, triangulate_pgs<TriangulatorType>, qual_param, &toprocess);

}

template <typename TriangulatorType>
E_Int NUGA::Splitter::triangulate_external_pgs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Int in_or_out, const transfo_t& qual_param, ngon_type& ngo)
{
  //detect exterior faces
  ngi.flag_external_pgs(INITIAL_SKIN);
  if (in_or_out == 0) ngon_type::discard_by_box(crd, ngi, false/*i.e we want to discard external contour*/); //INTERNAL ONLY
  else if (in_or_out == 1) ngon_type::discard_by_box(crd, ngi, true/*i.e we want to discard holes*/); //EXTERNAL ONLY

  E_Int nb_pgs = ngi.PGs.size(), count(0);
  Vector_t<bool> toprocess(nb_pgs, false);
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    toprocess[i] = (ngi.PGs._type[i] == INITIAL_SKIN);
    if (toprocess[i])++count;
  }
  
  //std::cout << "nb of remaining PG tp process : " << count << std::endl;

  // triangulate any specified PG
  ngon_unit split_graph;
  __split_pgs(crd, ngi, ngo, split_graph, triangulate_pgs<TriangulatorType>, qual_param, &toprocess);

  return 0;
}

///
E_Int NUGA::Splitter::split_pgs
(transfo_func func, const K_FLD::FloatArray& crd, const ngon_type& ngi, E_Float convexity_tol, ngon_type& ngo, const Vector_t<E_Int>* PHlist)
{
  //
  transfo_t p;
  p.convexity_tol = convexity_tol;
  p.improve_qual = true;

  Vector_t<bool> *PGs_to_process(0), tmp;
  if (PHlist)
  {
    ngi.flag_facets_of_elts(*PHlist, 1, tmp);
    PGs_to_process = &tmp;
  }

  ngon_unit split_graph;
  E_Int err = __split_pgs(crd, ngi, ngo, split_graph, func, p, PGs_to_process);

  return err;
}

///
template <typename TriangulatorType>
E_Int NUGA::Splitter::starify_pgs_at_chain_nodes
(const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orient, E_Float concave_threshold, E_Float convex_threshold, E_Float convexity_tol, ngon_type& ngo, const Vector_t<E_Int>* PHlist)
{
  // 
  transfo_t p;
  p.concavity_threshold = concave_threshold;
  p.convexity_threshold = convex_threshold;
  p.convexity_tol = convexity_tol;

  E_Int err = NUGA::Splitter::__get_chain_nodes<TriangulatorType>(crd, ngi, orient, p, PHlist);
  if (err)
    return err;

  Vector_t<bool> *PGs_to_process(0), tmp;
  if (PHlist)
  {
    ngi.flag_facets_of_elts(*PHlist, 1, tmp);
    PGs_to_process = &tmp;
  }

  ngon_unit split_graph;
  err = __split_pgs(crd, ngi, ngo, split_graph, starify_pgs_at_specified_nodes<TriangulatorType>, p, PGs_to_process);

  return err;
}

template <typename TriangulatorType>
E_Int NUGA::Splitter::prepareCellsSplit(const K_FLD::FloatArray& crd, ngon_type& ngi,
  transfo_t::ePHset phset, transfo_t::eSplitPolicy policy, E_Float PH_conc_threshold, E_Float PH_cvx_threshold,
  E_Float PG_cvx_threshold, ngon_type& ngo)
{
  // Create orientation info
  ngon_unit orient;
  E_Int err = ngon_type::build_orientation_ngu<TriangulatorType>(crd, ngi, orient);
  if (err) return 1;

  std::vector<E_Int> PHlist;
  if (phset == transfo_t::CONCAVES) // all concave cells
    err = ngon_type::detect_concave_PHs(crd, ngi, orient, PH_conc_threshold, PH_cvx_threshold, PHlist);
  else if (phset == transfo_t::NONSTARS) // non-centroid-star-shaped only
  {
    Vector_t<ngon_type::ePathoPH> PHtypes;
    ngon_type::detect_pathological_PHs<TriangulatorType>(crd, ngi, orient, PHtypes);
    for (size_t i=0; i < PHtypes.size(); ++i) if (PHtypes[i] == ngon_type::CONCAVITY_TO_SPLIT) PHlist.push_back(i);
  }

  if (err) return 2;

  if (policy == transfo_t::CONVEXIFY)      // convexify any concave polygon over PH set
    err = split_pgs(convexify_pgs<TriangulatorType>, crd, ngi, PG_cvx_threshold, ngo, &PHlist);
  else if (policy == transfo_t::STARIFY_CONCAVES) //at any reflex node on PH set's concave polygons: the worst if more than 1
    err = split_pgs(starify_pgs_on_reflex<TriangulatorType>, crd, ngi, PG_cvx_threshold, ngo, &PHlist);
  else if (policy == transfo_t::STARIFY_FROM_CHAIN_ENDS) // at concave-chain ending nodes over the PH set
    err = starify_pgs_at_chain_nodes<TriangulatorType>(crd, ngi, orient, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold, ngo, &PHlist);
  
  if (err) return 3;

  return 0;
}

///
template<typename TriangulatorType>
E_Int NUGA::Splitter::single_concave_split
(const K_FLD::FloatArray& crd, E_Int index_start, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient,
E_Float concave_threshold, E_Float convex_threshold, E_Float rtol, E_Float GRmin, E_Float Fluxmax, ngon_type& twoPHs, const E_Float** normals)
{

  // 1. DETECT CONCAVITIES /////////////////////////////////////////////////////////////////
  bool concave = false;
  edge_angle_t reflex_edges;
  eset_t convex_edges;
  E_Int err = K_MESH::Polyhedron<UNKNOWN>::is_concave(crd, PGS, first_pg, nb_pgs, false/*i.e. closed PH*/, orient, concave, reflex_edges, convex_edges, concave_threshold, convex_threshold, normals); // reflex edges are 0-based indices

  if (err /*|| !concave*/)
    return err;
  
  E_Float minA, maxA;
  E_Int maxAPG1, maxAPG2;
  err = K_MESH::Polyhedron<UNKNOWN>::min_max_angles(crd, PGS, first_pg, nb_pgs, false, orient, minA, maxA, maxAPG1, maxAPG2);
  

  // 2.  SIMPLIFIED VIEW : COLORS //////////////////////////////////////////////////
  std::vector<E_Int> simple_colors;
  ngon_unit lneighbors;
  {
    eset_t wire = convex_edges;
    for (edge_angle_t::const_iterator ii = reflex_edges.begin(); ii != reflex_edges.end(); ++ii)
      wire.insert(ii->first);
    // neighbor and color
    K_MESH::Polygon::build_pg_neighborhood(PGS, lneighbors, first_pg, nb_pgs, &wire);
    NUGA::EltAlgo<K_MESH::Polygon>::coloring(lneighbors, simple_colors);
  }

  //3. GET VALID CONCAVE CHAINS
  idlists_t chains;
  ivec_t chains_type;
  eset_t rflx_edges;
  for (edge_angle_t::const_iterator ii = reflex_edges.begin(); ii != reflex_edges.end(); ++ii)
    rflx_edges.insert(ii->first);

#ifdef DEBUG_PHSPLITTER
  //if (PHi == 160)
  {
    K_FLD::IntArray cT3;
    std::vector<E_Int> nT3_to_PG, scolorT3;
    TriangulatorType dt;
    E_Int err = K_MESH::Polyhedron<UNKNOWN>::triangulate(dt, PGS, first_pg, nb_pgs, crd, cT3, nT3_to_PG, true, false);
    scolorT3.resize(cT3.cols());
    for (E_Int i=0; i < cT3.cols(); ++i) scolorT3[i] = simple_colors[nT3_to_PG[i]];
    medith::write("simplified.mesh", crd, cT3, "TRI", 0, &scolorT3);

    K_FLD::IntArray cntE;
    for (std::set<K_MESH::NO_Edge>::iterator it = rflx_edges.begin(); it != rflx_edges.end(); ++it)
      cntE.pushBack(it->begin(), it->end());
    cntE.shift(-1);
    medith::write("reflex_edges.mesh", crd, cntE, "BAR");
    
    K_FLD::IntArray cntSk;
    for (std::set<K_MESH::NO_Edge>::iterator it = convex_edges.begin(); it != convex_edges.end(); ++it)
      cntSk.pushBack(it->begin(), it->end());
    cntSk.shift(-1);
    medith::write("skeleton.mesh", crd, cntSk, "BAR");
  }
#endif

  __get_valid_concave_chains(rflx_edges, chains, chains_type);

  if (chains.empty()) // not handled with the current strategy
    return 1;

  E_Float abstol = __get_abs_tol(crd, rtol, PGS, first_pg, nb_pgs);
  E_Float angle_max = std::max/*min*/(concave_threshold, convex_threshold) * NUGA::PI;

  size_t nb_chains = chains.size();

  E_Float v, G[3];
  TriangulatorType dt;
  K_MESH::Polyhedron<UNKNOWN>::metrics<TriangulatorType>(dt, crd, PGS, first_pg, nb_pgs, v, G);

  // get worst chain's angle : should create a cut better than that
  std::vector<E_Float> chain_angle(nb_chains, 0.);
  for (size_t c = 0; c < nb_chains; ++c)
  {
    size_t sz = chains[c].size();
    for (size_t n = 0; n < sz-1; ++n)
    {
      K_MESH::NO_Edge E(chains[c][n], chains[c][n + 1]);
      edge_angle_t::const_iterator it = reflex_edges.find(E);
      assert(it != reflex_edges.end());
      E_Float a = it->second;
      chain_angle[c] = std::max(chain_angle[c], ::fabs(NUGA::PI - a));//max deviation in the chain
    }
  }

  //4. IF VALID CLOSED PATH EXIST, DONE
  typedef std::vector< std::pair<E_Float, ngon_type> > palmares_t;
  palmares_t palma_split;
  ivec_t chain;
  for (size_t c = 0; c < nb_chains; ++c)
  {
    if (chains_type[c] == 1)
      continue;

    chain.clear();
    chain.insert(chain.end(), chains[c].begin(), chains[c].end()); //using &a[0] whan a is a deque is coorupted
    
    E_Float angle;
    if (!__valid_split<TriangulatorType>(crd, chain, index_start,
                       PGS, first_pg, nb_pgs, v, GRmin, Fluxmax, simple_colors, rflx_edges, angle_max, chain_angle[c], angle, minA, maxA, twoPHs))
      continue;

    if (angle < EPSILON) // perfectly planar cut
      return 0;

    palma_split.push_back(std::make_pair(angle, twoPHs));
  }

  //5. OVERALL NODES GRAPH
  eset_t all_edges;
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    E_Int PGi = *(first_pg + i) - 1;
    const E_Int* nodes = PGS.get_facets_ptr(PGi);
    E_Int nb_nodes = PGS.stride(PGi);
    for (E_Int n = 0; n < nb_nodes; ++n)
      all_edges.insert(K_MESH::NO_Edge(*(nodes + n), *(nodes + (n + 1) % nb_nodes)));
  }

  id_to_ids_t nodes_graph;
  __get_nodes_graph(all_edges, nodes_graph);

  //6. ENCLOSE OPEN CHAIN TO FORM PATH
  idlists_t paths;
  for (size_t c = 0; c < nb_chains; ++c)
  {
    if (chains_type[c] == 0)
      continue;

    E_Int sz = chains[c].size();

    iset_t bad_nodes;
    __get_bad_nodes(PGS, first_pg, nb_pgs, chains[c], simple_colors, bad_nodes);
    //except begin and end chain nodes
    bad_nodes.erase(chains[c][0]); //
    bad_nodes.erase(chains[c][sz - 1]);

    paths.clear();
    __append_all_topo_paths(crd, chains[c], index_start, nodes_graph, bad_nodes, paths);
    __append_all_geom_paths(crd, chains[c], index_start, nodes_graph, bad_nodes, all_edges, abstol, paths);

    size_t nb_paths = paths.size();
    if (nb_paths == 0)
      continue;

    //check validty
    for (size_t p = 0; p < nb_paths; ++p)
    {
      //for (size_t k = 0; k < paths[p].size(); ++k)
        //std::cout << paths[p][k] << std::endl;
      E_Float angle;

      chain.clear();
      chain.insert(chain.end(), paths[p].begin(), paths[p].end());//using &a[0] whan a is a deque is coorupted
      
      if (!__valid_split<TriangulatorType>(crd, chain, index_start,
                         PGS, first_pg, nb_pgs, v, GRmin, Fluxmax, simple_colors, rflx_edges, angle_max, chain_angle[c], angle, minA, maxA, twoPHs))
        continue;

      if (angle < EPSILON) // perfectly planar cut
        return 0;

      palma_split.push_back(std::make_pair(angle, twoPHs));
    }
  }

  if (palma_split.empty())
    return 1;

  // 7. CHOOSE THE BEST (MOST PLANAR)
  //std::sort(palma_split.begin(), palma_split.end());
  E_Int bestc = IDX_NONE;
  E_Float besta = NUGA::FLOAT_MAX;
  for (size_t i = 0; i < palma_split.size(); ++i)
  {
    if (palma_split[i].first < besta)
    {
      besta = palma_split[i].first;
      bestc = i;
    }
  }

  twoPHs = palma_split[bestc].second;

  return 0;
}

// PRIVATE ////////////////////////////////////////////////////////////////////

///
E_Int NUGA::Splitter::__split_pgs
(const K_FLD::FloatArray& crd, const ngon_type& ngi, ngon_type& ngo, ngon_unit& split_graph, NUGA::transfo_func transFunc, const NUGA::transfo_t& params, Vector_t<bool>* to_process)
{
  ngon_unit splitpgs/*list of split pgs*/;

  // do the transform and produce the split pgs (and their history).
  Vector_t<E_Int> oids;
  transFunc(ngi.PGs, crd, splitpgs, oids, params, to_process);

  E_Int nb_pgs = ngi.PGs.size();
  ngo = ngi; //init
  // append the split pgs and synchronize the graph.
  ngo.PGs.append(splitpgs);

  // now ngo contain both old and new pgs.

  // update oids to have it for all the polygons
  {
    std::vector<E_Int> new_oids;
    K_CONNECT::IdTool::init_inc(new_oids, ngo.PGs.size());
    for (size_t k = 0; k < oids.size(); ++k) {
      new_oids[k + nb_pgs] = oids[k];
      //force oids to be NONE when an entity has been split
      // so an entity cannot be self-referring anf having children referring to it
      if (oids[k] < nb_pgs)new_oids[oids[k]] = IDX_NONE;
    }

    oids = new_oids;
  }

  // create the split graph using the history.
  split_graph.clear();
  K_CONNECT::IdTool::reverse_indirection(nb_pgs, &oids[0], oids.size(), split_graph);

  // apply the modifications at the PHs level.
  __apply_pgs_splits(ngo.PHs, split_graph);

  // remove unreferenced (among olds) pgs.
  Vector_t<E_Int> pgnids, phnids;
  ngo.remove_unreferenced_pgs(pgnids, phnids);

  return 0;
}

/// same as above with specified pgs to replace (splitpgs + oids)
E_Int NUGA::Splitter::__split_pgs
(const K_FLD::FloatArray& crd, ngon_type& ngio, const ngon_unit& splitpgs, const Vector_t<E_Int>& oids)
{
  E_Int nb_pgs = ngio.PGs.size();
  
  ngio.PGs.append(splitpgs);

  // now ngo contain both old and new pgs.

  // update oids to have it for all the polygons
  std::vector<E_Int> new_oids;
  {
    K_CONNECT::IdTool::init_inc(new_oids, ngio.PGs.size());
    for (size_t k = 0; k < oids.size(); ++k) {
      new_oids[k + nb_pgs] = oids[k];
      //force oids to be NONE when an entity has been split
      // so an entity cannot be self-referring anf having children referring to it
      if (oids[k] < nb_pgs)new_oids[oids[k]] = IDX_NONE;
    }
  }

  // create the split graph using the history.
  ngon_unit split_graph;
  K_CONNECT::IdTool::reverse_indirection(nb_pgs, &new_oids[0], new_oids.size(), split_graph);

  // apply the modifications at the PHs level.
  __apply_pgs_splits(ngio.PHs, split_graph);

  // remove unreferenced (among olds) pgs.
  Vector_t<E_Int> pgnids, phnids;
  ngio.remove_unreferenced_pgs(pgnids, phnids);

  return 0;
}

///
void NUGA::Splitter::__apply_pgs_splits(ngon_unit& PHs, const ngon_unit& split_graph)
{
  // split_graph is sized as PGs BEFORE APPENDING, iE. what is referenced in PHs, but not what is stored in PGs
  // if a given PG is gone, molec is (1,IDX_NONE).
#ifdef DEBUG_NGON_T
  assert(PHs.get_facets_max_id() == split_graph.size());
#endif

  split_graph.updateFacets();
  PHs.updateFacets();

  E_Int nb_phs(PHs.size()), nnbe(0);
  ngon_unit new_phs;

  Vector_t<E_Int> molecPH;
  Vector_t<E_Int> nids(nb_phs, IDX_NONE);

  for (E_Int i = 0; i < nb_phs; ++i)
  {
    const E_Int* pPGi = PHs.get_facets_ptr(i);
    E_Int nb_pgs = PHs.stride(i);

    molecPH.clear();

    for (E_Int p = 0; p < nb_pgs; ++p)
    {
      E_Int PGi = *(pPGi + p) - 1;
      E_Int nb_split = split_graph.stride(PGi);
      const E_Int* new_facets = split_graph.get_facets_ptr(PGi);

      if (nb_split == 1 && *new_facets == IDX_NONE) // disappeared entity
        continue;
      else
      {
        for (E_Int k = 0; k < nb_split; ++k)
          molecPH.push_back(new_facets[k] + 1);
      }
    }

    E_Int sz = molecPH.size();
    if (sz == 0) // PH is gone
      continue;

    new_phs.add(molecPH.size(), &molecPH[0]);
    nids[i] = nnbe++;
  }

  //transfer attributes
  new_phs.compact_attributes(PHs, nids);

  PHs = new_phs;
  PHs.updateFacets();
}

///
void NUGA::Splitter::__get_valid_concave_chains
(const eset_t& reflex_edges, idlists_t& chains, ivec_t& chains_type)
{
  if (reflex_edges.empty()) return;

  //GRAPH "node to nodes" for REFLEX ONLY
  id_to_ids_t nodes_graph;
  __get_nodes_graph(reflex_edges, nodes_graph);
  
  //type : 0 for closed, 1 for open
  std::vector<E_Int> edgecol;
  BARSplitter::split_eset_into_manifold_chains<eset_t>(reflex_edges, chains, edgecol);

  size_t nb_chains = chains.size();
  //2.d remove any non-valid chain (having non-manifold node)
  {
    idlists_t tmp;
    id_to_ids_t::const_iterator it;
    for (size_t c = 0; c < nb_chains; ++c)
    {
      bool valid = true;
      E_Int sz = chains[c].size();
      for (E_Int n = 0; (n < sz) && valid; ++n)
      {
        const E_Int& N = chains[c][n];
        it = nodes_graph.find(N);
        assert(it != nodes_graph.end());
        valid &= (it->second.size() < 3); //arity less than 3
      }
       
      if (valid)
        tmp.push_back(chains[c]);
    }
    chains = tmp;
    nb_chains = chains.size();
  }
  
  //2.e. assign a type to each remaining valid chain
  chains_type.resize(nb_chains, 1);
  //bool has_closed = false;
  for (size_t c = 0; c < nb_chains; ++c)
  {
    E_Int sz = chains[c].size();
    if (chains[c][0] == chains[c][sz - 1])
    {
      chains_type[c] = 0;
      chains[c].pop_back(); // redundant node : would break the logic in some following algos
    }
  }
}

///
template <typename TriangulatorType>
E_Int NUGA::Splitter::__get_chain_nodes
(const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orient, transfo_t& params, const Vector_t<E_Int>* PHlist)
{

  params.nodes.clear();

  if (PHlist) //only over selected PHs
  {
    E_Int nb_phs = PHlist->size();
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      E_Int PHi = (*PHlist)[i];
      const E_Int* first_pg = ngi.PHs.get_facets_ptr(PHi);
      E_Int nb_pgs = ngi.PHs.stride(PHi);
      NUGA::Splitter::__append_chain_nodes(crd, 1, ngi.PGs, first_pg, nb_pgs, orient.get_facets_ptr(PHi), params.concavity_threshold, params.convexity_threshold, params.nodes);
    }
  }
  else // over all
  {
    E_Int nb_phs = ngi.PHs.size();
    for (E_Int PHi = 0; PHi < nb_phs; ++PHi)
    {
      const E_Int* first_pg = ngi.PHs.get_facets_ptr(PHi);
      E_Int nb_pgs = ngi.PHs.stride(PHi);
      NUGA::Splitter::__append_chain_nodes(crd, 1, ngi.PGs, first_pg, nb_pgs, orient.get_facets_ptr(PHi), params.concavity_threshold, params.convexity_threshold, params.nodes);
    }
  }

  return 0;
}

///
void NUGA::Splitter::__append_all_topo_paths
(const K_FLD::FloatArray& crd, const std::deque<E_Int>& chain, E_Int index_start, 
 const id_to_ids_t& nodes_graph, const iset_t& bad_nodes, idlists_t& paths)
{
  E_Int sz = chain.size();
  if (sz < 2)
    return;

  E_Int nbegin = chain[sz - 1];
  E_Int nend = chain[0];
  //E_Int nprev = chain[sz - 2];

  id_to_ids_t::const_iterator itb = nodes_graph.find(nbegin);
  assert(itb != nodes_graph.end());

  // FALSE LINK MISSING
  if (sz > 2)
  {
    for (iset_t::const_iterator itN = itb->second.begin(); (itN != itb->second.end()); ++itN)
    {
      if (*itN == nend) // in fact a type 0 : the missing edge exists but is convex
      {
        paths.push_back(chain);
        return;
      }
    }
  }

  // ONE LINK MISSING
  id_to_ids_t::const_iterator ite = nodes_graph.find(nend);
  assert(ite != nodes_graph.end());
  {
    std::vector<E_Int> commonN;
    std::set_intersection(itb->second.begin(), itb->second.end(), ite->second.begin(), ite->second.end(), std::back_inserter(commonN));
    
#ifdef DEBUG_SPLITTER
    for (E_Int i=0; i < commonN.size(); ++i) std::cout << commonN[i] << " ";
     std::cout << std::endl;
#endif

    iset_t final_commonN;
    for (size_t i = 0; i < commonN.size(); ++i)
      if ((bad_nodes.find(commonN[i]) == bad_nodes.end()))
        final_commonN.insert(commonN[i]);
    
#ifdef DEBUG_SPLITTER
    for (iset_t::iterator i=final_commonN.begin(); i != final_commonN.end(); ++i) std::cout << *i << " ";
     std::cout << std::endl;
#endif
    
    final_commonN.erase(nend);
    final_commonN.erase(nbegin);

    for (iset_t::const_iterator n = final_commonN.begin(); n != final_commonN.end(); ++n)
    {
      paths.push_back(chain);
      paths.back().push_back(*n); //close the chain
    }
  }

  // TWO LINKS MISSING
  for (iset_t::const_iterator nb = itb->second.begin(); nb != itb->second.end(); ++nb)
  {
    E_Int Ni = *nb;
    if (bad_nodes.find(Ni) != bad_nodes.end())
      continue;
    if (Ni == nend)
      continue;
    
    for (std::set<E_Int>::const_iterator ne = ite->second.begin(); ne != ite->second.end(); ++ne)
    {
      E_Int Nj = *ne;
      if (bad_nodes.find(Nj) != bad_nodes.end())
        continue;
      if (Nj == nbegin)
        continue;

      id_to_ids_t::const_iterator itNj = nodes_graph.find(Nj);
      assert(itNj != nodes_graph.end());

      for (std::set<E_Int>::const_iterator ne2 = itNj->second.begin(); ne2 != itNj->second.end(); ++ne2)
      {
        E_Int Nk = *ne2;
        if (Nk != Ni)
          continue;

        paths.push_back(chain);
        paths.back().push_back(Ni);
        paths.back().push_back(Nj); //close the chain
      }
    }
  }
  
  //if (!paths.empty()) return;
  
  // 3 LINKS MISSING
  for (iset_t::const_iterator nb = itb->second.begin(); nb != itb->second.end(); ++nb)
  {
    E_Int Ni = *nb;
    if (bad_nodes.find(Ni) != bad_nodes.end())
      continue;
    if (Ni == nend)
      continue;
    
    id_to_ids_t::const_iterator itb1 = nodes_graph.find(Ni);
    
    for (std::set<E_Int>::const_iterator ne = ite->second.begin(); ne != ite->second.end(); ++ne)
    {
      E_Int Nj = *ne;
      if (bad_nodes.find(Nj) != bad_nodes.end())
        continue;
      if (Nj == nbegin)
        continue;
      
      id_to_ids_t::const_iterator ite1 = nodes_graph.find(Nj);

      std::vector<E_Int> commonN;
      std::set_intersection(itb1->second.begin(), itb1->second.end(), ite1->second.begin(), ite1->second.end(), std::back_inserter(commonN));

      iset_t final_commonN;
      for (size_t i = 0; i < commonN.size(); ++i)
        if ((bad_nodes.find(commonN[i]) == bad_nodes.end()))
          final_commonN.insert(commonN[i]);
    
      final_commonN.erase(Ni);
      final_commonN.erase(Nj);

      for (iset_t::const_iterator n = final_commonN.begin(); n != final_commonN.end(); ++n)
      {
        paths.push_back(chain);
        paths.back().push_back(Ni); //close the chain
        paths.back().push_back(*n); //close the chain
        paths.back().push_back(Nj); //close the chain
      }
    }
  }
  
}

///
void NUGA::Splitter::__append_all_geom_paths
(const K_FLD::FloatArray& crd, const std::deque<E_Int>& chain, E_Int index_start, 
 const id_to_ids_t& nodes_graph, const iset_t& bad_nodes, 
 const eset_t all_edges, E_Float abstol,
 idlists_t& paths)
{
  E_Int sz = chain.size();
  if (sz < 2)
    return;

  eset_t split_edges;
  K_FLD::IntArray cntE;
  
  E_Int nbegin = chain[sz - 1];
  //E_Int nend = chain[0];
  E_Int nprev = chain[sz - 2];

  id_to_ids_t::const_iterator itb = nodes_graph.find(nbegin);
  assert(itb != nodes_graph.end());

  for (std::set<E_Int>::const_iterator nb = itb->second.begin(); nb != itb->second.end(); ++nb)
  {
    E_Int N = *nb;
    if (bad_nodes.find(N) != bad_nodes.end())
      continue;

    //cutting direction
    E_Float cutDir[3];
    K_MESH::Triangle::normal(crd.col(nprev - 1), crd.col(nbegin - 1), crd.col(N - 1), cutDir);

    // Get lying edges
    split_edges.clear();
    NUGA::MeshTool::get_edges_lying_on_plane(crd, 1, all_edges, nbegin - 1, cutDir, abstol, split_edges);

#ifdef DEBUG_PHSPLITTER
    /*{
      K_FLD::IntArray tmp;
      for (eset_t::iterator it = split_edges.begin(); it != split_edges.end(); ++it)
        tmp.pushBack(it->begin(), it->end());
      tmp.shift(-1);
      //std::cout << tmp << std::endl;
      std::ostringstream o;
      o << "lying_edges_" << PHi << ".mesh";
      medith::write(o.str().c_str(), crd, tmp, "BAR");
    }*/
#endif

    //nodes arity
    std::map<E_Int, E_Int> node_to_count;
    NUGA::MeshTool::build_node_arity(split_edges, node_to_count);

#ifdef DEBUG3_POLYHEDRON
    {
      K_FLD::IntArray tmp;
      for (eset_t::iterator it = split_edges.begin(); it != split_edges.end(); ++it)
        tmp.pushBack(it->begin(), it->end());
      tmp.shift(-1);
      //std::cout << tmp << std::endl;
      std::ostringstream o;
      o << "withadded_edges_" << PHi << ".mesh";
      medith::write(o.str().c_str(), crd, tmp, "BAR");
    }
#endif

    // Burn free branches
    NUGA::MeshTool::burn_free_branches(split_edges, node_to_count);
    if (split_edges.size() < 3)
      continue;

    // Remove any edge having a bad node
    eset_t tmp(split_edges);
    for (eset_t::iterator it = tmp.begin(); it != tmp.end(); ++it)
    {
      if (bad_nodes.find(it->node(0)) != bad_nodes.end())
        split_edges.erase(*it);
      else if (bad_nodes.find(it->node(1)) != bad_nodes.end())
        split_edges.erase(*it);
    }

    if (split_edges.size() < 3)
      continue;

    // Check if one single closed contour
    cntE.clear();
    for (eset_t::iterator it = split_edges.begin(); it != split_edges.end(); ++it)
      cntE.pushBack(it->begin(), it->end());
    //
    K_FLD::IntArray neighE;
    std::vector<E_Int> colorsE;
    NUGA::EltAlgo<K_MESH::Edge>::getManifoldNeighbours(cntE, neighE, false);
    NUGA::EltAlgo<K_MESH::Edge>::coloring(neighE, colorsE);
    E_Int nb_loops = 1 + *std::max_element(colorsE.begin(), colorsE.end());

    if (nb_loops == 1) //bingo
    {
      std::vector<E_Int> snodes;
      BARSplitter::getSortedNodes(cntE, snodes);
      std::deque<E_Int> pth;
      pth.insert(pth.end(), snodes.begin(), snodes.end());
      paths.push_back(pth);
    }
  }
}

///
E_Float NUGA::Splitter::__get_abs_tol
(const K_FLD::FloatArray& crd, E_Float rtol, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs)
{
  E_Float LRef = 1.;
  ivec_t nodes;
  K_MESH::Polyhedron<UNKNOWN>::unique_nodes(PGS, first_pg, nb_pgs, nodes);
  K_CONNECT::IdTool::shift(nodes, -1);
  //std::cout << nodes.size() << std::endl;
  K_SEARCH::BBox3D box;
  box.compute(crd, nodes);
  LRef = std::max(box.maxB[0] - box.minB[0], box.maxB[1] - box.minB[1]);
  LRef = std::max(box.maxB[2] - box.minB[2], LRef);

  return LRef * rtol;
}

///
void NUGA::Splitter::__get_nodes_graph(const eset_t& edges, id_to_ids_t& nodes_graph)
{
  for (eset_t::const_iterator e = edges.begin(); e != edges.end(); ++e)
  {
    E_Int Ni = e->node(0);
    E_Int Nj = e->node(1);
    nodes_graph[Ni].insert(Nj);
    nodes_graph[Nj].insert(Ni);
  }
}

///
void NUGA::Splitter::__get_bad_nodes
(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const std::deque<E_Int>& chain, const ivec_t& simple_colors, iset_t& bad_nodes)
{
  bad_nodes.clear();
  eset_t chain_edges;
  size_t sz = chain.size();
  for (size_t i = 0; i < sz; ++i)
    chain_edges.insert(K_MESH::NO_Edge(chain[i], chain[(i + 1) % sz]));

  // first set bad colors for this chain
  iset_t bad_colors;
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    E_Int PGi = *(first_pg + i) - 1;
    const E_Int* nodes = PGS.get_facets_ptr(PGi);
    E_Int nb_nodes = PGS.stride(PGi);

    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      //E_Int N = *(nodes + n);
      //E_Int np1 = *(nodes + (n + 1) % nb_nodes);
      K_MESH::NO_Edge E(*(nodes + n), *(nodes + (n + 1) % nb_nodes));
      if (chain_edges.find(E) == chain_edges.end())
        continue;

      bad_colors.insert(simple_colors[i]);
    }
  }
  // bad nodes are those belonging to bad colors 
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    E_Int PGi = *(first_pg + i) - 1;

    if (bad_colors.find(simple_colors[i]) == bad_colors.end())//fixme : change i insetad PGi
      continue;

    const E_Int* nodes = PGS.get_facets_ptr(PGi);
    E_Int nb_nodes = PGS.stride(PGi);

    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      //E_Int N = *(nodes + n);
      bad_nodes.insert(*(nodes + n));
    }
  }
#ifdef DEBUG_PHSPLITTER
  for (iset_t::const_iterator i = bad_nodes.begin(); i != bad_nodes.end(); ++i) std::cout << *i << " ";
  std::cout << std::endl;
#endif
}

///
E_Int NUGA::Splitter::__append_chain_nodes
(const K_FLD::FloatArray& crd, E_Int index_start, const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, const E_Int* orient,
 E_Float concave_threshold, E_Float convex_threshold, iset_t& chain_nodes, const E_Float** normals)
{

  // 1. DETECT CONCAVITIES /////////////////////////////////////////////////////////////////
  bool concave = false;
  edge_angle_t reflex_edges;
  eset_t convex_edges;
  E_Int err = K_MESH::Polyhedron<UNKNOWN>::is_concave(crd, PGS, first_pg, nb_pgs, false/* i.e. closed PH*/, orient, concave, reflex_edges, convex_edges, concave_threshold, convex_threshold, normals); // reflex edges are 0-based indices

  if (err || !concave)
    return err;  

  //4. GET VALID CONCAVE CHAINS
  idlists_t chains;
  ivec_t chains_type;
  eset_t rflx_edges;
  for (edge_angle_t::const_iterator ii = reflex_edges.begin(); ii != reflex_edges.end(); ++ii)
    rflx_edges.insert(ii->first);
  
  __get_valid_concave_chains(rflx_edges, chains, chains_type);

  for (size_t c = 0; c < chains.size(); ++c)
  {
    chain_nodes.insert(chains[c][0]);
    chain_nodes.insert(chains[c][chains[c].size()-1]);
  }

  return 0;
}

///
template <typename TriangulatorType>
bool NUGA::Splitter::__valid_split
(const K_FLD::FloatArray& crd, const std::vector<E_Int>& chain, E_Int index_start,
const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs, E_Float V, E_Float GRmin, E_Float Fluxmax,
const ivec_t& simple_colors, const eset_t& reflex_edges, E_Float angle_max, E_Float chain_angle, E_Float & angle, E_Float minA, E_Float maxA, ngon_type& twoPH)
{
  angle = 0;
  twoPH.clear();

  E_Int nb_ch_nodes = chain.size();

#ifdef DEBUG_PHSPLITTER
if (PHi == faultyPH)
{
  K_FLD::IntArray cntE;
  for (size_t n = 0; n < nb_ch_nodes; ++n)
  {
    E_Int e[] = { chain[n] - index_start, chain[(n + 1) % nb_ch_nodes] - index_start };
    //std::cout << e[0] << "-" << e[1] << std::endl;
    cntE.pushBack(e, e + 2);
  }
  //std::cout << cntE << std::endl;
  medith::write("path.mesh", crd, cntE, "BAR");
}
#endif

  // check chain validity using normal : 
  // 3 consecutive nodes in the chain belonging to a same face should not make a normal aligned with that face
  {
    std::map<int, std::set<int>> node_to_faces;
    for (E_Int f = 0; f < nb_pgs; ++f)
    {
      E_Int PGi = first_pg[f] - 1;
      const E_Int* pn = PGS.get_facets_ptr(PGi);
      int nnodes = PGS.stride(PGi);

      for (int n = 0; n < nnodes; ++n)
      {
        E_Int Ni = pn[n];
        node_to_faces[Ni].insert(PGi);
      }
    }

    size_t nc = chain.size();
    std::vector<int> common_faces1, common_faces;
    for (size_t n = 0; n < nc; ++n)
    {
      E_Int N1 = chain[n];
      E_Int N2 = chain[(n + 1) % nc];
      E_Int N3 = chain[(n + 2) % nc];

      common_faces1.clear();
      std::set_intersection(ALL(node_to_faces[N1]), ALL(node_to_faces[N2]), std::back_inserter(common_faces1));
      if (common_faces1.empty()) continue;
      std::sort(ALL(common_faces1));
      common_faces.clear();
      std::set_intersection(ALL(common_faces1), ALL(node_to_faces[N3]), std::back_inserter(common_faces));
      if (common_faces.empty()) continue;

      const E_Float* P0 = crd.col(N1-1);
      const E_Float* P1 = crd.col(N2-1);
      const E_Float* P2 = crd.col(N3-1);

      E_Float N[3];
      K_MESH::Triangle::normal(P0, P1, P2, N);

      E_Float nPGi[3];
      for (auto PGi : common_faces)
      {
        const E_Int* pn = PGS.get_facets_ptr(PGi);
        int nnodes = PGS.stride(PGi);

        K_MESH::Polygon::normal<K_FLD::FloatArray, 3>(crd, pn, nnodes, 1, nPGi);
        E_Float ps = ::fabs(NUGA::dot<3>(N, nPGi));
        if (ps > 0.99) return false;
      }
    }
  }

  // try the cut
  eset_t split_edges;
  for (E_Int n = 0; n < nb_ch_nodes; ++n)
    split_edges.insert(K_MESH::NO_Edge(chain[n], chain[(n + 1) % nb_ch_nodes]));

  ngon_unit lneighbors;
  K_MESH::Polygon::build_pg_neighborhood(PGS, lneighbors, first_pg, nb_pgs, &split_edges);
  ivec_t bits_colors;
  NUGA::EltAlgo<K_MESH::Polygon>::coloring(lneighbors, bits_colors);
  E_Int nb_bits = 1 + *std::max_element(bits_colors.begin(), bits_colors.end());

  if (nb_bits != 2) // not a closed PG cut
    return false;

  TriangulatorType dt;
  //static int count = 0;
  K_FLD::IntArray cT3;
  bool prefer_decomp = false;
  if (nb_ch_nodes > 3)
  {
    // check if the PG cut will produce smaller concavities than the input
    K_FLD::IntArray neighT3;
    E_Int err = K_MESH::Polygon::triangulate(dt, crd, &chain[0], nb_ch_nodes, index_start, cT3, neighT3);
    if (err)
      return false;

    angle = NUGA::MeshTool::get_max_deviation(crd, cT3, neighT3);

    if (angle >= chain_angle) // the cut PG produce worst concavities
      return false;

    if (angle >= angle_max) // prefer T3 decomp : fixme : criterion
    {
      //std::cout << "decomp " << ++count << std::endl;
      prefer_decomp = true;
    }
  }
  else // the cut PG is a T3
    angle = 0.;

  // check the the 2 bits have non-null volume (relative to input)

  // create them
  twoPH.PGs.clear();
  twoPH.PGs.append(PGS, first_pg, nb_pgs);

  E_Int pgbits = 1;
  if (!prefer_decomp)
    twoPH.PGs.add(nb_ch_nodes, &chain[0]);//fixme : assume index_start is 1
  else
  {
    pgbits = cT3.cols();
    ngon_unit ngu;
    ngon_unit::convert_fixed_stride_to_ngon_unit(cT3, 1, ngu);
    twoPH.PGs.append(ngu);
  }

  if (!twoPH.PGs._type.empty())
    for (E_Int b = 0; b < pgbits; ++b) twoPH.PGs._type.push_back(INNER);
  E_Int none[] = { IDX_NONE, IDX_NONE };
  if (twoPH.PGs._ancEs.cols() != 0)
    for (E_Int b = 0; b < pgbits; ++b) twoPH.PGs._ancEs.pushBack(none, none + 2);

  ivec_t molecPH1, molecPH2;
  for (E_Int p = 0; p < nb_pgs; ++p)
  {
    if (bits_colors[p] == 0)
      molecPH1.push_back(p + 1);
    else
      molecPH2.push_back(p + 1);
  }

  for (E_Int b = 0; b < pgbits; ++b)
  {
    molecPH1.push_back(nb_pgs + b + 1); // add the cut PG bit
    molecPH2.push_back(nb_pgs + b + 1); // add the cut PG bit
  }

  if (molecPH1.size() - (pgbits - 1) < 4) // invalid PH bit
    return false;

  if (molecPH2.size() - (pgbits - 1) < 4) // invalid bit
    return false;

  twoPH.PHs.add(molecPH1.size(), &molecPH1[0]);
  twoPH.PHs.add(molecPH2.size(), &molecPH2[0]);

  twoPH.PHs._type.resize(2, IDX_NONE);// unknwon : a bit might be SKIN, the orhter one INNER
  
  twoPH.PHs._dirty=twoPH.PGs._dirty=true;
  twoPH.PGs.updateFacets();
  twoPH.PHs.updateFacets();

  /*if (ng.PHs._ancEs.cols() > PHi) FIXME
  {
    E_Int* anc = twoPH.PHs._ancEs.col(PHi);
    twoPH.PHs._ancEs.pushBack(anc, anc + 2);
    twoPH.PHs._ancEs.pushBack(anc, anc + 2);
  }*/
  
  /*for (size_t k=0; k < twoPH.PGs._NGON.size(); ++k)
    std::cout << twoPH.PGs._NGON[k] << std::endl;*/
    
#ifdef DEBUG_SPLITTER
  static int count = 0;
  ++count;
  std::ostringstream o;
  
  {
    o << "D:\\slandier\\DATA\\tmp\\split\\twoPH_" << count << "a.tp";
    ngon_type ngo;
    ngo.PHs.add(twoPH.PHs.stride(0), twoPH.PHs.get_facets_ptr(0));
    ngo.PGs = twoPH.PGs;
    ngo.PHs.updateFacets();
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    tp::write(o.str().c_str(), crd, cnto, "NGON");
  }

  {
    o.str("");
    o << "D:\\slandier\\DATA\\tmp\\split\\twoPH_" << count << "b.tp";
    ngon_type ngo;
    ngo.PHs.add(twoPH.PHs.stride(1), twoPH.PHs.get_facets_ptr(1));
    ngo.PGs = twoPH.PGs;
    ngo.PHs.updateFacets();
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    tp::write(o.str().c_str(), crd, cnto, "NGON");
  }
#endif

  // check the the 2 bits have non-null volume (relative to input)
  E_Float G[3], v1, v2;
  E_Int nb_pgs0 = twoPH.PHs.stride(0);
  E_Int nb_pgs1 = twoPH.PHs.stride(1);
  const E_Int * pgs0 = twoPH.PHs.get_facets_ptr(0);
  const E_Int * pgs1 = twoPH.PHs.get_facets_ptr(1);
  
  K_MESH::Polyhedron<UNKNOWN>::metrics2(dt, crd, twoPH.PGs, pgs0, nb_pgs0, v1, G, false/*not all cvx*/);
  K_MESH::Polyhedron<UNKNOWN>::metrics2(dt, crd, twoPH.PGs, pgs1, nb_pgs1, v2, G, false/*not all cvx*/);

  if (::fabs(v1) < V * GRmin) //null or two small (doesn't worth a cut)
    return false;
  if (::fabs(v2) < V* GRmin) //null or too small (doewnt worth a cut)
    return false;
  
  // check also they are not patholical bits
  ngon_unit orient;
  E_Int err = ngon_type::build_orientation_ngu<TriangulatorType>(crd, twoPH, orient);//fixme hpc : should be deduced from the input PH orientation
  E_Int res = K_MESH::Polyhedron<UNKNOWN>::is_pathological(dt, crd, twoPH.PGs, pgs0, nb_pgs0, orient.get_facets_ptr(0));
  if (res != dCONCAVITY_TO_SPLIT && res != dPATHO_PH_NONE)
    return false;

  res = K_MESH::Polyhedron<UNKNOWN>::is_pathological(dt, crd, twoPH.PGs, pgs1, nb_pgs1, orient.get_facets_ptr(1));
  if (res != dCONCAVITY_TO_SPLIT && res != dPATHO_PH_NONE)
    return false;

  E_Float THRESHOLD = 5.e-2;
  
  E_Float minA1, maxA1;
  E_Int maxAPG11, maxAPG12;
  err = K_MESH::Polyhedron<UNKNOWN>::min_max_angles(crd, twoPH.PGs, pgs0, nb_pgs0, false/*i.e open cell is error*/, orient.get_facets_ptr(0), minA1, maxA1, maxAPG11, maxAPG12);
  if (err) return false; //open cell
  if ((minA1 < minA) && ::fabs(minA1) < THRESHOLD) return false; // degen (first cond is there to allow splitting cells that have bad minA upon entry
  if ((maxA1 > maxA) && ::fabs(2.*NUGA::PI - maxA1) < THRESHOLD) return false; // same comment
  
  E_Float minA2, maxA2;
  E_Int maxAPG21, maxAPG22;
  err = K_MESH::Polyhedron<UNKNOWN>::min_max_angles(crd, twoPH.PGs, pgs1, nb_pgs1, false/*i.e open cell is error*/, orient.get_facets_ptr(1), minA2, maxA2, maxAPG21, maxAPG22);
  if (err) return false; //open cell
  if ((minA2 < minA) && ::fabs(minA2) < THRESHOLD) return false; // degen (first cond is there to allow splitting cells that have bad minA upon entry
  if ((maxA2 > maxA) && ::fabs(2.*NUGA::PI - maxA2) < THRESHOLD) return false; // same comment

  K_MESH::Polyhedron<0> PH0(twoPH, 0);
  E_Float fluxvec[3];
  PH0.flux(crd, orient.get_facets_ptr(0), fluxvec);
  E_Float f = sqrt(NUGA::sqrNorm<3>(fluxvec));
  E_Float s = PH0.surface(crd);
  f /= s;
  if (f > Fluxmax) 
  {
    std::cout << "rejected by flux" << std::endl;
    return false;
  }

  K_MESH::Polyhedron<0> PH1(twoPH, 1);
  PH1.flux(crd, orient.get_facets_ptr(1), fluxvec);
  f = sqrt(NUGA::sqrNorm<3>(fluxvec));
  s = PH0.surface(crd);
  f /= s;
  if (f > Fluxmax) 
  {
    std::cout << "rejected by flux" << std::endl;
    return false;
  }
  
  
  // count edges sharing number
  
//  for (E_Int i=0; i < 2; ++i)
//  {
//    const E_Int* pgs = twoPH.PHs.get_facets_ptr(i);
//    E_Int nb_pg = twoPH.PHs.stride(i);
//
//    std::map<K_MESH::NO_Edge, E_Int> edge_to_count;
//    std::map<K_MESH::NO_Edge, E_Int>::iterator it;
//    K_MESH::NO_Edge E;
//    
//    for (E_Int j = 0; j < nb_pg; ++j)
//    {
//      E_Int PGi = *(pgs + j) - 1;
//    
//      E_Int* pN = twoPH.PGs.get_facets_ptr(PGi);
//      E_Int nb_nodes = twoPH.PGs.stride(PGi);
//      for (E_Int n = 0; n < nb_nodes; ++n)
//      {
//        E_Int ni = *(pN + n);
//        E_Int nj = *(pN + (n + 1) % nb_nodes);
//        E.setNodes(ni, nj);
//        it = edge_to_count.find(E);
//        if (it == edge_to_count.end())
//          edge_to_count.insert(std::make_pair(E, 1));
//        else
//        {
//          ++it->second;
//          if (it->second > 2) return false;
//        }
//      }
//    }
//  }
  
  //
  

    return true;
  }

#endif
