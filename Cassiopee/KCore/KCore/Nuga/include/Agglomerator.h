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

#ifndef NUGA_AGGLOMERATOR_H
#define NUGA_AGGLOMERATOR_H

#include "Nuga/include/DynArray.h"
#include "Nuga/include/ngon_t.hxx"
#include "Nuga/include/Splitter.h"

#ifdef DEBUG_AGGLOMERATOR
#include "Nuga/include/medit.hxx"
#endif

#define NONEVAL -2
#define UNCHANGED -1

namespace NUGA
{

  class Agglomerator
  {

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
    
  public: //WRAPPERS
   
    ///
    template<typename TriangulatorType>
    inline static void agglomerate_small_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Float vmin, E_Float vratio, ngon_type& ngo, E_Int& nb_aggs, E_Float angle_threshold=1.e-12, int method=0);
    ///
    template<typename TriangulatorType>
    inline static void agglomerate_non_star_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, ngon_type& ngo, E_Int& nb_aggs, double angle_threshold = 1.e-12);

    //
    template<typename TriangulatorType>
    inline static void shell_agglomerate_small_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Float vmin, E_Float vratio, ngon_type& ngo, E_Int& nb_aggs);
  
    ///
    inline static void agglomerate_phs_having_pgs(const K_FLD::FloatArray& crd, ngon_type& ngi, const E_Int* PGlist, E_Int sz, std::vector<E_Int>& pgnids);
    
    template<typename TriangulatorType>
    inline static void collapse_uncomputable_pgs(K_FLD::FloatArray& crd, ngon_type& ngio);
    
    inline static E_Int collapse_pgs(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids);
    inline static E_Int collapse_pgs2(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids);

    // PROTO
    template<typename TriangulatorType>
    inline static E_Int collapse_small_tetras(K_FLD::FloatArray& crd, ngon_type& ngio, double vmin, double vratio);
    template<typename TriangulatorType>
    inline static E_Int collapse_small_tetras2(K_FLD::FloatArray& crd, ngon_type& ngio, double vmin, double vratio);

  public:
    /// agglomerate superfluous polygons (multiply-shared by the same polyhedra. within the flatness tolerance only for area-computable polygons)
    inline static void simplify_phs
      (const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, bool process_externals,
       ngon_type& ngo, ngon_unit& oriento, ngon_unit& phneighborso, const Vector_t<E_Int>* PHlist = nullptr, const Vector_t<E_Int>* skipPGlist = nullptr);

    /// 
    template<typename TriangulatorType>
    inline static void agglomerate_phs_NEW(const K_FLD::FloatArray& crd,
                                           const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
                                           ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, double angle_threshold, int enforce_reflex_criteria_and_or_badagglo_allowance);

    /// 
    template<typename TriangulatorType>
    inline static void agglomerate_phs_OLD(const K_FLD::FloatArray& crd,
                                           const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
                                           ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, bool force, double angle_threshold);

    ///
    template<typename TriangulatorType>
    inline static void agglomerate_phs (const K_FLD::FloatArray& crd, 
                                        const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
                                        ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, double angle_threshold, int enforce_reflex_mode);

    /// TRYING TO FOCUSE MORE ON REFLEX SITUATION
    template<typename TriangulatorType>
    inline static void agglomerate_phs2(const K_FLD::FloatArray& crd,
      const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
      ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, double angle_threshold, int enforce_reflex_mode);
  
  
  private:
      
    inline static void __simplify_phs
      (const K_FLD::FloatArray& crd, const ngon_type& ngi, E_Int PHi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, std::vector<bool> PGprocess,
      ngon_unit& gagg_pgs, std::vector<E_Int>& nids, ngon_unit& wagg_pgs, std::map<E_Int, std::vector<E_Int> >& wneigh_to_faces);

  };

  ///
  void NUGA::Agglomerator::simplify_phs
  (const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, bool process_externals,
   ngon_type& ngo, ngon_unit& oriento, ngon_unit& phneighborso, const Vector_t<E_Int>* PHlist, const Vector_t<E_Int>* skipPGlist)
  {
    typedef std::map<E_Int, std::vector<E_Int> > map_t;
    map_t neigh_to_faces;
    map_t::const_iterator it;

    std::vector<E_Int> nids;
    ngon_unit lagg_pgs, gagg_pgs;

    nids.clear();
    nids.resize(ngi.PGs.size(), NONEVAL);

    std::vector<bool> toprocess(ngi.PGs.size(), true);
    if (!process_externals)
    {
      for (E_Int i = 0; i < ngi.PHs.size(); ++i)
      {
        E_Int nb_neighs = phneighborsi.stride(i);
        const E_Int* neighs = phneighborsi.get_facets_ptr(i);
        const E_Int* pgs = ngi.PHs.get_facets_ptr(i);

        for (E_Int n = 0; n < nb_neighs; ++n)
        {
          E_Int PHn = *(neighs + n);
          if (PHn != IDX_NONE) continue;

          E_Int PGi = *(pgs + n) - 1;
          toprocess[PGi] = false;
        }
      }
    }

    if (skipPGlist != nullptr)
    {
      for (size_t i = 0; i < skipPGlist->size(); ++i)
      {
        E_Int PGi = (*skipPGlist)[i];
        toprocess[PGi] = false;
      }
    }
    
#ifdef DEBUG_AGGLOMERATOR   
    std::cout << "simplify_phs : initial nb of pgs : " << ngi.PGs.size() << std::endl;
#endif

    if (PHlist)
    {
      size_t nb_phs = PHlist->size();
      for (size_t i = 0; i < nb_phs; ++i)
      {
        E_Int PHi = (*PHlist)[i];
        NUGA::Agglomerator::__simplify_phs
          (crd, ngi, PHi, orienti, phneighborsi, angular_max, toprocess, gagg_pgs, nids, lagg_pgs, neigh_to_faces);
      }
    }
    else
    {
      E_Int nb_phs = ngi.PHs.size();
      for (E_Int PHi = 0; PHi < nb_phs; ++PHi)
      {
        NUGA::Agglomerator::__simplify_phs
          (crd, ngi, PHi, orienti, phneighborsi, angular_max, toprocess, gagg_pgs, nids, lagg_pgs, neigh_to_faces);
      }
    }
    
#ifdef DEBUG_AGGLOMERATOR
  {
    E_Int idmax=-1;
    for (size_t i=0; i < nids.size(); ++i) if (nids[i] != IDX_NONE) idmax = std::max(nids[i], idmax);
    assert (idmax == gagg_pgs.size()-1);
    K_FLD::IntArray cnto;
    //ngon_type ng(gagg_pgs);
    //ng.export_to_array(cnto);
    //medith::write("agg.plt", crd, cnto, "NGON");
    medith::write("agged", crd, gagg_pgs);
  }
#endif
    
    ngo = ngi;
    
    //if (gagg_pgs.size() == 0)
      //return;

    // append new pgs
    E_Int nb_pgs_init = ngo.PGs.size();
    ngo.PGs.append(gagg_pgs);
    
#ifdef DEBUG_AGGLOMERATOR
  {
    K_FLD::IntArray cnto;
    ngon_type ng(ngo.PGs);
    ng.export_to_array(cnto);
    //medith::write("aggandinit.plt", crd, cnto, "NGON");
    medith::write("aggandinit", crd, ngo.PGs);
  }
#endif
   
    // sync nids
    size_t nsz = nids.size();
    for (size_t i=0; i < nsz;++i)
    {
      E_Int& ni = nids[i];
      ni = (ni == IDX_NONE) ? IDX_NONE : (ni <0/*UNCHANGED OR NONEVAL*/) ? i : ni + nb_pgs_init;
    }
    
    //NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH("ph220.plt", crd, ngi, 220);

    // sync PHs
    //E_Int nb_phs; // = ngo.PHs.size();
    std::vector<E_Int> dummy;
    ngo.PHs.remove_facets(nids, dummy);
    //nb_phs = ngo.PHs.size();
    
    Vector_t<E_Int> pgnids, phnids;
    ngo.remove_unreferenced_pgs(pgnids, phnids);

#ifdef DEBUG_AGGLOMERATOR
    std::cout << "simplify_phs : final nb of pgs : " << ngo.PGs.size() << std::endl;
    bool ok = ngo.is_consistent(crd.cols());
#endif

    // clean superfluous nodes
    ngon_type::simplify_pgs(ngo, crd, process_externals);
  }

  ///
  template<typename TriangulatorType>
  void NUGA::Agglomerator::agglomerate_phs2
  (const K_FLD::FloatArray& crd, 
   const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
   ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, double angle_threshold, int enforce_reflex_criteria)
  {
    ngo.clear();
    oriento.clear();
    nb_aggs = 0;

    E_Float concave_threshold = angle_threshold;
    E_Float convex_threshold = angle_threshold;
    
    // Fast return
    if (PHlist.size() * ngi.PHs.size() * ngi.PGs.size() == 0) return;
  
    TriangulatorType dt;
  
    std::vector<bool> frozen(ngi.PHs.size(), false); //process one agglo at a time
    std::vector<bool> bad(ngi.PHs.size(), false);
    for (size_t k = 0; k < PHlist.size(); ++k) bad[PHlist[k]] = true;
    
    std::vector<E_Int> shared_pgs;
    std::map<K_MESH::NO_Edge, E_Float> reflex_edges;
    std::set<K_MESH::NO_Edge> convex_edges;

    ngon_unit all_agg, cur_agg, all_ori, cur_ori, best_agg, best_ori;

    E_Int nb_phs = PHlist.size();
    //std::cout << "NB OF BAD CELLS TO AGGLOMERATE : " << nb_phs << std::endl;
    //E_Int agg_id(0);
    for (E_Int ii = 0; ii < nb_phs; ++ii)
    {
      E_Int i = PHlist[ii];

      if (frozen[i]) continue;

      E_Int bestn = IDX_NONE;
      //E_Float worst_reflex_a = -1.;

#ifdef DEBUG_AGGLOMERATOR 
      //size_t nbfmax=0;
      //E_Float smax=0.;
#endif
      E_Float qmax=0.;
      E_Int nb_reflex_edges_1(0);
      E_Int min_delta_reflex = 1000;
      
      E_Int nb_neighs = ngi.PHs.stride(i);
      const E_Int* neighs = neighborsi.get_facets_ptr(i);
      const E_Int* pgsi = ngi.PHs.get_facets_ptr(i);
      
      bool conc1;
      K_MESH::Polyhedron<UNKNOWN>::is_concave
                (crd, ngi.PGs, pgsi, nb_neighs, false/*deal with closed PH*/, orienti.get_facets_ptr(i), conc1, nb_reflex_edges_1, concave_threshold, convex_threshold);

      E_Float vi;
      K_MESH::Polyhedron<0>::volume<DELAUNAY::Triangulator>(crd, ngi.PGs, ngi.PHs.get_facets_ptr(i), ngi.PHs.stride(i), vi, true);

      // the best is the one sharing the most number of faces
      for (E_Int n = 0; (n < nb_neighs); ++n)
      {
        E_Int j = *(neighs + n);
        if (j == IDX_NONE)
          continue;

        if (bad[j]) continue;

        if (frozen[j]) 
        {
          // continue; fixme : hack to have the targeted logic (but bad perfo) : 
          // a cell has to be aggregated with its best mate, i.e the one with most surface in common, not the best available one.
          bestn = IDX_NONE;
          break;
        }

        const E_Int* pgsj = ngi.PHs.get_facets_ptr(j);
        E_Int nb_pgsj = ngi.PHs.stride(j);

        cur_agg.clear();
        cur_ori.clear();

        K_MESH::Polyhedron<UNKNOWN>::merge_two_phs(crd, ngi.PGs, pgsi, nb_neighs, pgsj, nb_pgsj, orienti.get_facets_ptr(i), orienti.get_facets_ptr(j), cur_agg, cur_ori);

        reflex_edges.clear();
        convex_edges.clear();

        E_Int res = K_MESH::Polyhedron<UNKNOWN>::is_pathological
                (dt, crd, ngi.PGs, cur_agg.get_facets_ptr(0), cur_agg.stride(0), cur_ori.get_facets_ptr(0), reflex_edges, convex_edges, concave_threshold, convex_threshold);
        
        //if (res == dOPEN_PHS) continue; // just wrong
        if (res != dPATHO_PH_NONE) continue;

        // IS IT BETTER ?
        E_Float vj;
        K_MESH::Polyhedron<0>::volume<DELAUNAY::Triangulator>(crd, ngi.PGs, pgsj, nb_pgsj, vj, true);

        // -1 number of reflex edges MUST not increase
       
        //compuTe now for j
        bool conc2;
        E_Int nb_reflex_edges_2(0);
        K_MESH::Polyhedron<UNKNOWN>::is_concave
                  (crd, ngi.PGs, pgsj, nb_pgsj, false/*deal with closed PH*/, orienti.get_facets_ptr(j), conc2, nb_reflex_edges_2, concave_threshold, convex_threshold);
        
        E_Int nb_reflex_new = (E_Int)reflex_edges.size();
        
        // -2 number of shared faces and shared surface SHOULD be max
    
        ngon_type::shared_faces_list(ngi, i, j, shared_pgs);
        size_t nbf = shared_pgs.size();
        // compute shared surface
        E_Float s=0.;
        for (size_t f = 0; f < nbf; ++f)
        {
          E_Int PGi = shared_pgs[f]-1;
          s += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, ngi.PGs.get_facets_ptr(PGi), ngi.PGs.stride(PGi), 1);
        }
        
        //fixme proto : compute the total surface of PHi
        E_Float stot=0.;
        for (E_Int f = 0; f < nb_neighs; ++f)
        {
          E_Int PGi = *(pgsi+f)-1;
          stot += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, ngi.PGs.get_facets_ptr(PGi), ngi.PGs.stride(PGi), 1);
        }
        
        E_Float face_ratio = E_Float(nbf) / E_Float (nb_neighs);
        E_Float surface_ratio = s /stot;
        E_Float reflex_ratio = (1. + E_Float(nb_reflex_edges_1 + nb_reflex_edges_2)) / (1. + E_Float(nb_reflex_new));
        //E_Float volume_ratio = 1. / (1. + ::fabs(vi / vj));
        E_Int delta_reflex = nb_reflex_new - (nb_reflex_edges_1 + nb_reflex_edges_2);
        
        // -3 worst reflex angle SHOULD decrease
        
        E_Float reflex_a = NUGA::PI;

        if (!reflex_edges.empty())
        {
          for (std::map<K_MESH::NO_Edge, E_Float>::const_iterator it = reflex_edges.begin(); it != reflex_edges.end(); ++it) 
            reflex_a = std::min(reflex_a, it->second);
        }

        reflex_a /= NUGA::PI;
        E_Int N = 3;
        if (enforce_reflex_criteria == 1)
        {
          // prioritize contributions
          reflex_ratio *= 2;// sqrt(reflex_ratio); // increase impact
          N++;
          reflex_a *= 2;// sqrt(worst_reflex_a);   // increase impact
          N++;
          //volume_ratio *= volume_ratio;            // decrease impact
          //N++
        }

        bool is_better = false;
        E_Float q = (face_ratio + surface_ratio + reflex_a) /  (E_Float)N;// +reflex_ratio); // *volume_ratio;

        if (delta_reflex < min_delta_reflex)
          is_better = true;
        else if (delta_reflex == min_delta_reflex)
          is_better = (q > qmax);

        if (is_better)
        {
          bestn = n;

          //smax = s;
          //nbfmax = nbf;
          best_agg = cur_agg;
          best_ori = cur_ori;
          //worst_reflex_a = reflex_a;
          min_delta_reflex = delta_reflex;
          qmax = q;
          
#ifdef DEBUG_AGGLOMERATOR
         //if (ii == 5 || ii == 6)
          /*{
           std::ostringstream o;
           o << "D://slandier//DATA//tmp//agglo//agg_" << ii << "_with_neigh_" << n << ".tp";
           K_FLD::IntArray cnto;
           ngon_type ng(ngi.PGs, best_agg);
           ng.export_to_array(cnto);
           //medith::write(o.str().c_str(), crd, cnto, "NGON");
           tp::write(o.str().c_str(), crd, cnto, "NGON");
         }*/
#endif
        } 
      }

      if (bestn == IDX_NONE) continue;
      
      E_Int j = *(neighs+bestn);
      frozen[i] = frozen[j] = true;

      //std::cout << " agglomerate " << agg_id++ << " was " << i << " and " << j << std::endl;
      
#ifdef DEBUG_AGGLOMERATOR
      std::cout << "AGGLOMERATION : " << i << " with " << j << std::endl;
      //if (ii==5 || ii ==6)
      {
        std::ostringstream o;
        o << "/visu/slandier/ATOM/TURBINE_REFROIDIE/ROTOR2/RXXX/TEST_10/best_agg_" << ii << ".tp";//";//D://slandier//DATA//tmp//agglo//best_agg_" << i << ".tp";
        K_FLD::IntArray cnto;
        ngon_type ng(ngi.PGs, best_agg);
        ng.export_to_array(cnto);
        medith::write(o.str().c_str(), crd, cnto, "NGON");
        //tp::write(o.str().c_str(), crd, cnto, "NGON");
      }
#endif

      all_agg.append(best_agg);
      all_ori.append(best_ori);
    }

    ngo.PGs = ngi.PGs;
    ngo.PHs = all_agg;
    
    //std::cout << ngo.PHs.size() << " small cells have been agglomerated." << std::endl;

    // now add untouched ones
    for (size_t i = 0; i < frozen.size(); ++i)
    {
      if (!frozen[i])
      {
        ngo.PHs.add(ngi.PHs.stride(i), ngi.PHs.get_facets_ptr(i));
        if (!ngi.PHs._type.empty()) ngo.PHs._type.push_back(ngi.PHs._type[i]);
      }
    }
    
    ngo.PGs.updateFacets();
    ngo.PHs.updateFacets();
    
    //std::cout << "effective nb of phs in agglomerated ngo : " << ngo.PHs.size() << std::endl;

    std::vector<E_Int> pgnids, phnids;
    ngo.remove_unreferenced_pgs(pgnids, phnids);

  }
  
  ///
  template<typename TriangulatorType>
  void NUGA::Agglomerator::agglomerate_phs
  (const K_FLD::FloatArray& crd, 
   const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
   ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, double angle_threshold, int enforce_reflex_criteria)
  {
    ngo.clear();
    oriento.clear();
    nb_aggs = 0;

    E_Float concave_threshold = angle_threshold;
    E_Float convex_threshold = angle_threshold;
    
    // Fast return
    if (PHlist.size() * ngi.PHs.size() * ngi.PGs.size() == 0) return;
  
    TriangulatorType dt;
  
    std::vector<bool> frozen(ngi.PHs.size(), false); //process one agglo at a time
    std::vector<bool> bad(ngi.PHs.size(), false);
    for (size_t k = 0; k < PHlist.size(); ++k) bad[PHlist[k]] = true;
    
    std::vector<E_Int> shared_pgs;
    std::map<K_MESH::NO_Edge, E_Float> reflex_edges;
    std::set<K_MESH::NO_Edge> convex_edges;

    ngon_unit all_agg, cur_agg, all_ori, cur_ori, best_agg, best_ori;

    E_Int nb_phs = PHlist.size();
    //std::cout << "NB OF BAD CELLS TO AGGLOMERATE : " << nb_phs << std::endl;
    //E_Int agg_id(0);
    for (E_Int ii = 0; ii < nb_phs; ++ii)
    {
      E_Int i = PHlist[ii];

      if (frozen[i])
        continue;

      E_Int bestn = IDX_NONE;

#ifdef DEBUG_AGGLOMERATOR 
      //size_t nbfmax=0;
      //E_Float smax=0.;
#endif
      E_Float qmax=0.;
      E_Int nb_reflex_edges_1(0);
      
      E_Int nb_neighs = ngi.PHs.stride(i);
      const E_Int* neighs = neighborsi.get_facets_ptr(i);
      const E_Int* pgsi = ngi.PHs.get_facets_ptr(i);
      
      bool conc1;
      K_MESH::Polyhedron<UNKNOWN>::is_concave
                (crd, ngi.PGs, pgsi, nb_neighs, false/*deal with closed PH*/, orienti.get_facets_ptr(i), conc1, nb_reflex_edges_1, concave_threshold, convex_threshold);

      E_Float vi;
      K_MESH::Polyhedron<0>::volume<DELAUNAY::Triangulator>(crd, ngi.PGs, ngi.PHs.get_facets_ptr(i), ngi.PHs.stride(i), vi, true);

      // the best is the one sharing the most number of faces
      for (E_Int n = 0; (n < nb_neighs); ++n)
      {
        E_Int j = *(neighs + n);
        if (j == IDX_NONE)
          continue;

        if (bad[j]) continue;

        if (frozen[j]) 
        {
          // continue; fixme : hack to have the targeted logic (but bad perfo) : 
          // a cell has to be aggregated with its best mate, i.e the one with most surface in common, not the best available one.
          bestn = IDX_NONE;
          break;
        }

        const E_Int* pgsj = ngi.PHs.get_facets_ptr(j);
        E_Int nb_pgsj = ngi.PHs.stride(j);

        cur_agg.clear();
        cur_ori.clear();

        K_MESH::Polyhedron<UNKNOWN>::merge_two_phs(crd, ngi.PGs, pgsi, nb_neighs, pgsj, nb_pgsj, orienti.get_facets_ptr(i), orienti.get_facets_ptr(j), cur_agg, cur_ori);

        reflex_edges.clear();
        convex_edges.clear();

        E_Int res = K_MESH::Polyhedron<UNKNOWN>::is_pathological
                (dt, crd, ngi.PGs, cur_agg.get_facets_ptr(0), cur_agg.stride(0), cur_ori.get_facets_ptr(0), reflex_edges, convex_edges, concave_threshold, convex_threshold);
        
        //if (res == dOPEN_PHS) continue; // just wrong
        if (res != dPATHO_PH_NONE) continue;

        // IS IT BETTER ?
        E_Float vj;
        K_MESH::Polyhedron<0>::volume<DELAUNAY::Triangulator>(crd, ngi.PGs, pgsj, nb_pgsj, vj, true);

        // -1 number of reflex edges MUST not increase
       
        //compuTe now for j
        bool conc2;
        E_Int nb_reflex_edges_2(0);
        K_MESH::Polyhedron<UNKNOWN>::is_concave
                  (crd, ngi.PGs, pgsj, nb_pgsj, false/*deal with closed PH*/, orienti.get_facets_ptr(j), conc2, nb_reflex_edges_2, concave_threshold, convex_threshold);
        
        E_Int nb_reflex_new = (E_Int)reflex_edges.size();
        
        // -2 number of shared faces and shared surface SHOULD be max
    
        ngon_type::shared_faces_list(ngi, i, j, shared_pgs);
        size_t nbf = shared_pgs.size();
        // compute shared surface
        E_Float s=0.;
        for (size_t f = 0; f < nbf; ++f)
        {
          E_Int PGi = shared_pgs[f]-1;
          s += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, ngi.PGs.get_facets_ptr(PGi), ngi.PGs.stride(PGi), 1);
        }
        
        //fixme proto : compute the total surface of PHi
        E_Float stot=0.;
        for (E_Int f = 0; f < nb_neighs; ++f)
        {
          E_Int PGi = *(pgsi+f)-1;
          stot += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, ngi.PGs.get_facets_ptr(PGi), ngi.PGs.stride(PGi), 1);
        }
        
        E_Float face_ratio = E_Float(nbf) / E_Float (nb_neighs);
        E_Float surface_ratio = s /stot;
        E_Float reflex_ratio = (1. + E_Float(nb_reflex_edges_1 + nb_reflex_edges_2)) / (1. + E_Float(nb_reflex_new));
        //E_Float volume_ratio = 1. / (1. + ::fabs(vi / vj));
        
        // -3 worst reflex angle SHOULD decrease
        
        E_Float worst_reflex_a = NUGA::PI;

        if (!reflex_edges.empty())
        {
          for (std::map<K_MESH::NO_Edge, E_Float>::const_iterator it = reflex_edges.begin(); it != reflex_edges.end(); ++it) 
            worst_reflex_a = std::min(worst_reflex_a, it->second);

          //std::cout << "worst angle : " << worst_reflex_a << std::endl;
        }

        worst_reflex_a /= NUGA::PI;

        E_Float q{ 0 };
        
        if (enforce_reflex_criteria == 1)
        {
          // prioritize contributions
          reflex_ratio = sqrt(reflex_ratio);       // increase impact
          worst_reflex_a = sqrt(worst_reflex_a);   // increase impact
          //volume_ratio *= volume_ratio;            // decrease impact

          q = face_ratio * surface_ratio * worst_reflex_a * reflex_ratio; // *volume_ratio;
        }
        else if ((enforce_reflex_criteria == 0) || (enforce_reflex_criteria == 2))
        {
          E_Float F = 1.;
          E_Float N = 4.;
          if (enforce_reflex_criteria == 2) {F = 2; N = 6;}

          q = (face_ratio + surface_ratio + F * (worst_reflex_a + reflex_ratio) )  / N;
        }

      
        //maximum surface is best. if equality, count the number of shared faces.
        //bool is_better = (s > smax) || ((::fabs(s-smax) < EPSILON) && (nbf > nbfmax));
        bool is_better = (q > qmax);
      
        if (is_better)
        {
          bestn = n;

          //smax = s;
          //nbfmax = nbf;
          best_agg = cur_agg;
          best_ori = cur_ori;
          
          qmax = q;
          
#ifdef DEBUG_AGGLOMERATOR
         {
           std::ostringstream o;
           o << "D://slandier//DATA//tmp//agglo//agg_" << ii << "_with_neigh_" << n << ".tp";
           K_FLD::IntArray cnto;
           ngon_type ng(ngi.PGs, best_agg);
           ng.export_to_array(cnto);
           //medith::write(o.str().c_str(), crd, cnto, "NGON");
           tp::write(o.str().c_str(), crd, cnto, "NGON");
         }
#endif
        } 
      }

      if (bestn == IDX_NONE) continue;
      
      E_Int j = *(neighs+bestn);
      frozen[i] = frozen[j] = true;

      //std::cout << " agglomerate " << agg_id++ << " was " << i << " and " << j << std::endl;
      
#ifdef DEBUG_AGGLOMERATOR
      //std::cout << "AGGLOMERATION : " << i << " with " << j << std::endl;
      //if (ii==73)
      {
        std::ostringstream o;
        o << "D://slandier//DATA//tmp//agglo//best_agg_" << ii << ".tp";
        K_FLD::IntArray cnto;
        ngon_type ng(ngi.PGs, best_agg);
        ng.export_to_array(cnto);
        //medith::write(o.str().c_str(), crd, cnto, "NGON");
        tp::write(o.str().c_str(), crd, cnto, "NGON");
      }
#endif

      all_agg.append(best_agg);
      all_ori.append(best_ori);
    }

    ngo.PGs = ngi.PGs;
    ngo.PHs = all_agg;
    
    //std::cout << ngo.PHs.size() << " small cells have been agglomerated." << std::endl;

    // now add untouched ones
    for (size_t i = 0; i < frozen.size(); ++i)
    {
      if (!frozen[i])
      {
        ngo.PHs.add(ngi.PHs.stride(i), ngi.PHs.get_facets_ptr(i));
        if (!ngi.PHs._type.empty()) ngo.PHs._type.push_back(ngi.PHs._type[i]);
      }
    }
    
    ngo.PGs.updateFacets();
    ngo.PHs.updateFacets();
    
    //std::cout << "effective nb of phs in agglomerated ngo : " << ngo.PHs.size() << std::endl;

    std::vector<E_Int> pgnids, phnids;
    ngo.remove_unreferenced_pgs(pgnids, phnids);

  }

  ///
  template<typename TriangulatorType>
  void NUGA::Agglomerator::agglomerate_phs_NEW
  (const K_FLD::FloatArray& crd,
    const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
    ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, double angle_threshold, int enforce_reflex_criteria_and_or_badagglo_allowance)
  {
    // enforce_reflex_criteria_and_or_badagglo_allowance == 0 => do not enforce + forbid bad agglo
    // enforce_reflex_criteria_and_or_badagglo_allowance == 1 =>        enforce + forbid bad agglo
    // enforce_reflex_criteria_and_or_badagglo_allowance == 2 => do not enforce + allow bad agglo
    // enforce_reflex_criteria_and_or_badagglo_allowance == 3 =>        enforce + allow bad agglo

    ngo.clear();
    oriento.clear();
    nb_aggs = 0;

    E_Float concave_threshold = angle_threshold;
    E_Float convex_threshold = angle_threshold;

    // Fast return
    if (PHlist.size() * ngi.PHs.size() * ngi.PGs.size() == 0) return;

    TriangulatorType dt;

    std::vector<bool> frozen(ngi.PHs.size(), false); //process one agglo at a time
    std::vector<bool> bad(ngi.PHs.size(), false);
    for (size_t k = 0; k < PHlist.size(); ++k) bad[PHlist[k]] = true;

    std::vector<E_Int> shared_pgs;
    std::map<K_MESH::NO_Edge, E_Float> reflex_edges, reflex_edges_i, reflex_edges_j;
    std::set<K_MESH::NO_Edge> convex_edges, convex_edges_i, convex_edges_j;

    ngon_unit all_agg, cur_agg, all_ori, cur_ori, best_agg, best_ori;

    E_Int nb_phs = PHlist.size();
    //std::cout << "NB OF BAD CELLS TO AGGLOMERATE : " << nb_phs << std::endl;
    //E_Int agg_id(0);
    for (E_Int ii = 0; ii < nb_phs; ++ii)
    {
      E_Int i = PHlist[ii];

      if (frozen[i])
        continue;

      E_Int bestn = IDX_NONE;
      //E_Float worst_reflex_a = -1.;

#ifdef DEBUG_AGGLOMERATOR 
      //size_t nbfmax=0;
      //E_Float smax=0.;
#endif
      E_Float qmax = 0.;
      E_Int nb_reflex_edges_1(0);
      E_Int min_delta_reflex = 1000;
      
      E_Int nb_neighs = ngi.PHs.stride(i);
      const E_Int* neighs = neighborsi.get_facets_ptr(i);
      const E_Int* pgsi = ngi.PHs.get_facets_ptr(i);

      reflex_edges_i.clear();
      convex_edges_i.clear();

      E_Int resi = K_MESH::Polyhedron<UNKNOWN>::is_pathological
      (dt, crd, ngi.PGs, pgsi, nb_neighs, orienti.get_facets_ptr(i), reflex_edges_i, convex_edges_i, concave_threshold, convex_threshold);

      nb_reflex_edges_1 = convex_edges_i.size();

      E_Float vi;
      K_MESH::Polyhedron<0>::volume<DELAUNAY::Triangulator>(crd, ngi.PGs, ngi.PHs.get_facets_ptr(i), ngi.PHs.stride(i), vi, true);

      // the best is the one sharing the most number of faces
      for (E_Int n = 0; (n < nb_neighs); ++n)
      {
        E_Int j = *(neighs + n);
        if (j == IDX_NONE) continue;

        if ((bad[j]) && (enforce_reflex_criteria_and_or_badagglo_allowance == 0 || enforce_reflex_criteria_and_or_badagglo_allowance == 1))
          continue;
       
        if (frozen[j])
        {
          // continue; fixme : hack to have the targeted logic (but bad perfo) : 
          // a cell has to be aggregated with its best mate, i.e the one with most surface in common, not the best available one.
          bestn = IDX_NONE;
          break;
        }

        const E_Int* pgsj = ngi.PHs.get_facets_ptr(j);
        E_Int nb_pgsj = ngi.PHs.stride(j);

        cur_agg.clear();
        cur_ori.clear();

        reflex_edges.clear();
        convex_edges.clear();

        //compuTe now for j
        E_Int nb_reflex_edges_2(0);

        reflex_edges_j.clear();
        convex_edges_j.clear();

        E_Int resj = K_MESH::Polyhedron<UNKNOWN>::is_pathological
        (dt, crd, ngi.PGs, pgsj, nb_pgsj, orienti.get_facets_ptr(j), reflex_edges_j, convex_edges_j, concave_threshold, convex_threshold);

        K_MESH::Polyhedron<UNKNOWN>::merge_two_phs(crd, ngi.PGs, pgsi, nb_neighs, pgsj, nb_pgsj, orienti.get_facets_ptr(i), orienti.get_facets_ptr(j), cur_agg, cur_ori);

#ifdef DEBUG_AGGLOMERATOR
        {
          std::ostringstream o;
          o << "D://slandier//DATA//tmp//agglo//curagg_" << ii << "_with_neigh_" << n << ".tp";
          K_FLD::IntArray cnto;
          ngon_type ng(ngi.PGs, cur_agg);
          ng.export_to_array(cnto);
          //medith::write(o.str().c_str(), crd, cnto, "NGON");
          tp::write(o.str().c_str(), crd, cnto, "NGON");
        }
#endif


        E_Int res = K_MESH::Polyhedron<UNKNOWN>::is_pathological
        (dt, crd, ngi.PGs, cur_agg.get_facets_ptr(0), cur_agg.stride(0), cur_ori.get_facets_ptr(0), reflex_edges, convex_edges, concave_threshold, convex_threshold);

        if (res == dOPEN_PHS) continue; // just wrong

        bool stop = true;

        if (res != dPATHO_PH_NONE)
        {
          if (dPATHO_PH_NONE == resi && dPATHO_PH_NONE == resj) // clean comps, dirty aggregate
          {
            stop = true;
          }
          else if (res == resi || res == resj)  // same pathology as one of the component : the added part is not necessarily making things worst
          {
            E_Float worst_reflex_a_i{ NUGA::PI }, worst_reflex_a_j{ NUGA::PI };
            for (auto it = reflex_edges_i.begin(); it != reflex_edges_i.end(); ++it)
              worst_reflex_a_i = std::min(worst_reflex_a_i, it->second);
            for (auto it = reflex_edges_j.begin(); it != reflex_edges_j.end(); ++it)
              worst_reflex_a_j = std::min(worst_reflex_a_j, it->second);
            E_Float worst_ref_a{ NUGA::PI };
            for (auto it = reflex_edges.begin(); it != reflex_edges.end(); ++it)
              worst_ref_a = std::min(worst_ref_a, it->second);

            if (worst_ref_a >= std::min(worst_reflex_a_i, worst_reflex_a_j)) // the gluing edges does not increase the pathology compared to separated components => OK 
              stop = false;
          }
        }
        else
          stop = false;


        if (stop == true) continue;

        // IS IT BETTER ?
        E_Float vj;
        K_MESH::Polyhedron<0>::volume<DELAUNAY::Triangulator>(crd, ngi.PGs, pgsj, nb_pgsj, vj, true);

        // -1 number of reflex edges MUST not increase

        E_Int nb_reflex_new = (E_Int)reflex_edges.size();

        // -2 number of shared faces and shared surface SHOULD be max

        ngon_type::shared_faces_list(ngi, i, j, shared_pgs);
        size_t nbf = shared_pgs.size();
        // compute shared surface
        E_Float s = 0.;
        for (size_t f = 0; f < nbf; ++f)
        {
          E_Int PGi = shared_pgs[f] - 1;
          s += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, ngi.PGs.get_facets_ptr(PGi), ngi.PGs.stride(PGi), 1);
        }

        //fixme proto : compute the total surface of PHi
        E_Float stot = 0.;
        for (E_Int f = 0; f < nb_neighs; ++f)
        {
          E_Int PGi = *(pgsi + f) - 1;
          stot += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, ngi.PGs.get_facets_ptr(PGi), ngi.PGs.stride(PGi), 1);
        }

        E_Float face_ratio = E_Float(nbf) / E_Float(nb_neighs);
        E_Float surface_ratio = s / stot;
        E_Float reflex_ratio = (1. + E_Float(nb_reflex_edges_1 + nb_reflex_edges_2)) / (1. + E_Float(nb_reflex_new));
        //E_Float volume_ratio = 1. / (1. + ::fabs(vi / vj));
        E_Int delta_reflex = nb_reflex_new - (nb_reflex_edges_1 + nb_reflex_edges_2);
        
        // -3 worst reflex angle SHOULD decrease

        E_Float worst_reflex_a = NUGA::PI;

        if (!reflex_edges.empty())
        {
          for (std::map<K_MESH::NO_Edge, E_Float>::const_iterator it = reflex_edges.begin(); it != reflex_edges.end(); ++it)
            worst_reflex_a = std::min(worst_reflex_a, it->second);

          //std::cout << "worst angle : " << worst_reflex_a << std::endl;
        }

        worst_reflex_a /= NUGA::PI;
        E_Int N = 3;
       
        if (enforce_reflex_criteria_and_or_badagglo_allowance == 1 || enforce_reflex_criteria_and_or_badagglo_allowance == 3)
        {
          // prioritize contributions
          reflex_ratio *= 2;// sqrt(reflex_ratio); // increase impact
          N++;
          worst_reflex_a *= 2;// sqrt(worst_reflex_a);   // increase impact
          N++;
          //volume_ratio *= volume_ratio;            // decrease impact
          //N++
        }

        bool is_better = false;
        E_Float q = (face_ratio + surface_ratio + worst_reflex_a) /  (E_Float)N;// +reflex_ratio); // *volume_ratio;

        if (delta_reflex < min_delta_reflex)
          is_better = true;
        else if (delta_reflex == min_delta_reflex)
          is_better = (q > qmax);

        if (is_better)
        {
          bestn = n;

          //smax = s;
          //nbfmax = nbf;
          best_agg = cur_agg;
          best_ori = cur_ori;
          //worst_reflex_a = reflex_a;
          min_delta_reflex = delta_reflex;
          qmax = q;
          
#ifdef DEBUG_AGGLOMERATOR
         //if (ii == 5 || ii == 6)
          /*{
           std::ostringstream o;
           o << "D://slandier//DATA//tmp//agglo//agg_" << ii << "_with_neigh_" << n << ".tp";
           K_FLD::IntArray cnto;
           ngon_type ng(ngi.PGs, best_agg);
           ng.export_to_array(cnto);
           //medith::write(o.str().c_str(), crd, cnto, "NGON");
           tp::write(o.str().c_str(), crd, cnto, "NGON");
         }*/
#endif
        } 
      }

      if (bestn == IDX_NONE) continue;
      
      E_Int j = *(neighs+bestn);
      frozen[i] = frozen[j] = true;

      //std::cout << " agglomerate " << agg_id++ << " was " << i << " and " << j << std::endl;
      
#ifdef DEBUG_AGGLOMERATOR
      std::cout << "AGGLOMERATION : " << i << " with " << j << std::endl;
      //if (ii==5 || ii ==6)
      {
        std::ostringstream o;
        o << "D://slandier//DATA//tmp//agglo//best_agg_" << ii << ".tp";//";//D://slandier//DATA//tmp//agglo//best_agg_" << i << ".tp";
        K_FLD::IntArray cnto;
        ngon_type ng(ngi.PGs, best_agg);
        ng.export_to_array(cnto);
        medith::write(o.str().c_str(), crd, cnto, "NGON");
        //tp::write(o.str().c_str(), crd, cnto, "NGON");
      }
#endif

      all_agg.append(best_agg);
      all_ori.append(best_ori);
    }

    ngo.PGs = ngi.PGs;
    ngo.PHs = all_agg;
    
    //std::cout << ngo.PHs.size() << " small cells have been agglomerated." << std::endl;

    // now add untouched ones
    for (size_t i = 0; i < frozen.size(); ++i)
    {
      if (!frozen[i])
      {
        ngo.PHs.add(ngi.PHs.stride(i), ngi.PHs.get_facets_ptr(i));
        if (!ngi.PHs._type.empty()) ngo.PHs._type.push_back(ngi.PHs._type[i]);
      }
    }
    
    ngo.PGs.updateFacets();
    ngo.PHs.updateFacets();
    
    //std::cout << "effective nb of phs in agglomerated ngo : " << ngo.PHs.size() << std::endl;

    std::vector<E_Int> pgnids, phnids;
    ngo.remove_unreferenced_pgs(pgnids, phnids);
  }

  ///
  template<typename TriangulatorType>
  void NUGA::Agglomerator::agglomerate_phs_OLD(const K_FLD::FloatArray& crd,
    const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
    ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, bool force, double angle_threshold)
  {
    ngo.clear();
    oriento.clear();
    nb_aggs = 0;

    E_Float concave_threshold = angle_threshold;
    E_Float convex_threshold = angle_threshold;
    
    //std::cout << "ngi : initial nb of phs : " << ngi.PHs.size() << std::endl;
  
    // Fast return
    if (PHlist.size() * ngi.PHs.size() * ngi.PGs.size() == 0) return;
  
    TriangulatorType dt;
  
    std::vector<bool> frozen(ngi.PHs.size(), false); //process one agglo at a time
    
    std::vector<E_Int> shared_pgs;
    std::map<K_MESH::NO_Edge, E_Float> reflex_edges;
    std::set<K_MESH::NO_Edge> convex_edges;

    ngon_unit all_agg, cur_agg, all_ori, cur_ori, best_agg, best_ori;

    E_Int nb_phs = PHlist.size();
    //E_Int agg_id(0);
    for (E_Int ii = 0; ii < nb_phs; ++ii)
    {
      E_Int i = PHlist[ii];

      if (frozen[i])
        continue;

      E_Int bestn = IDX_NONE;

#ifdef DEBUG_AGGLOMERATOR 
      //size_t nbfmax=0;
      //E_Float smax=0.;
#endif
      E_Float qmax=0.;
      E_Int nb_reflex_edges_1(0);
      
      E_Int nb_neighs = ngi.PHs.stride(i);
      const E_Int* neighs = neighborsi.get_facets_ptr(i);
      const E_Int* pgsi = ngi.PHs.get_facets_ptr(i);
      
      bool conc1;
      K_MESH::Polyhedron<UNKNOWN>::is_concave
                (crd, ngi.PGs, pgsi, nb_neighs, false/*deal with closed PH*/, orienti.get_facets_ptr(i), conc1, nb_reflex_edges_1, concave_threshold, convex_threshold);

      // the best is the one sharing the most number of faces
      for (E_Int n = 0; (n < nb_neighs); ++n)
      {
        E_Int j = *(neighs + n);
        if (j == IDX_NONE)
          continue;

        if (frozen[j]) 
        {
          // continue; fixme : hack to have the targeted logic (but bad perfo) : 
          // a cell has to be aggregated with its best mate, i.e the one with most surface in common, not the best available one.
          bestn = IDX_NONE;
          break;
        }

        const E_Int* pgsj = ngi.PHs.get_facets_ptr(j);
        E_Int nb_pgsj = ngi.PHs.stride(j);

        cur_agg.clear();
        cur_ori.clear();

        K_MESH::Polyhedron<UNKNOWN>::merge_two_phs(crd, ngi.PGs, pgsi, nb_neighs, pgsj, nb_pgsj, orienti.get_facets_ptr(i), orienti.get_facets_ptr(j), cur_agg, cur_ori);

        reflex_edges.clear();
        convex_edges.clear();

        E_Int res = K_MESH::Polyhedron<UNKNOWN>::is_pathological
                (dt, crd, ngi.PGs, cur_agg.get_facets_ptr(0), cur_agg.stride(0), cur_ori.get_facets_ptr(0), reflex_edges, convex_edges, concave_threshold, convex_threshold);
        
        if (res == dOPEN_PHS) continue; // just wrong

        E_Float Q0;
        if (res == dPATHO_PH_NONE) Q0 = 10.; //give more credit to non-patho result
        else Q0 = 0.;
               
        if (!force && res != dPATHO_PH_NONE) continue;

        // IS IT BETTER ?

        // -1 number of reflex edges MUST not increase
       
        //compuTe now for j
        bool conc2;
        E_Int nb_reflex_edges_2(0);
        K_MESH::Polyhedron<UNKNOWN>::is_concave
                  (crd, ngi.PGs, pgsj, nb_pgsj, false/*deal with closed PH*/, orienti.get_facets_ptr(j), conc2, nb_reflex_edges_2, concave_threshold, convex_threshold);
        
        E_Int nb_reflex_new = (E_Int)reflex_edges.size();

        if ( (  nb_reflex_new > nb_reflex_edges_1 + nb_reflex_edges_2) ) continue;
        
        // -2 number of shared faces and shared surface SHOULD be max
    
        ngon_type::shared_faces_list(ngi, i, j, shared_pgs);
        size_t nbf = shared_pgs.size();
        // compute shared surface
        E_Float s=0.;
        for (size_t f = 0; f < nbf; ++f)
        {
          E_Int PGi = shared_pgs[f]-1;
          s += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, ngi.PGs.get_facets_ptr(PGi), ngi.PGs.stride(PGi), 1);
        }
        
        //fixme proto : compute the total surface of PHi
        E_Float stot=0.;
        for (E_Int f = 0; f < nb_neighs; ++f)
        {
          E_Int PGi = *(pgsi+f)-1;
          stot += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, ngi.PGs.get_facets_ptr(PGi), ngi.PGs.stride(PGi), 1);
        }
        
        E_Float face_ratio = E_Float(nbf) / E_Float (nb_neighs);
        E_Float surface_ratio = s /stot;
        
        // -3 worst reflex angle SHOULD decrease
        
        E_Float worst_reflex_a=2.*NUGA::PI;

        if (!reflex_edges.empty())
        {
          for (std::map<K_MESH::NO_Edge, E_Float>::const_iterator it = reflex_edges.begin(); it != reflex_edges.end(); ++it) 
            worst_reflex_a = std::min(worst_reflex_a, it->second);
        
          //std::cout << "worst angle : " << worst_reflex_a << std::endl;
        }
        
        E_Float q = Q0 + face_ratio * surface_ratio * worst_reflex_a;
      
        //maximum surface is best. if equality, count the number of shared faces.
        //bool is_better = (s > smax) || ((::fabs(s-smax) < EPSILON) && (nbf > nbfmax));
        bool is_better = (q > qmax);
      
        if (is_better)
        {
          bestn = n;

          //smax = s;
          //nbfmax = nbf;
          best_agg = cur_agg;
          best_ori = cur_ori;
          
          qmax = q;
          
#ifdef DEBUG_AGGLOMERATOR
          //if (ii==73)
          /*{
            K_FLD::IntArray cnto;
            ngon_type ng(ngi.PGs, best_agg);
            ng.export_to_array(cnto);
            std::ostringstream o;
            o << "agg_" << n << ".tp";
            medith::write(o.str().c_str(), crd, cnto, "NGON"); 
            
            E_Float isoG[3];
            typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
            acrd_t acrd(crd);
            std::vector<E_Int> phnodes;
            K_MESH::Polyhedron<UNKNOWN>::unique_nodes(ng.PGs, ng.PHs.get_facets_ptr(0), ng.PHs.stride(0), phnodes);
            K_MESH::Polyhedron<UNKNOWN>::iso_barycenter(acrd, &phnodes[0], phnodes.size(), 1, isoG);
            E_Float v, centroid[3];
            K_MESH::Polyhedron<UNKNOWN>::metrics2<DELAUNAY::Triangulator>(dt, crd, ng.PGs, ng.PHs.get_facets_ptr(0), ng.PHs.stride(0), v, centroid);
        
            K_FLD::IntArray cn(2,1,0);
            cn(1,0)=1;
            K_FLD::FloatArray crd;
            crd.pushBack(isoG, isoG+3);
            crd.pushBack(centroid, centroid+3);

            medith::write("baryTocentroid.mesh", crd, cn, "BAR"); 
            
          }*/
#endif
        } 
      }

      if (bestn == IDX_NONE) continue;
      
      E_Int j = *(neighs+bestn);
      frozen[i] = frozen[j] = true;

      //std::cout << " agglomerate " << agg_id++ << " was " << i << " and " << j << std::endl;
      
#ifdef DEBUG_AGGLOMERATOR
      //std::cout << "AGGLOMERATION : " << i << " with " << j << std::endl;
      //if (ii==73)
      {
        std::ostringstream o;
        o << "best_agg_" << ii << ".tp";
        K_FLD::IntArray cnto;
        ngon_type ng(ngi.PGs, best_agg);
        ng.export_to_array(cnto);
        medith::write(o.str().c_str(), crd, cnto, "NGON");
      }
#endif

      all_agg.append(best_agg);
      all_ori.append(best_ori);
    }

    ngo.PGs = ngi.PGs;
    ngo.PHs = all_agg;
    
    //std::cout << ngo.PHs.size() << " small cells have been agglomerated." << std::endl;

    // now add untouched ones
    for (size_t i = 0; i < frozen.size(); ++i)
    {
      if (!frozen[i])
      {
        ngo.PHs.add(ngi.PHs.stride(i), ngi.PHs.get_facets_ptr(i));
        if (!ngi.PHs._type.empty()) ngo.PHs._type.push_back(ngi.PHs._type[i]);
      }
    }
    
    ngo.PGs.updateFacets();
    ngo.PHs.updateFacets();
    
    //std::cout << "effective nb of phs in agglomerated ngo : " << ngo.PHs.size() << std::endl;

    std::vector<E_Int> pgnids, phnids;
    ngo.remove_unreferenced_pgs(pgnids, phnids);
  }
  
  ///
  template<typename TriangulatorType>
  void NUGA::Agglomerator::agglomerate_small_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Float vmin, E_Float vratio, ngon_type& ngo, E_Int& nb_aggs, double angle_threshold, int method)
  {
    ngon_unit neighborsi;
    ngi.build_ph_neighborhood(neighborsi);
    
    ngon_unit orienti;
    ngon_type::build_orientation_ngu<TriangulatorType>(crd, ngi, orienti); //WARNING : ngi previous types are lost
    
    std::vector<E_Int> PHlist;
    ngon_type::detect_bad_volumes<TriangulatorType>(crd, ngi, neighborsi, vmin, vratio, PHlist);
    
    if (PHlist.empty()) 
    {
      ngo = ngi;
      return;
    }

    ngon_unit oriento;
    if (method == 0) // OLD METHOD : bad cells can agglomerate together / reflex edges CANNOT increase
      NUGA::Agglomerator::agglomerate_phs_OLD<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, false/*force*/, angle_threshold);
    else if (method == 1) // method 0 + force mode allows pathological aggregates
      NUGA::Agglomerator::agglomerate_phs_OLD<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, true/*force*/, angle_threshold);
    else if (method == 2) // NEW METHOD : bad cells CANNOT agglomerate together / reflex edges CAN increase => reflex criterion / volume criterion
      NUGA::Agglomerator::agglomerate_phs<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, angle_threshold, 0/*enforce_reflex_criteria*/);
    else if (method == 3) // method 2 + putting more weight to reflex criteria than volume criterion
      NUGA::Agglomerator::agglomerate_phs<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, angle_threshold, 1/*enforce_reflex_criteria*/);
    else if (method == 4) // method 2 + putting more weight to reflex criteria than volume criterion
      NUGA::Agglomerator::agglomerate_phs2<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, angle_threshold, 0/*enforce_reflex_criteria*/);
    else if (method == 5) // method 2 + putting more weight to reflex criteria than volume criterion
      NUGA::Agglomerator::agglomerate_phs2<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, angle_threshold, 1/*enforce_reflex_criteria*/);
    else if (method == 6) // method 2 + putting more weight to reflex criteria than volume criterion
      NUGA::Agglomerator::agglomerate_phs<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, angle_threshold, 2/*enforce_reflex_criteria*/);
    else if (method == 7) // NEW METHOD :
      NUGA::Agglomerator::agglomerate_phs_NEW<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, angle_threshold, 0/*do not enforce_reflex_criteria*/);
    else if (method == 8) // method 2 + putting more weight to reflex criteria than volume criterion
      NUGA::Agglomerator::agglomerate_phs_NEW<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, angle_threshold, 1/*enforce_reflex_criteria*/);
    else if (method == 9) // method 2 + putting more weight to reflex criteria than volume criterion
      NUGA::Agglomerator::agglomerate_phs_NEW<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, angle_threshold, 2/*do not enforce_reflex_criteria and allow between-bad agglo*/);
    else if (method == 10) // method 2 + putting more weight to reflex criteria than volume criterion
      NUGA::Agglomerator::agglomerate_phs_NEW<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, angle_threshold, 3/*enforce_reflex_criteria and allow between-bad agglo*/);
  }
  
  ///
  template<typename TriangulatorType>
  void NUGA::Agglomerator::agglomerate_non_star_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, ngon_type& ngo, E_Int& nb_aggs, double angle_threshold)
  {
    ngon_unit orienti;
    ngon_type::build_orientation_ngu<TriangulatorType>(crd, ngi, orienti);
    
    std::vector<ngon_type::ePathoPH> PHtypes;
    ngon_type::detect_pathological_PHs<TriangulatorType>(crd, ngi, orienti, PHtypes);
    
    std::vector<E_Int> PHlist;
    for (size_t i=0; i< PHtypes.size(); ++i) if (PHtypes[i] == ngon_type::CONCAVITY_TO_SPLIT)PHlist.push_back(i);
    
    if (PHlist.empty()) 
    {
      ngo = ngi;
      return;
    }
    
    ngon_unit neighborsi;
    ngi.build_ph_neighborhood(neighborsi);

    ngon_unit oriento;
    NUGA::Agglomerator::agglomerate_phs<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, angle_threshold, 1/*give more weight to reflex criteria*/); //we choose here to agglomerate non star iff it gets better,ie. non-star (any pathological aggragte is rejected)
  }

  ///
  template<typename TriangulatorType>
  void NUGA::Agglomerator::shell_agglomerate_small_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Float vmin, E_Float vratio, ngon_type& ngo, E_Int& nb_aggs)
  {
    ngon_unit neighborsi;
    ngi.build_ph_neighborhood(neighborsi);
    
    ngon_unit orienti;
    ngon_type::build_orientation_ngu<TriangulatorType>(crd, ngi, orienti); //WARNING : ngi previous types are lost
    
    std::vector<E_Int> PHlist;
    ngon_type::detect_bad_volumes<TriangulatorType>(crd, ngi, neighborsi, vmin, vratio, PHlist);
    
    ngo = ngi;

    if (PHlist.empty())
      return;

    Vector_t<bool> in_shell(ngi.PGs.size(), false), wprocessed/*externalized to not reallocate it each time*/;
    Vector_t<E_Int> to_remove;
    ngon_unit aggs;

    for (size_t i=0; i < PHlist.size(); ++i)
    {
      //std::cout << "compute shell for " << PHlist[i] << std::endl;
      Vector_t<E_Int> shellPHs, boundPGs;
      ngon_type::ph_shell(ngi, PHlist[i], neighborsi, shellPHs, boundPGs, wprocessed);

      // ensure current shell is not colliding with another shell
      bool valid=true;
      for (size_t k=0; (k < shellPHs.size()) && valid; ++k)
        valid &= (in_shell[shellPHs[k]] == false); //check if not involved in a shell

      if (!valid) continue;

      for (size_t k=0; (k < shellPHs.size()) && valid; ++k)
      {
        in_shell[shellPHs[k]] = true;

      }

      to_remove.insert(to_remove.end(), ALL(shellPHs));
      aggs.add(boundPGs.size(), &boundPGs[0]);
    }

    aggs.updateFacets();
 
    Vector_t<E_Int> phnids;
    ngo.PHs.remove_entities(to_remove, phnids);
    ngo.PHs.append(aggs);

  }
  
  ///  
  void NUGA::Agglomerator::agglomerate_phs_having_pgs
  (const K_FLD::FloatArray& crd, ngon_type& ngi, const E_Int* PGlist, E_Int sz, std::vector<E_Int>& pgnids)
  {
    K_FLD::IntArray F2E;
    ngi.build_noF2E(F2E);
    
    ngon_unit mergedPH;
    
    E_Int nb_phs = ngi.PHs.size();
    
    Vector_t<E_Int> newPGlist, PHtoremove;
    const E_Int* pPGlist(&PGlist[0]);
    
    bool carry_on(true);
    while (carry_on)
    {
      std::vector<E_Int> phnids;
      nb_phs = ngi.PHs.size(); //it changes
      K_CONNECT::IdTool::init_inc(phnids, nb_phs);
      
      newPGlist.clear();
      
      for (E_Int i=0; i < sz; ++i)
      {
        E_Int PGi = pPGlist[i];

        if (PGi >= F2E.size())
          continue;

        E_Int lPH = F2E(0, PGi);
        E_Int rPH = F2E(1, PGi);
        
        if (lPH == rPH) continue; //the PG is gone already in a merge
        if (lPH == IDX_NONE || rPH == IDX_NONE) continue; // discard skin PGs for aggregattion
        
        if ( (phnids[lPH] != lPH) || (phnids[rPH] != rPH) ) //at least one ph has been already processed in that iteration
        {
          //so keep that face for the nex iteration
          newPGlist.push_back(PGi);
          continue;
        }

        const E_Int* pgs1 = ngi.PHs.get_facets_ptr(lPH);
        E_Int nb_pgs1 = ngi.PHs.stride(lPH);
        
        const E_Int* pgs2 = ngi.PHs.get_facets_ptr(rPH);
        E_Int nb_pgs2 = ngi.PHs.stride(rPH);
        
        K_MESH::Polyhedron<UNKNOWN>::merge_two_phs(crd, ngi.PGs, pgs1, nb_pgs1, pgs2, nb_pgs2, mergedPH);
        
        //assert (mergedPH.size());

        E_Int newPHi = ngi.PHs.size();
        ngi.PHs.append(mergedPH);
        
        phnids[rPH] = phnids[lPH] = newPHi;
        PHtoremove.push_back(lPH);
        PHtoremove.push_back(rPH);
        //phnids.push_back(newPHi);
      }
      
      carry_on = !newPGlist.empty();
      if (carry_on)
      {
        pPGlist = &newPGlist[0];
        sz = newPGlist.size();
        ngi.PHs.updateFacets();
        K_FLD::IntArray::changeIndices(F2E, phnids);

#ifdef DEBUG_AGGLOMERATOR
//        ngon_unit PGs;
//        for (size_t i=0; i< sz; ++i)
//          PGs.add(ngi.PGs.stride(newPGlist[i]), ngi.PGs.get_facets_ptr(newPGlist[i]));
//        PGs.updateFacets();
//        NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGs("new_ps", crd, PGs);
//        
//        ngon_type ngo = ngi;
//        
//        ngo.PHs.updateFacets();
//    
//        Vector_t<E_Int> pgnids, phnids;
//    
//        ngo.PHs.remove_entities(PHtoremove, phnids);
//        ngo.remove_unreferenced_pgs(pgnids, phnids); 
//    
//        K_FLD::IntArray cnto;
//        ngo.export_to_array(cnto);
//        medith::write("current.plt", crd, cnto, "NGON");
//
//        for (size_t i=0; i < F2E.cols(); ++i)
//        {
//          if (F2E(0,i) == F2E(1,i)){
//            std::cout << "baffle" << std::endl;
//            NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH("PHbaffle.plt", crd, ngi, F2E(0,i));
//          }
//        }
#endif   
      }
    }
    
    ngi.PHs.updateFacets();
    
    Vector_t<E_Int> phnids;
    pgnids.clear();
    
    ngi.PHs.remove_entities(PHtoremove, phnids);
    ngi.remove_unreferenced_pgs(pgnids, phnids); 
  }
 
  ///
  template<typename TriangulatorType>
  void NUGA::Agglomerator::collapse_uncomputable_pgs(K_FLD::FloatArray& crd, ngon_type& ngio)
  {
    std::vector<E_Int> PGlist;
    std::vector<ngon_type::ePathoPG> flagPG;
    ngon_type::detect_uncomputable_pgs<TriangulatorType>(crd, ngio.PGs, flagPG);

#ifdef DEBUG_AGGLOMERATOR
    std::vector<E_Int> badPGs;
    for (size_t i=0; i < flagPG.size(); ++i) if (flagPG[i] != ngon_type::NONE/*== ngon_type::SPIKE*/) badPGs.push_back(i);
    
    NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGs(crd, ngio.PGs, badPGs);
    
    ngon_type tmp = ngio;
    tmp.keep_PHs_having_PGs(badPGs);
    K_FLD::IntArray ctmp;
    tmp.export_to_array(ctmp);
    medith::write("bad_phs.plt", crd, ctmp, "NGON");
#endif
    
    ///  COLLAPSING THE SMALLEST EDGE => NUGA::Agglomerator::collapse_pgs2
    E_Int nb_pathos(0);
    for (size_t i=0; i < flagPG.size(); ++i)
    {
      if (flagPG[i] != ngon_type::PATHO_PG_NONE)
      {
        ++nb_pathos;
        PGlist.push_back(i);
      }
    }
    
    if (nb_pathos == 0) 
      std::cout << "OK : COMPUTABLE MESH" << std::endl;
    else
      std::cout << PGlist.size() << " uncomputable pgs will be fixed over " <<  nb_pathos << " pathologies." << std::endl;
    
    //do the healing
    NUGA::Agglomerator::collapse_pgs2(crd, ngio, PGlist);
    
  }
  
  /// collapse to the barycenter
  E_Int NUGA::Agglomerator::collapse_pgs(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids)
  {
    if (pgids.empty())
      return 0;
    
    E_Float G[3];
    std::vector<E_Int> nid(crd.cols());
    for (size_t i=0; i < nid.size(); ++i) nid[i]=i;
    //E_Int ng_pgs = pgids.size();
  
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray > acrd_t;
    acrd_t acrd(crd);
  
    for (size_t i=0; i < pgids.size(); ++i)
    {
      E_Int PGi = pgids[i];
        
      K_MESH::Polygon::iso_barycenter<acrd_t, 3>(acrd, ng.PGs.get_facets_ptr(PGi), ng.PGs.stride(PGi), 1, G);
    
      E_Int Ni = ng.PGs.get_facet(PGi, 0)-1;
        for (E_Int j=0; j < ng.PGs.stride(PGi); ++j)
          nid[ng.PGs.get_facet(PGi, j)-1]=Ni;
    
      crd(0,Ni)=G[0]; crd(1,Ni)=G[1]; crd(2,Ni)=G[2];
    }
  
    ng.PGs.change_indices(nid);
  
    ngon_type::clean_connectivity(ng, crd, -1/*ngon_dim*/, EPSILON/*tolerance*/, false/*remove_dup_phs*/, true/*do_omp*/);
  
  return 0;
}
  /// collaspe the smallest edge keeping the oldest node.
  E_Int NUGA::Agglomerator::collapse_pgs2(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids)
  {
    if (pgids.empty())
      return 0;
    
    //E_Float G[3];
    std::vector<E_Int> nids(crd.cols());
    for (size_t i=0; i < nids.size(); ++i) nids[i]=i;
    //E_Int ng_pgs = pgids.size();
  
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray > acrd_t;
    acrd_t acrd(crd);
  
    for (size_t i=0; i<pgids.size(); ++i)
    {
      E_Int PGi = pgids[i];
      const E_Int* nodes = ng.PGs.get_facets_ptr(PGi);
      E_Int nb_nodes = ng.PGs.stride(PGi);
      E_Float d2min = NUGA::FLOAT_MAX;
      E_Int Nsource(IDX_NONE), Ntarget(IDX_NONE);
      
      // for the smallest edge of the current polygon, assign the non-moved and biggest id to the smallest one
      for (E_Int n=0; n<nb_nodes; ++n)
      {
        E_Int Ni = *(nodes + n) - 1;
        E_Int Nj = *(nodes + (n+1) % nb_nodes) - 1;
        E_Float d2 = NUGA::sqrDistance(crd.col(Ni), crd.col(Nj), 3);
        
        if (d2 < d2min)
        {
          Ntarget = (nids[Ni] == Ni && nids[Nj] == Nj) ? MIN(Ni, Nj) : (nids[Ni] == Ni) ? Ni : Nj;
          d2min=d2;
          Nsource = (Ntarget == Ni) ? Nj : Ni;//the other one
        }
      }
      
      if (nids[Nsource] != Nsource || nids[Ntarget] != Ntarget) continue; //both have moved so this polygon is discarded
      //degenerate an edge of this polygon the pg by setting all the nodes to 
      nids[Nsource]= Ntarget;
    }
  
    ng.PGs.change_indices(nids);
  
    ngon_type::clean_connectivity(ng, crd, 3/*ngon_dim*/, 0./*tolerance*/, false/*remove_dup_phs*/, true/*do_omp*/);
  
  return 0;
}

  ///
  template<typename TriangulatorType>
  E_Int NUGA::Agglomerator::collapse_small_tetras(K_FLD::FloatArray& crd, ngon_type& ngio, double vmin, double vratio)
  {
    ngon_unit neighborsi;
    ngio.build_ph_neighborhood(neighborsi);

    ngon_unit orienti;
    ngon_type::build_orientation_ngu<TriangulatorType>(crd, ngio, orienti); //WARNING : ngi previous types are lost

    std::vector<E_Int> PHlist;
    ngon_type::detect_bad_volumes<TriangulatorType>(crd, ngio, neighborsi, vmin, vratio, PHlist);

    if (PHlist.empty())
    {
      return 0;
    }

    std::vector<E_Int> nids;
    K_CONNECT::IdTool::init_inc(nids, crd.cols());

    for (size_t i = 0; i < PHlist.size(); ++i)
    {
      E_Int PHi = PHlist[i];

      const E_Int* pf = ngio.PHs.get_facets_ptr(PHi);
      int npf = ngio.PHs.stride(PHi);
      
      if (!K_MESH::Polyhedron<0>::is_of_type<K_MESH::Tetrahedron>(ngio.PGs, pf, npf)) continue;

      //std::cout << "collapsing a tetra " << std::endl;

      //find tha max id : more chances to be an X point, so more chance to lie on the X interface :
      // it it better if entities crash on the interface to not deform it
      E_Int Ntarget=-1; 
      for (size_t p = 0; p < 2; ++p) // only needed to cross 2 faces to pass through the four nodes
      {
        E_Int* pnodes = ngio.PGs.get_facets_ptr(pf[p]-1);
        for (size_t k=0; k < 3; ++k)
        {
          E_Int Ni = pnodes[k]-1;
          Ntarget = (Ntarget < Ni) ? Ni : Ntarget;
        }
      }

      // now collapse the tetra on Ntrget
      for (size_t p = 0; p < 2; ++p) // only needed to cross 2 faces to pass through the four nodes
      {
        E_Int* pnodes = ngio.PGs.get_facets_ptr(pf[p]-1);
        for (size_t k=0; k < 3; ++k)
          nids[pnodes[k]-1]=Ntarget;
      }
    }

    // update the pointers to point to the leaves
    for (size_t i =0; i < nids.size(); ++i)
    {
      E_Int Fi=nids[i];
      while (Fi != nids[Fi])Fi=nids[Fi];
      nids[i]=Fi;
    }

    ngio.PGs.change_indices(nids);

    ngon_type::clean_connectivity(ngio, crd, -1/*ngon_dim*/, EPSILON/*tolerance*/, false/*remove_dup_phs*/, false/*do_omp*/);

    return 0;

  }

  ///
  template<typename TriangulatorType>
  E_Int NUGA::Agglomerator::collapse_small_tetras2(K_FLD::FloatArray& crd, ngon_type& ngio, double vmin, double vratio)
  {
    ngon_unit neighborsi;
    ngio.build_ph_neighborhood(neighborsi);

    ngon_unit orienti;
    ngon_type::build_orientation_ngu<TriangulatorType>(crd, ngio, orienti); //WARNING : ngi previous types are lost

    std::vector<E_Int> PHlist;
    ngon_type::detect_bad_volumes<TriangulatorType>(crd, ngio, neighborsi, vmin, vratio, PHlist);
    //K_CONNECT::IdTool::init_inc(PHlist, 38);

    if (PHlist.empty())
    {
      return 0;
    }

    std::vector<E_Int> nids;
    K_CONNECT::IdTool::init_inc(nids, crd.cols());

    std::map<K_MESH::NO_Edge, std::deque<E_Int> > edge_to_refine_nodes;
    ngon_unit splitpgs;
    std::vector<E_Int> splitoids;
    ///
    for (size_t i = 0; i < PHlist.size(); ++i)
    {
      E_Int PHi = PHlist[i];

      const E_Int* pf = ngio.PHs.get_facets_ptr(PHi);
      int npf = ngio.PHs.stride(PHi);

      if (!K_MESH::Polyhedron<0>::is_of_type<K_MESH::Tetrahedron>(ngio.PGs, pf, npf)) continue;

      K_MESH::Tetrahedron th4(ngio.PGs, pf);
      K_MESH::Tetrahedron::eShapeType shape = th4.shape_type(crd);

#ifdef DEBUG_AGGLOMERATOR
      std::ostringstream o;
      o << "shape_" << PHi;
      if (shape == K_MESH::Tetrahedron::eShapeType::REGULAR) o << "_REG.mesh";
      else if (shape == K_MESH::Tetrahedron::eShapeType::SLICE1) o << "_SLICE1.mesh";
      else if (shape == K_MESH::Tetrahedron::eShapeType::SLICE2) o << "_SLICE2.mesh";
      else if (shape == K_MESH::Tetrahedron::eShapeType::SPIKE) o << "_SPIKE.mesh";
      else if (shape == K_MESH::Tetrahedron::eShapeType::DELTA) o << "_DELTA.mesh";
      else if (shape == K_MESH::Tetrahedron::eShapeType::KNIFE1) o << "_KNIFE1.mesh";
      else if (shape == K_MESH::Tetrahedron::eShapeType::KNIFE2) o << "_KNIFE2.mesh";
      
      medith::write(o.str().c_str(), crd, ngio, PHi);
      medith::write("check.mesh",crd, th4.nodes(), 4, 0);
#endif

      if (shape == K_MESH::Tetrahedron::REGULAR)
      {
        //break;
        //find tha max id : more chances to be an X point, so more chance to lie on the X interface :
        // it it better if entities crash on the interface to not deform it
        //E_Int Ntarget = -1;
        /*for (size_t p = 0; p < 2; ++p) // only needed to cross 2 faces to pass through the four nodes
        {
          E_Int* pnodes = ngio.PGs.get_facets_ptr(pf[p] - 1);
          for (size_t k = 0; k < 3; ++k)
          {
            E_Int Ni = pnodes[k] - 1;
            Ntarget = (Ntarget < Ni) ? Ni : Ntarget;
          }
        }

        // now collapse the tetra on Ntrget
        for (size_t p = 0; p < 2; ++p) // only needed to cross 2 faces to pass through the four nodes
        {
          E_Int* pnodes = ngio.PGs.get_facets_ptr(pf[p] - 1);
          for (size_t k = 0; k < 3; ++k)
            nids[pnodes[k] - 1] = Ntarget;
        }*/
        E_Float Lmin2 = NUGA::FLOAT_MAX;
        E_Int* pf = ngio.PHs.get_facets_ptr(PHi);
        E_Int Nmin=0, Nmax=0;
        for (size_t f = 0; f < 4; ++f)
        {
          E_Int* pnodes = ngio.PGs.get_facets_ptr(pf[f] - 1);

          for (size_t n = 0; n < 3; ++n)
          {
            E_Int Ni = pnodes[n] - 1;
            E_Int Nj = pnodes[(n+1)%3] - 1;

            E_Float d2 = NUGA::sqrDistance(crd.col(Ni), crd.col(Nj), 3);

            if (d2 < Lmin2)
            {
              Nmin = Ni;
              Nmax = Nj;
              Lmin2 = d2;
            }
          }
        }

        if (Nmin > Nmax) std::swap(Nmin, Nmax);
        nids[Nmin] = Nmax;

      }
      else if (shape == K_MESH::Tetrahedron::SPIKE)
      {
        // collapse the smallest side (first 3 nodes)
        E_Int Nmax = std::max(th4.node(0), th4.node(1));
        Nmax = std::max(Nmax, th4.node(2));

        nids[th4.node(0)] = Nmax;
        nids[th4.node(1)] = Nmax;
        nids[th4.node(2)] = Nmax;
      }
      else if (shape == K_MESH::Tetrahedron::SLICE1)
      {
        // collapse the smallest side (first 2 nodes)
        E_Int Nmax = std::max(th4.node(0), th4.node(1));
     
        nids[th4.node(0)] = Nmax;
        nids[th4.node(1)] = Nmax;
      }
      else if (shape == K_MESH::Tetrahedron::SLICE2)
      {
        // collapse separately each pair of smmalest edge
        E_Int Nmax = std::max(th4.node(0), th4.node(1));

        nids[th4.node(0)] = Nmax;
        nids[th4.node(1)] = Nmax;

        Nmax = std::max(th4.node(2), th4.node(3));

        nids[th4.node(2)] = Nmax;
        nids[th4.node(3)] = Nmax;
      }
      else if (shape == K_MESH::Tetrahedron::DELTA)
      {
        // create N1's sibling on edge to split [N2, N3]
        E_Int N1 = th4.node(0);
        E_Int N2 = th4.node(1);
        E_Int N3 = th4.node(2);
        E_Int N4 = th4.node(3);

#ifdef DEBUG_AGGLOMERATOR
        medith::write("delta0.mesh", crd, ngio, PHi);
        medith::write("check0.mesh", crd, th4.nodes(), 4, 0);
#endif

        E_Float lambda;
        K_MESH::Edge::edgePointMinDistance2<3>(crd.col(N2), crd.col(N3), crd.col(N1), lambda);
        assert(lambda > 0. && lambda < 1.);
        E_Float N2N3[3], Psibling[3];
        NUGA::diff<3>(crd.col(N3), crd.col(N2), N2N3);
        NUGA::sum<3>(lambda, N2N3, crd.col(N2), Psibling);
        crd.pushBack(Psibling, Psibling + 3);
        E_Int Nsibling = crd.cols() - 1;

        edge_to_refine_nodes[K_MESH::NO_Edge(N2+1, N3+1)].push_back(Nsibling+1);
        nids.push_back(N1);

        // append splitpgs : N2N3N4 has to be replaced by N2NsiblingN4 & NsiblingN3N4
        
        // find the face id in ngon to split
        E_Int Fopp = K_MESH::Tetrahedron::get_opposite_face_to_node(ngio.PGs, pf, N1);
        assert(Fopp != IDX_NONE);

        E_Int molecule[3];
        molecule[0] = N2 + 1;
        molecule[1] = Nsibling + 1;
        molecule[2] = N4 + 1;
        splitpgs.add(3, molecule);
        splitoids.push_back(Fopp);
        splitpgs._type.push_back(ngio.PGs._type[Fopp]);

        molecule[0] = Nsibling + 1;
        molecule[1] = N3 + 1;
        molecule[2] = N4 + 1;
        splitpgs.add(3, molecule);
        splitoids.push_back(Fopp);
        splitpgs._type.push_back(ngio.PGs._type[Fopp]);
      }
      else if (shape == K_MESH::Tetrahedron::KNIFE1)
      {
      }
      else if (shape == K_MESH::Tetrahedron::KNIFE2)
      {
        E_Int Nmax = std::max(th4.node(0), th4.node(1));
        nids[th4.node(0)] = Nmax;
        nids[th4.node(1)] = Nmax;

        //create Nmax's sibling on edge to split [node(2), node(3)]
        // and ling Nsibling to Nmax (this way is easier to revert back in case of error
        E_Int N2 = th4.node(2);
        E_Int N3 = th4.node(3);
        E_Float lambda;
        K_MESH::Edge::edgePointMinDistance2<3>(crd.col(N2), crd.col(N3), crd.col(Nmax), lambda);
        
        assert(lambda > 0. && lambda < 1.);
        
        E_Float N2N3[3], Psibling[3];
        NUGA::diff<3>(crd.col(N3), crd.col(N2), N2N3);
        NUGA::sum<3>(lambda, N2N3, crd.col(N2), Psibling);
        crd.pushBack(Psibling, Psibling + 3);
        E_Int Nsibling = crd.cols() - 1;
        
        edge_to_refine_nodes[K_MESH::NO_Edge(th4.node(2)+1, th4.node(3)+1)].push_back(Nsibling+1);
        nids.push_back(Nmax) ; // nids[Nsibling]=Nmax;
      }
    }

    // update the pointers to point to the leaves
    for (size_t i = 0; i < nids.size(); ++i)
    {
      E_Int Fi = nids[i];
      while (Fi != nids[Fi])Fi = nids[Fi];
      nids[i] = Fi;
    }

    //medith::write<ngon_type>("before", crd, ngio, 1);*/

    //////////////// SPLIT EDGE STAGE /////////////////////////////////
    if (true)
    {
      // propagate new ids and 1-start => NOT NEEDED BECAUSE change_indices is now called after moves validation
      /*std::map<K_MESH::NO_Edge, std::deque<int> > new_edge_to_refine_nodes;
      for (auto it : edge_to_refine_nodes)
      {
        auto & ie = it.first;

        K_MESH::NO_Edge e(nids[ie.node(0)]+1, nids[ie.node(1)]+1);

        for (size_t u = 0; u < it.second.size(); ++u)
          new_edge_to_refine_nodes[e].push_back(nids[it.second[u]] + 1);
      }
      edge_to_refine_nodes = new_edge_to_refine_nodes;*/

      // complete refinement with ends and eventually sort in between if more than 1 point
      std::vector < std::pair<E_Float, E_Int>> sorter;
      for (auto& e : edge_to_refine_nodes)
      {
        e.second.push_front(e.first.node(0));
        if (e.second.size() > 2)
          NUGA::MeshTool::reorder_nodes_on_edge<std::deque<E_Int>, 3>(crd, e.second, 1, sorter);
        e.second.push_back(e.first.node(1));
      }

      // Refine the PGs
      Vector_t<E_Int> pg_molec;
      ngon_unit refinedPGs;

      for (E_Int PGi = 0; PGi < ngio.PGs.size(); ++PGi)
      {
        ngon_type::refine_pg(ngio.PGs.get_facets_ptr(PGi), ngio.PGs.stride(PGi), edge_to_refine_nodes, pg_molec);
        refinedPGs.add(pg_molec);
      }

      // update PGs ngon unit
      refinedPGs._type = ngio.PGs._type;  // hack to preserve flags (externality)
      refinedPGs._ancEs = ngio.PGs._ancEs;// hack
      ngio.PGs = refinedPGs;
      ngio.PGs.updateFacets();
    }

    //////////////// REPLACE FACE STAGE (DELTA) /////////////////////////////////
    if (splitpgs.size() > 0)
    {
      //medith::write("delta1.mesh", crd, ngio, PHlist[0]);
      NUGA::Splitter::__split_pgs(crd, ngio, splitpgs, splitoids);
      //medith::write("delta2.mesh", crd, ngio, PHlist[0]);
      //medith::write("splitpgs.mesh", crd, splitpgs);
    }

    // NIDS VALIDATION : now validate/invalidates mooves by group for each cell : looking at the ph shell flux evolution
    ngon_type::validate_moves_by_fluxes<DELAUNAY::Triangulator>(nids, crd, ngio, neighborsi, PHlist);

    // now apply the validated moves globally
    ngio.PGs.change_indices(nids);
    ngon_type::clean_connectivity(ngio, crd, 3/*ngon_dim*/, 0./*tolerance*/, false/*remove_dup_phs*/, true/*do_omp*/);

    return 0;

  }

/// PRIVATE ///

  ///
  void NUGA::Agglomerator::__simplify_phs
  (const K_FLD::FloatArray& crd, const ngon_type& ngi, E_Int PHi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, std::vector<bool> PGprocess,
   ngon_unit& gagg_pgs, std::vector<E_Int>& nids, ngon_unit& wagg_pgs, std::map<E_Int, std::vector<E_Int> >& wneigh_to_faces)
  {
    // nids[i] == -1         =>  UNCHANGED
    // nids[i] == IDX_NONE   =>  DELETED
    // nids[i] == j          =>  j-th PG in agg_pgs
    
    typedef std::map<E_Int, std::vector<E_Int> > map_t;
    map_t::const_iterator it;

    std::vector<E_Int> lnids, ori, pgsi;

    wneigh_to_faces.clear();
 
    E_Int nb_neighs = phneighborsi.stride(PHi);
    const E_Int* neighs = phneighborsi.get_facets_ptr(PHi);
    const E_Int* pgs = ngi.PHs.get_facets_ptr(PHi);
    const E_Int* orient = orienti.get_facets_ptr(PHi);

    for (E_Int n = 0; n < nb_neighs; ++n)
    {
      E_Int PGi = *(pgs + n) - 1;
      if (!PGprocess[PGi]) continue; // include external choice

      E_Int PHn = *(neighs + n); // can be IDX_NONE in case of externals

      wneigh_to_faces[PHn].push_back(n);
    }

    if ((E_Int)wneigh_to_faces.size() == nb_neighs) //no possible agglo
      return;

    for (it = wneigh_to_faces.begin(); it != wneigh_to_faces.end(); ++it)
    {
      //const E_Int& PHn = it->first;
      const std::vector<E_Int>& common_pg_pos = it->second;
      size_t sz = common_pg_pos.size();
      if (sz == 1)
        continue;

      ori.clear();  ori.resize(sz);
      pgsi.clear(); pgsi.resize(sz);

      size_t already_done = 0;
      for (size_t k = 0; k < sz; ++k)
      {
        E_Int ni = common_pg_pos[k];
        ori[k] = orient[ni];
        pgsi[k] = pgs[ni];
        if (nids[pgsi[k] - 1] != NONEVAL) ++already_done;
        
       // std::cout << "Face  : " << pgsi[k] << " / id : " << nids[pgsi[k] - 1] << " . done already ? : " << (nids[pgsi[k] - 1] != -1) << std::endl;
      }
      
#ifdef DEBUG_AGGLOMERATOR
      assert(already_done == 0 || already_done == sz);
#endif
      
      if (already_done)
        continue;

      K_MESH::Polygon::full_agglomerate(crd, ngi.PGs, &pgsi[0], E_Int(sz), angular_max, &ori[0], wagg_pgs, lnids /*,normals*/);
      
#ifdef DEBUG_AGGLOMERATOR
      assert (lnids.size() == sz);
      const E_Int& PHn = it->first;
      if (PHn == 31923)
      {
        std::vector<E_Int> tmppg = pgsi;
        K_CONNECT::IdTool::shift(tmppg, -1);
        medith::write("PGs", crd, ngi.PGs, &tmppg, 1);
        medith::write("PHi", crd, ngi, PHi);
        medith::write("PHn", crd, ngi, PHn);
      }
#endif

      //local to global
      E_Int shft = E_Int(gagg_pgs.size());
      E_Int nb_com_pgs=E_Int(sz);

#ifdef DEBUG_AGGLOMERATOR
      assert (lnids.size() == sz);
#endif

      for (E_Int k=0; k < nb_com_pgs; ++k)
      {
        E_Int PGi = pgsi[k] - 1;
        E_Int nid = lnids[k];

        nids[PGi] = (nid == IDX_NONE) ? IDX_NONE : (nid >= nb_com_pgs) ? nid -nb_com_pgs +shft : UNCHANGED; // deleted ? agglomerated ? unchanged
        
#ifdef DEBUG_AGGLOMERATOR  
        std::string tmp = (nids[PGi] == IDX_NONE) ? "DELETED" : (nids[PGi] == UNCHANGED) ? "UNCHANGED" : "NEW";
        if (tmp != "NEW") std::cout << PGi << " status : " <<  tmp << std::endl;
        else std::cout << PGi << " status : " << nid -nb_com_pgs +shft << std::endl;
#endif
        
      }

      if (wagg_pgs.size())
      {
        gagg_pgs.append(wagg_pgs);
        gagg_pgs._type.resize(gagg_pgs.size(), INNER);
        gagg_pgs._ancEs.resize(2, gagg_pgs.size(), IDX_NONE);
      }
    }
    
#ifdef DEBUG_AGGLOMERATOR
    assert (nids.size() == ngi.PGs.size());
#endif
  }
}

#endif
