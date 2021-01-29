#ifndef NUGA_AGGLOMERATOR_H
#define NUGA_AGGLOMERATOR_H

#include "Nuga/include/DynArray.h"
#include "Nuga/include/ngon_t.hxx"
#ifdef DEBUG_AGGLOMERATOR
#include "NGON_debug.h"
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
    inline static void agglomerate_small_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Float vmin, E_Float vratio, ngon_type& ngo, E_Int& nb_aggs, bool force, double angle_threshold=1.e-12);
    ///
    template<typename TriangulatorType>
    inline static void agglomerate_non_star_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, ngon_type& ngo, E_Int& nb_aggs, double angle_threshold = 1.e-12);
  
//    template<typename TriangulatorType>
//    inline static void agglomerate_uncomputable_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, ngon_type& ngo);
    
    ///
    inline static void agglomerate_phs_having_pgs(const K_FLD::FloatArray& crd, ngon_type& ngi, const E_Int* PGlist, E_Int sz, std::vector<E_Int>& pgnids);
    
    template<typename TriangulatorType>
    inline static void collapse_uncomputable_pgs(K_FLD::FloatArray& crd, ngon_type& ngio);
    
    inline static E_Int collapse_pgs(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids);
    inline static E_Int collapse_pgs2(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids);

  public:
    /// agglomerate superfluous polygons (multiply-shared by the same polyhedra. within the flatness tolerance only for area-computable polygons)
    inline static void simplify_phs
      (const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, bool process_externals,
       ngon_type& ngo, ngon_unit& oriento, ngon_unit& phneighborso, const Vector_t<E_Int>* PHlist = nullptr, const Vector_t<E_Int>* skipPGlist = nullptr);

    ///
//    template<typename TriangulatorType>
//    inline static void agglomerate_phs(const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& neighborsi, const Vector_t<E_Int>& PHlist, ngon_type& ngo);

    ///
    template<typename TriangulatorType>
    inline static void agglomerate_phs (const K_FLD::FloatArray& crd, 
                                        const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
                                        ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, bool force, double angle_threshold);
  
  
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
    //MIO::write("agg.plt", crd, cnto, "NGON");
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
    //MIO::write("aggandinit.plt", crd, cnto, "NGON");
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

//  template<typename TriangulatorType>
//  void NUGA::Agglomerator::agglomerate_phs
//  (const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& neighborsi, const Vector_t<E_Int>& PHlist, ngon_type& ngo)
//  {
//    ngo.clear();
//    
//    // Fast return
//    if (PHlist.size() * ngi.PHs.size() * ngi.PGs.size() == 0) 
//    {
//      return;
//    }
//    
//    //std::cout << "ngi : initial nb of phs : " << ngi.PHs.size() << std::endl;
//    
//    std::vector<bool> frozen(ngi.PHs.size(), false); //process one agglo at a time
//    std::vector<E_Int> shared_pgs;
//
//    ngon_unit all_agg_phs, agg_phs;
//
//    E_Int nb_phs = PHlist.size();
//    for (E_Int ii = 0; ii < nb_phs; ++ii)
//    {
//      E_Int i = PHlist[ii];
//
//      if (frozen[i])
//        continue;
//
//      E_Int bestn = IDX_NONE;
//      size_t maxf = 0;
//
//      E_Int nb_neighs = ngi.PHs.stride(i);
//      const E_Int* neighs = neighborsi.get_facets_ptr(i);
//
//      // the best is the one sharing the most number of faces
//      for (E_Int n = 0; (n < nb_neighs); ++n)
//      {
//        E_Int j = *(neighs + n);
//        if (j == IDX_NONE)
//          continue;
//
//        if (frozen[j])
//          continue;
//
//        ngon_type::shared_faces_list(ngi, i, j, shared_pgs);
//        size_t nbf = shared_pgs.size();
//
//        if (nbf < maxf)
//          continue;
//
//        maxf = nbf;
//        bestn = n;
//      }
//
//      if (bestn == IDX_NONE) continue;
//
//      E_Int j = *(neighs + bestn);
//      frozen[i] = frozen[j] = true;
//
//      K_MESH::Polyhedron<UNKNOWN>::merge_two_phs(crd, ngi.PGs, 
//                                                 ngi.PHs.get_facets_ptr(i), ngi.PHs.stride(i),
//                                                 ngi.PHs.get_facets_ptr(j), ngi.PHs.stride(j), 
//                                                 agg_phs);
//
//      all_agg_phs.append(agg_phs);
//
//    }
//
//    ngo.PGs = ngi.PGs;
//    ngo.PHs = all_agg_phs;
//    
//    //std::cout << ngo.PHs.size()/2 << " uncomputable cells have been agglomerated." << std::endl;
//
//    // now add untouched ones
//    for (size_t i = 0; i < frozen.size(); ++i)
//    {
//      if (!frozen[i])
//      {
//        ngo.PHs.add(ngi.PHs.stride(i), ngi.PHs.get_facets_ptr(i));
//        if (!ngi.PHs._type.empty()) ngo.PHs._type.push_back(ngi.PHs._type[i]);
//      }
//    }
//    
//    //std::cout << "effective nb of phs in agglomerated ngo : " << ngo.PHs.size() << std::endl;
//
//    std::vector<E_Int> pgnids, phnids;
//    ngo.remove_unreferenced_pgs(pgnids, phnids);
//
//  }
  
  template<typename TriangulatorType>
  void NUGA::Agglomerator::agglomerate_phs
  (const K_FLD::FloatArray& crd, 
   const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
   ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, bool force, double angle_threshold)
  {
    ngo.clear();
    oriento.clear();
    nb_aggs = 0;

    double concave_threshold = angle_threshold;
    double convex_threshold = angle_threshold;
    
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
      //size_t nbfmax=0;
      //E_Float smax=0.;
      //E_Float worst_reflex_angle=0.;
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
         if (worst_reflex_angle > worst_reflex_a)
            std::cout << "better but generating more concavity !!" << std::endl;
#endif          
          //worst_reflex_angle = worst_reflex_a;
          
#ifdef DEBUG_AGGLOMERATOR
          //if (ii==73)
          /*{
            K_FLD::IntArray cnto;
            ngon_type ng(ngi.PGs, best_agg);
            ng.export_to_array(cnto);
            std::ostringstream o;
            o << "agg_" << n << ".tp";
            MIO::write(o.str().c_str(), crd, cnto, "NGON"); 
            
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

            MIO::write("baryTocentroid.mesh", crd, cn, "BAR"); 
            
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
        MIO::write(o.str().c_str(), crd, cnto, "NGON");
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
  void NUGA::Agglomerator::agglomerate_small_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Float vmin, E_Float vratio, ngon_type& ngo, E_Int& nb_aggs, bool force, double angle_threshold)
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
    NUGA::Agglomerator::agglomerate_phs<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, force, angle_threshold);
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
    NUGA::Agglomerator::agglomerate_phs<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs, false, angle_threshold); //we choose here to agglomerate non star iff it gets better,ie. non-star
  }
  
//  ///
//  template<typename TriangulatorType>
//  void NUGA::Agglomerator::agglomerate_uncomputable_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, ngon_type& ngo)
//  {
//    std::vector<E_Int> PHlist;
//    ngon_type::detect_uncomputable_phs<TriangulatorType>(crd, ngi, PHlist);
//    
//    if (PHlist.empty()) 
//    {
//      std::cout << "OK : There are no uncomputable phs" << std::endl;
//      ngo = ngi;
//      return;
//    }
//    
//    std::cout << PHlist.size() << " uncomputable phs were detected" << std::endl;
//    
//    ngon_unit neighborsi;
//    ngi.build_ph_neighborhood(neighborsi);
//
//    NUGA::Agglomerator::agglomerate_phs<TriangulatorType>(crd, ngi, neighborsi, PHlist, ngo);
//  }
  
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
//        MIO::write("current.plt", crd, cnto, "NGON");
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
    MIO::write("bad_phs.plt", crd, ctmp, "NGON");
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
  
    ngon_type::clean_connectivity(ng, crd);
  
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
  
    ngon_type::clean_connectivity(ng, crd);
  
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
      if (PHn == 31923)
      {
        std::vector<E_Int> tmppg = pgsi;
        K_CONNECT::IdTool::shift(tmppg, -1);
        NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGs(crd, ngi.PGs, tmppg);
        NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH("PHi.tp", crd, ngi, PHi);
        NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH("PHn.tp", crd, ngi, PHn);
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
