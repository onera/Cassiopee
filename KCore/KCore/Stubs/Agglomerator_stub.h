#ifndef NUGA_AGGLOMERATOR_H
#define NUGA_AGGLOMERATOR_H

#include "Fld/DynArray.h"
#include "Fld/ngon_t.hxx"
#ifdef DEBUG_AGGLOMERATOR
#include "NGON_debug.h"
#endif


//#include  "Nuga/Boolean/NGON_debug.h"

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
    inline static void agglomerate_small_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Float vmin, E_Float vratio, ngon_type& ngo, E_Int& nb_aggs, bool force){}
    ///
    template<typename TriangulatorType>
    inline static void agglomerate_non_star_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, ngon_type& ngo, E_Int& nb_aggs){}
  
//    template<typename TriangulatorType>
//    inline static void agglomerate_uncomputable_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, ngon_type& ngo);
    
    template<typename TriangulatorType>
    inline static void collapse_uncomputable_pgs(K_FLD::FloatArray& crd, ngon_type& ngio){}
    
    inline static E_Int collapse_pgs(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids){}
    inline static E_Int collapse_pgs2(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids){}

  public:
    /// agglomerate superfluous polygons (multiply-shared by the same polyhedra. within the flatness tolerance only for area-computable polygons)
    inline static void simplify_phs
      (const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, bool process_externals,
       ngon_type& ngo, ngon_unit& oriento, ngon_unit& phneighborso, const Vector_t<E_Int>* PHlist = 0){}

    ///
//    template<typename TriangulatorType>
//    inline static void agglomerate_phs(const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& neighborsi, const Vector_t<E_Int>& PHlist, ngon_type& ngo);

    ///
    template<typename TriangulatorType>
    inline static void agglomerate_phs (const K_FLD::FloatArray& crd, 
                                        const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
                                        ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs, bool force){}
  
  
  private:
      
    inline static void __simplify_phs
      (const K_FLD::FloatArray& crd, const ngon_type& ngi, E_Int PHi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, bool process_externals,
      ngon_unit& gagg_pgs, std::vector<E_Int>& nids, ngon_unit& wagg_pgs, std::map<E_Int, std::vector<E_Int> >& wneigh_to_faces){}

  };

#endif
