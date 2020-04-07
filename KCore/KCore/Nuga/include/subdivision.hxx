/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_SUBDIVISION_HXX
#define NUGA_SUBDIVISION_HXX

#include<vector>
#include "Nuga/include/macros.h"


namespace NUGA
{
  enum eSUBDIV_TYPE { ISO = 0, ISO2/*iso metric field : spheres*/, DIR, ANISO/*aniso medtric field : ellipses*/ };
  enum eDIR { NONE = 0, X, Y, XY, /*XZ, YZ*/XYZ };

  
  // Isotropic subdivision data : Vector_t<E_Int>

  // Directional subdivision data
  struct dir_incr_type
  {
    Vector_t<E_Int>      _adap_incr;
    Vector_t<NUGA::eDIR> _ph_dir, _pg_dir;
    
    E_Int& operator[](E_Int i) {return _adap_incr[i]; };
    
    void clear() {_adap_incr.clear();}
    void resize(E_Int sz, E_Int val) {_adap_incr.resize(sz, val);}
    Vector_t<E_Int>::iterator begin() {return _adap_incr.begin();}
    Vector_t<E_Int>::iterator end() {return _adap_incr.end();}
  };

  
  template <typename mesh_t, typename adap_incr_t>
  static void discard_disabledand_unhandled(mesh_t& hmesh, adap_incr_t& adap_incr)
  {
    E_Int nb_phs = hmesh._ng->PHs.size();
    // prevent to adapt on unhandled elements
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      if (adap_incr[i] == 0) continue;

      E_Int nb_faces = hmesh._ng->PHs.stride(i);
      E_Int* faces = hmesh._ng->PHs.get_facets_ptr(i);
      bool admissible_elt = K_MESH::Polyhedron<0>::is_basic(hmesh._ng->PGs, faces, nb_faces);

      if (!admissible_elt || !hmesh._PHtree.is_enabled(i)) // fixme : why the latter can happen ? due to smoothing neighbor traversal ?
        adap_incr[i] = 0;
    }
  }

  //
  template <typename mesh_t>
  static void fix_adap_incr(mesh_t& hmesh, Vector_t<E_Int>& adap_incr)
  {
    discard_disabledand_unhandled(hmesh, adap_incr);
  }

  //
  template <typename mesh_t>
  static void fix_adap_incr(mesh_t& hmesh, dir_incr_type& adap_incr)
  {
    discard_disabledand_unhandled(hmesh, adap_incr);

    //solve inconsistencies
    //alexis : todo

    // premiere version : desactiver l'adaptation dans les cellules XYZ (iso) qui sont connectées à une cellule "layer" par un QUAD lateral
    adap_incr._ph_dir.clear();
    E_Int nb_phs = hmesh._ng->PHs.size();
    adap_incr._ph_dir.resize(nb_phs, XYZ);

    for (E_Int i = 0; i < nb_phs; ++i)
    {
      // if type is layer => adap_incr._ph_dir[i] = XY
    }

    E_Int nb_pgs = hmesh._ng->PGs.size();
    adap_incr._pg_dir.clear();
    adap_incr._pg_dir.resize(nb_pgs, NONE);
    
    // boucle sur les layer

    // appel à get_local

    // remplissage approprié de adap_incr._pg_dir pour les 6 faces : X, Y, ou XY

    // boucle sur le reste : 

    // appel à get_local

    // remplissage de adap_incr._pg_dir ou desactivation de adap_incr._ph_dir

    // si NONE => remplissage
    // sinon, si valeur différente => desactivation
  }

  };
#endif
