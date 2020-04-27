/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef NUGA_PHQ4_HXX
#define NUGA_PHQ4_HXX

#include "splitting_t.hxx"
#include "MeshElement/Polyhedron.h"

namespace NUGA
{
  template <>
  class splitting_t<K_MESH::Polyhedron<0>, ISO_HEX, 1> : public splitting_base_t
  {

  public:

    ///
    template <typename arr_t>
    splitting_t(K_FLD::FloatArray& crd, const ngon_type& ng, E_Int PHi, E_Int centroidId, 
        const K_FLD::IntArray & F2E, const tree<arr_t> & PGtree)
    {
      //todo JP : calcul et stockage du centroid de PHi dans crd (see HX27 class)
    }

    ///
    template <typename arr_t>
    void split(ngon_type& ng,E_Int PHi,  tree<arr_t>& PHtree, tree<arr_t>& PGtree, K_FLD::IntArray& F2E,
               E_Int firstIntPG, E_Int firstPHchild)
    {
      //todo JP : proposition d'une approche

      //1. dispatcher les quads découpés (decoupage des faces de PHi a l'etape refine_PGs ) par sommet

      // 2. pour chaque paquet de quad :
      
        //a. determiner la ligne polygonale frontière du paquet

        //b. "tisser" les quads internes (gestion des doublons ?) avec le centroid et la ligne

      // 3. MAJ de PHtree: PHtree.set_children(PHi, firstPHchild, nbc);
 
      //MAJ F2E
    }        
   };

   using PHQ4 = splitting_t<K_MESH::Polyhedron<0>, ISO_HEX, 1>;
}

#endif
