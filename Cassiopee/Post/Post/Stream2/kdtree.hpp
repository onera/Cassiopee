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
#ifndef _POST_STREAM_KDTREE_HPP_
#define _POST_STREAM_KDTREE_HPP_
#include <vector>
#include <memory>
#include <array>
#include "point3d.hpp"
#include "Memory/vector_view.hpp"
#include "kcore.h"

namespace K_POST
{
    /**
     * @brief      Cette classe décrit un kd-tree (ici en fait un 3d-tree)
     * 
     * Un arbre k-d (ou k-d tree, pour k-dimensional tree) est une structure de données de partition de 
     * l'espace permettant de stocker des points, et de faire des recherches (recherche par plage, plus proche voisin, etc.) 
     * plus rapidement qu'en parcourant linéairement le tableau de points. 
     * Les arbres k-d sont des cas particuliers d'arbres BSP (binary space partition trees). 
     * 
     * Cette structure a été proposée par Jon Louis Bentley de l'Université Stanford en 1975. (source wikipedia)
     * 
     * Ici, on ne stocke pas les points mais des indices qui correspondent à des points stockés dans des tableaux
     * (un par coordonnée). 
     * 
     */
    class kdtree
    {
    public:
        static constexpr const int dimensions = 3;
        using direction_type = short;
        /// Constructeurs et destructeur
        //@{
        /// Constructeur par défaut
        kdtree() = default;
        /**
         * @brief      Construit un arbre k-d à partir d'un ensemble de sommets
         *
         * Les sommets sont donnés par trois tableaux de coordonnées. Le stride indique le saut à effectuer pour
         * accéder au prochain sommet dans les tableaux (dans le cas où ces tableaux sont intermêlés)
         *
         * @param[in]  xCoords  Les abcisses des sommets
         * @param[in]  yCoords  Les hauteurs des sommets
         * @param[in]  zCoords  Les profondeurs des sommets
         * @param[in]  min_nodes_per_leaf Arête de construire l'arbre k-d tree dès qu'une feuille contient moins que cette valeur
         */
        kdtree( const K_MEMORY::vector_view<const double>& xCoords, const K_MEMORY::vector_view<const double>& yCoords,
                const K_MEMORY::vector_view<const double>& zCoords, unsigned min_nodes_per_leaf = 10 );
        /**
         * @brief      Construit un arbre k-d à partir d'un ensemble de sommets
         *
         * Les sommets sont donnés par un tableau de trois tableaux de coordonnées. Le stride indique le saut à effectuer pour
         * accéder au prochain sommet dans les tableaux (dans le cas où ces tableaux sont intermêlés)
         *
         * @param[in]  coords   Le tableau des trois tableaux de sommets
         * @param[in]  min_nodes_per_leaf Arête de construire l'arbre k-d tree dès qu'une feuille contient moins que cette valeur
         */
        kdtree( const std::array<K_MEMORY::vector_view<const double>,3>& coords, unsigned min_nodes_per_leaf = 10 );
        /// Constructeur de copie (interdit)
        kdtree( const kdtree& ) = delete;
        /// Constructeur de deplacement (par defaut)
        kdtree( kdtree&& )      = default;
        /// Destructeur
        ~kdtree() = default;
        //@}

        /// Opérateurs
        //@{
        /// Opérateur de copie (interdit)
        kdtree& operator = ( const kdtree& ) = delete;
        /// Opérateur de déplacement ( par défaut )
        kdtree& operator = ( kdtree&& )      = default;
        //@}

        /// Méthodes
        //@{
        /**
         * @brief      Retourne l'indice du sommet le plus proche d'un point donné
         *
         * @param[in]  pt  Le point dont il faut trouver le sommet le plus proche
         *
         * @return     L'indice du sommet le plus proche
         */
        std::pair<E_Int,double> nearest(const point3d& pt) const;
        //@}
    private:
        struct node;
        std::shared_ptr<node> make_tree(direction_type, unsigned min_nodes,
                                        K_MEMORY::vector_view<E_Int>& indices);
        std::array<K_MEMORY::vector_view<const double>,3> m_coords;
        std::shared_ptr<node> m_root;
        std::vector<E_Int> m_indices_nodes;
    };
}

#endif
