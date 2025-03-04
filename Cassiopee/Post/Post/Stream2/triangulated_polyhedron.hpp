/*
    Copyright 2013-2025 Onera.

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
#ifndef _POST_STREAM_TRIANGULATED_POLYHEDRON_HPP_
#define _POST_STREAM_TRIANGULATED_POLYHEDRON_HPP_
#include <array>
#include <vector>
#include "kcore.h"
#include "point3d.hpp"
#include "Memory/vector_view.hpp"
using K_MEMORY::vector_view;

namespace K_POST
{
    /**
     * @brief      Classe décrivant un polyèdre quelconque dont les faces sont triangularisées
     * 
     * Cette classe permet en particuliers, de façon robuste, de détecter si un point est contenu
     * ou non par un polyèdre.
     * 
     * Pour le cas où les faces sont non triangularisées, allez voir la classe polyhedron générale
     * 
     */
    class triangulated_polyhedron
    {
    public:
        using triangle_type = std::array<E_Int,3>;

        /// Constructeurs et destructeur
        //@{
        
        /**
         * @brief      Construit un polyèdre dont toutes les faces sont triangularisées
         *
         * @param      triangular_faces  L'ensemble des indices des sommets des faces triangulaires du polyèdre (indices globale provenant du maillage)
         * @param      coords            Les coordonnées de l'ensemble des sommets du maillage dont est extrait le polygone
         */
        triangulated_polyhedron( const std::vector<triangle_type>& triangular_faces, 
                                 const std::array<vector_view<const double>,3>& coords );
        
        /// Constructeur de copie
        triangulated_polyhedron( const triangulated_polyhedron& ) = default;

        /// Constructeur de déplacement
        triangulated_polyhedron( triangulated_polyhedron&& ) = default;

        /// Destructeur
        ~triangulated_polyhedron() = default;
        //@}

        std::array<point3d,3> triangle_position( const triangle_type& triangle_indices) const
        {
            return {point3d{coords[0][triangle_indices[0]],coords[1][triangle_indices[0]], coords[2][triangle_indices[0]]},
                           {coords[0][triangle_indices[1]],coords[1][triangle_indices[1]], coords[2][triangle_indices[1]]},
                           {coords[0][triangle_indices[2]],coords[1][triangle_indices[2]], coords[2][triangle_indices[2]]} };
        }

        /**
         * @brief      Retourne le volume de ce polynôme dans le cas triangulaire
         *
         * Utilise les volumes signés en additionnant les volumes signées des prismes créés par chaque triangle avec sa projection sur
         * le plan XY
         *
         * @return     { description_of_the_return_value }
         */
        double volume() const;

        int winding_number(const point3d& pt) const;

        /**
         * @brief      Détermine si le point pt est dans le polyèdre
         *
         * @param[in]  pt    Le point à tester
         *
         * Renvoie une exception de type std::invalid_argument si le point pt est sur la surface du polyèdre.
         *
         * @return     Vrai si le point pt est dedans, faux sinon.
         */
        bool is_containing(const point3d& pt) const
        {
            int w = winding_number(pt);
            return (w==1);
        }
    private:
        std::vector<triangle_type> triangular_faces;
        const std::array<vector_view<const double>,3> coords;
    };
}

#endif
