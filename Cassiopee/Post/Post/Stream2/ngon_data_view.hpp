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
#ifndef _POST_NGON_DATA_VIEW_HPP_
#define _POST_NGON_DATA_VIEW_HPP_
#include <array>
#include <vector>
#include "Numpy/vector.hpp"
#include "Memory/vector_view.hpp"
#include "../post.h"
#include "zone_data_view.hpp"

namespace K_POST
{
    /**
     * @brief      Vue sur les données contenues par une zone non structurée de type n-gon
     */
    class ngon_data_view : public zone_data_view
    {
    public:
        template<typename K>
        using view_type          = zone_data_view::view_type<K>;
        template<typename K>
        using const_view_type    = zone_data_view::const_view_type<K>;
        using coordinates_type   = zone_data_view::coordinates_type; // Coordonnées par composantes
        using fields_type        = zone_data_view::fields_type;
        using const_fields_type  = zone_data_view::const_fields_type;
        using impl_pointer_type  = zone_data_view::impl_pointer_type;

        ngon_data_view( E_Int nb_faces, const const_view_type<E_Int>& face2verts, 
                        E_Int nb_elts,  const const_view_type<E_Int>& elt2faces, const fields_type& fields,
                        const coordinates_npos& pos_coords, const coordinates_npos& pos_velocity);

        /// Getters and setters
        //@{
        /**
         * @brief Retourne le nombre d'éléments contenus par le maillage de la zone non structurée.
         */
        E_Int number_of_elements() const;

        /**
         * @brief      Retourne la connectivité facette vers sommets du maillage de la zone non structurée
         * 
         * Retourne la connectivité décrivant pour chaque facette les indices des sommets la constituant.
         *
         * @return     la connectivité facette vers sommets
         */
        std::pair<const std::vector<E_Int>&, const std::vector<E_Int>&> get_face_to_vertices() const;
        /**
         * @brief      Retourne la connectivité facette vers élements
         *
         * Retourne la connectivité décrivant pour chaque facette les indices des deux éléments contenant cette face.
         * Si la facette est à la frontière du domaine, le deuxième indice est mis à -1.
         * 
         * @return     La connectivité facette vers éléments
         */
        const std::vector<std::pair<E_Int,E_Int>>& get_face_to_elements() const;

        /**
         * @brief      Retourne la connectivité élément vers facettes
         *
         *  Retourne la connectivité décrivant pour chaque facette les indices des faces la constituant. Le tableau contient nf
         *  valeurs par élément (nf étant le nombre de facettes par élément), rangés par éléments de sorte que l'indice de la première
         *  facette d'un élément e est à la position e*nf dans le tableau retourné.
         *
         * @return     La connectivité élément vers facettes
         */
        std::pair<const std::vector<E_Int>&, const std::vector<E_Int>&> get_element_to_faces() const;

        //@}
    private:
        struct Implementation;
    };
}

#endif
