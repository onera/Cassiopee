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
#ifndef _POST_STRUCTURED_DATA_VIEW_HPP_
#define _POST_STRUCTURED_DATA_VIEW_HPP_
#include <array>
#include <vector>
#include "Numpy/vector.hpp"
#include "Memory/vector_view.hpp"
#include "../post.h"
#include "zone_data_view.hpp"

namespace K_POST
{
    /**
     * @brief      Vue sur les données contenues par une zone structurée.
     */
    class structured_data_view : public zone_data_view
    {
    public:
        template<typename K>
        using view_type          = zone_data_view::view_type<K>;
        template<typename K>
        using const_view_type    = zone_data_view::const_view_type<K>;
        using dimension_type     = std::array<E_Int,3>;
        using coordinates_type   = zone_data_view::coordinates_type; // Coordonnées par composantes
        using fields_type        = zone_data_view::fields_type;
        using const_fields_type  = zone_data_view::const_fields_type;
        using interpdata_pointer = K_INTERP::InterpData*;
        using impl_pointer_type  = zone_data_view::impl_pointer_type;

        structured_data_view( const dimension_type& dim, const fields_type& fields, const coordinates_npos& pos_coords, 
                              const coordinates_npos& pos_velocity, E_Int cellN);

        structured_data_view( const structured_data_view& ) = default;
        structured_data_view( structured_data_view&& ) = default;

        ~structured_data_view() = default;

        
        /// Getters and setters
        //@{
        /**
         * @brief Retourne la dimension dans chaque direction du maillage structuré
         */
        const dimension_type& dimension() const;        
        
        const_view_type<E_Float> get_cellN() const
        {
            return getField(get_position_of_cellN());
        }

        /**
         * @brief      Retourne le numéro de champs des cellules masquées dans les champs de la zone
         *
         * @return     Le numéro du champs correspondant à cellN
         */
        E_Int get_position_of_cellN() const;
        //@}
    private:
        struct Implementation;
    };
}

#endif
