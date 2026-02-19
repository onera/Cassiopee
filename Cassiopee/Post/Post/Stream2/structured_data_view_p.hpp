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
#include "structured_data_view.hpp"
#include "zone_data_view_p.hpp"

namespace K_POST
{
    struct structured_data_view::Implementation : public zone_data_view::Implementation
    {
        using zone_data_view::Implementation::fields;
        dimension_type dimensions;
        E_Int  pos_cellN;

        Implementation(const dimension_type& dim, const coordinates_npos& pos_coords, const fields_type& f,
                       const coordinates_npos& pos_vel, E_Int cN) :
            zone_data_view::Implementation(zone_data_view::STRUCTURED, pos_coords, f, pos_vel),
            dimensions(dim), pos_cellN(cN)
        {
            //std::cout << "Construction de l'implementation pour zone structurée à l'adresse " << (void*)this << std::endl;
        }

        ~Implementation() = default;
        //{
        //    std::cout << "Destruction de l'implementation structurée situé à l'adresse " << (void*)this << std::endl;
        //}
        /**
         * @brief      Test si la cellule d'indice indices contient le point pt :
         *
         * @param[in]  indices   le triplet d'indices de la cellule
         * @param[in]  pt        Le point à tester
         *
         * @return     True if containing, False otherwise.
         */
        bool is_containing( const std::array<E_Int,3>& indices, const point3d& pt) const;

        std::array<E_Int,3> get_indices_of_cell_from_unique_index(E_Int index) const;
        E_Int get_unique_index_from_indices_of_cell(const std::array<E_Int,3>& indices) const;
        E_Int get_unique_index_from_indices_of_vertex(const std::array<E_Int,3>& indices) const;
        std::array<E_Int,3> get_indices_of_vertex_from_unique_index(E_Int index) const;

        virtual std::vector<face> get_faces_of_element( E_Int number, E_Int no_zone ) const override;
        virtual std::vector<E_Int> get_indices_of_vertices(E_Int icell) const override;        
        virtual E_Int get_interpolation_cell( const point3d& point ) const override;
        virtual void compute_interpolated_field( const point3d& pt, E_Int ind_cell, 
                                                 E_Int ipos, FldArrayF& interpolatedField ) const override;
        virtual vector3d compute_rotational_in_cell( E_Int ind_cell ) const override;
        virtual E_Float compute_volume_of_cell( E_Int ind_cell ) const override;
    };
}
