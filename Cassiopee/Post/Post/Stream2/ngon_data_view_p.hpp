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
#include <utility>
#include "ngon_data_view.hpp"
#include "zone_data_view_p.hpp"

namespace K_POST
{
    struct ngon_data_view::Implementation : public zone_data_view::Implementation
    {
        using zone_data_view::Implementation::fields;

        std::vector<E_Int> m_beg_elt2verts, m_elt2verts;
        std::vector<E_Int> m_beg_vert2elts, m_vert2elts;
        std::vector<E_Int> m_beg_face2verts, m_face2verts;
        std::vector<E_Int> m_beg_elt2faces, m_elt2faces;
        // m_face2elts : Ici, si l'index de l'élément est positif, sens direct, si négatif, sens indirect 
        // (on prendre index commençant à 1 pour les éléments dans ce tableau uniquement ?)
        std::vector<std::pair<E_Int,E_Int>> m_face2elts;

        Implementation( E_Int nb_faces, const const_view_type<E_Int>& face2verts, 
                        E_Int nb_elts,  const const_view_type<E_Int>& elt2faces,
                        const fields_type& fields,
                        const coordinates_npos& pos_coords, const coordinates_npos& pos_velocity);

        bool is_containing(E_Int ind_elt, const point3d& pt) const;

        virtual std::vector<face> get_faces_of_element( E_Int number, E_Int no_zone ) const override;
        virtual std::vector<E_Int> get_indices_of_vertices(E_Int icell) const override;        
        virtual E_Int get_interpolation_cell( const point3d& point ) const override;
        virtual void compute_interpolated_field( const point3d& pt, E_Int ind_cell, 
                                                 E_Int ipos, FldArrayF& interpolatedField ) const override;
        virtual vector3d compute_rotational_in_cell( E_Int ind_cell ) const override;
        virtual E_Float compute_volume_of_cell( E_Int ind_cell ) const override;
    };
}
