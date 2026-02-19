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
#include "unstructured_data_view.hpp"
#include "zone_data_view_p.hpp"
#include "aligned_axis_bounding_box.hpp"
#include "kdtree.hpp"

namespace K_POST
{
    struct unstructured_data_view::Implementation : public zone_data_view::Implementation
    {
        using zone_data_view::Implementation::fields;

        enum element_type 
        {
            tetraedre = 0,
            pyramide      ,
            pentaedre    ,
            hexaedre     ,
            number_of_element_types
        };
        
        static constexpr const std::array<unsigned char,number_of_element_types> number_of_vertices_per_element = {
           {4, 5, 6, 8}
        };

        static constexpr const std::array<unsigned char,number_of_element_types> number_of_vertices_for_polyhedron = {
            {4, 6, 9, 14}
        };

        static constexpr const std::array<unsigned char,number_of_element_types> number_of_faces_per_element = {
            {4, 5, 5, 6}
        };

        static constexpr const std::array<unsigned char,number_of_element_types> total_nb_of_vertices_for_faces_per_element = {
            {12, 16, 18, 24}
        };

        static constexpr const std::array<unsigned char,number_of_element_types> number_of_triangles_per_element = {
            {4, 8, 14, 24}
        };

        static const std::array<std::vector<unsigned char>,number_of_element_types> number_of_vertices_per_face_per_element;

        static const std::array<std::vector<std::vector<unsigned char>>,number_of_element_types> vertices_per_face_per_element;

        static element_type str_type_of_element_to_element_type(const char* eltType)
        {
            if (K_STRING::cmp(eltType, "TETRA") == 0 || K_STRING::cmp(eltType, "TETRA*") == 0)
                return tetraedre;
            else if (K_STRING::cmp(eltType, "PYRA") == 0 || K_STRING::cmp(eltType, "PYRA*") == 0)
                return pyramide;
            else if (K_STRING::cmp(eltType, "PENTA") == 0 || K_STRING::cmp(eltType, "PENTA*") == 0)
                return pentaedre;
            else if (K_STRING::cmp(eltType, "HEXA") == 0 || K_STRING::cmp(eltType, "HEXA*") == 0)
                return hexaedre;
            else
                throw(std::runtime_error("Bad type of element"));
        }

        element_type m_type_of_element;
        const pt_connect_type m_elt2verts;
        std::vector<E_Int> m_beg_face2verts, m_face2verts;
        std::vector<E_Int> m_elt2faces;
        std::vector<std::pair<E_Int,E_Int>> m_face2elts;
        std::vector<E_Int> m_beg_vert2elts, m_vert2elts;
        //const_view_type<E_Int> m_beg_vert2elts, m_vert2elts;

        Implementation( const char* str_elt_type, const pt_connect_type elt2verts, const fields_type& fields,
                        const coordinates_npos& pos_coords, const coordinates_npos& pos_velocity);

        bool is_containing(E_Int ind_elt, const point3d& pt) const;

        virtual std::vector<face> get_faces_of_element( E_Int number, E_Int no_zone ) const override;
        virtual std::vector<E_Int> get_indices_of_vertices(E_Int icell) const override;        
        virtual E_Int get_interpolation_cell( const point3d& point ) const override;
        virtual void compute_interpolated_field( const point3d& pt, E_Int ind_cell, 
                                                 E_Int ipos, FldArrayF& interpolatedField ) const override;
        virtual vector3d compute_rotational_in_cell( E_Int ind_cell ) const override;
        virtual E_Float compute_volume_of_cell( E_Int ind_cell ) const override;
    private:
        void compute_faces_connectivity();
    };
}
