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
#include <cassert>
#include "zone_data_view_p.hpp"

namespace K_POST
{
    auto zone_data_view::Implementation::getField(int no) const -> const_view_type<E_Float>
    {
        assert(no>=0);
        assert(no<=this->fields->getNfld());
        if (no == 0) return const_view_type<E_Float>(this->fields->begin(1), this->fields->begin(1));
        E_Int size = this->fields->getSize();
        return const_view_type<E_Float>(this->fields->begin(no), this->fields->begin(no) + size);
    }
    // #############################################################################################
    zone_data_view::zone_data_view( Implementation* impl ) :
        implementation(impl)
    {
        assert(implementation != nullptr);
    }
    // ================================================================================
    auto zone_data_view::get_type() const -> kind_of_zone
    {
        assert(implementation != nullptr);
        return implementation->kind;
    }
    // ================================================================================
    auto zone_data_view::get_fields() const -> const_fields_type
    {
        assert(implementation != nullptr);
        return implementation->fields;
    }
    // ================================================================================
    auto zone_data_view::getField(int no) const -> const_view_type<E_Float>
    {
        assert(implementation != nullptr);
        return implementation->getField(no);
    }
    // ================================================================================ 
    auto zone_data_view::get_position_of_coordinates() const -> coordinates_npos 
    {
        assert(implementation != nullptr);
        return this->implementation->pos_coordinates;
    }
    // ================================================================================ 
    std::vector<E_Int> zone_data_view::get_indices_of_vertices(E_Int icell) const
    {
        assert(implementation != nullptr);
        return this->implementation->get_indices_of_vertices(icell);
    } 
    //_ ______________________________________________________________________________
    auto zone_data_view::get_position_of_velocity() const -> coordinates_npos
    {
        assert(implementation != nullptr);
        return this->implementation->pos_velocity;
    }
    //_ ______________________________________________________________________________
    E_Int zone_data_view::get_interpolation_cell(const point3d &point) const
    {
        assert(implementation != nullptr);
        return implementation->get_interpolation_cell(point);
    }
    //_ ______________________________________________________________________________
    void zone_data_view::compute_interpolated_field( const point3d &pt, E_Int ind_cell, 
                                                     E_Int ipos, FldArrayF &interpolatedField ) const
    {
        assert(implementation != nullptr);
        implementation->compute_interpolated_field(pt, ind_cell, ipos, interpolatedField);
    }
    //_ ______________________________________________________________________________
    std::vector<face> zone_data_view::get_faces_of_element( E_Int number, E_Int no_zone ) const
    {
        assert(implementation != nullptr);
        return implementation->get_faces_of_element(number, no_zone);
    }
    //_ _______________________ Calcul du rotationnel dans une cellule _______________
    vector3d zone_data_view::compute_rotational_in_cell(E_Int ind_cell) const
    {
        assert(implementation != nullptr);
        return implementation->compute_rotational_in_cell(ind_cell);
    }
    //_ ___________________ Calcul du volume d'une cellule ___________________________
    double zone_data_view::compute_volume_of_cell(E_Int ind_cell) const
    {
        assert(implementation != nullptr);
        return implementation->compute_volume_of_cell(ind_cell);        
    }
}
