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
#include "zone_data_view.hpp"
#include "aligned_axis_bounding_box.hpp"
#include "kdtree.hpp"

namespace K_POST
{
    struct zone_data_view::Implementation
    {
        kind_of_zone kind;
        coordinates_npos pos_coordinates;
        coordinates_npos pos_velocity;
        fields_type fields;
        //interpdata_pointer interpdata;
        aligned_axis_bounding_box aabbox;
        kdtree tree;

        Implementation( kind_of_zone type_zone, const coordinates_npos& pos_coords, const fields_type& f, const coordinates_npos& vel) :
            kind(type_zone), pos_coordinates(pos_coords), pos_velocity(vel), fields(f)//, interpdata(interp)
        {
            const auto& coords = this->getCoordinates();
            tree = std::move(kdtree(coords, 100u));
            aabbox = aligned_axis_bounding_box(coords[0], coords[1], coords[2]);
        }

        Implementation(const Implementation& impl)   = delete;
        Implementation(const Implementation&& impl)  = delete;

        virtual ~Implementation() {
            //std::cout << "Destruction de l'implementation de base située à l'adresse " << (void*)this << std::endl;
        }

        Implementation& operator = ( const Implementation& ) = delete;
        Implementation& operator = ( Implementation&& )      = delete;

        const_view_type<E_Float> getField(int no) const;
        const_coordinates_type getCoordinates() const
        {
            return {getField(pos_coordinates[0]), getField(pos_coordinates[1]), getField(pos_coordinates[2])};
        }

        virtual std::vector<face> get_faces_of_element( E_Int number, E_Int no_zone ) const = 0;
        virtual std::vector<E_Int> get_indices_of_vertices(E_Int icell) const = 0;
        virtual E_Int get_interpolation_cell( const point3d& point ) const = 0;
        virtual void compute_interpolated_field( const point3d& pt, E_Int ind_cell, 
                                                 E_Int ipos, FldArrayF& interpolatedField ) const = 0;
        virtual vector3d compute_rotational_in_cell( E_Int ind_cell ) const = 0;
        virtual double compute_volume_of_cell( E_Int ind_cell ) const = 0;
    };

}
