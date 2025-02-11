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
#include "connectivity.hpp"

namespace K_POST
{
    std::pair<std::vector<E_Int>,std::vector<E_Int>> 
    compute_vertex_to_elements( E_Int number_of_vertices, const FldArrayI& element_to_vertices )
    {
        // Attention, les indices de element_to_vertices commencent à un
        E_Int stride = element_to_vertices.getStride();

        std::vector<E_Int> nb_elts_containing_vertex(number_of_vertices, 0);
        for ( E_Int ivert = 0; ivert < element_to_vertices.getNfld(); ++ivert )
        {
            const E_Int* e2v = element_to_vertices.begin(ivert+1);
            for ( E_Int ielt = 0; ielt < element_to_vertices.getSize(); ++ielt )
            {
                nb_elts_containing_vertex[e2v[ielt*stride]-1] += 1;
            }
        }
        std::vector<E_Int> beg_vert_to_elt(number_of_vertices+1);
        beg_vert_to_elt[0] = 0;
        for ( E_Int ivert = 0; ivert < number_of_vertices; ++ivert )
        {
            beg_vert_to_elt[ivert+1] = beg_vert_to_elt[ivert] + nb_elts_containing_vertex[ivert];
        }
        nb_elts_containing_vertex.assign(number_of_vertices, 0);
        std::vector<E_Int> vert_to_elt(beg_vert_to_elt[number_of_vertices]);

        for ( E_Int ivert = 0; ivert < element_to_vertices.getNfld(); ++ivert )
        {
            const E_Int* e2v = element_to_vertices.begin(ivert+1);
            for ( E_Int ielt = 0; ielt < element_to_vertices.getSize(); ++ielt )
            {
                E_Int ind_vert = e2v[ielt*stride]-1;
                vert_to_elt[beg_vert_to_elt[ind_vert]+nb_elts_containing_vertex[ind_vert]] = ielt;
                nb_elts_containing_vertex[ind_vert] += 1;
                assert(nb_elts_containing_vertex[ind_vert] <= beg_vert_to_elt[ind_vert+1] - beg_vert_to_elt[ind_vert]);
            }
        }
        return {beg_vert_to_elt, vert_to_elt};
    }
}
