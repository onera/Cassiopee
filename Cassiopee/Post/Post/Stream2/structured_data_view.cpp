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
#include <cmath>
#include "vector3d.hpp"
#include "structured_data_view_p.hpp"
#include "triangulated_polyhedron.hpp"
#include "linear_algebra.hpp"
#include "rotational.hpp"
#include "volume.hpp"
//#define VERBOSE_DEBUG

namespace K_POST
{
    //# ###################################################################################################
    //_ _                              Méthodes pour la classe privée                                     _
    //# ###################################################################################################
    std::array<E_Int,3> 
    structured_data_view::Implementation::get_indices_of_cell_from_unique_index(E_Int index) const
    {
        dimension_type dim = this->dimensions;
        E_Int nb_cells_i  = dim[0]-1;
        E_Int nb_cells_j  = dim[1]-1;
        E_Int nb_cells_ij = nb_cells_i * nb_cells_j;
        E_Int kcell = index/nb_cells_ij;
        E_Int jcell = (index - kcell*nb_cells_ij)/nb_cells_i;
        E_Int icell =  index - kcell*nb_cells_ij - jcell*nb_cells_i;
        return {icell,jcell,kcell};
    }
    // -----------------------------------------------------------------------------------------------------
    E_Int 
    structured_data_view::Implementation::get_unique_index_from_indices_of_cell(const std::array<E_Int,3>& indices) const
    {
        dimension_type dim = this->dimensions;
        E_Int nb_cells_i  = dim[0]-1;
        E_Int nb_cells_j  = dim[1]-1;
        return indices[0] + (indices[1] + indices[2]*nb_cells_j)*nb_cells_i;
    }
    // -----------------------------------------------------------------------------------------------------
    std::array<E_Int,3> 
    structured_data_view::Implementation::get_indices_of_vertex_from_unique_index(E_Int index) const
    {
        dimension_type dim = this->dimensions;
        E_Int nb_verts_i  = dim[0];
        E_Int nb_verts_j  = dim[1];
        E_Int nb_verts_ij = nb_verts_i * nb_verts_j;
        E_Int kvert = index/nb_verts_ij;
        E_Int jvert = (index - kvert*nb_verts_ij)/nb_verts_i;
        E_Int ivert =  index - kvert*nb_verts_ij - jvert*nb_verts_i;
        return {ivert,jvert,kvert};
    }
    // -----------------------------------------------------------------------------------------------------
    E_Int 
    structured_data_view::Implementation::get_unique_index_from_indices_of_vertex(const std::array<E_Int,3>& indices) const
    {
        dimension_type dim = this->dimensions;
        E_Int nb_vertices_i  = dim[0];
        E_Int nb_vertices_j  = dim[1];
        return indices[0] + (indices[1] + indices[2]*nb_vertices_j)*nb_vertices_i;
    }
    // -----------------------------------------------------------------------------------------------------
    std::vector<face> 
    structured_data_view::Implementation::get_faces_of_element( E_Int number, E_Int no_zone ) const
    {
#       if defined(VERBOSE_DEBUG)
        std::cout << "number of element : " << number << std::endl;
#       endif
        dimension_type dim = this->dimensions;
#       if defined(VERBOSE_DEBUG)
        std::cout << "dimension du maillage structuré : " << dim[0] << "," << dim[1] << "," << dim[2] << std::flush << std::endl;
#       endif
        std::array<E_Int,3> nbCells_per_dir{ dim[0]-1, dim[1]-1, dim[2]-1};
        std::array<E_Int,3> indices_of_cell;
        indices_of_cell[2] = number/(nbCells_per_dir[1]*nbCells_per_dir[0]);
        indices_of_cell[1] = (number - indices_of_cell[2]*nbCells_per_dir[1]*nbCells_per_dir[0])/nbCells_per_dir[0];
        indices_of_cell[0] = number - indices_of_cell[2]*nbCells_per_dir[1]*nbCells_per_dir[0] 
                                    - indices_of_cell[1]*nbCells_per_dir[0];

#       if defined(VERBOSE_DEBUG)
        std::cout << "indices of cell : " << indices_of_cell[0] << "; " << indices_of_cell[1] << "; "
                  << indices_of_cell[2] << std::endl;
#       endif
        std::array<E_Int,3> indices_min_of_vertices(indices_of_cell);
        std::vector<face> faces; faces.reserve(6);
        std::vector<E_Int> coords(4);

        E_Int iv = indices_min_of_vertices[0];
        E_Int jv = indices_min_of_vertices[1];
        E_Int kv = indices_min_of_vertices[2];
        // Première face :
        // -------------
        E_Int vert_ind = iv + jv * dim[0] + kv * dim[0]*dim[1]; // sommet (i,j,k)
        coords[0] = vert_ind;                     // Sommet (i,j,k)
        coords[1] = vert_ind + dim[0];            // Sommet (i, j+1, k)
        coords[2] = vert_ind + 1 + dim[0];        // Sommet (i+1,j+1,k)
        coords[3] = vert_ind + 1;                 // Sommet (i+1,j,k)

        E_Int op_elt = (kv == 0 ? -1 : number - nbCells_per_dir[1]*nbCells_per_dir[0]);

        faces.emplace_back(coords, std::pair<E_Int,E_Int>{number, op_elt}, this->getCoordinates());

        // Deuxième face :
        // -------------
        vert_ind = iv + jv * dim[0] + (kv+1) * dim[0]*dim[1]; // sommet (i,j,k+1)
        coords[0] = vert_ind;                     // Sommet (i,j,k+1)
        coords[1] = vert_ind + 1;                 // Sommet (i+1,j,k+1)
        coords[2] = vert_ind + dim[0] + 1;        // Sommet (i+1,j+1,k+1)
        coords[3] = vert_ind + dim[0];            // Sommet (i,j+1,k+1)

        op_elt = (kv == nbCells_per_dir[2]-1 ? -1 : number + nbCells_per_dir[1]*nbCells_per_dir[0]);

        faces.emplace_back(coords, std::pair<E_Int,E_Int>{number, op_elt}, this->getCoordinates());
        
        // Troisième face :
        // --------------
        vert_ind = iv + jv * dim[0] + kv * dim[0]*dim[1]; // sommet (i,j,k)
        coords[0] = vert_ind;                             // Sommet (i,j,k)
        coords[1] = vert_ind + 1;                         // Sommet (i+1,j,k)
        coords[2] = vert_ind + dim[0]*dim[1] + 1;         // Sommet (i+1,j,k+1)
        coords[3] = vert_ind + dim[0]*dim[1];             // Sommet (i,j,k+1)

        op_elt = (jv == 0 ? -1 : number - nbCells_per_dir[0]);

        faces.emplace_back(coords, std::pair<E_Int,E_Int>{number, op_elt}, this->getCoordinates());

        // Quatrième face :
        // --------------
        vert_ind = iv + (jv+1) * dim[0] + kv * dim[0]*dim[1]; 
        coords[0] = vert_ind;                                 // sommet (i,j+1,k)
        coords[1] = vert_ind + dim[0]*dim[1];                 // Sommet (i,j+1,k+1)
        coords[2] = vert_ind + 1 + dim[0]*dim[1];             // sommet (i+1,j+1,k+1)
        coords[3] = vert_ind + 1;                             // sommet (i+1,j+1,k)

        op_elt = (jv == nbCells_per_dir[1]-1 ? -1 : number + nbCells_per_dir[0]);

        faces.emplace_back(coords, std::pair<E_Int,E_Int>{number, op_elt}, this->getCoordinates());

        // Cinquième face :
        // --------------
        vert_ind = iv + jv * dim[0] + kv * dim[0]*dim[1]; 
        coords[0] = vert_ind;                             // sommet (i,j,k)
        coords[1] = vert_ind + dim[0]*dim[1];             // sommet (i,j,k+1)
        coords[2] = vert_ind + dim[0]*dim[1] + dim[0];    // sommet (i,j+1,k+1)
        coords[3] = vert_ind + dim[0];                    // Sommet (i,j+1,k)

        op_elt = (iv == 0 ? -1 : number - 1);

        faces.emplace_back(coords, std::pair<E_Int,E_Int>{number, op_elt}, this->getCoordinates());

        // Sixième face :
        // ------------
        vert_ind = iv + 1 + jv * dim[0] + kv * dim[0]*dim[1]; 
        coords[0] = vert_ind;                                 // sommet (i+1,j,k)
        coords[1] = vert_ind + dim[0];                        // sommet (i+1,j+1,k)
        coords[2] = vert_ind + dim[0] + dim[0]*dim[1];        // sommet (i+1,j+1,k+1)
        coords[3] = vert_ind + dim[0]*dim[1];                 // Sommet (i+1,j,k+1)

        op_elt = (iv == nbCells_per_dir[0]-1 ? -1 : number + 1);

        faces.emplace_back(coords, std::pair<E_Int,E_Int>{number, op_elt}, this->getCoordinates());

        return faces;
    }
    // --------------------------------------------------------------------------------------------
    std::vector<E_Int> 
    structured_data_view::Implementation::get_indices_of_vertices(E_Int icell) const
    {
        std::vector<E_Int> indices;
        indices.reserve(8);
        auto indices_cell = this->get_indices_of_cell_from_unique_index(icell);
        for (const std::array<E_Int,3>& inds : std::vector<std::array<E_Int,3>>{
                            indices_cell,                                                // i  ,j  ,k
                            {indices_cell[0]+1, indices_cell[1]  , indices_cell[2]  },   // i+1,j  ,k
                            {indices_cell[0]+1, indices_cell[1]+1, indices_cell[2]  },   // i+1,j+1,k
                            {indices_cell[0]  , indices_cell[1]+1, indices_cell[2]  },   // i  ,j+1,k
                            {indices_cell[0]  , indices_cell[1]  , indices_cell[2]+1},   // i  ,j  ,k+1
                            {indices_cell[0]+1, indices_cell[1]  , indices_cell[2]+1},   // i+1,j  ,k+1
                            {indices_cell[0]+1, indices_cell[1]+1, indices_cell[2]+1},   // i+1,j+1,k+1
                            {indices_cell[0]  , indices_cell[1]+1, indices_cell[2]+1} }) // i  ,j+1,k+1
        {
            indices.push_back(this->get_unique_index_from_indices_of_vertex(inds));
        }
        return indices;
    }
    // --------------------------------------------------------------------------------------------
    bool 
    structured_data_view::Implementation::is_containing( const std::array<E_Int,3>& indices_cell, const point3d& pt) const
    {
        const auto& crds = this->getCoordinates();
        std::array<std::vector<double>,3> coords; // 8 points cellules + 6 points barycentres
        coords[0].reserve(14); coords[1].reserve(14); coords[2].reserve(14);
        // On extrait tous les points de la cellule :
        for (const std::array<E_Int,3>& indices : std::vector<std::array<E_Int,3>>{
        //for (const auto& indices : {
                            indices_cell,                                                // i  ,j  ,k
                            {indices_cell[0]+1, indices_cell[1]  , indices_cell[2]  },   // i+1,j  ,k
                            {indices_cell[0]+1, indices_cell[1]+1, indices_cell[2]  },   // i+1,j+1,k
                            {indices_cell[0]  , indices_cell[1]+1, indices_cell[2]  },   // i  ,j+1,k
                            {indices_cell[0]  , indices_cell[1]  , indices_cell[2]+1},   // i  ,j  ,k+1
                            {indices_cell[0]+1, indices_cell[1]  , indices_cell[2]+1},   // i+1,j  ,k+1
                            {indices_cell[0]+1, indices_cell[1]+1, indices_cell[2]+1},   // i+1,j+1,k+1
                            {indices_cell[0]  , indices_cell[1]+1, indices_cell[2]+1} }) // i  ,j+1,k+1
        {
            E_Int indv = this->get_unique_index_from_indices_of_vertex(indices);
            coords[0].push_back(crds[0][indv]);
            coords[1].push_back(crds[1][indv]);
            coords[2].push_back(crds[2][indv]);
        }
        // On doit créer 24 triangles pour définir de façon non ambivalente si le point appartient ou non à cette cellule
        // (quatre triangles par face en prenant le barycentre à chaque fois :
        using triangle_type = triangulated_polyhedron::triangle_type;
        std::vector<triangle_type> faces; faces.reserve(24);
        // Première face : (p1,p0,p3,p2) :+ barycentre p8
        //    => calcul barycentre :
        coords[0].push_back(0.25*(coords[0][0]+coords[0][1]+coords[0][2]+coords[0][3]));
        coords[1].push_back(0.25*(coords[1][0]+coords[1][1]+coords[1][2]+coords[1][3]));
        coords[2].push_back(0.25*(coords[2][0]+coords[2][1]+coords[2][2]+coords[2][3]));
        //    => Et les quatre triangles constituant la face :
        faces.emplace_back(triangle_type{1,0,8});
        faces.emplace_back(triangle_type{0,3,8});
        faces.emplace_back(triangle_type{3,2,8});
        faces.emplace_back(triangle_type{1,8,2});
        // Deuxième face : (5,1,2,6) + barycentre p9
        //    => calcul barycentre :
        coords[0].push_back(0.25*(coords[0][5]+coords[0][1]+coords[0][2]+coords[0][6]));
        coords[1].push_back(0.25*(coords[1][5]+coords[1][1]+coords[1][2]+coords[1][6]));
        coords[2].push_back(0.25*(coords[2][5]+coords[2][1]+coords[2][2]+coords[2][6]));
        //    => Et les quatre triangles constituant la face :
        faces.emplace_back(triangle_type{5,1,9});
        faces.emplace_back(triangle_type{1,2,9});
        faces.emplace_back(triangle_type{2,6,9});
        faces.emplace_back(triangle_type{6,5,9});
        // Troisième face : (7,6,2,3) + barycentre p10
        //    => calcul barycentre :
        coords[0].push_back(0.25*(coords[0][7]+coords[0][6]+coords[0][2]+coords[0][3]));
        coords[1].push_back(0.25*(coords[1][7]+coords[1][6]+coords[1][2]+coords[1][3]));
        coords[2].push_back(0.25*(coords[2][7]+coords[2][6]+coords[2][2]+coords[2][3]));
        //    => Et les quatre triangles constituant la face :
        faces.emplace_back(triangle_type{7,6,10});
        faces.emplace_back(triangle_type{6,2,10});
        faces.emplace_back(triangle_type{2,3,10});
        faces.emplace_back(triangle_type{3,7,10});
        // Quatrième face : (0,4,7,3) + barycentre p11
        //    => calcul barycentre :
        coords[0].push_back(0.25*(coords[0][0]+coords[0][4]+coords[0][7]+coords[0][3]));
        coords[1].push_back(0.25*(coords[1][0]+coords[1][4]+coords[1][7]+coords[1][3]));
        coords[2].push_back(0.25*(coords[2][0]+coords[2][4]+coords[2][7]+coords[2][3]));
        //    => Et les quatre triangles constituant la face :
        faces.emplace_back(triangle_type{0,4,11});
        faces.emplace_back(triangle_type{4,7,11});
        faces.emplace_back(triangle_type{7,3,11});
        faces.emplace_back(triangle_type{3,0,11});
        // Cinquième face : (0,1,5,4) + barycentre p12
        //    => calcul barycentre :
        coords[0].push_back(0.25*(coords[0][0]+coords[0][1]+coords[0][5]+coords[0][4]));
        coords[1].push_back(0.25*(coords[1][0]+coords[1][1]+coords[1][5]+coords[1][4]));
        coords[2].push_back(0.25*(coords[2][0]+coords[2][1]+coords[2][5]+coords[2][4]));
        //    => Et les quatre triangles constituant la face :
        faces.emplace_back(triangle_type{0,1,12});
        faces.emplace_back(triangle_type{1,5,12});
        faces.emplace_back(triangle_type{5,4,12});
        faces.emplace_back(triangle_type{4,0,12});
        // Sixième face : (4,5,6,7) + barycentre p13
        coords[0].push_back(0.25*(coords[0][4]+coords[0][5]+coords[0][6]+coords[0][7]));
        coords[1].push_back(0.25*(coords[1][4]+coords[1][5]+coords[1][6]+coords[1][7]));
        coords[2].push_back(0.25*(coords[2][4]+coords[2][5]+coords[2][6]+coords[2][7]));
        //    => Et les quatre triangles constituant la face :
        faces.emplace_back(triangle_type{4,5,13});
        faces.emplace_back(triangle_type{5,6,13});
        faces.emplace_back(triangle_type{6,7,13});
        faces.emplace_back(triangle_type{7,4,13});
#       if defined(VERBOSE_DEBUG)
        std::cout << "=====================================================================" << std::endl;
        std::cout << "pt : " << std::string(pt) << std::endl;
        std::cout << "faces : ";
        for ( const auto& f : faces )
            std::cout << "(" << f[0] << "," << f[1] << "," << f[2] << ")";
        std::cout << std::endl << "coords : ";
        for ( E_Int ip = 0; ip < coords[0].size(); ++ip )
            std::cout << "{" << coords[0][ip] << "," << coords[1][ip] << "," << coords[2][ip] << "}";
        std::cout << std::endl << std::flush;
#       endif
        // Construction du polyèdre :
        triangulated_polyhedron hexaedre( faces, {
                                     K_MEMORY::vector_view<const double>(coords[0].begin(),coords[0].end()),
                                     K_MEMORY::vector_view<const double>(coords[1].begin(),coords[1].end()),
                                     K_MEMORY::vector_view<const double>(coords[2].begin(),coords[2].end())
                                                 } );
        bool is_inside;
        try {
            is_inside = hexaedre.is_containing(pt);
        } catch(std::invalid_argument& err)
        {
            // On est sur la frontiere de l'element :
            //std::cerr << "Warning: streamLine2: interpolated point is on interface. Possibility to have two points in same location in the stream line"
            //          << std::flush << std::endl;
            is_inside = true; // Dans ce cas, on considère qu'on est à l'intérieur (on prend l'élément comme un fermé topologique)
        }
        return is_inside;
    }
    // --------------------------------------------------------------------------------------------
    E_Int 
    structured_data_view::Implementation::get_interpolation_cell( const point3d& point ) const
    {
#       if defined(VERBOSE_DEBUG)
        std::cout << __PRETTY_FUNCTION__ << " : " << std::string(point) << std::endl;    
        std::cout << "Boîte englobante : " << std::string(this->aabbox) << std::endl;    
#       endif
        if (not this->aabbox.contains(point)) return -1;
#       if defined(VERBOSE_DEBUG)        
        std::cout << "ok, bien dans la boite englobante..." << std::endl;
#       endif
        E_Int ind_nearest_vertex;
        double dist_nearest_vertex;
        std::tie(ind_nearest_vertex, dist_nearest_vertex) = this->tree.nearest(point);
        std::array<E_Int,3> indices_nearest_vertex = this->get_indices_of_vertex_from_unique_index(ind_nearest_vertex);
#       if defined(VERBOSE_DEBUG)
        std::cout << "Cellule potentielle : " << indices_nearest_vertex[0] << ", " << indices_nearest_vertex[1]
                  << ", " << indices_nearest_vertex[2] << std::endl;
#       endif
        // On recherche sur les huit cellules potentielles contenant ce sommet :
        dimension_type dim = this->dimensions;
        E_Int ind_cell = -1;
        for (const std::array<E_Int,3>& neighbour_cell_indices : std::vector<std::array<E_Int,3>>{
        //for ( const auto& neighbour_cell_indices : {
                        {indices_nearest_vertex[0]-1,indices_nearest_vertex[1]-1,indices_nearest_vertex[2]-1},
                        {indices_nearest_vertex[0]  ,indices_nearest_vertex[1]-1,indices_nearest_vertex[2]-1},
                        {indices_nearest_vertex[0]  ,indices_nearest_vertex[1]  ,indices_nearest_vertex[2]-1},
                        {indices_nearest_vertex[0]-1,indices_nearest_vertex[1]  ,indices_nearest_vertex[2]-1},
                        {indices_nearest_vertex[0]-1,indices_nearest_vertex[1]-1,indices_nearest_vertex[2]  },
                        {indices_nearest_vertex[0]  ,indices_nearest_vertex[1]-1,indices_nearest_vertex[2]  },
                        indices_nearest_vertex,
                        {indices_nearest_vertex[0]-1,indices_nearest_vertex[1]  ,indices_nearest_vertex[2]  } })
        {
            if ( (neighbour_cell_indices[0]>=0) && (neighbour_cell_indices[1]>=0) && (neighbour_cell_indices[2]>=0) &&
                 (neighbour_cell_indices[0]<dim[0]-1) && (neighbour_cell_indices[1]<dim[1]-1) &&
                 (neighbour_cell_indices[2]<dim[2]-1) )
            {
#       if defined(VERBOSE_DEBUG)                
                std::cout << "Recherche si appartient à cellule : " << neighbour_cell_indices[0] << ", "
                          << neighbour_cell_indices[1] << ", " << neighbour_cell_indices[2] << std::endl;
#       endif
                bool found = this->is_containing(neighbour_cell_indices, point);
                if (found)
                {
                    ind_cell = this->get_unique_index_from_indices_of_cell(neighbour_cell_indices);
                    break;
                }
            }
        }
        return ind_cell;
    }
    // ===================================================================================================
    void structured_data_view::Implementation::compute_interpolated_field( const point3d& pt, E_Int ind_cell, 
                                                                           E_Int ipos, FldArrayF& interpolatedField ) const
    {
        // Soit notre hexaèdre représenté par huit sommets : p₀₀₀, p₁₀₀, p₀₁₀, p₀₀₁, p₁₁₀, p₁₀₁, p₀₁₁, p₁₁₁
        // tels que le sommet pᵢⱼₖ a pour coordonnées {xᵢⱼₖ, yᵢⱼₖ, zᵢⱼₖ}.
        //         →                             →                              →
        // On note e₁₀₀ le vecteur (p₀₀₀, p₁₀₀), e₀₁₀ le vecteur (p₀₀₀,p₀₁₀) et e₀₀₁ le vecteur (p₀₀₀,p₀₀₁)
        // 
        // On va utiliser une interpolation trilinéaire :
        // On prend pour polynôme d'interpolation f({𝛼,𝛽,𝛾}) = a₀ + a₁.𝛼 + a₂.𝛽 + a₃.𝛾 + a₄.𝛼𝛽 + a₅.𝛼𝛾 + a₆.𝛽𝛾 + a₇.𝛼𝛽𝛾
        // avec {𝛼,𝛽,𝛾} ∈ [0;1]³               →        →        →
        // Il est clair que :  p₀₀₀ = p₀₀₀ + 0.e₁₀₀ + 0.e₀₁₀ + 0.e₀₀₁ soit {𝛼,𝛽,𝛾} = {0,0,0}
        //                                     →        →        →
        //                     p₁₀₀ = p₀₀₀ + 1.e₁₀₀ + 0.e₀₁₀ + 0.e₀₀₁ soit {𝛼,𝛽,𝛾} = {1,0,0}
        //                                     →        →        →
        //                     p₀₁₀ = p₀₀₀ + 0.e₁₀₀ + 1.e₀₁₀ + 0.e₀₀₁ soit {𝛼,𝛽,𝛾} = {0,1,0}
        //                                     →        →        →
        //                     p₀₀₁ = p₀₀₀ + 0.e₁₀₀ + 0.e₀₁₀ + 1.e₀₀₁ soit {𝛼,𝛽,𝛾} = {0,0,1}
        //                                                              →     →       →
        // Pour trouver les autres sommets de l'hexahèdre par rapport à e₁₀₀, e₀₁₀ et e₀₀₁, il faut résoudre le système
        // linéaire suivant (pour le sommet pᵢⱼₖ) :
        //   →      →      →     ⎛𝛼⎞
        //  (e₁₀₀ | e₀₁₀ | e₀₀₁) ⎜𝛽⎟ = pᵢⱼₖ - p₀₀₀
        //                       ⎝𝛾⎠                                                                       →     →     →
        // On obtient alors un triplet (𝛼ᵢⱼₖ,𝛽ᵢⱼₖ,𝛾ᵢⱼₖ) permettant de représenter le sommet dans le repère e₁₀₀, e₀₁₀, e₀₀₁
        // Puisqu'on connaît la valeur du champs aux sommets de l'élément, c'est à dire que
        //  f(p₀₀₀), f(p₁₀₀), f(p₀₁₀), f(p₀₀₁), f(p₁₁₀), f(p₁₀₁), f(p₀₁₁) et f(p₁₁₁) sont connus.
        //  D'après le choix du polynôme d'interpolation qu'on a fait plus haut :
        //      f(p₀₀₀) = a₀
        //      f(p₁₀₀) = a₀ + a₁ ⇒ a₁ = f(p₁₀₀) - a₀
        //      f(p₀₁₀) = a₀ + a₂ ⇒ a₂ = f(p₀₁₀) - a₀
        //      f(p₀₀₁) = a₀ + a₃ ⇒ a₃ = f(p₀₀₁) - a₀
        //  et pour les autres sommets, on aura donc :
        //      f(pᵢⱼₖ) = a₀ + a₁.𝛼ᵢⱼₖ + a₂.𝛽ᵢⱼₖ + a₃.𝛾ᵢⱼₖ + a₄.𝛼ᵢⱼₖ𝛽ᵢⱼₖ + a₅.𝛼ᵢⱼₖ𝛾ᵢⱼₖ + a₆.𝛽ᵢⱼₖ𝛾ᵢⱼₖ + a₇.𝛼ᵢⱼₖ𝛽ᵢⱼₖ𝛾ᵢⱼₖ
        //  ce qui nous donne un système linéaire de dimension quatre à résoudre :
        //      a₄.𝛼ᵢⱼₖ𝛽ᵢⱼₖ + a₅.𝛼ᵢⱼₖ𝛾ᵢⱼₖ + a₆.𝛽ᵢⱼₖ𝛾ᵢⱼₖ + a₇.𝛼ᵢⱼₖ𝛽ᵢⱼₖ𝛾ᵢⱼₖ = f(pᵢⱼₖ) - a₀ - a₁.𝛼ᵢⱼₖ - a₂.𝛽ᵢⱼₖ - a₃.𝛾ᵢⱼₖ
        //   avec i+j+k=2 ou 3 (i,j et k valant 0 ou 1).
        //   Il faudra ensuite déterminer les coordonnées (𝛼,𝛽,𝛾) de notre point à interpoler, puis appliquer la fonction
        //   polynomiale ainsi calculée.
        //   
        const auto& crds = this->getCoordinates();
        auto indices_cell = this->get_indices_of_cell_from_unique_index(ind_cell);
#       if defined(VERBOSE_DEBUG)
        std::cout << "indices_cell : " << indices_cell[0] << " " << indices_cell[1] << " " << indices_cell[2] << std::endl;
#       endif
        std::array<point3d,8> coords; // 8 points/cellules
        // On extrait tous les points de la cellule ainsi que les champs à interpoler :
        FldArrayF* fld = this->fields;
        E_Int nfld = fld->getNfld();
        std::vector<std::array<double,8>> values;
        values.reserve(nfld);

        E_Int ivert = 0;
        for (const std::array<E_Int,3>& indices : std::vector<std::array<E_Int,3>>{
        //for (const auto& indices : {
                            indices_cell,                                                // i  ,j  ,k
                            {indices_cell[0]+1, indices_cell[1]  , indices_cell[2]  },   // i+1,j  ,k
                            {indices_cell[0]+1, indices_cell[1]+1, indices_cell[2]  },   // i+1,j+1,k
                            {indices_cell[0]  , indices_cell[1]+1, indices_cell[2]  },   // i  ,j+1,k
                            {indices_cell[0]  , indices_cell[1]  , indices_cell[2]+1},   // i  ,j  ,k+1
                            {indices_cell[0]+1, indices_cell[1]  , indices_cell[2]+1},   // i+1,j  ,k+1
                            {indices_cell[0]+1, indices_cell[1]+1, indices_cell[2]+1},   // i+1,j+1,k+1
                            {indices_cell[0]  , indices_cell[1]+1, indices_cell[2]+1} }) // i  ,j+1,k+1
        {
            E_Int indv = this->get_unique_index_from_indices_of_vertex(indices);
#           if defined(VERBOSE_DEBUG)
            std::cout << "indice sommet interpolation " << indv << std::endl;
#           endif
            coords[ivert] = {crds[0][indv],crds[1][indv],crds[2][indv]};
            for ( E_Int ifld = 1; ifld <= nfld; ++ifld )
            {
                values[ifld-1][ivert] = (*fld)(indv,ifld);
#               if defined(VERBOSE_DEBUG)
                std::cout << "field(" << indv << "," << ifld << ") = " << values[ifld-1][ivert] << std::endl;
#               endif
            }
            ivert += 1;
        }
        vector3d e100(coords[0], coords[1]);
        vector3d e010(coords[0], coords[3]);
        vector3d e001(coords[0], coords[4]);
#       if defined(VERBOSE_DEBUG)
        std::cout << "e100 = " << std::string(e100) << ", e010 = " << std::string(e010) << ", e001 = "
                  << std::string(e001) << std::endl;
#       endif
        matrix_3x3_type A{std::array<double,3>{ e100.x ,e010.x, e001.x},
                                              { e100.y ,e010.y, e001.y},
                                              { e100.z ,e010.z, e001.z}
                         };
        matrix_3x3_type invA = inverse(A);
        vector3d bar110 = invA*vector3d(coords[0],coords[2]);
        vector3d bar101 = invA*vector3d(coords[0],coords[5]);
        vector3d bar111 = invA*vector3d(coords[0],coords[6]);
        vector3d bar011 = invA*vector3d(coords[0],coords[7]);

        vector3d bary_pt= invA*vector3d(coords[0], pt);
#       if defined(VERBOSE_DEBUG)
        std::cout << "bar110 : " << std::string(bar110) << ", bar101 : " << std::string(bar101)
                  << "bar111 : " << std::string(bar111) << ", bar011 : " << std::string(bar011)
                  << ", bary pt : " << std::string(bary_pt) << std::endl;
#       endif
        //      a₄.𝛼ᵢⱼₖ𝛽ᵢⱼₖ + a₅.𝛼ᵢⱼₖ𝛾ᵢⱼₖ + a₆.𝛽ᵢⱼₖ𝛾ᵢⱼₖ + a₇.𝛼ᵢⱼₖ𝛽ᵢⱼₖ𝛾ᵢⱼₖ = f(pᵢⱼₖ) - a₀ - a₁.𝛼ᵢⱼₖ - a₂.𝛽ᵢⱼₖ - a₃.𝛾ᵢⱼₖ
        matrix_4x4_type B{ std::array<double,4>{bar110.x*bar110.y,bar110.x*bar110.z,bar110.y*bar110.z,bar110.x*bar110.y*bar110.z},
                                               {bar101.x*bar101.y,bar101.x*bar101.z,bar101.y*bar101.z,bar101.x*bar101.y*bar101.z},
                                               {bar111.x*bar111.y,bar111.x*bar111.z,bar111.y*bar111.z,bar111.x*bar111.y*bar111.z},
                                               {bar011.x*bar011.y,bar011.x*bar011.z,bar011.y*bar011.z,bar011.x*bar011.y*bar011.z}
                         };
#       if defined(VERBOSE_DEBUG)                         
        std::cout << "B : " << std::endl << to_string(B) << std::endl;
#       endif
        auto LUB = factorize(B);
#       if defined(VERBOSE_DEBUG)
        std::cout << "champs interpolé ";
#       endif
        for ( E_Int ifld = 0; ifld < nfld; ++ifld)
        {
#           if defined(VERBOSE_DEBUG)
            std::cout << "Pour le champs n°" << ifld << std::endl;
#           endif
            double a0 = values[ifld][0];        // a₀ = f(p₀₀₀)
            double a1 = values[ifld][1] - a0;   // a₁ = f(p₁₀₀) - a₀
            double a2 = values[ifld][3] - a0;   // a₂ = f(p₀₁₀) - a₀
            double a3 = values[ifld][4] - a0;   // a₃ = f(p₀₀₁) - a₀
            vector4d b{
                values[ifld][2] - a0 -a1*bar110.x -a2*bar110.y - a3*bar110.z, 
                values[ifld][5] - a0 -a1*bar101.x -a2*bar101.y - a3*bar101.z, 
                values[ifld][6] - a0 -a1*bar111.x -a2*bar111.y - a3*bar111.z, 
                values[ifld][7] - a0 -a1*bar011.x -a2*bar011.y - a3*bar011.z
                      };
            auto x = inverse_linear_system(LUB, b);
            double a4 = x[0];
            double a5 = x[1];
            double a6 = x[2];
            double a7 = x[3];
#           if defined(VERBOSE_DEBUG)
            std::cout << "a0 : " << a0 << ", a1 : " << a1 << ", a2 : " << a2 << ", a3 : " << a3
                      << ", a4 : " << a4 << ", a5 : " << a5 << ", a6 : " << a6 << ", a7 : " << a7 << std::endl;
#           endif
            interpolatedField(ipos,ifld+1) = a0 + a1*bary_pt.x + a2*bary_pt.y + a3*bary_pt.z +
                                             a4*bary_pt.x*bary_pt.y + a5*bary_pt.x*bary_pt.z + a6*bary_pt.y*bary_pt.z +
                                             a7*bary_pt.x*bary_pt.y*bary_pt.z;
#           if defined(VERBOSE_DEBUG)
            std::cout << interpolatedField(ipos,ifld+1) << " " << std::endl;
#           endif
        }
    }
    // _ _________________________ Calcul du rotationel dans une cellule d'un maillage structuré __________
    vector3d structured_data_view::Implementation::compute_rotational_in_cell( E_Int ind_cell ) const
    {
        const auto& crds = this->getCoordinates();
        auto indices_cell = this->get_indices_of_cell_from_unique_index(ind_cell);
        std::array<point3d,8> coords; // 8 points/cellules
        // On extrait tous les points de la cellule ainsi que le champs vitesse :
        FldArrayF* fld = this->fields;
        std::array<vector3d,8> vel;
        //coordinates_npos pos_vel = this->pos_velocity;

        E_Int ivert = 0;
        for (const std::array<E_Int,3>& indices : std::vector<std::array<E_Int,3>>{
                            indices_cell,                                                // i  ,j  ,k
                            {indices_cell[0]+1, indices_cell[1]  , indices_cell[2]  },   // i+1,j  ,k
                            {indices_cell[0]+1, indices_cell[1]+1, indices_cell[2]  },   // i+1,j+1,k
                            {indices_cell[0]  , indices_cell[1]+1, indices_cell[2]  },   // i  ,j+1,k
                            {indices_cell[0]  , indices_cell[1]  , indices_cell[2]+1},   // i  ,j  ,k+1
                            {indices_cell[0]+1, indices_cell[1]  , indices_cell[2]+1},   // i+1,j  ,k+1
                            {indices_cell[0]+1, indices_cell[1]+1, indices_cell[2]+1},   // i+1,j+1,k+1
                            {indices_cell[0]  , indices_cell[1]+1, indices_cell[2]+1} }) // i  ,j+1,k+1
        {
            E_Int indv = this->get_unique_index_from_indices_of_vertex(indices);
#           if defined(VERBOSE_DEBUG)
            std::cout << "indice sommet interpolation " << indv << std::endl;
#           endif
            coords[ivert] = {crds[0][indv],crds[1][indv],crds[2][indv]};
            vel[ivert] = { (*fld)(indv, this->pos_velocity[0]),
                           (*fld)(indv, this->pos_velocity[1]),
                           (*fld)(indv, this->pos_velocity[2]) };
            //std::cout << "coords(rot," << ivert << ") : " << std::string(coords[ivert]) << ", vel(" << ivert << ") = " 
            //          << std::string(vel[ivert]) << " ";
            ivert += 1;
        }
        return compute_rotational_on_hexaedra(coords[0], coords[1], coords[2], coords[3], 
                                              coords[4], coords[5], coords[6], coords[7], 
                                              vel   [0], vel   [1], vel   [2], vel   [3], 
                                              vel   [4], vel   [5], vel   [6], vel   [7]);

    }
    //_ ___________________ Calcul du volume d'une cellule du maillage structuré __________________________
    double structured_data_view::Implementation::compute_volume_of_cell(E_Int ind_cell) const
    {
        const auto& crds = this->getCoordinates();
        auto indices_cell = this->get_indices_of_cell_from_unique_index(ind_cell);
        std::array<point3d,8> coords; // 8 points/cellules
        E_Int ivert = 0;
        for (const std::array<E_Int,3>& indices : std::vector<std::array<E_Int,3>>{
                            indices_cell,                                                // i  ,j  ,k
                            {indices_cell[0]+1, indices_cell[1]  , indices_cell[2]  },   // i+1,j  ,k
                            {indices_cell[0]+1, indices_cell[1]+1, indices_cell[2]  },   // i+1,j+1,k
                            {indices_cell[0]  , indices_cell[1]+1, indices_cell[2]  },   // i  ,j+1,k
                            {indices_cell[0]  , indices_cell[1]  , indices_cell[2]+1},   // i  ,j  ,k+1
                            {indices_cell[0]+1, indices_cell[1]  , indices_cell[2]+1},   // i+1,j  ,k+1
                            {indices_cell[0]+1, indices_cell[1]+1, indices_cell[2]+1},   // i+1,j+1,k+1
                            {indices_cell[0]  , indices_cell[1]+1, indices_cell[2]+1} }) // i  ,j+1,k+1
        {
            E_Int indv = this->get_unique_index_from_indices_of_vertex(indices);
            coords[ivert] = {crds[0][indv],crds[1][indv],crds[2][indv]};
            ivert += 1;
        }
        return compute_volume_hexaedra(coords[0], coords[1], coords[2], coords[3], 
                                       coords[4], coords[5], coords[6], coords[7] );
    }
    //# ###################################################################################################
    //_ _                              Méthodes pour la classe publique                                   _
    //# ###################################################################################################
    structured_data_view::structured_data_view( const dimension_type& dim, const fields_type& fields, const coordinates_npos& pos_coords, 
                                                const coordinates_npos& pos_velocity, E_Int cellN) :
        zone_data_view(new Implementation(dim, pos_coords, fields, pos_velocity, cellN))
    {}
    //===================================================================================================
    auto structured_data_view::dimension() const -> const dimension_type&
    {
        assert(implementation != nullptr);
        assert(implementation->kind == zone_data_view::STRUCTURED);
        return static_cast<Implementation&>(*implementation).dimensions;
    }
    // ===================================================================================================
    E_Int structured_data_view::get_position_of_cellN() const
    {
        assert(implementation != nullptr);
        assert(implementation->kind == zone_data_view::STRUCTURED);
        return static_cast<Implementation&>(*implementation).pos_cellN;
    }
}
