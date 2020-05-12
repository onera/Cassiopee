/*
    Copyright 2013-2020 Onera.

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
#ifndef _POST_STREAM_FACE_HPP_
#define _POST_STREAM_FACE_HPP_
#include "../post.h"
#include "Memory/vector_view.hpp"
#include "point3d.hpp"
#include "vector3d.hpp"

namespace K_POST
{
    /**
     * @brief      Classe définissant une face d'un maillage orienté par rapport à un élément (direct, normale sortante)
     * 
     * Pour les méthodes proposées dans cette classe, il est fondamental que la facette soit orientée avec la normale
     * sortante par rapport à l'élément à laquelle elle réfère (par l'indice de la première cellule passée dans la
     * paire d'entier donnée par ind_cells dans le constructeur). Cette orientation est de la responsabilité du
     * programme qui construit la facette, à l'aide des indices des coordonnées des sommets de la facette, dont l'ordre
     * nous donne l'orientation de la face.
     */
    class face
    {
    public:
        template<typename K> using const_view_type    = K_MEMORY::vector_view<const K>;
        using const_coordinates_type = std::array<const_view_type<E_Float>,3>;
        using intersection_data = std::pair<point3d, E_Int>;

        /// Indices globaux des sommets constituant la face
        std::vector<E_Int>      indices_vertices;
        /// Indices globaux des éléments contenant la face
        std::pair<E_Int,E_Int> indices_cells;
        /// Vue sur les coordonnées des sommets du maillage contenant la face
        const_coordinates_type  coordinates_of_zone;

        /**
         * @brief      Construit une face orienté direct par rapport à un élément
         *
         *  L'ordre des indices donné dans ind_coords oriente la facette ! 
         *  ind_cells est une paire d'entier : le premier entier donne l'index de l'élément contenant la facette
         *  telle que la facette soit sortante. Le deuxième entier donne l'index de l'élément opposé qui est tel
         *  que la facette contenu également par cet élément, est orienté rentrant pour l'élément opposé.
         *
         * @param[in]  ind_coords   Les indices des sommets constituant la facette. 
         * @param[in]  ind_cells    Paire d'entier donnant l'élément courant  et l'élément opposé.
         * @param[in]  zone_coords  Les coordonnées des points de tout le maillage.
         */
        face( const std::vector<E_Int>& ind_coords, std::pair<E_Int,E_Int> ind_cells, const const_coordinates_type& zone_coords ) :
              indices_vertices(ind_coords), indices_cells(ind_cells), coordinates_of_zone(zone_coords)
        {
            barycenter.x = std::numeric_limits<E_Float>::quiet_NaN();
            normal.x     = std::numeric_limits<E_Float>::quiet_NaN();
        }

        /// Retourne le nombre de sommets constituant la face
        E_Int number_of_vertices() const { return this->indices_vertices.size(); }

        /**
         * @brief      Renvoie le barycentre de la face
         *
         * Renvoie le barycentre de la face. Au premier appel, le barycentre est calculé, puis
         * il est simplement retourné (en tant que proxy).
         *
         * @return     Un point, barycentre de la face.
         */
        const point3d& get_barycenter() const;

        /**
         * @brief      Retourne la normale à la face (supposée plane)
         *
         *  La direction de la normale dépendra également de l'orientation de la face.
         *
         * @return     La normale à la face
         */
        const vector3d& get_normal() const;

        /**
         * @brief      Détermine si la face est intersectée par un rayon et si le rayon est rentrant ou sortant
         *             selon l'orientation de la face.
         *
         * Cette méthode détecte rapidement si la face est intersectée par un rayon défini par une origine
         * et une direction. Cette méthode ne calcule en aucun cas l'intersection, elle ne fait que détecter
         * si l'intersection a lieu ou non.
         * 
         * Le rayon est dit rentrant si il est orienté dans la même direction que la normale à la face (le produit scalaire
         * entre la direction du rayon et la normale est positif), sortant sinon.
         * 
         * @param[in]  origin     L'origine du rayon
         * @param[in]  direction  La direction du rayon
         *
         * @return     Première valeur vrai si l'intersection a lieu, faux sinon, deuxième valeur vraie si le rayon est rentrant, faux si il est sortant
         */
        std::pair<bool,bool> is_intersecting_ray( const point3d& origin, const vector3d& direction, bool is_with_perturbation=false ) const;

        /**
         * @brief      Calcule l'intersection entre la face et un rayon
         * 
         * Recherche les coordonnées de l'intersection de la face avec un rayon défini par son origine et sa direction.
         * La méthode retourne les coordonnées de l'intersection et un entier qui correspond à un numéro de triangle
         * correspondant à un triangle issu de la tesselation de la face. Cet entier ne sert pas à priori à un utilisateur
         * de cette méthode mais à la méthode compute_interpolation qui s'en servira pour effectuer l'interpolation d'un
         * champs en ce point d'intersection.
         * 
         * Si le rayon n'intersecte pas la face, la méthode renvoie une exception de type domain_error.
         *
         * @param[in]  origin     L'origine du rayon
         * @param[in]  direction  La direction du rayon
         *
         * @return     Les coordonnées de l'intersection de la face avec le rayon et l'indice du triangle issu de la tessellation de
         *             la face qui s'intersecte avec le rayon.
         */
        intersection_data compute_intersection( const point3d& origin, const vector3d& direction ) const;

        /**
         * @brief      Calcul le champs field interpolé au point d'intersection calculé par compute_intersection
         *
         * @param[in]  intersection  Les données de l'intersectin utiles à l'interpolation
         * @param[in]  field         Le champs à interpoler
         *
         * @return     Le champs interpolé.
         */
        std::vector<E_Float> compute_interpolated_field( const intersection_data& intersection, const FldArrayF& field ) const;

        private:
            using triangle_indices_type = std::array<E_Int,3>;
            mutable std::vector<triangle_indices_type> triangles;
            void compute_tessellation() const;
            mutable point3d barycenter;
            mutable vector3d normal;

        };

}

#endif
