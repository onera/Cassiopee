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
#ifndef _POST_ZONE_DATA_VIEW_HPP_
#define _POST_ZONE_DATA_VIEW_HPP_
#include <array>
#include <vector>
#include <limits>
#include <memory>
#include "Numpy/vector.hpp"
#include "Memory/vector_view.hpp"
#include "../post.h"
#include "face.hpp"

namespace K_POST
{
    /**
     * @brief      Vue sur les données contenues par une zone.
     */
    class zone_data_view
    {
    public:
        struct Implementation;
        std::shared_ptr<Implementation> implementation;

    protected:
        zone_data_view( Implementation* impl );
    public:
        zone_data_view(const zone_data_view&) = default;
        zone_data_view(zone_data_view&&) = default;
        ~zone_data_view() = default;
        
        template<typename K>
        using view_type          = K_MEMORY::vector_view<K>;
        template<typename K>
        using const_view_type    = K_MEMORY::vector_view<const K>;
        using coordinates_type   = std::array<view_type<E_Float>,3>; // Coordonnées par composantes
        using const_coordinates_type = std::array<const_view_type<E_Float>,3>;
        using coordinates_npos   = std::array<E_Int,3>; // numéro du champs correspondant à chaque composante des coordonnées dans la zone 
        using fields_type        = FldArrayF*;          // Ce qui contient tous les champs de la zone ( dont les coordonnées )
        using const_fields_type  = const fields_type;   // Idem mais constant
        //using interpdata_pointer = K_INTERP::InterpData*;// Les données pour l'interpolation
        using impl_pointer_type  = std::unique_ptr<Implementation>; // Le  pointeur sur l'implémentation

        enum kind_of_zone
        {
            STRUCTURED=0, ELEMENT_TYPE, NGON, END_KIND
        };

        /// Getters and setters
        //@{

        /**
         * @brief      Retourne le type de zone (structuré, par éléments, ngon)
         *
         * @return     The type de la zone
         */
        kind_of_zone get_type() const;

        /**
         * @brief      Retourne tous les champs de la zone
         *
         * @return     Les champs de la zone
         */
        const_fields_type get_fields() const;

        /**
         * @brief      Retourne un champs stocké sur la zone
         *
         * @param[in]  no    Numéro du champs à retourner
         *
         * @return     Le champs numéro no stocké dans la zone.
         */
        const_view_type<E_Float> getField(int no) const;

        /**
         * @brief      Retourne les coordonnées des sommets contenus dans la zone
         *
         * @return     les coordonnées sous forme d'un tableau à trois entrées ( contenant chacun une composante des coordonnées )
         */
        const_coordinates_type getCoordinates() const
        {
            auto npos = get_position_of_coordinates();
            return { getField(npos[0]), getField(npos[1]), getField(npos[2])};
        }

        /**
         * @brief      Retourne la position des composantes X,Y et Z des coordonnées dans le champs associé à la zone
         *
         * @return     Le numéro du champs pour les composantes X,Y,Z dans l'ensemble des champs contenus dans la zone
         */
        coordinates_npos get_position_of_coordinates() const;

        /**
         * @brief Retourne les indices dans la zone des sommets de la cellule d'index icell.
         * 
         * @param icell Le n° de la cellule
         * @return Retourne dans la convention CGNS les n° des sommets de la cellule d'index icell
         */
        std::vector<E_Int> get_indices_of_vertices(E_Int icell) const;

        /**
         * @brief      Retourne la position des composantes u,v et w du vecteur vitesse du champs associé à la zone
         *
         * @return     Les numéros de poisition des trois composantes de la vitesse du champs
         */
        coordinates_npos get_position_of_velocity() const;

        /**
         * @brief      Retourne le champs de vecteur vitesse de la zone courante
         *
         * @return     Le champs de vecteur vitesse
         */
        const_coordinates_type get_velocity() const
        {
            auto npos = get_position_of_velocity();
            return { getField(npos[0]), getField(npos[1]), getField(npos[2])};
        }

        /**
         * @brief      Retourne les faces appartenant à un élément de la zone repéré par son indice
         *
         * @param[in]  number  L'indice de l'élément dans la zone
         *
         * @return     Une liste des faces de l'éléments avec les données permettant de passer à l'élément opposé
         */
        std::vector<face> get_faces_of_element( E_Int number, E_Int no_zone ) const;

        /**
         * @brief      Retourne l'indice de la cellule de la zone contenant un point donné
         *
         *  Recherche l'élément du maillage contenant le point donné en paramètre. Si un tel élément n'existe pas,
         *  retourne l'indice -1.
         *
         * @param[in]  point  Le point à localiser
         *
         * @return     L'indice de la cellule contenant le point (-1 si une telle cellule n'existe pas)
         */
        E_Int get_interpolation_cell( const point3d& point ) const;
        /**
         * @brief      Calcule le champs interpolé au point pt situé dans la cellule d'indice ind_cell et le 
         * stocke à la ième position dans le champs interpolatedField.
         *
         * Calcule le champs interpolé au point pt situé dans la cellule d'indice ind_cell, indice donné par la
         * méthode @get_interpolation_cell. Le champs interpolé est ensuite stocké à la ipos ème position dans le tableau
         * de champs interpolatedField.
         *
         * @param[in]  pt                 Le point où interpoler le champs
         * @param[in]  ind_cell           L'index de la cellule contenant le point pt
         * @param[in]  ipos               La position dans interpolatedField où stocker le champs interpolé            
         * @param      interpolatedField  Le tableau de champs où stocker le champ interpolé.
         */
        void compute_interpolated_field( const point3d& pt, E_Int ind_cell, 
                                         E_Int ipos, FldArrayF& interpolatedField) const;

        /**
         * @brief Calcule le rotationel d'une cellule donnée de la zone
         * @details Calcule le rotationnel donné d'une cellule selon son type géométrique
         * 
         * @param ind_cell Le numéro de la cellule
         * @return Retourne le vecteur rotationnel de la cellule donnée
         */
        vector3d compute_rotational_in_cell( E_Int ind_cell ) const;

        /**
         * @brief Calcule le volume d'une cellule de la zone
         * @details Calcule le volume d'une cellule de la zone selon son type
         * 
         * @param ind_cell Numéro de la cellule
         * @return Un réel en virgule flottante donnant le volume de la cellule
         */
        double compute_volume_of_cell( E_Int ind_cell ) const;
        //@}
    };
}

#endif
