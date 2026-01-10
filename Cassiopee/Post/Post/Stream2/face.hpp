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
#ifndef _POST_STREAM_FACE_HPP_
#define _POST_STREAM_FACE_HPP_
#include "../post.h"
#include "Memory/vector_view.hpp"
#include "point3d.hpp"
#include "vector3d.hpp"
#include <limits>
#include <array>

namespace K_POST
{
    //@brief      Définit une face d'un maillage orienté                       
    //@details Pour les méthodes proposées dans cette classe,il est fondamental
    //-        que la face soit orientée avec la normale sortante par rapport à
    //-        l'élément à laquelle elle réfère(par l'indice de la 1ère cellule
    //-        passée dans la paire d'entier ind_cells donnée au constructeur).
    //-        Cette  orientation  est  la  responsabilité  de  l'appelant  qui
    //-        construit la face, dont l'ordre des indices des coordonnées  des
    //-        sommets de la face nous donne l'orientation de la face.         
    class face
    {
        //~                          Partie publique                           
    public:
        //_________________ Définition de types (traits) ______________________
        template<typename K> 
        using const_view_type = K_MEMORY::vector_view<const K>;
        using const_coordinates_type = std::array<const_view_type<E_Float>,3>;
        using intersection_data = std::pair<point3d, E_Int>;

        //________________________ Membres publiques __________________________
        //@brief Indices globaux des sommets constituant la face               
        std::vector<E_Int>      indices_vertices;
        //@brief Indices globaux des elements contenant la face                
        std::pair<E_Int,E_Int> indices_cells;
        //@brief Vue des coordonnees des sommets du maillage contenant la face 
        const_coordinates_type  coordinates_of_zone;

        //_________________ Constructeurs et destructeurs _____________________
        //@brief   Construit une face orienté direct par rapport à un élément  
        //@details L'ordre des indices donné dans ind_coords oriente la face ! 
        //- ind_cells est une paire d'entier : le premier entier donne  l'index
        //- de l'élément contenant la face telle que la face soit sortante.  Le
        //- deuxième entier donne l'index de l'élément opposé qui est  tel  que
        //- la face contenue également par cet élément, est  orientée rentrante
        //- pour l'élément opposé.                                             
        //-                                                                    
        //@param[in]  ind_coords   Indices des sommets constituant la face.    
        //@param[in]  ind_cells    Paire d'entier donnant les éléments  courant
        //-                        et opposés.                                 
        //@param[in]  zone_coords  Coordonnées des points de tout le maillage. 
        face( const std::vector<E_Int>& ind_coords, 
              std::pair<E_Int,E_Int> ind_cells, 
              const const_coordinates_type& zone_coords )
                : indices_vertices(ind_coords), 
                  indices_cells(ind_cells), 
                  coordinates_of_zone(zone_coords)
        {
            barycenter.x = std::numeric_limits<E_Float>::quiet_NaN();
            normal.x     = std::numeric_limits<E_Float>::quiet_NaN();
        }

        //_____________________ Accesseurs et modifieurs ______________________
        //@brief Retourne le nombre de sommets constituant la face             
        E_Int number_of_vertices() const 
        { return this->indices_vertices.size(); }

        //@brief Retourne les coordonnées du ième sommet de la face            
        point3d get_vertex(E_Int ivert) const
        {
            assert(ivert >=0);
            assert(ivert < this->number_of_vertices());
            E_Int ind_vert = this->indices_vertices[ivert];
            return point3d{this->coordinates_of_zone[0][ind_vert],
                           this->coordinates_of_zone[1][ind_vert],
                           this->coordinates_of_zone[2][ind_vert]};
        }

        //@brief      Renvoie le barycentre de la face                         
        //@details Renvoie le barycentre de la face.Au 1er appel, le barycentre
        //-        est calculé, puis il est simplement retourné ( en  tant  que
        //-        proxy).                                                     
        //@return  Un point, barycentre de la face.                            
        const point3d& get_barycenter() const;

        //@brief      Retourne la normale à la face (supposée plane)           
        //@details La direction de la normale dépendra de l'orientation  de  la
        //-        face.                                                       
        //@return  La normale à la face                                        
        const vector3d& get_normal() const;

        //_______________________ Méthodes publiques __________________________
        //@brief Détecte intersection  face par un rayon rentrant ou sortant   
        //@details Cette méthode détecte rapidement si la face  est intersectée
        //-        par un rayon défini par une origine et une direction.  Cette
        //-        méthode ne calcule en aucun cas l'intersection, elle ne fait
        //-        que détecter si l'intersection a lieu ou non.               
        //-                                                                    
        //-        Un rayon orienté dans la même direction que la  normale  est
        //-        dit rentrant ( produit scalaire entre la direction du  rayon
        //-        et la normale est positif), sortant sinon.                  
        //@param[in]  origin     L'origine du rayon                            
        //@param[in]  direction  La direction du rayon                         
        //@return     Première valeur vrai si intersection a lieu,  faux sinon,
        //-           deuxième valeur vraie si rayon rentrant, faux si sortant 
        std::pair<bool,bool> 
        is_intersecting_ray( const point3d& origin, const vector3d& direction, 
                             bool is_with_perturbation=false ) const;

        //@brief      Calcule l'intersection entre la face et un rayon         
        //@details Recherche les coordonnées de l'intersection de la face  avec
        //-        un rayon défini par son origine et sa direction. La  méthode
        //-        retourne les coordonnées de l'intersection et un entier  qui
        //-        donne un numéro de triangle correspondant à un triangle issu
        //-        de la tesselation de la face.Cet entier ne sert pas a priori
        //-        à  un  utilisateur  de  cette  méthode  mais  à  la  méthode
        //-        compute_interpolation qui s'en  servira pour  interpoler  un
        //-        champ en ce point d'intersection.                           
        //-                                                                    
        //-        Si le rayon n'intersecte pas la face, la méthode renvoie une
        //-        exception de type domain_error.                             
        //@param[in]  origin     L'origine du rayon                            
        //@param[in]  direction  La direction du rayon                         
        //@return     Coordonnées de l'intersection de la face avec le rayon et
        //-           l'indice du triangle issu de la  tessellation  de la face
        //-           qui s'intersecte avec le rayon.
        intersection_data 
        compute_intersection( const point3d& origin, 
                              const vector3d& direction ) const;

        //@brief    Calcul le champ interpolé au point d'intersection          
        //@details  Calcul le champ interpolé au point d'intersection,utilisant
        //-         le résultat retourné par compute_intersection              
        //@param[in]  intersection  Données retournées par compute_intersection
        //@param[in]  field         Champ à interpoler                         
        //@return     Le champ interpolé.                                      
        std::vector<E_Float> 
        compute_interpolated_field( const intersection_data& intersection, 
                                    const FldArrayF& field ) const;
        //~                           Partie privée                            
        private:
            using triangle_indices_type = std::array<E_Int,3>;
            mutable std::vector<triangle_indices_type> triangles{};
            void compute_tessellation() const;
            mutable point3d barycenter;
            mutable vector3d normal;

        };

}

#endif
