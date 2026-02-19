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
#ifndef _POST_STREAM_ALIGNED_AXIS_BOUNDING_BOX_HPP_
#define _POST_STREAM_ALIGNED_AXIS_BOUNDING_BOX_HPP_
#include <algorithm>
#include "Memory/vector_view.hpp"
using K_MEMORY::vector_view;
#include "point3d.hpp"

namespace K_POST
{
    //@brief    Description d'une boîte englobante parallèle aux axes du repère
    //@details Par convention, on prendra pour repère le repère suivant        
    //-        ( i : gauche-droite, j : bas-haut, k : derrière-devant )        
    //-                                                                        
    //-              j                                                         
    //-            ­  ⬆                                                         
    //-              O ⮕ i                                                    
    //-            ⬋                                                           
    //-           k                                                            
    //-                                                                        
    //- de sorte que le sommet de la boîte ayant  les  plus petites coordonnées
    //- soit le plus à gauche, le plus derrière et le plus bas.                
    static constexpr const E_Float tolerance_bbox = 1.E-14;
    struct aligned_axis_bounding_box
    {
        point3d p_min; //_ Le point le plus à gauche, derrière et en bas
        point3d p_max; //_ Le point le plus à droite, devant et en haut 

        //_______________ Constructeurs et destructeur _________________

        //@brief  Constructeur par défaut                               
        aligned_axis_bounding_box() = default;
        //@brief Construit la boîte englobante avec  un nuage de  points
        aligned_axis_bounding_box( 
            const vector_view<const E_Float>& xcoords, 
            const vector_view<const E_Float>& ycoords,
            const vector_view<const E_Float>& zcoords )
        {
            auto minmax_x = std::minmax_element(xcoords.begin(), xcoords.end());
            auto minmax_y = std::minmax_element(ycoords.begin(), ycoords.end());
            auto minmax_z = std::minmax_element(zcoords.begin(), zcoords.end());
            p_min = {(*minmax_x.first)-tolerance_bbox, 
                     (*minmax_y.first)-tolerance_bbox, 
                     (*minmax_z.first)-tolerance_bbox};
            p_max = {(*minmax_x.second)+tolerance_bbox, 
                     (*minmax_y.second)+tolerance_bbox, 
                     (*minmax_z.second)+tolerance_bbox};
        }
        //@brief Constructeur de copie                                  
        aligned_axis_bounding_box( const aligned_axis_bounding_box& ) = default;
        //@brief Constructeur de déplacement                            
        aligned_axis_bounding_box( aligned_axis_bounding_box&& ) = default;
        //@brief Destructeur                                            
        ~aligned_axis_bounding_box() = default;

        //________________________ Opérateurs __________________________

        //@brief Opérateur de copie                                     
        aligned_axis_bounding_box& operator = ( const aligned_axis_bounding_box& ) = default;
        //@brief Opérateur de déplacement                               
        aligned_axis_bounding_box& operator = ( aligned_axis_bounding_box&& ) = default;
        //@brief      Opérateur de conversion en chaîne de caractère    
        //@details Transforme la boîte englobante en chaîne de caractère
        //-        facile à lire par un humain                          
        explicit operator std::string() const
        {
            std::ostringstream sout;
            sout << "[ " << std::string(this->p_min) << ", " 
                 << std::string(this->p_max) << " ]";
            return sout.str();
        }

        //_______________________ Autres méthodes_______________________
        //@brief   Teste l'appartenance d'un point à la boîte englobante
        //@param[in]  pt    Le point à tester                           
        //@return     Renvoie  true  si  le  point  est  dans  la  boîte
        //-           englobante, faux sinon                            
        bool contains( const point3d& pt ) const
        {
            return ( (pt.x >= p_min.x) && (pt.x <= p_max.x) &&
                     (pt.y >= p_min.y) && (pt.y <= p_max.y) &&
                     (pt.z >= p_min.z) && (pt.z <= p_max.z) );
        }
        //@brief      Retourne la dimension de la boîte                 
        //@return Retourne la longueur, la hauteur et la  profondeur  de
        //-       la boîte                                              
        std::array<E_Float,3> size() const
        {
            return { p_max.x-p_min.x, p_max.y-p_min.y, p_max.z-p_min.z };
        }
    };
}

#endif
