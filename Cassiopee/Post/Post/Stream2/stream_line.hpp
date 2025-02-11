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
#ifndef _POST_STREAM_LINE_HPP_
#define _POST_STREAM_LINE_HPP_
#include "../post.h"
#include <array>
#include "zone_data_view.hpp"
#include "point3d.hpp"

namespace K_POST
{
    //# ################################################################
    //# #              Définition de la classe streamline              #
    //# ################################################################

    //@brief      This class describes a streamline.                    
    class streamline
    {
    public:
        //_______________ Constructeurs et destructeur _________________

        // @brief     Construit une nouvelle instance                   
        //
        // @param[in]  init_pos          The initialize position        
        // @param[in]  zones             The zones                      
        // @param[in]  is_bidirectional  Indicates if bidirectional     
        // 
        // @details Construit une streamline à partir d'un point initial
        //-         donné dans  l'espace. Si  l'option  is_bidirectional
        //-         est vrai, on construit la streamline dans les  deux 
        //-         sens, en  aval  et  en  amont  du  flux, sinon on ne
        //-         construit que dans le sens amont.                   
        streamline( 
            const point3d& init_pos, 
            const std::vector<zone_data_view>&  zones, 
            E_Int max_vertices_for_streamline, bool is_bidirectional = false);

        //@brief Destructeur                                            
        //@details Détruit l'implémentation du streamline...            
        ~streamline();

        //__________________ Accesseurs et modifieurs __________________

        //@brief      Retourne les segments constituant la streamline   
        //
        //@return     liste des points constituant la streamline        
        std::vector<point3d> segments() const;

        //@brief Retourne le champs correspondant à la streamline       
        // 
        // @details Le champs contient également les coordonnées        
        const FldArrayF& field() const;

        //@brief Retourne le champs correspondant à la streamline       
        // 
        //@details Le champ contient également les coordonnées          
        FldArrayF& field();
        //~                  Partie privée de la classe                 
    private:
        struct Implementation;
        Implementation* ptr_impl;
    };
}

#endif
