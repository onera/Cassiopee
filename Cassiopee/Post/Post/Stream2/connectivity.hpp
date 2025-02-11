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
#ifndef _POST_STREAM_CONNECTIVITY_HPP_
#define _POST_STREAM_CONNECTIVITY_HPP_
#include "../post.h"

namespace K_POST
{
    //@brief      Calcule la connectivité sommets vers elements                
    //@details Calcule la  connectivité  sommets vers  éléments à  partir de la
    //-        connectivité  éléments  vers  sommets, sachant  que  pour  cette
    //-        fonction, on suppose que le nombre de sommets par  éléments  est
    //-        unique. Ainsi cette routine ne marche que pour des maillages non
    //-        structurés avec le même type d'éléments.                        
    //-                                                                        
    //-        La fonction retourne deux  tableaux, le  premier  pointant  pour
    //-        chaque sommet, le début de la liste des éléments incidents  dans
    //-        le deuxième tableau, et le deuxième tableau, pour chaque sommet,
    //-        la liste des éléments incidents. Le format des deux tableaux est
    //-        du même type que celui utilisé pour gérer les  matrices  creuses
    //-        en algèbre lineaire (format CSR ou CSC ou encore stockage morse)
    //-                                                                        
    //-        Exemple :                                                       
    //-        --------                                                        
    //-        Soit le maillage suivant :                                      
    //-                0───1───2                                               
    //-                │   │   │                                               
    //-                3───4───5                                               
    //-        Le tableau décrivant les sommets de chaque  élément  du maillage
    //-        (avec quatre sommets par élément ) est :                        
    //-        [ 1, 0, 3, 4, 2, 1, 4, 5] ⇦ deux éléments 0 et 1               
    //-                                                                        
    //-        La connectivité inverse donnant pour chaque sommet les  éléments
    //-        le contenant sera donné par les deux tableaux :                 
    //-                                                                        
    //-        [ 0, 1,    3, 4, 5,  , 7, 8] : dernier élt ⇨ taille 2ᵉ tableau 
    //-          ↓  ↓     ↓  ↓  ↓     ↓  ↓  : début des sommets dans 2ᵉ tableau
    //-        [ 0, 0, 1, 1, 0, 0, 1, 1 ]°  : indice éléments incidents/sommet 
    //-                                                                        
    //@param[in]  number_of_vertices   Nombre de sommets du maillage           
    //@param[in]  element_to_vertices  Connectivité éléments vers sommets      
    //@return     Connectivité sommets vers elements (deux tableaux)           
    std::pair<std::vector<E_Int>,std::vector<E_Int>> 
    compute_vertex_to_elements( E_Int number_of_vertices, 
                                const FldArrayI& element_to_vertices );
}

#endif
