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
#ifndef _POST_STREAM_POINT3D_HPP_
#define _POST_STREAM_POINT3D_HPP_
#include <cassert>
#include <cmath> 
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace K_POST
{
    //# ################################################################
    //# #   Définition d'un point en 3 dimensions via une structure.   #
    //# ################################################################

    //@brief Coordonnées d'un point en trois dimensions                 
    struct point3d
    {
        double x,y,z;//_ Coordonnées du point __________________________

        //______________ Constructeurs et destructeur __________________

        //@brief Constructeur par défaut.                               
        point3d() = default;
        //@brief      Constructeur initialisation un point              
        //
        //@param[in]  _x    La valeur de l'abscisse                     
        //@param[in]  _y    La valeur de l'ordonnée                     
        //@param[in]  _z    La valeur de la hauteur                     
        point3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
        {}
        //@brief Constructeur par copie                                 
        point3d(const point3d& ) = default;
        //@brief Constructeur de déplacement                            
        point3d(point3d&&) = default;
        //@brief Destructeur                                            
        ~point3d() = default;

        //_________________ Opérateurs sur les points __________________

        //@brief Opérateur de copie                                     
        point3d& operator = ( const point3d& p ) = default;
        //@brief Opérateur de déplacement                               
        point3d& operator = ( point3d&& p )      = default;
        //@brief Opérateur d'accès à la ième coordonnée en lecture seule
        double operator[] (unsigned i) const
        {
            assert(i<3);
            return (i==0 ? this->x : i == 1 ? this->y : this->z);
        }
        //@brief      Conversion d'un point en  chaîne  de caractère    
        //@details    Cet opérateur transforme le point  en une  chaîne 
        //-           de caractère agréable à lire pour un être humain. 
        explicit operator std::string() const
        {
            std::ostringstream sout;
            sout << std::setprecision(14) << "{" << this->x << ", " << this->y << ", " << this->z
                 << "}";
            return sout.str();
        }
    };

    //@brief      Calcul la distance au carré entre deux points         
    //@param[in]  p1    Le premier point                                
    //@param[in]  p2    Le deuxième point                               
    //@return     Renvoie la distance au carré entre p1 et p2           
    inline double square_distance( const point3d& p1, const point3d& p2)
    {
        double dx = p2.x-p1.x, dy = p2.y-p1.y, dz = p2.z - p1.z;
        return dx*dx + dy*dy + dz*dz;
    }

    //@brief      Calcul la distance entre deux points                  
    //@param[in]  p1    Le premier point                                
    //@param[in]  p2    Le deuxième point                               
    //@return     Renvoie la distance entre p1 et p2                    
    inline double distance( const point3d& p1, const point3d& p2 )
    { return std::sqrt(square_distance(p1, p2)); }
}

//@brief      Opérateur de sortie dans un flux pour le point            
//@details Cet  opérateur  écrit  un  point dans  un flux de sortie, de 
//-        sorte que le point écrit soit facilement  relisable  par  le 
//-         flux d'entrée, mais pas forcément pour un être humain.      
//-                                                                     
//@param      out   Le flux de sortie                                   
//@param[in]  p     Le point à écrire dans le flux de sortie            
//@return     Le nouvel état du flux de sortie                          
inline std::ostream& 
operator << ( std::ostream& out, const K_POST::point3d& p)
{
    out << p.x << " " << p.y << " " << p.z;
    return out;
}

//@brief      Opérateur d'entrée d'un flux pour lire un point           
//@details Cet opérateur va lire un point à partir  d'un flux  d'entrée 
//-         et écrire les coordonnées de ce point dans le point donné  à
//-         droite de l'opérateur >>.                                   
//-                                                                     
//@param      inp   Le flux d'entrée                                    
//@param      p     Le point recevant les données du flux d'entrée      
//@return     Le nouvel état du flux d'entrée                           
inline std::istream& 
operator >> ( std::istream& inp, K_POST::point3d& p)
{
    inp >> p.x >> p.y >> p.z;
    return inp;
}

#endif
