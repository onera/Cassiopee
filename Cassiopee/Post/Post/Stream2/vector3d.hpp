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
#ifndef _POST_STREAM_VECTOR3D_HPP_
#define _POST_STREAM_VECTOR3D_HPP_
#include <cassert>
#include <cmath>
#include "point3d.hpp"

namespace K_POST
{
    /**
     * @brief      Structure décrivant un vecteur géométrique en trois dimensions
     */
    struct vector3d
    {
        double x,y,z; /// Coordonnées du vecteur

        /// Constructeurs et destructeur
        //@{
        /// Constructeur par défaut
        vector3d() = default;
        /**
         * @brief      Constructeur initialisant un vecteur à partir de ses valeurs
         *
         * @param[in]  _x    La valeur de l'abscisse
         * @param[in]  _y    La valeur de l'ordonnée
         * @param[in]  _z    La valeur de la hauteur
         */
        vector3d( double _x, double _y, double _z ) : x(_x), y(_y), z(_z)
        {}
        /**
         * @brief      Constructeur initialisant un vecteur à partir de deux points
         *
         * @param[in]  p1    Le point origine du vecteur
         * @param[in]  p2    Le point extrémité du vecteur
         */
        vector3d( const point3d& p1, const point3d& p2) : x(p2.x-p1.x), y(p2.y-p1.y), z(p2.z-p1.z)
        {}
        /// Constructeur par copie
        vector3d( const vector3d& ) = default;
        /// Constructeur par déplacement
        vector3d( vector3d&& ) = default;
        /// Destructeur
        ~vector3d() = default;
        //@}

        /// Opérateurs
        //@{
        /// Opérateur de copie
        vector3d& operator = ( const vector3d& ) = default;
        /// Opérateur de déplacement
        vector3d& operator = ( vector3d&& ) = default;
        /// Addition inplace
        vector3d& operator += (const vector3d& u)
        {
            this->x += u.x; this->y += u.y; this->z += u.z;
            return *this;
        }
        /// Addition non inplace
        vector3d operator + ( const vector3d& u ) const
        {
            return { this->x + u.x, this->y + u.y, this->z + u.z };
        }
        /// Soustraction inplace
        vector3d& operator -= (const vector3d& u)
        {
            this->x -= u.x; this->y -= u.y; this->z -= u.z;
            return *this;
        }
        /// Soustraction non inplace
        vector3d operator - ( const vector3d& u ) const
        {
            return { this->x - u.x, this->y - u.y, this->z - u.z };
        }
        /// Inversion du vecteur
        vector3d operator -() const
        {
            return { -this->x, -this->y, -this->z };
        }
        /// Produit vectoriel
        vector3d operator ^( const vector3d& u ) const
        {
            return { this->y * u.z - this->z * u.y,
                     this->z * u.x - this->x * u.z,
                     this->x * u.y - this->y * u.x };
        }
        /// Produit scalaire
        double operator | ( const vector3d& u ) const
        {
            return this->x*u.x + this->y*u.y + this->z*u.z;
        }
        /// Accesseur par indice aux coefficients du vecteur
        const double& operator [] ( unsigned i ) const
        {
            assert(i<3);
            return (i==0 ? this->x : i == 1 ? this->y : this->z );
        }
        /// Accesseur par indice aux coefficients du vecteur
        double& operator [] ( unsigned i )
        {
            assert(i<3);
            return (i==0 ? this->x : i == 1 ? this->y : this->z );
        }
        /**
         * @brief      Opérateur de conversion d'un vecteur en chaîne de caractère
         * 
         * Cet opérateur transforme le vecteur en une chaîne de caractère agréable à lire pour
         * un être humain.
         */
        explicit operator std::string() const
        {
            std::ostringstream sout;
            sout << "(" << this->x << ", " << this->y << ", " << this->z << ")";
            return sout.str();
        }
        //@}
    };
    /**
     * @brief      Opérateur d'homothétie
     *
     * @param[in]  alpha  Le coefficient d'homothétie
     * @param[in]  u      Le vecteur sur lequel on applique l'homothétie
     *
     * @return     Le résultat de l'homothétie
     */
    inline vector3d operator * ( double alpha, const vector3d& u )
    {
        return {alpha*u.x, alpha*u.y, alpha*u.z};
    }
    /**
     * @brief      Translation d'un point par un vecteur
     *
     * @param[in]  p     Le point à translater
     * @param[in]  t     Le vecteur de translation
     *
     * @return     Le résultat de la translation
     */
    inline point3d operator + ( const point3d& p, const vector3d& t )
    {
        return {p.x + t.x, p.y + t.y, p.z + t.z};
    }
    /**
    * @brief      Calcul la norme L2 au carré du vecteur u
    *
    *  Calcul la norme L2 au carré du vecteur. Pour avoir la norme L2, utiliser la méthode abs
    *  qui sera plus chère à calculer (à cause de la racine carrée).
    *  Remarque : l'interface choisie ici suit la règle des moindres surprises. En C++, abs sert pour calculer
    *             la norme tandis que norm sert pour calculer la norme au carré...
    *
    * @return     Retourne ||u||_{L2}^{2}
    */
    inline double norm( const vector3d& u )
    {
        return (u|u);
    }
    /**
     * @brief      Retourne la norme L2 (cartésienne) du vecteur u
     *
     * @param[in]  u     Le vecteur u
     *
     * @return     ||u||_{L2}
     */
    inline double abs( const vector3d& u )
    {
        return std::sqrt(norm(u));
    }

    /**
     * @brief Effectue une rotation selon un axe et un angle theta d'un vecteur u
     * @details Effectue une rotation selon un axe a (vecteur normalisé) et un angle theta d'un vecteur u
     * 
     * On applique la formule suivante pour effectuer la rotation :
     * 
     # ⎛vₓ⎞     ⎛uₓ⎞          ⎧ cos²(θ/2) + sin²(θ/2)(2aᵢ²-1) si i = j 
     # ⎜vᵧ⎟ = R.⎜uᵧ⎟ où Rᵢⱼ = ⎨                                        
     # ⎝vz⎠     ⎝uz⎠          ⎩ 2aᵢaⱼsin²(θ/2) - εᵢⱼₖaₖ sin(θ) si i ≠ j 
     * 
     * @param axis Le vecteur a donnant l'axe de la rotation
     * @param theta L'angle de la rotation
     * @param u Le vecteur sur lequel on applique la rotation
     * 
     * @return Le vecteur résultat après rotation
     */
    inline vector3d rotate( const vector3d& axis, double theta, const vector3d& u)
    {
        assert(std::abs(abs(axis) - 1.) < 1.E-14);
        double ct = std::cos(theta);
        double st = std::sin(theta);
        double umct = 1. - ct;
        return {
            (ct + axis.x*axis.x*umct)*u.x + (axis.x*axis.y*umct-axis.z*st)*u.y + (axis.x*axis.z*umct+axis.y*st)*u.z,
            (axis.y*axis.x*umct+axis.z*st)*u.x + (ct + axis.y*axis.y*umct)*u.y + (axis.y*axis.z*umct-axis.x*st)*u.z,
            (axis.x*axis.z*umct-axis.y*st)*u.x + (axis.z*axis.y*umct+axis.x*st)*u.y + (ct + axis.z*axis.z*umct)*u.z
        };
    }
}
/**
 * @brief      Opérateur de sortie dans un flux pour le vecteur
 *
 * Cet opérateur écrit un vecteur dans un flux de sortie, de sorte que le vecteur écrit
 * soit facilement relisable par le flux d'entrée, mais pas forcément pour un être humain.
 *
 * @param      out   Le flux de sortie
 * @param[in]  p     Le vecteur à écrire dans le flux de sortie
 *
 * @return     Le nouvel état du flux de sortie
 */
inline std::ostream& 
operator << ( std::ostream& out, const K_POST::vector3d& p)
{
    out << p.x << " " << p.y << " " << p.z;
    return out;
}

/**
 * @brief      Opérateur d'entrée d'un flux pour lire un vecteur
 *
 * Cet opérateur va lire un vecteur à partir d'un flux d'entrée et écrire les coordonnées
 * de ce vecteur dans le vecteur donné à droite de l'opérateur >>.
 *
 * @param      inp   Le flux d'entrée
 * @param      p     Le vecteur recevant les données du flux d'entrée
 *
 * @return     Le nouvel état du flux d'entrée
 */
inline std::istream& 
operator >> ( std::istream& inp, K_POST::vector3d& p)
{
    inp >> p.x >> p.y >> p.z;
    return inp;
}

#endif
