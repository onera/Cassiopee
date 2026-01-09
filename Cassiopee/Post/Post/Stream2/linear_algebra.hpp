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
#ifndef _POST_STREAM_LINEAR_ALGEBRA_HPP_
#define _POST_STREAM_LINEAR_ALGEBRA_HPP_
#include <array>
#include "vector3d.hpp"

namespace K_POST
{
    /// Définition d'un vecteur de ℝ²
    using vector2d        = std::array<double,2>;
    /// Matrice définie de ℝ² dans ℝ²
    using matrix_2x2_type = std::array<std::array<double,2>,2>;
    /// Affichage d'un vecteur en deux dimensions (lisible par l'humain)
    std::ostream& operator << ( std::ostream& out, const vector2d& u);
    /// Multiplication matrice ℝ²xℝ² par une vecteur ℝ²
    vector2d operator * ( const matrix_2x2_type& A, const vector2d& u);
    /// Calcul du déterminant d'une matrice ℝ²xℝ²
    double det(const matrix_2x2_type& A);
    /// Calcul de la matrice inverse d'une matrice ℝ²xℝ². Si la matrice est non inversible, retourne une exception de type std::invalid_argument
    
    /**
     * @brief      Retourne l'inverse d'une matrice ℝ²xℝ²
     * 
     * Retourne la matrice inverse. Si la matrice A est non inversible, renvoie une exception 
     * de type std::invalid_argument.
     *
     * @param[in]  A     La matrice à inverser
     *
     * @return     Une matrice ℝ²xℝ², inverse de A
     */
    matrix_2x2_type inverse(const matrix_2x2_type& A);
    /**
     * @brief      Résoud le système linéaire A.v = u
     *
     * Si le système est non inversible, renvoie une exception de type std::invalid_argument
     * 
     * @param[in]  A     La matrice du système linéaire (lhs)
     * @param[in]  u     Le vecteur du second membre    (rhs)
     *
     * @return     Un vecteur v de ℝ² solution du système linéaire A.v=u
     */
    vector2d inverse_linear_system(const matrix_2x2_type& A, const vector2d& u);

    /**
     * @brief      Retourne une chaîne de caractère représentant la matrice ℝ²xℝ²
     *
     * Retourne une chaîne de caractère lisible par un humain (moins par un ordinateur)
     * représentant la matrice A passée en argument
     *
     * @param[in]  A     La matrice à convertir en chaîne de caractère.
     *
     * @return     La chaîne de caractère représentant la matrice A.
     */
    std::string to_string(const matrix_2x2_type& A);


    /// Matrice définie de ℝ³ dans ℝ³
    using matrix_3x3_type = std::array<std::array<double,3>,3>;
    /// Produit matrice-vecteur d'une matrice ℝ³xℝ³ avec un vecteur de ℝ³
    vector3d operator * ( const matrix_3x3_type& A, const vector3d& u);
    /// Calcul le déterminant d'une matrice ℝ³xℝ³
    double det(const matrix_3x3_type& A);
    /**
     * @brief      Retourne l'inverse d'une matrice ℝ³xℝ³
     *
     *  Calcule et retourne la matrice inverse d'une matrice ℝ³xℝ³.
     *  Si la matrice n'est pas inversible, renvoie une exception de type std::invalid_argument.
     *  L'intérêt de calcul directement la matrice inverse pour une matrice ℝ³xℝ³, c'est que cela
     *  coûte globalement moins cher de calculer la matrice inverse puis faire des produits
     *  matrice inverse avec un vecteur que de résoudre à chaque fois un système linéaire avec
     *  la même matrice mais un second membre différent. De plus, la factorisation de Gauss d'une
     *  telle matrice est aussi chère que le calcul direct de l'inverse.
     *
     * @param[in]  A     La matrice à inverser
     *
     * @return     Une matrice ℝ³xℝ³ inverse de la matrice A
     */
    matrix_3x3_type inverse(const matrix_3x3_type& A);
    /**
     * @brief      Résoud un système linéaire de dimension 3
     *
     * Résoud le système linéaire A.v=u où v est un vecteur de ℝ³ contenant les inconnues.
     * Si la matrice A n'est pas inversible, renvoie une exception de type std::invalid_argument.
     *
     * @param[in]  A     La matrice ℝ³xℝ³ du système linéaire (lhs)
     * @param[in]  u     Le second membre dans ℝ³ du système linéaire (rhs)
     *
     * @return     Le vecteur v dans ℝ³ solution du système linéaire A.v  = u
     */
    vector3d inverse_linear_system(const matrix_3x3_type& A, const vector3d& u);
    /**
     * @brief      Retourne une chaîne de caractère représentant la matrice ℝ³xℝ³
     *
     * Retourne une chaîne de caractère lisible par un humain (moins par un ordinateur)
     * représentant la matrice A passée en argument
     *
     * @param[in]  A     La matrice à convertir en chaîne de caractère.
     *
     * @return     La chaîne de caractère représentant la matrice A.
     */    
    std::string to_string(const matrix_3x3_type& A);

    /// Définition d'un vecteur de ℝ⁴
    using vector4d        = std::array<double,4>;
    /// Définition d'une matrice de ℝ⁴ vers ℝ⁴
    using matrix_4x4_type = std::array<std::array<double,4>,4>;
    /// Définition d'une paire contenant une matrice factorisée et un tableau des permutations effectuées lors de la factorisation PLU
    using factorized_matrix_4x4_type = std::pair<matrix_4x4_type,std::array<int,4>>;
    /// Produit d'une matrice ℝ⁴xℝ⁴ avec une vecteur de ℝ⁴
    vector4d operator * ( const matrix_4x4_type& A, const vector4d& u);
    /**
     * @brief      Factorisation PLU d'une matrice de ℝ⁴xℝ⁴
     *
     * Factorise une matrice ℝ⁴xℝ⁴. Si la matrice est non inversible, retourne une exception de type std::invalid_argument
     * sinon retourne dans une paire la matrice factorisée et le tableau des permutations effectuées lors de la recherche
     * partielle de pivots durant la factorisation PLU.
     *
     * @param[in]  A     La matrice à factoriser
     *
     * @return     La matrice factorisée et le tableau des permutations
     */
    factorized_matrix_4x4_type factorize(const matrix_4x4_type& A);
    /**
     * @brief      Résoud un système linéaire à partir de la forme factorisée d'une matrice ℝ⁴xℝ⁴ et d'un vecteur de ℝ⁴ second membre.
     *
     * Résoud le système linéaire A.v = P.L.U.v = u avec P le tableau de permutation issu de la factorisation d'une matrice A = P.L.U,
     * L la matrice triangulaire inférieure et U la matrice supérieure issues de la factorisation.
     *
     * @param[in]  PLU   Une matrice factorisée retournée par la fonction factorize.
     * @param[in]  u     Le second membre
     *
     * @return     La solution dans ℝ⁴ du système linéaire A.v=u
     */
    vector4d inverse_linear_system(const factorized_matrix_4x4_type& PLU, const vector4d& u);
    /// Affiche dans un format lisible par un humain un vecteur de ℝ⁴
    std::ostream& operator << ( std::ostream& out, const vector4d& u);
    /// Convertit en chaîne de caractère lisible par un humain une matrice ℝ⁴xℝ⁴
    std::string to_string(const matrix_4x4_type& A);
    /// Convertit en chaîne de caractère lisisble par un humain une matrice ℝ⁴xℝ⁴ factorisée en P.L.U.
    std::string to_string(const factorized_matrix_4x4_type& A);
}


#endif
