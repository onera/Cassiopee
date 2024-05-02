#ifndef _DIST2WALLS_EIKONAL_QUADRATICSOLVER_H_
# define _DIST2WALLS_EIKONAL_QUADRATICSOLVER_H_
#include "kcore.h"

namespace {
   /**
     * @brief Résolution de l'équation quadratique
     * @details Résolution de l'équation quadratique issue de la discrétisation décentrée du gradient de l'équation
     *          Eikonale.
     *          L'équation quadratique est la suivante : \f$\begin{array}{lclcl}
     *   \max\left(D^{x}_{-}u_{i,j,k},-D^{x}_{+}u_{i,j,k},0\right)^{2} & + &
     *   \max\left(D^{y}_{-}u_{i,j,k},-D_{+}^{y}u_{i,j,k},0\right)^{2} &  & \\
     *   & + & \max\left(D^{z}_{-}u_{i,j,k},-D^{z}_{+}u_{i,j,k},0\right)^{2} & = &
     *   \frac{1}{c_{i,j,k}^{2}}
     *   \end{array}\f$
     *   La discrétisation des opérateurs différentiels en décentré est :
     *   \f$\left\{\begin{array}{ll}
     *   D^{x}_{-}u_{i,j[,k]} = \frac{u_{i,j[,k]}-u_{i-1,j[,k]}}{h} & D^{x}_{+}u_{i,j[,k]} =
     * \frac{u_{i+1,j[,k]}-u_{i,j[,k]}}{h} \\
     *   D^{y}_{-}u_{i,j[,k]} = \frac{u_{i,j[,k]}-u_{i,j-1[,k]}}{h} & D^{y}_{+}u_{i,j[,k]} =
     * \frac{u_{i,j+1[,k]}-u_{i,j[,k]}}{h} \\
     *   D^{z}_{-}u_{i,j,k} = \frac{u_{i,j,k}-u_{i,j,k-1}}{h} & D^{z}_{+}u_{i,j,k} = \frac{u_{i,j,k+1}-u_{i,j,k}}{h} \\
     *   \end{array}
     *   \right.\f$
     *   Puis on joue sur le fait que max(u-u1,0)=max(u-u2,0)=0 si et seulement si u <= u1
     *   et max(u-u2,0) = 0 si et seulement si u <= u2
     *
     * @param u0 Plus petite valeurs des trois valeurs connues passées à l'équation
     * @param u1 Valeur intermédiaire parmi les trois valeurs connues passées à l'équation quadratique
     * @param u2 Plus grande valeur parmis les trois valeurs connues passées à l'équation
     * @param inv_f Valeur inverse de f au sommet où on veut calculer la solution de l'équation quadratique
     * @return La solution de l'équation quadratique
     */
    E_Float solveQuadraticEquation( E_Float u0, E_Float u1, E_Float u2, E_Float inv_f ) {
        // On suppose que max(u-u1,0)=max(u-u2,0)=0. L'équation est simple à résoudre :
        E_Float u = u0 + inv_f;
        if ( u > u1 ) {
            // L'hypothèse précédente était fausse, donc on suppose seulement que max(u-u2,0)=0
            u = 0.5 * ( u0 + u1 + std::sqrt( -u0 * u0 - u1 * u1 + 2 * u0 * u1 + 2. * inv_f * inv_f ) );
        }
        if ( u > u2 ) {
            E_Float s = u0 + u1 + u2;
            // L'hypothèse précédente était fausse, donc on résoud l'équation complète :
            u = ( 2 * s + std::sqrt( 4 * s * s - 12 * ( u0 * u0 + u1 * u1 + u2 * u2 - inv_f * inv_f ) ) ) / 6.;
        }
        return u;
    }

/**
 * @brief Résoud l'hamiltonien au point courant à l'ordre un
 * @details Résoud l'équation hamiltonienne \f$\left|\nabla \phi(x)\right|^{2}-\frac{1}{f(x)}=0, \forall x\in\Omega\f$
 *          au point courant $(i,j,k)$ à l'aide d'un schéma décentré au premier ordre de Godunov.
 * 
 * @param i Index i du point courant
 * @param j Index j du point courant
 * @param k Index k du point courant
 * @param ind Index global du point courant dans le tableau
 * @param ni Nombre de points dans la direction OI
 * @param nj Nombre de points dans la direction OJ
 * @param nk Nombre de points dans la direction OK
 * @param nbcells Nombre total de cellules
 * @param stride_i Saut en mémoire pour avoir la valeur précédente ou suivante en OI
 * @param stride_j Saut en mémoire pour avoir la valeur précédente ou suivante en OJ
 * @param stride_k Saut en mémoire pour avoir la valeur précédente ou suivante en OK
 * @param current_sol L'état de la solution avant calcul de la valeur au point courant
 * @param speed La fonction vitesse ( ou également appelée la fonction de poids pour la distance )
 * @return Retourne la solution du problème hamiltonien à l'ordre un en fonction de la valeur de points voisins.
 */
    E_Float solve_hamiltonian_order_1( size_t i, size_t j, size_t k, size_t ind, 
                                       size_t ni, size_t nj, size_t nk, size_t nbcells,
                                      size_t stride_i, size_t stride_j, size_t stride_k,
                                      const E_Float* sol, E_Float speed)
    {
        size_t ind_left  = ind-stride_i;
        size_t ind_right = ind+stride_i;
        size_t ind_top   = ind-stride_j;
        size_t ind_bottom= ind+stride_j;
        size_t ind_front = ind-stride_k;
        size_t ind_back  = ind+stride_k;

        E_Float u_left  = ( i > 0 ? sol[ind_left] : std::numeric_limits<E_Float>::max());
        E_Float u_right = ( i < ni-1 ? sol[ind_right] : std::numeric_limits<E_Float>::max() );
        E_Float ulr = std::min(u_left,u_right);

        E_Float u_top = ( j > 0 ? sol[ind_top] : std::numeric_limits<E_Float>::max());
        E_Float u_bot = ( j < nj-1 ? sol[ind_bottom] : std::numeric_limits<E_Float>::max());
        E_Float utb = std::min(u_top,u_bot);

        E_Float u_front = ( k > 0 ? sol[ind_front] : std::numeric_limits<E_Float>::max());
        E_Float u_back  = ( k < nk-1 ? sol[ind_back] : std::numeric_limits<E_Float>::max());
        E_Float ufb = std::min(u_front,u_back);

        if (ulr > utb) std::swap(ulr,utb);
        if (utb > ufb) std::swap(utb,ufb);
        if (ulr > utb) std::swap(ulr,utb);

        return std::min(sol[ind], solveQuadraticEquation(ulr, utb, ufb, 1./speed));
    }
}

#endif