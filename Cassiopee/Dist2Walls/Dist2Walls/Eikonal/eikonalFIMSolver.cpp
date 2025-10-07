#include "eikonalFIMSolver.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "kcore.h"
#include "quadratic_solver.h"

namespace {
    struct fim_t {
        struct indice_t {
            // Pour des raisons d'optimisation, on stocke les indices (i,j,k) et l'indice global de la vue a plat
            size_t i, j, k, ind;
            indice_t( ) {}
            indice_t( size_t ii, size_t ij, size_t ik, size_t iind ) : i( ii ), j( ij ), k( ik ), ind( iind ) {}
        };
        typedef std::vector< indice_t > active_list_t;
        enum states { source = 0, active = 1, computed, far };
        std::vector< states > state_cells;
        active_list_t         L;

        /**
         * @brief Scanne le champs solution initial pour localiser les points sources ( != infini )
         * @details [long description]
         *
         * @param ni Nombre de noeuds dans la direction Oi
         * @param nj Nombre de noeuds dans la direction Oj
         * @param nk Nombre de noeuds dans la direction Ok
         * @param sol La solution initiale
         */
        fim_t( unsigned ni, unsigned nj, unsigned nk, E_Float* sol, E_Float max_float = std::numeric_limits< E_Float >::max( ) )
            : state_cells( ni * nj * nk, far ), L( ), m_dimensions( ni, nj, nk, ni * nj * nk ) {
            const unsigned stride_i = 1;
            unsigned       stride_j = ni;
            unsigned       stride_k = ni * nj;
            /* Préalloue de la place pour la liste des noeuds actifs */
            L.reserve( ni );
            /* On prealloue le deuxieme tableau de liste des noeuds actifs qui s'echangera a chaque iteration
               avec le tableau L pour pouvoir inserer et supprimer sans beaucoup d'efforts et d'allocations/desallocations
               inutiles les nouveaux et anciens noeuds actifs )
               */
            m_active_swapper.reserve( ni );
            /* Pour chaque point de la grille */
            for ( size_t k = 0; k < nk; ++k ) {
                size_t ind_k = k * stride_k;
                for ( size_t j = 0; j < nj; ++j ) {
                    size_t ind_j = ind_k + j * stride_j;
                    for ( size_t i = 0; i < ni; ++i ) {
                        size_t ind_i = ind_j + i;
                        /* Si la valeur de la solution initiale n'est pas "infini", alors c'est un point source
                        */
                        if ( sol[ind_i] < max_float ) {
                            // On a trouvé un point source, on le marque en tant que tel.
                            state_cells[ind_i] = source;
                            // On rajoute ses voisins dans la liste des points actifs :
                            active_neighbours( i, j, k, ind_i, stride_i, stride_j, stride_k );
                        }
                    }
                }
            }
        }

        /**
         * @brief Recherche les points voisins d'un point et regarde si il faut les rajouter aux points actifs
         * @details [long description]
         *
         * @param i L'indice i du point courant
         * @param j L'indice j du point courant
         * @param k L'indice k du point courant
         * @param ind L'indice global de la vue à plat
         * @param stride_i Saut en mémoire pour avoir la valeur gauche ou droite ( en nombre de valeurs )
         * @param stride_j Saut en mémoire pour avoir la valeur dessus ou dessous ( en nombre de valeurs )
         * @param stride_k Saut en mémoire pour avoir la valeur devant ou derrière ( en nombre de valeurs )
         */
        void active_neighbours( size_t i, size_t j, size_t k, size_t ind, size_t stride_i, size_t stride_j,
                                size_t stride_k ) {
            if ( i > 0 ) {
                size_t  ind_m_i = ind - stride_i;
                states& st      = state_cells[ind_m_i];
                if ( st == far ) {
                    st = active;
                    L.push_back( indice_t( i - 1, j, k, ind_m_i ) );
                }
            }
            if ( i < m_dimensions.i - 1 ) {
                size_t  ind_p_i = ind + stride_i;
                states& st      = state_cells[ind_p_i];
                if ( st == far ) {
                    st = active;
                    L.push_back( indice_t( i + 1, j, k, ind_p_i ) );
                }
            }
            if ( j > 0 ) {
                size_t  ind_m_j = ind - stride_j;
                states& st      = state_cells[ind_m_j];
                if ( st == far ) {
                    st = active;
                    L.push_back( indice_t( i, j - 1, k, ind_m_j ) );
                }
            }
            if ( j < m_dimensions.j - 1 ) {
                size_t  ind_p_j = ind + stride_j;
                states& st      = state_cells[ind_p_j];
                if ( st == far ) {
                    st = active;
                    L.push_back( indice_t( i, j + 1, k, ind_p_j ) );
                }
            }
            if ( k > 0 ) {
                size_t  ind_m_k = ind - stride_k;
                states& st      = state_cells[ind_m_k];
                if ( st == far ) {
                    st = active;
                    L.push_back( indice_t( i, j, k - 1, ind_m_k ) );
                }
            }
            if ( k < m_dimensions.k - 1 ) {
                size_t  ind_p_k = ind + stride_k;
                states& st      = state_cells[ind_p_k];
                if ( st == far ) {
                    st = active;
                    L.push_back( indice_t( i, j, k + 1, ind_p_k ) );
                }
            }
        }
        /**
         * @brief Fait avancer le front d'onde selon l'algorithme FIM avec mise à jour éventuel des noeuds déjà
         * calculés.
         * @details Algorithme de mise à jour pour faire avancer un front d'onde en mettant à jour éventuellement
         *          des noeuds déjà calculés. Un critère permet de savoir si un noeud doit être conservé ou non
         *          dans le front d'onde et on rajoute de nouveaux noeuds avec l'avancement du front.
         *          La gestion du front se fait à l'aide de deux vecteurs permettant la mise à jour des noeuds sans
         *          allocations/désallocations inutiles.
         *
         */
        void update( E_Float* sol, const E_Float* speed, const E_Float& eps = 1.E-6 ) {
            size_t ninj = m_dimensions.i * m_dimensions.j;
            m_active_swapper.resize( 0 );  // Garde la même taille réservée pour optimisation !
//#pragma omp parallel for schedule(static)// <--- plus lent avec que sans !?????
            for ( size_t i = 0; i < L.size( ); ++i ) {
                E_Float p = sol[L[i].ind];
                // Calculer sa valeur :
                E_Float q = solve_hamiltonian_order_1( L[i].i, L[i].j, L[i].k, L[i].ind, m_dimensions.i, m_dimensions.j,
                                                       m_dimensions.k, m_dimensions.ind, 1, m_dimensions.i, ninj, sol,
                                                       ( speed == NULL ? 1. : speed[L[i].ind] ) );
                sol[L[i].ind] = std::min(q,p);
                if ( std::abs( p - q ) < eps ) {
                    // Indices du point courant et nombre de points courants :
                    indice_t neighbours_ind[6];
                    size_t   nb_neighbours = 0;
                    if ( L[i].i > 0 )
                        neighbours_ind[nb_neighbours++] = indice_t( L[i].i - 1, L[i].j, L[i].k, L[i].ind - 1 );
                    if ( L[i].i < m_dimensions.i - 1 )
                        neighbours_ind[nb_neighbours++] = indice_t( L[i].i + 1, L[i].j, L[i].k, L[i].ind + 1 );
                    if ( L[i].j > 0 )
                        neighbours_ind[nb_neighbours++] =
                            indice_t( L[i].i, L[i].j - 1, L[i].k, L[i].ind - m_dimensions.i );
                    if ( L[i].j < m_dimensions.j - 1 )
                        neighbours_ind[nb_neighbours++] =
                            indice_t( L[i].i, L[i].j + 1, L[i].k, L[i].ind + m_dimensions.i );
                    if ( L[i].k > 0 )
                        neighbours_ind[nb_neighbours++] = indice_t( L[i].i, L[i].j, L[i].k - 1, L[i].ind - ninj );
                    if ( L[i].k < m_dimensions.k - 1 )
                        neighbours_ind[nb_neighbours++] = indice_t( L[i].i, L[i].j, L[i].k + 1, L[i].ind + ninj );
                    // On parcourt chaque voisin du point courant :
                    for ( unsigned n = 0; n < nb_neighbours; ++n ) {
                        if ( state_cells[neighbours_ind[n].ind] == far ) {
                            size_t ii, ij, ik, iind;
                            ii        = neighbours_ind[n].i;
                            ij        = neighbours_ind[n].j;
                            ik        = neighbours_ind[n].k;
                            iind      = neighbours_ind[n].ind;
                            E_Float p = sol[iind];
                            E_Float q = solve_hamiltonian_order_1( ii, ij, ik, iind, m_dimensions.i, m_dimensions.j,
                                                                   m_dimensions.k, m_dimensions.ind, 1, m_dimensions.i,
                                                                   ninj, sol, ( speed == NULL ? 1. : speed[iind] ) );
//#pragma omp critical
                            if ( p > q ) {
                                sol[iind]         = q;
                                state_cells[iind] = active;
                                m_active_swapper.push_back( neighbours_ind[n] );
                            }
                        }
                    }
                    state_cells[L[i].ind] = computed;
                } else {
//#pragma omp critical
                    m_active_swapper.push_back( L[i] );
                }
            }  // End for (size_t i=0; ...)
            // On échange le nouveau front avec l'ancien.
            m_active_swapper.swap( L );
        }

    private:
        indice_t m_dimensions;
        /*
         * L'idée ici est de gérer la mise à jour de la liste active en switchant entre une nouvelle liste crée à partir
         * de l'ancienne
         * et l'ancienne liste elle-même. Cela permet :
         * 1. De paralléliser la boucle pour la fonction update
         * 2. D'éviter trop d'allocation/désallocation
         * Par contre cela double la mémoire prise par la liste de noeuds actifs.
         */
        active_list_t m_active_swapper;
    };

/*    struct multi_blk_fim_t {
    private:

    };*/
}

namespace Eikonal {
    namespace FIM {
        /**
         * @brief Résolution de l'équation Eikonal à l'aide de l'algorithme FIM
         * @details Résolution de l'équation Eikonal à l'aide de l'algorithme FIM séquentiel.
         *
         * @param ni Nombre d'inconnus dans la direction i
         * @param nj Nombre d'inconnus dans la direction j
         * @param nk Nombre d'inconnus dans la direction k
         * @param lbx Abcisse du premier point
         * @param lby [description]
         * @param lbz [description]
         * @param h [description]
         * @param sol [description]
         * @param speed [description]
         * @param max_float [description]
         */
        void solveOnIsotropGrid( unsigned ni, unsigned nj, unsigned nk, E_Float lbx, E_Float lby, E_Float lbz,
                                 E_Float h, E_Float* sol, E_Float* speed,
                                 E_Float max_float ) {
            // Initialisation de l'algorithme FIM :
            fim_t fim( ni, nj, nk, sol, max_float );
            // Mise à jour des points de L :
            while (fim.L.empty() == false)// Tant qu'il y a des noeuds actifs
            {
                fim.update(sol, speed);
            }
        }
    }
}
