#include "eikonalFMMSolver.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include "kcore.h"

#include "quadratic_solver.h"
// Espace de nommage anonyme pour caché les fonctions et attributs de l'extérieur.
namespace {
    // Structure permettant de stocker dans la queue des priorités une valeur u avec l'indice (i,j,k) associé
    struct queue_item {
        queue_item( ) {}
        queue_item( E_Float* v, int ii, int ij, int ik ) : pt_val( v ), i( ii ), j( ij ), k( ik ) {}
        E_Float* pt_val;
        unsigned i, j, k;
    };

    class queue_item_comp {
    public:
        queue_item_comp( ) {}
        bool operator( )( const queue_item& it1, const queue_item& it2 ) {
            //:std::cout << "Comparing" << it1.val << " with " << it2.val << std::endl;
            if ( *( it1.pt_val ) > *( it2.pt_val ) ) return true;
            return false;
        }
    };
}

namespace Eikonal {
    namespace FMM {
        void solveOnIsotropGrid( unsigned ni, unsigned nj, unsigned nk, E_Float lbx, E_Float lby, E_Float lbz,
                                 E_Float h, E_Float* sol, E_Float* speed ) {
            std::priority_queue< queue_item, std::vector< queue_item >, queue_item_comp > tests_nodes;
            std::size_t         nb_to_compute = ni * nj * nk;
            std::size_t         ninj          = ni * nj;
            std::vector< char > mask( nb_to_compute, char( 0 ) );
            // Initialisation de l'algorithme : on recherche la plus petite valeur dans la solution initiale
            // ( les valeurs non connues sont supposées mises à std::numeric_limits< E_Float >::max() )
            int     im, jm, km;
            im = 0; jm = 0; km = 0;
            E_Float v = std::numeric_limits< E_Float >::max( );
            for ( unsigned k = 0; k < nk; ++k )
                for ( unsigned j = 0; j < nj; ++j )
                    for ( unsigned i = 0; i < ni; ++i ) {
                        if ( sol[k * ninj + j * ni + i] < v ) {
                            v  = sol[k * ninj + j * ni + i];
                            im = i;
                            jm = j;
                            km = k;
                        }
                    }
            tests_nodes.push( queue_item( sol + im + jm * ni + km * ninj, im, jm, km ) );
            while ( nb_to_compute > 0 ) {
                queue_item a = tests_nodes.top( );
                tests_nodes.pop( );
                size_t indA = a.i + a.j * ni + a.k * ninj;
                mask[indA]  = char( 1 );
                nb_to_compute -= 1;
                int ind_left   = indA - 1;
                int ind_right  = indA + 1;
                int ind_top    = indA - ni;
                int ind_bottom = indA + ni;
                int ind_front  = indA - ninj;
                int ind_back   = indA + ninj;
                if ( ( a.i > 0 ) and ( mask[ind_left] != 1 ) )  // Si indice valable et valeur non connue
                {
                    // Calculer sa valeur :
                    sol[ind_left] =
                        std::min( sol[ind_left],
                                  solve_hamiltonian_order_1( a.i - 1, a.j, a.k, ind_left, ni, nj, nk, ninj * nk, 1, ni,
                                                             ninj, sol, ( speed == NULL ? 1. : speed[ind_left] ) ) );
                    if ( mask[ind_left] == 0 ) tests_nodes.push( queue_item( sol + ind_left, a.i - 1, a.j, a.k ) );
                    mask[ind_left] = 2;
                }
                if ( ( a.i < ni - 1 ) and ( mask[ind_right] != 1 ) )  // Si indice valable et valeur non connue
                {
                    // Calculer sa valeur :
                    sol[ind_right] =
                        std::min( sol[ind_right],
                                  solve_hamiltonian_order_1( a.i + 1, a.j, a.k, ind_right, ni, nj, nk, ninj * nk, 1, ni,
                                                             ninj, sol, ( speed == NULL ? 1. : speed[ind_right] ) ) );
                    if ( mask[ind_right] == 0 ) tests_nodes.push( queue_item( sol + ind_right, a.i + 1, a.j, a.k ) );
                    mask[ind_right] = 2;
                }
                if ( ( a.j > 0 ) and ( mask[ind_top] != 1 ) )  // Si indice valable et valeur non connue
                {
                    // Calculer sa valeur :
                    sol[ind_top] = std::min( sol[ind_top], solve_hamiltonian_order_1(
                                                               a.i, a.j - 1, a.k, ind_top, ni, nj, nk, ninj * nk, 1, ni,
                                                               ninj, sol, ( speed == NULL ? 1. : speed[ind_top] ) ) );
                    if ( mask[ind_top] == 0 ) tests_nodes.push( queue_item( sol + ind_top, a.i, a.j - 1, a.k ) );
                    mask[ind_top] = 2;
                }
                if ( ( a.j < nj - 1 ) and ( mask[ind_bottom] != 1 ) )  // Si indice valable et valeur non connue
                {
                    // Calculer sa valeur :
                    sol[ind_bottom] = std::min(
                        sol[ind_bottom],
                        solve_hamiltonian_order_1( a.i, a.j + 1, a.k, ind_bottom, ni, nj, nk, ninj * nk, 1, ni, ninj,
                                                   sol, ( speed == NULL ? 1. : speed[ind_bottom] ) ) );
                    if ( mask[ind_bottom] == 0 ) tests_nodes.push( queue_item( sol + ind_bottom, a.i, a.j + 1, a.k ) );
                    mask[ind_bottom] = 2;
                }
                if ( ( a.k > 0 ) and ( mask[ind_front] != 1 ) )  // Si indice valable et valeur non connue
                {
                    // Calculer sa valeur :
                    sol[ind_front] =
                        std::min( sol[ind_front],
                                  solve_hamiltonian_order_1( a.i, a.j, a.k - 1, ind_front, ni, nj, nk, ninj * nk, 1, ni,
                                                             ninj, sol, ( speed == NULL ? 1. : speed[ind_front] ) ) );
                    if ( mask[ind_front] == 0 ) tests_nodes.push( queue_item( sol + ind_front, a.i, a.j, a.k - 1 ) );
                    mask[ind_front] = 2;
                }
                if ( ( a.k < nk - 1 ) and ( mask[ind_back] != 1 ) )  // Si indice valable et valeur non connue
                {
                    // Calculer sa valeur :
                    sol[ind_back] =
                        std::min( sol[ind_back],
                                  solve_hamiltonian_order_1( a.i, a.j, a.k + 1, ind_back, ni, nj, nk, ninj * nk, 1, ni,
                                                             ninj, sol, ( speed == NULL ? 1. : speed[ind_back] ) ) );
                    if ( mask[ind_back] == 0 ) tests_nodes.push( queue_item( sol + ind_back, a.i, a.j, a.k + 1 ) );
                    mask[ind_back] = 2;
                }
            }  // While nb_to_compute
        }      // End Function solveOnIsotropGrid
    }          // End namespace FMM
}  // End namespace Eikonal