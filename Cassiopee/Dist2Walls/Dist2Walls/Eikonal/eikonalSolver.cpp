// Solveur Eikonal sur une grille cartesienne

#include "eikonalSolver.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <list>
#include <vector>

// Resoud
static E_Float solveQuadratic( E_Float a, E_Float b, E_Float c, E_Float evalf ) {
    E_Float u = c + 1. / evalf;
    if ( u <= b ) return u;
    u = ( b + c + sqrt( -b * b - c * c + 2 * b * c + 2. / ( evalf * evalf ) ) ) * 0.5;
    if ( u <= a ) return u;
    u = ( 2 * ( a + b + c ) +
          sqrt( 4 * ( a + b + c ) * ( a + b + c ) - 12 * ( a * a + b * b + c * c - 1 / ( evalf * evalf ) ) ) ) /
        6.;
    return u;
}

static E_Float updateSolOn( unsigned i, unsigned j, unsigned k, unsigned stride_i, unsigned stride_j, unsigned stride_k,
                            unsigned ni, unsigned nj, unsigned nk, E_Float cx, E_Float cy, E_Float cz,
                            // speedFctPtr speed,
                            E_Float* speed, const E_Float* sol ) {
    // std::cerr << "Traitement point " << i << ", " << j << ", " << k << std::endl;
    unsigned ind = i * stride_i + j * stride_j + k * stride_k;
    E_Float  abc[3];
    abc[0] = abc[1] = abc[2] = std::numeric_limits< E_Float >::max( );
    if ( i > 0 ) abc[0]      = sol[ind - stride_i];
    if ( i < ni - 1 ) abc[0] = std::min( abc[0], sol[ind + stride_i] );
    if ( j > 0 ) abc[1]      = sol[ind - stride_j];
    if ( j < nj - 1 ) abc[1] = std::min( abc[1], sol[ind + stride_j] );
    if ( k > 0 ) abc[2]      = sol[ind - stride_k];
    if ( k < nk - 1 ) abc[2] = std::min( abc[2], sol[ind + stride_k] );
    if ( abc[0] > abc[1] ) std::swap( abc[0], abc[1] );
    if ( abc[1] > abc[2] ) std::swap( abc[1], abc[2] );
    if ( abc[0] > abc[1] ) std::swap( abc[0], abc[1] );
    // std::cerr << "a,b,c = " << abc[0] << " : " << abc[1] << " : " << abc[2] << std::endl;
    // E_Float f = speed(cx, cy, cz);
    E_Float f = speed[ind];
    E_Float q = solveQuadratic( abc[2], abc[1], abc[0], f );
    // std::cerr << "q = " << q << std::endl;
    return q;
}

// Modifie sol et met max pour les pts non source
static void initSolutionFIM( unsigned ni, unsigned nj, unsigned nk, E_Float lbx, E_Float lby, E_Float lbz,
                             // E_Float h,
                             E_Float* sol ) {
    for ( unsigned ind = 0; ind < ni * nj * nk; ind++ ) {
        if ( sol[ind] > 1.e10 ) sol[ind] = std::numeric_limits< E_Float >::max( );
    }
}

static void ind2Blk( unsigned szBlki, unsigned szBlkj, unsigned szBlkk, unsigned i, unsigned j, unsigned k,
                     unsigned& iBlk, unsigned& jBlk, unsigned& kBlk ) {
    kBlk = k / szBlkk;
    jBlk = j / szBlkj;
    iBlk = i / szBlki;
}

// initialise le front, on suppose que sol vaut 0 ou max
static void initActiveList( unsigned ni, unsigned nj, unsigned nk, unsigned stride_i, unsigned stride_j,
                            unsigned stride_k, const E_Float* sol, std::list< unsigned >& L, unsigned char* Lmask,
                            unsigned stride_bi = 0, unsigned stride_bj = 0, unsigned stride_bk = 0, unsigned szBlki = 0,
                            unsigned szBlkj = 0, unsigned szBlkk = 0, unsigned char* Cbb = NULL ) {
    unsigned iBlk, jBlk, kBlk;
    size_t   ind = 0;
    for ( unsigned k = 0; k < nk; ++k ) {
        for ( unsigned j = 0; j < nj; ++j ) {
            for ( unsigned i = 0; i < ni; ++i ) {
                if ( sol[ind] == std::numeric_limits< E_Float >::max( ) ) {  // Le point courant n'est pas une source :
                    // bool isNeighbourOfSource = false;
                    if ( ( i > 0 ) && ( sol[ind - stride_i] < std::numeric_limits< E_Float >::max( ) ) ) {
                        L.push_back( ind );
                        Lmask[ind] = 1;
                        if ( Cbb ) {
                            // Quel bloc ?
                            ind2Blk( szBlki, szBlkj, szBlkk, i, j, k, iBlk, jBlk, kBlk );
                            Cbb[iBlk * stride_bi + jBlk * stride_bj + kBlk * stride_bk] = 1;
                        }
                    } else if ( ( i < ni - 1 ) && ( sol[ind + stride_i] < std::numeric_limits< E_Float >::max( ) ) ) {
                        L.push_back( ind );
                        Lmask[ind] = 1;
                        if ( Cbb ) {
                            // Quel bloc ?
                            ind2Blk( szBlki, szBlkj, szBlkk, i, j, k, iBlk, jBlk, kBlk );
                            Cbb[iBlk * stride_bi + jBlk * stride_bj + kBlk * stride_bk] = 1;
                        }
                    } else if ( ( j > 0 ) && ( sol[ind - stride_j] < std::numeric_limits< E_Float >::max( ) ) ) {
                        L.push_back( ind );
                        Lmask[ind] = 1;
                        if ( Cbb ) {
                            // Quel bloc ?
                            ind2Blk( szBlki, szBlkj, szBlkk, i, j, k, iBlk, jBlk, kBlk );
                            Cbb[iBlk * stride_bi + jBlk * stride_bj + kBlk * stride_bk] = 1;
                        }
                    } else if ( ( j < nj - 1 ) && ( sol[ind + stride_j] < std::numeric_limits< E_Float >::max( ) ) ) {
                        L.push_back( ind );
                        Lmask[ind] = 1;
                        if ( Cbb ) {
                            // Quel bloc ?
                            ind2Blk( szBlki, szBlkj, szBlkk, i, j, k, iBlk, jBlk, kBlk );
                            Cbb[iBlk * stride_bi + jBlk * stride_bj + kBlk * stride_bk] = 1;
                        }
                    } else if ( ( k > 0 ) && ( sol[ind - stride_k] < std::numeric_limits< E_Float >::max( ) ) ) {
                        L.push_back( ind );
                        Lmask[ind] = 1;
                        if ( Cbb ) {
                            // Quel bloc ?
                            ind2Blk( szBlki, szBlkj, szBlkk, i, j, k, iBlk, jBlk, kBlk );
                            Cbb[iBlk * stride_bi + jBlk * stride_bj + kBlk * stride_bk] = 1;
                        }
                    } else if ( ( k < nk - 1 ) && ( sol[ind + stride_k] < std::numeric_limits< E_Float >::max( ) ) ) {
                        L.push_back( ind );
                        Lmask[ind] = 1;
                        if ( Cbb ) {
                            // Quel bloc ?
                            ind2Blk( szBlki, szBlkj, szBlkk, i, j, k, iBlk, jBlk, kBlk );
                            Cbb[iBlk * stride_bi + jBlk * stride_bj + kBlk * stride_bk] = 1;
                        }
                    }
                }
                ind++;
            }  // Fin (for unsigned i ... )
        }      // Fin ( for unsigned j ... )
    }          // Fin ( for unsigned k ... )
}

// IN: ni,nj,nk: nbre de pts de la grille cartesienne
// IN: lbx,lby,lbz: coord du premier pt de la grille cartesienne
// IN: h: pas de la grille cartesienne (identique dans toutes les directions)
// IN: speed: vitesse (aux noeuds) dans l'equation eikonale
// IN: sol: solution (aux noeuds)
void solveEikonalOnIsotropGrid( unsigned ni, unsigned nj, unsigned nk, E_Float lbx, E_Float lby, E_Float lbz, E_Float h,
                                // isSourcePtr isSrc,
                                // speedFctPtr speed,
                                // bool* isSrc,
                                E_Float* speed, E_Float* sol ) {
    unsigned       stride_i = 1;
    unsigned       stride_j = ni;
    unsigned       stride_k = ni * nj;
    unsigned char* Lmask    = new unsigned char[ni * nj * nk];
    std::fill_n( Lmask, ni * nj * nk, (unsigned char)( 0 ) );
    std::list< unsigned > L;  //<---- activeList;
    // Initialisation de la solution pour la methode FIM
    initSolutionFIM( ni, nj, nk, lbx, lby, lbz, sol );
    // ........................................................................
    // Initialisation de la activeList pour la methode FIM :
    initActiveList( ni, nj, nk, stride_i, stride_j, stride_k, sol, L, Lmask );
    // ........................................................................
    // Iterations avec maj de l'activeList : les valeurs de sol doivent converger vers la solution
    while ( !L.empty( ) ) {
        // std::cerr << "Size of L " << L.size() << std::endl;
        std::list< unsigned >::iterator itL = L.begin( );
        std::list< unsigned >           append;
        while ( itL != L.end( ) ) {
            std::list< unsigned >::iterator it = itL;
            it++;  // it pointe sur l'element suivant...
            unsigned ind = *itL;
            unsigned k   = ind / stride_k;
            unsigned j   = ( ind - k * stride_k ) / stride_j;
            unsigned i   = ( ind - k * stride_k - j * stride_j );
            E_Float  cx = lbx + i * h, cy = lby + j * h, cz = lbz + k * h;
            E_Float  p = sol[ind];
            E_Float  q = updateSolOn( i, j, k, stride_i, stride_j, stride_k, ni, nj, nk, cx, cy, cz, speed, sol );
            sol[ind]   = std::min( q, p );
            if ( std::fabs( p - q ) < 1.E-6 ) {  // Valeurs d'epsilon a changer en parametre ( externe ? ) )
                // Pour chaque voisin xnb de x:
                if ( ( i > 0 ) &&
                     ( Lmask[ind - stride_i] == 0 ) ) {  // Si i-1 existe et n'est pas dans la liste active :
                    p = sol[ind - stride_i];
                    q = updateSolOn( i - 1, j, k, stride_i, stride_j, stride_k, ni, nj, nk, cx - h, cy, cz, speed,
                                     sol );
                    if ( p > q ) {
                        sol[ind - stride_i] = q;
                        append.push_back( ind - stride_i );
                        Lmask[ind - stride_i] = 1;
                    }
                }

                if ( ( i < ni - 1 ) &&
                     ( Lmask[ind + stride_i] == 0 ) ) {  // Si i-1 existe et n'est pas dans la liste active :
                    p = sol[ind + stride_i];
                    q = updateSolOn( i + 1, j, k, stride_i, stride_j, stride_k, ni, nj, nk, cx + h, cy, cz, speed,
                                     sol );
                    if ( p > q ) {
                        sol[ind + stride_i] = q;
                        append.push_back( ind + stride_i );
                        Lmask[ind + stride_i] = 1;
                    }
                }

                if ( ( j > 0 ) &&
                     ( Lmask[ind - stride_j] == 0 ) ) {  // Si j-1 existe et n'est pas dans la liste active :
                    p = sol[ind - stride_j];
                    q = updateSolOn( i, j - 1, k, stride_i, stride_j, stride_k, ni, nj, nk, cx, cy - h, cz, speed,
                                     sol );
                    if ( p > q ) {
                        sol[ind - stride_j] = q;
                        append.push_back( ind - stride_j );
                        Lmask[ind - stride_j] = 1;
                    }
                }

                if ( ( j < nj - 1 ) &&
                     ( Lmask[ind + stride_j] == 0 ) ) {  // Si j+1 existe et n'est pas dans la liste active :
                    p = sol[ind + stride_j];
                    q = updateSolOn( i, j + 1, k, stride_i, stride_j, stride_k, ni, nj, nk, cx, cy + h, cz, speed,
                                     sol );
                    if ( p > q ) {
                        sol[ind + stride_j] = q;
                        append.push_back( ind + stride_j );
                        Lmask[ind + stride_j] = 1;
                    }
                }

                if ( ( k > 0 ) &&
                     ( Lmask[ind - stride_k] == 0 ) ) {  // Si k-1 existe et n'est pas dans la liste active :
                    p = sol[ind - stride_k];
                    q = updateSolOn( i, j, k - 1, stride_i, stride_j, stride_k, ni, nj, nk, cx, cy, cz - h, speed,
                                     sol );
                    if ( p > q ) {
                        sol[ind - stride_k] = q;
                        append.push_back( ind - stride_k );
                        Lmask[ind - stride_k] = 1;
                    }
                }

                if ( ( k < nk - 1 ) &&
                     ( Lmask[ind + stride_k] == 0 ) ) {  // Si k+1 existe et n'est pas dans la liste active :
                    p = sol[ind + stride_k];
                    q = updateSolOn( i, j, k + 1, stride_i, stride_j, stride_k, ni, nj, nk, cx, cy, cz + h, speed,
                                     sol );
                    if ( p > q ) {
                        sol[ind + stride_k] = q;
                        append.push_back( ind + stride_k );
                        Lmask[ind + stride_k] = 1;
                    }
                }
                // On en a fini avec les voisins.
                // On enleve x de la liste active :
                L.erase( itL );
                // printf("inside %d\n", L.size());
                Lmask[ind] = 0;
            }          // Fin if ( abs(p-q) < epsilon
            itL = it;  // itL pointe sur l'element "suivant"
        }              // ( while (itL != end() );
        L.splice( L.end( ), append );
    }  // Fin while (!L.empty())
    delete[] Lmask;
}
// ========================================================================
static void solveBlock( unsigned iBlk, unsigned jBlk, unsigned kBlk, unsigned nbiBlk, unsigned nbjBlk, unsigned nbkBlk,
                        unsigned ni, unsigned nj, unsigned nk, E_Float lbx, E_Float lby, E_Float lbz, E_Float h,
                        unsigned stride_i, unsigned stride_j, unsigned stride_k,
                        // speedFctPtr speed,
                        E_Float* speed, unsigned char* Cbn, E_Float* sol ) {
    unsigned startk = kBlk * nbkBlk, startj = jBlk * nbjBlk, starti = iBlk * nbiBlk;
    unsigned endk = std::min( startk + nbkBlk, nk );
    unsigned endj = std::min( startj + nbjBlk, nj );
    unsigned endi = std::min( starti + nbiBlk, ni );
    for ( unsigned k = startk; k < endk; ++k ) {
        unsigned indk = k * stride_k;
        for ( unsigned j = startj; j < endj; ++j ) {
            unsigned indj = indk + j * stride_j;
            for ( unsigned i = starti; i < endi; ++i ) {
                unsigned ind = indj + i;
                if ( Cbn[ind] == 1 ) {
                    E_Float cx, cy, cz;
                    cx        = lbx + i * h;
                    cy        = lby + j * h;
                    cz        = lbz + k * h;
                    E_Float p = sol[ind];
                    E_Float q =
                        updateSolOn( i, j, k, stride_i, stride_j, stride_k, ni, nj, nk, cx, cy, cz, speed, sol );
                    sol[ind] = std::min( q, p );
                    if ( std::fabs( p - q ) < 1.E-6 )  // Valeurs d'epsilon à changer en parametre ( externe ? ) )
                    {                                  // Pour chaque voisin xnb de x :
                        if ( ( i > 0 ) &&
                             ( Cbn[ind - stride_i] == 0 ) )  // Si i-1 existe et n est pas dans la liste active :
                        {
                            p = sol[ind - stride_i];
                            q = updateSolOn( i - 1, j, k, stride_i, stride_j, stride_k, ni, nj, nk, cx - h, cy, cz,
                                             speed, sol );
                            if ( p > q ) {
                                sol[ind - stride_i] = q;
                                Cbn[ind - stride_i] = 1;
                            }
                        }
                        if ( ( i < ni - 1 ) &&
                             ( Cbn[ind + stride_i] == 0 ) )  // Si i-1 existe et n est pas dans la liste active :
                        {
                            p = sol[ind + stride_i];
                            q = updateSolOn( i + 1, j, k, stride_i, stride_j, stride_k, ni, nj, nk, cx + h, cy, cz,
                                             speed, sol );
                            if ( p > q ) {
                                sol[ind + stride_i] = q;
                                Cbn[ind + stride_i] = 1;
                            }
                        }
                        if ( ( j > 0 ) &&
                             ( Cbn[ind - stride_j] == 0 ) )  // Si j-1 existe et n est pas dans la liste active :
                        {
                            p = sol[ind - stride_j];
                            q = updateSolOn( i, j - 1, k, stride_i, stride_j, stride_k, ni, nj, nk, cx, cy - h, cz,
                                             speed, sol );
                            if ( p > q ) {
                                sol[ind - stride_j] = q;
                                Cbn[ind - stride_j] = 1;
                            }
                        }
                        if ( ( j < nj - 1 ) &&
                             ( Cbn[ind + stride_j] == 0 ) )  // Si j+1 existe et n est pas dans la liste active :
                        {
                            p = sol[ind + stride_j];
                            q = updateSolOn( i, j + 1, k, stride_i, stride_j, stride_k, ni, nj, nk, cx, cy + h, cz,
                                             speed, sol );
                            if ( p > q ) {
                                sol[ind + stride_j] = q;
                                Cbn[ind + stride_j] = 1;
                            }
                        }
                        if ( ( k > 0 ) &&
                             ( Cbn[ind - stride_k] == 0 ) )  // Si k-1 existe et n'est pas dans la liste active :
                        {
                            p = sol[ind - stride_k];
                            q = updateSolOn( i, j, k - 1, stride_i, stride_j, stride_k, ni, nj, nk, cx, cy, cz - h,
                                             speed, sol );
                            if ( p > q ) {
                                sol[ind - stride_k] = q;
                                Cbn[ind - stride_k] = 1;
                            }
                        }
                        if ( ( k < nk - 1 ) &&
                             ( Cbn[ind + stride_k] == 0 ) )  // Si k+1 existe et n est pas dans la liste active :
                        {
                            p = sol[ind + stride_k];
                            q = updateSolOn( i, j, k + 1, stride_i, stride_j, stride_k, ni, nj, nk, cx, cy, cz + h,
                                             speed, sol );
                            if ( p > q ) {
                                sol[ind + stride_k] = q;
                                Cbn[ind + stride_k] = 1;
                            }
                        }
                        Cbn[ind] = 0;
                    }  // abs(p-q)<tol
                }      // test cbn = 1
            }          //  i
        }              // j
    }                  // k
}
// ------------------------------------------------------------------------
unsigned char convReduction( unsigned iBlk, unsigned jBlk, unsigned kBlk, unsigned nbi, unsigned nbj, unsigned nbk,
                             unsigned nbiBlk, unsigned nbjBlk, unsigned nbkBlk, unsigned ni, unsigned nj, unsigned nk,
                             unsigned stride_i, unsigned stride_j, unsigned stride_k, const unsigned char* Cbn ) {
    // unsigned indBlk = iBlk + jBlk*nbiBlk + kBlk*nbiBlk*nbjBlk;
    unsigned char Cbb    = 0;
    unsigned      startk = kBlk * nbkBlk, startj = jBlk * nbjBlk, starti = iBlk * nbiBlk;
    unsigned      endk = std::min( startk + nbkBlk, nk );
    unsigned      endj = std::min( startj + nbjBlk, nj );
    unsigned      endi = std::min( starti + nbiBlk, ni );
    for ( unsigned k = startk; k < endk; ++k ) {
        unsigned indk = k * stride_k;
        for ( unsigned j = startj; j < endj; ++j ) {
            unsigned indj = indk + j * stride_j;
            for ( unsigned i = starti; i < endi; ++i ) {
                unsigned ind = indj + i * stride_i;
                Cbb |= Cbn[ind];
            }
        }
    }
    return Cbb;
}

// ------------------------------------------------------------------------
void blockFIM( unsigned ni, unsigned nj, unsigned nk, E_Float lbx, E_Float lby, E_Float lbz, E_Float h, unsigned niBlk,
               unsigned njBlk, unsigned nkBlk, unsigned nbSubIter, E_Float* speed, E_Float* sol ) {
    initSolutionFIM( ni, nj, nk, lbx, lby, lbz, sol );
    unsigned       stride_i = 1, stride_j = ni, stride_k = ni * nj;
    unsigned char* Lmask = new unsigned char[ni * nj * nk];
    std::fill_n( Lmask, ni * nj * nk, (unsigned char)( 0 ) );
    std::list< unsigned > LNd;  //<---- activeList;
    // Calcul du nombre de blocs ( par direction ) :
    unsigned nbi, nbj, nbk;
    nbi                  = ( ni + niBlk - 1 ) / niBlk;
    nbj                  = ( nj + njBlk - 1 ) / njBlk;
    nbk                  = ( nk + nkBlk - 1 ) / nkBlk;
    unsigned       nbBlk = nbi * nbj * nbk;
    unsigned char* Cbb   = new unsigned char[nbBlk];
    std::fill_n( Cbb, nbBlk, (unsigned char)( 0 ) );

    initActiveList( ni, nj, nk, stride_i, stride_j, stride_k, sol, LNd, Lmask, 1, nbi, nbi * nbj, niBlk, njBlk, nkBlk,
                    Cbb );
    std::vector< unsigned > L( nbBlk );
    unsigned                nbL = 0;
    for ( unsigned ib = 0; ib < nbBlk; ++ib ) {
        if ( Cbb[ib] == 1 ) {
            L[nbL] = ib;
            nbL++;
        }
    }
    // Update blocks b in active list L :
    while ( nbL > 0 ) {
// Update active blocks :
#pragma omp parallel for
        for ( unsigned ib = 0; ib < nbL; ++ib ) {
            unsigned indBlk = L[ib];
            unsigned kBlk   = indBlk / ( nbi * nbj );
            unsigned jBlk   = ( indBlk - kBlk * nbi * nbj ) / nbi;
            unsigned iBlk   = ( indBlk - nbi * ( kBlk * nbj + jBlk ) );
            for ( unsigned i = 0; i < nbSubIter; ++i ) {
                solveBlock( iBlk, jBlk, kBlk, niBlk, njBlk, nkBlk, ni, nj, nk, lbx, lby, lbz, h, stride_i, stride_j,
                            stride_k, speed, Lmask, sol );
            }

            // Reduction Cbn -> Cbb
            Cbb[indBlk] = convReduction( iBlk, jBlk, kBlk, nbi, nbj, nbk, niBlk, njBlk, nkBlk, ni, nj, nk, stride_i,
                                         stride_j, stride_k, Lmask );
        }  // for unsigned ib parallel
           // Step 2 : Check neighbours blocks :
#pragma omp parallel for
        for ( unsigned ib = 0; ib < nbL; ++ib ) {
            unsigned indBlk = L[ib];
            if ( Cbb[indBlk] == 0 ) {
                unsigned kBlk = indBlk / ( nbi * nbj );
                unsigned jBlk = ( indBlk - kBlk * nbi * nbj ) / nbi;
                unsigned iBlk = ( indBlk - nbi * ( kBlk * nbj + jBlk ) );
                if ( iBlk > 0 ) {
                    solveBlock( iBlk - 1, jBlk, kBlk, niBlk, njBlk, nkBlk, ni, nj, nk, lbx, lby, lbz, h, stride_i,
                                stride_j, stride_k, speed, Lmask, sol );
                    Cbb[indBlk - 1] = convReduction( iBlk - 1, jBlk, kBlk, nbi, nbj, nbk, niBlk, njBlk, nkBlk, ni, nj,
                                                     nk, stride_i, stride_j, stride_k, Lmask );
                }
                if ( iBlk < nbi - 1 ) {
                    solveBlock( iBlk + 1, jBlk, kBlk, niBlk, njBlk, nkBlk, ni, nj, nk, lbx, lby, lbz, h, stride_i,
                                stride_j, stride_k, speed, Lmask, sol );
                    Cbb[indBlk + 1] = convReduction( iBlk + 1, jBlk, kBlk, nbi, nbj, nbk, niBlk, njBlk, nkBlk, ni, nj,
                                                     nk, stride_i, stride_j, stride_k, Lmask );
                }
                // ........................................................................
                if ( jBlk > 0 ) {
                    solveBlock( iBlk, jBlk - 1, kBlk, niBlk, njBlk, nkBlk, ni, nj, nk, lbx, lby, lbz, h, stride_i,
                                stride_j, stride_k, speed, Lmask, sol );
                    Cbb[indBlk - nbi] = convReduction( iBlk, jBlk - 1, kBlk, nbi, nbj, nbk, niBlk, njBlk, nkBlk, ni, nj,
                                                       nk, stride_i, stride_j, stride_k, Lmask );
                }
                if ( jBlk < nbj - 1 ) {
                    solveBlock( iBlk, jBlk + 1, kBlk, niBlk, njBlk, nkBlk, ni, nj, nk, lbx, lby, lbz, h, stride_i,
                                stride_j, stride_k, speed, Lmask, sol );
                    Cbb[indBlk + nbi] = convReduction( iBlk, jBlk + 1, kBlk, nbi, nbj, nbk, niBlk, njBlk, nkBlk, ni, nj,
                                                       nk, stride_i, stride_j, stride_k, Lmask );
                }
                // ........................................................................
                if ( kBlk > 0 ) {
                    solveBlock( iBlk, jBlk, kBlk - 1, niBlk, njBlk, nkBlk, ni, nj, nk, lbx, lby, lbz, h, stride_i,
                                stride_j, stride_k, speed, Lmask, sol );
                    Cbb[indBlk - nbi * nbj] = convReduction( iBlk, jBlk, kBlk - 1, nbi, nbj, nbk, niBlk, njBlk, nkBlk,
                                                             ni, nj, nk, stride_i, stride_j, stride_k, Lmask );
                }
                if ( kBlk < nbk - 1 ) {
                    solveBlock( iBlk, jBlk, kBlk + 1, niBlk, njBlk, nkBlk, ni, nj, nk, lbx, lby, lbz, h, stride_i,
                                stride_j, stride_k, speed, Lmask, sol );
                    Cbb[indBlk + nbi * nbj] = convReduction( iBlk, jBlk, kBlk + 1, nbi, nbj, nbk, niBlk, njBlk, nkBlk,
                                                             ni, nj, nk, stride_i, stride_j, stride_k, Lmask );
                }
            }  // if ( Cbb[
        }      // for unsigned ib parallel

        // Step 3 : Update active list L:
        nbL = 0;
        for ( unsigned ib = 0; ib < nbBlk; ++ib ) {
            if ( Cbb[ib] == 1 ) {
                L[nbL] = ib;
                nbL++;
            }
        }
    }  // while L not empty
}
