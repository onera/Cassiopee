#ifndef _EIKONAL_HPP_
# define _EIKONAL_HPP_

//using isSourcePtr = bool (*) (double x, double y, double z);
//using speedFctPtr = double (*) (double x, double y, double z);

#include "kcore.h"

/*
  Definition grille:
     (lbx, lby, lbz): Coordonnee coin gauche de la grille
     h: Pas d'espace (le meme en x, y et z)
     (ni, nj, nk): Nombre de noeuds dans chaque direction

  Definition de l'equation Eikonal:
     isSrc: Dit si un noeud de la grille est un noeud source ou non
     speed: Fonction de vitesse de l'equation (= 1 si uniquement distance cartesienne)

  Sortie:
     sol: solution sous la forme: grid[i+ni*(j+k*nk)]
*/
void solveEikonalOnIsotropGrid(unsigned ni, unsigned nj, unsigned nk,
                               E_Float lbx, E_Float lby, E_Float lbz,
                               E_Float h,
                               E_Float* speed,
                               E_Float* sol);
/*
  Same algo as above but with block approach and OpenMP instructions
  niBlk, njBlk and nkBlk give the size of a block for cache optimization
*/
void blockFIM(unsigned ni, unsigned nj, unsigned nk,
              E_Float lbx, E_Float lby, E_Float lbz,
              E_Float h,
              unsigned niBlk, unsigned njBlk, unsigned nkBlk,
              unsigned nbSubIter,
              E_Float* speed,
              E_Float* sol);
#endif
