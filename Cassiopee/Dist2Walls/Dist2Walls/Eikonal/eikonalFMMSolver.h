#ifndef _EIKONALFMMSOLVER_HPP_
#define _EIKONALFMMSOLVER_HPP_

#include "kcore.h"

namespace Eikonal {
    namespace FMM {
      // Algorithme FMM pour résoudre une équation Eikonale sur une grille cartésienne régulière
      // où la vitesse en chaque point peut être définie. Si elle n'est pas définie, la valeur de la vitesse
      // sera mise à un, et on calculera alors un champs de distance

      /**
       * @brief Résoud une équation Eikonale
       * @details Résoud une équation Eikonale à l'aide de la méthode Fast Marching Method (FMM) qui consiste
       *          à prendre un point source où la valeur de la solution est connue et faire avancer un front
       *          d'onde pour calculer la valeur sur les autres sommets de la grille. Un champs de vitesse peut être
       *          défini traduisant des différences de milieu où l'onde se propage à des vitesses différentes. Si ce
       *          champs est omis en paramètre, il vaudra un et la résolution de l'équation Eikonale reviendra à
       *          résoudre un champs de distance.
       *
       * @param ni Nombre de sommets dans la direction Oi
       * @param nj Nombre de sommets dans la direction Oj
       * @param nk Nombre de sommets dans la direction Ok
       * @param lbx Abcisse du sommet d'indice (0,0,0)
       * @param lby Ordonnée du sommet d'indice (0,0,0)
       * @param lbz Profondeur du sommet d'indice (0,0,0)
       * @param h Le pas du maillage ( le même dans toutes les directions )
       * @param sol Le vecteur solution. Il doit être rempli avec une valeur de type
       *            std::numeric_limits< E_Float >::max( ) sur les points dont on ne connaît pas
       *            initialement la valeur de la solution et la valeur de la solution ( fixée )
       *            sur les points "connus" initialement.
       * @param speed Le champs de vitesse en chacun des sommets du maillage. Si le pointeur est NULL,
       *              on suppose que le champs de vitesse est constant et vaut un sur chaque sommet du maillage.
       */
        void solveOnIsotropGrid( unsigned ni, unsigned nj, unsigned nk, E_Float lbx, E_Float lby, E_Float lbz,
                                 E_Float h, E_Float* sol, E_Float* speed = NULL );
    }
}
#endif
