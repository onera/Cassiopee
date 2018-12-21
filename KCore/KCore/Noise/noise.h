/*    
    Copyright 2013-2019 Onera.

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

#ifndef _KCORE_NOISE_NOISE_H_
#define _KCORE_NOISE_NOISE_H_

#include "Def/DefTypes.h"

namespace K_NOISE
{
  /* Minimum standard random generator
     Retourne un reel entre 0 et 1.
     Tirages correles. A initialiser avec un entier negatif. */
  double stdRand(E_LONG *idum);

  /* Minimum shuffle random generator
     Retourne un reel entre 0 et 1.
     Non correle pour un nbre de tirages < a 10e8.
     A initialiser avec un entier negatif. */
  double shuffleRand(E_LONG *idum);

  /* Minimum combined random generator
     Retourne un reel entre 0 et 1.
     Tres longue periode.
     A initialiser avec un entier negatif. */
  double longRand(E_LONG *idum);

  /* Structure pour stocker les donnees du perlin noise. */
#define MAXB 0x100
  typedef struct
  {
      int p[MAXB + MAXB + 2];
      double g3[MAXB + MAXB + 2][3];
      double g2[MAXB + MAXB + 2][2];
      double g1[MAXB + MAXB + 2];
      int start;
      int B;
      int BM;
  } PDS;
    
  // Fonctions internes pour le perlin noise
  void normalize2(double* v);
  void normalize3(double* v);
  void initNoise(PDS& data);
  double noise1(double arg, PDS& data);
  double noise2(double vec[2], PDS& data);
  double noise3(double vec[3], PDS& data);

  /* Initialise le perlin noise. A appeler avant tout autre appel. */
  void initPerlinNoise(int frequency, PDS& data);

  /* Rend une valeur de bruit coherent entre -1 et 1 quand x varie.
     alpha est le poids dans la somme des bruits unitaires.
     Generalement, alpha=2, plus alpha approche 1 plus la fonction est bruitee.
     Beta est le pas, generalement 2. */
  double perlinNoise1D(double x, double alpha, double beta, int n, PDS& data);
  double perlinNoise2D(double x, double y, double alpha, double beta, int n, 
                       PDS& data);
  double perlinNoise3D(double x, double y, double z, double alpha,
                       double beta, int n, PDS& data);
}
#endif
