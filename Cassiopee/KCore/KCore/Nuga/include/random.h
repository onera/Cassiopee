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
//Authors : Sam Landier (sam.landier@onera.fr)

#include <random>

#ifndef NUGA_RANDOM_H
#define NUGA_RANDOM_H

#define NUGA_MERSENNE_TWISTER

namespace NUGA
{
#ifndef NUGA_MERSENNE_TWISTER
  struct random
  {
    inline void srand(int seed)
    {
      std::srand(seed);
    }

    inline unsigned int rand()
    {
      return std::rand();
    }
  };
#else
  struct random
  {
    inline void srand(int sd)
    {
      gen.seed(sd);
    }

    inline unsigned int rand()
    {
      return gen();
    }

    std::mt19937 gen;
};
#endif
}

#endif
