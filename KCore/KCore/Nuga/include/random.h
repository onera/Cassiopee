/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

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