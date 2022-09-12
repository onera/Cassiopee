/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __NUGA_MACROS_H__
#define __NUGA_MACROS_H__

#include <memory>
#include <cassert>

#define ALL(c)  (c).begin(), (c).end()

#define IS_IN(c,v) (c).find(v) != (c).end()

#define zSIGN(a, TOL) ( (a < -TOL)? -1 : (a > TOL) ? 1 : 0 )

#define STACK_ARRAY(T, n, name) std::unique_ptr<T[]> name(new T[n]);

#define NEIGHBOR(PHi, F2E, PGi) ( (F2E(0,PGi) == PHi) ? F2E(1,PGi) : F2E(0,PGi) )

#define ASSERT_IN_VECRANGE(arr, i) assert (((i)> -1) && ((i) < (int)arr.size()));

#define ASSERT_IN_DYNARANGE(arr, i) assert (((i)> -1) && ((i) < arr.cols()));

#endif
