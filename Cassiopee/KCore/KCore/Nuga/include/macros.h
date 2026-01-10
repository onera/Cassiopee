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
