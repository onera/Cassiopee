/*    
    Copyright 2013-2018 Onera.

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

#ifndef __DELAUNAY_INTERPOLATOR_H__
#define __DELAUNAY_INTERPOLATOR_H__

#include "Def/DefTypes.h"

template <typename T>
class Interpolator
{
public:

  Interpolator(void)
  {
  }

  virtual ~Interpolator(void)
  {
  }

  virtual T interpolate(const T& H0, const T& H1, E_Float u) const = 0;
};

template <typename T>
class LinearInterpolator : public Interpolator<T>
{
public:
  LinearInterpolator(void)
  {
  }

  ~LinearInterpolator(void)
  {
  }

  inline T interpolate(const T& H0, const T& H1, E_Float u) const {return H0*(1-u) + H1*u;}

};

#endif
