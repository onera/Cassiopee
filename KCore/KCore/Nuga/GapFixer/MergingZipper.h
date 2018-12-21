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

#pragma once
#include "Zipper.h"

class MergingZipper : public Zipper
{

public:

  ///
  MergingZipper(bool reorient);

  ///
  virtual ~MergingZipper(void);

  ///
  virtual void setMates(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                        const K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates);

  /// Merges a nodal mesh. Does nothing otherwise.
  virtual void merge(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                     K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates);

private:
  E_Float __computeTolerance(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                             const std::vector<E_Int>& nodes);

  E_Int __removeDuplicates(K_FLD::IntArray& connect);

};
