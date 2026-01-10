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

#pragma once
#include "NodeAssociator.h"
#include "Zipper.h"


class PostNodeAssociator : public NodeAssociator
{
public:
  typedef NodeAssociator parent_type;

public:
  ///
  PostNodeAssociator(bool nodal);

  ///
  virtual ~PostNodeAssociator(void);

  ///
  void make_pairs(const K_FLD::FloatArray& pos, const std::vector< K_FLD::IntArray* >& components,
                  std::vector<E_Int> &mates, K_FLD::IntArray& OneSurface);

private:

  ///
  void __setZipMates(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components,
                     std::vector<E_Int>& nmates);

private:
  Zipper* _zipper;
};


