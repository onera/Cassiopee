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

#ifndef __GENERATOR_GAP_FIXER_H__
#define __GENERATOR_GAP_FIXER_H__

#include "Def/DefContainers.h"

class GapFixer
{
public:
  typedef K_CONT_DEF::int_vector_type int_vector_type;
public:
  ///
  GapFixer();

  ///
  ~GapFixer(void);

  /// Main function to fix a gap.
  static E_Int run(const K_FLD::FloatArray& posC,
                   E_Int nj,
                   const K_FLD::FloatArray& posB0,
                   const K_FLD::IntArray& connectB0,
                   K_FLD::FloatArray& posG,
                   K_FLD::IntArray& connectG,
                   E_Bool refine,
                   const K_FLD::FloatArray* coordHP = 0); // hard points

private:
  /// Converts the parametrization to the quads global one.
  static void __convertToGlobalQInterp(const K_FLD::FloatArray& posB0, E_Int nj,
                                       const K_CONT_DEF::int_vector_type& cell_indices,
                                       K_FLD::FloatArray& posUV);

};

#endif
