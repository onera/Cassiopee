/*    
    Copyright 2013-2025 Onera.

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
#include "Nuga/include/DynArray.h"
#include "Nuga/include/defs.h"

class Zipper
{
public:
  enum eMateType {FREE = -9, OVERLAP = -1, UNKNO = IDX_NONE};

public:

  ///
  Zipper(bool reorient);

  ///
  virtual ~Zipper(void);

  ///
  virtual void setMates(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                        const K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates);

  /// Merges a nodal mesh. Does nothing otherwise.
  virtual void merge(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                     K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates){}

  ///
  static void zip(const K_FLD::IntArray& connectE2, const std::vector<E_Int>& nmates, K_FLD::IntArray& zipped);

protected:
  ///
  static void __setNodeMates(const K_FLD::FloatArray& vpos, const K_FLD::IntArray& vE2, const std::vector<E_Int>& realID,
                             std::vector<E_Int>& nmates);

private:
  ///
  static void __setBottoms(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectT3, const K_FLD::IntArray& connectE2, K_FLD::IntArray& bottomE2);

  ///
  static void __computeVirtualEdges(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectE2,
                                     const K_FLD::IntArray& bottomE2,
                                     K_FLD::FloatArray& vpos, K_FLD::IntArray& vE2,
                                     std::vector<E_Int>& realID);

/*
  ///
  static void __computeRadius(const K_FLD::FloatArray& vpos, const K_FLD::IntArray& vE2, const K_FLD::FloatArray& pos,
                              const std::vector<E_Int>& realID, std::vector<E_Float>& R2);
*/
protected:
  bool _reorient;

};
