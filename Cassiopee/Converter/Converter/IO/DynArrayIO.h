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

#ifndef __DYNARRAYIO_H__
#define	__DYNARRAYIO_H__

#include <vector>
# include "Nuga/include/DynArray.h"
#include <string>

namespace K_CONVERTER
{

class DynArrayIO {
  
public: // Common interface with meshIO

  static std::string rdir, wdir;

  // Reads one single zone file
  template <typename Coordinate_t, typename Cnt_t>
  static E_Int read (const char* filename, Coordinate_t& crd, Cnt_t& cnt);

  // Writes a single zone into a file (element type must be specified to distinguish QUAD/TETRA, it is not ambiguous otherwise.)
  static E_Int write(const char* filename, const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt,
	                 const char* elt_type=0, const std::vector<bool>* mask=0, const std::vector<E_Int>* colors=0);

  // Writes a point cloud into a file.
  static E_Int write(const char* filename, const K_FLD::FloatArray& coord);

public : // Interface specific to DynArrayIO for reading/writing multiple zone files at once.
  
  static E_Int read
  (const char* fileName, std::vector<K_FLD::FloatArray>& coords, std::vector<K_FLD::IntArray>& connects);
  
  static E_Int read
  (const char* fileName, std::vector<K_FLD::FldArrayF>& coords, std::vector<K_FLD::FldArrayI>& connects);

  static E_Int write
  (const char* fileName, const std::vector<K_FLD::FloatArray>& coords,
   const std::vector<K_FLD::IntArray>& connects,
   const std::vector<std::string>& elt_type);
  
  static E_Int write
  (const char* fileName, const K_FLD::FldArrayF& coords, const K_FLD::FldArrayI& connects,
   const char* elt_type=0, const std::vector<bool>* mask=0, const std::vector<E_Int>* colors=0);
  
private:
  static E_Int getElementTypeId(const char* eltType);
  static const char* getElementType(E_Int id);
  
  static const char* get_fmt(const char* fname);
  
private:
  DynArrayIO();
  DynArrayIO(const DynArrayIO& orig);
  virtual ~DynArrayIO();
};

// Reads one single zone file
template <typename Crd_t, typename Cnt_t>
E_Int DynArrayIO::read(const char* filename, Crd_t& crd, Cnt_t& cnt)
{
  std::vector<Crd_t> coords;
  std::vector<Cnt_t> connects;
  E_Int ret = read(filename, coords, connects);
  if (ret) return ret;
  crd = coords[0];
  cnt = connects[0];
  return 0;
}

}

#endif	/* __DYNARRAYIO_H__ */

