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
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef IMPORT_OCC_CAD_WRAPPER_H
#define	IMPORT_OCC_CAD_WRAPPER_H

#include "Def/DefTypes.h"
#include <vector>
# include "Nuga/include/DynArray.h"

namespace K_OCC
{
  
class import_OCC_CAD_wrapper {
public:
  static E_Int import_cad
  (
    const char* fname, const char*format, 
    std::vector<K_FLD::FloatArray> & crds,
    std::vector<K_FLD::IntArray>& connectMs,
    E_Float h=0., E_Float chordal_err=0.,  E_Float gr = 0. /*growth ratio*/,
    bool aniso = false, bool do_join = true
  );

private:
  import_OCC_CAD_wrapper();
  import_OCC_CAD_wrapper(const import_OCC_CAD_wrapper& orig);
  virtual ~import_OCC_CAD_wrapper();

};
}

#endif	/* IMPORT_OCC_CAD_WRAPPER_H */

