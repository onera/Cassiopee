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


# include <string>
# include <sstream> 
# include "intersector.h"
# include "Fld/ngon_t.hxx"
# include "Nuga/Delaunay/Triangulator.h"
# include "Nuga/Boolean/SelfX.h"

//#include <iostream>

using namespace std;
using namespace NUGA;


//=============================================================================
/* Detect any self intersection in a NGON volume mesh */
//=============================================================================
PyObject* K_INTERSECTOR::selfX(PyObject* self, PyObject* args)
{
  PyObject *arr;
  //E_Float vmin(0.), vratio(1000.);

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);
  
  std::vector<E_Int> xlist;
  NUGA::selfX<DELAUNAY::Triangulator>(crd, ngi, xlist);
  
  std::vector<bool> keep(ngi.PHs.size(), false);
  
  for (size_t i=0; i < xlist.size(); ++i)
    keep[xlist[i]]=true;
  
  ngon_t<cnt_t> ngo;
  std::vector<E_Int> npgids;
  ngi.select_phs(ngi, keep, npgids, ngo);

  PyObject* tpl = 0; 
  if (ngo.PHs.size())
  {
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
  
    tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);
  }
  //fixme : return itself to shutdown the python layer in case of emptyness
  else tpl = K_ARRAY::buildArray(crd, varString, cnt, -1, eltType, false);
  
  delete f; delete cn;
  return tpl;
}

//=======================  Intersector/PolyMeshTools/selfX.cpp ====================
