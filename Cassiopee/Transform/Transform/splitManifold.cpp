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

# include "transform.h"
# include "Nuga/include/MeshTool.h"
# include "Nuga/include/EltAlgo.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   splitTRI: 
   Splits an unstructured mesh (currently only TRI) into several manifold 
   pieces.
*/
//=============================================================================
PyObject* K_TRANSFORM::splitManifold(PyObject* self, PyObject* args)
{
  PyObject *array;
  if (!PYPARSETUPLE_(args, O_, &array))
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FloatArray* f; IntArray* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType); 

  if (res == 1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "splitManifold: can not be used on a structured array.");
    return NULL;
  }
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitManifold: unknown type of array.");
    return NULL;
  }
  
  if ((strcmp(eltType, "TRI") != 0) && (strcmp(eltType, "BAR") != 0))
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "splitManifold: currently only supported for TRI and BAR arrays.");
    return NULL;
  }
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
   
  if (posx == -1 || posy == -1)
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "splitManifold: can't find coordinates in array.");
    return NULL;
  }
  
  K_FLD::IntArray neighbors;
  std::vector<E_Int> colors;
  E_Int NB_NODES=3;
  std::vector<E_Int> dupIds;
  NUGA::MeshTool::removeDuplicated(*cn, dupIds, false);
  
  if (strcmp(eltType, "TRI") == 0)
  {
    NUGA::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(*cn, neighbors, false);
    NUGA::EltAlgo<K_MESH::Triangle>::coloring (neighbors, colors);
  }
  else if (strcmp(eltType, "BAR") == 0)
  {
    NUGA::EltAlgo<K_MESH::Edge>::getManifoldNeighbours(*cn, neighbors, false);
    NUGA::EltAlgo<K_MESH::Edge>::coloring(neighbors, colors);
    NB_NODES = 2;
  }
  
  size_t nb_bits = *std::max_element(colors.begin(), colors.end())+1;
  
  // Formation des array de sortie
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  
  if (nb_bits == 1)
  {
    tpl = K_ARRAY::buildArray(*f, varString, *cn, -1, eltType);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  else
  {
    K_FLD::FloatArray* fs = new K_FLD::FloatArray[nb_bits];
    K_FLD::IntArray*   cs = new K_FLD::IntArray[nb_bits];
    K_FLD::IntArray::const_iterator pS;
    
    for (size_t i = 0; i < colors.size(); ++i)
    {
      pS = cn->col(i);
      cs[colors[i]].pushBack(pS, pS+NB_NODES);
    }
    
    std::vector<E_Int> nids;
    for (size_t i=0; i < nb_bits; ++i)
    {
      fs[i] = *f;
      NUGA::MeshTool::compact_to_mesh(fs[i], cs[i], nids);
      tpl = K_ARRAY::buildArray(fs[i], varString, cs[i], -1, eltType);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    delete [] fs; delete [] cs;
  }
  
  delete f; delete cn;
  return l;
}
