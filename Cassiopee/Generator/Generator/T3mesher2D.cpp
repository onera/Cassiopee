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

// Mesh (or triangulate only) a 2D contour.

# include "generator.h"
# include "Nuga/include/T3Mesher.h"
using namespace std;

//=========================================================================
/* Maillage de type Delaunay a partir d'un ensemble de points */
//=========================================================================
PyObject* K_GENERATOR::T3mesher2D(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int triangulateOnly(0), metric_interp_type(0);
  E_Float grading(1.2);

  if (!PYPARSETUPLE_(args, O_ R_ II_, 
                    &array, &grading, &triangulateOnly, &metric_interp_type)) { return NULL; }

  // Check array
  E_Int ni, nj, nk;
  K_FLD::FloatArray* f;
  char* varString; char* eltType;
  K_FLD::IntArray* cn;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk,
                                    cn, eltType);
  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "T3mesher2D: invalid array.");
    return NULL;
  }
  
  if (res == 1) 
  {
    if (ni != 1 && nj != 1 && nk != 1)
    {
      delete f; 
      PyErr_SetString(PyExc_TypeError,
                      "T3mesher2D: array must define a plane.");
      return NULL;
    }
  }
  if (res == 2) 
  {
    if (strcmp(eltType, "TRI") != 0 && 
        strcmp(eltType, "QUAD") != 0 &&
        strcmp(eltType, "NODE") != 0 &&
        strcmp(eltType, "BAR") != 0)
    {
      delete f; delete cn;
      PyErr_SetString(PyExc_TypeError,
                      "T3mesher2D: array must define a plane.");
      return NULL;
    }
  }

  // Coordinates
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);

  if (posx == -1 || posy == -1)
  {
    delete f;
    if (res == 2) delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "T3mesher2D: can't find coordinates in array.");
    return NULL;
  }
 
  K_FLD::FloatArray& POS = *f;
  K_FLD::IntArray&   CONNECT = *cn;

  POS.resize(2, POS.cols());

  DELAUNAY::MeshData dT3(POS, CONNECT);

  // User metric
  E_Int posm11 = K_ARRAY::isNamePresent("m11", varString);
  E_Float zero = 0., sz = POS.cols();
  
  if (posm11 != -1) // Isotropic metric (so far).
  {  
    dT3.metrics.resize(1, POS.cols(), &zero);

    for (E_Int j = 0; j < sz; ++j)
      dT3.metrics(0,j) = POS(posm11-1, j);
  }

  E_Int posm12 = K_ARRAY::isNamePresent("m12", varString);
  E_Int posm22 = K_ARRAY::isNamePresent("m22", varString);

  if ((posm12 != -1) && (posm22 != -1)) // Anisotropic metric
  {
    dT3.metrics.resize(3, POS.cols(), &zero);
      
    for (E_Int j = 0; j < sz; ++j)
    {
      dT3.metrics(1,j) = POS(posm12-1, j);
      dT3.metrics(2,j) = POS(posm22-1, j);
    }
  }

  // Mesher mode.
  DELAUNAY::MesherMode mode;
  if (triangulateOnly)
    mode.mesh_mode = DELAUNAY::MesherMode::TRIANGULATION_MODE;
  else if (grading != 1.) mode.symmetrize = true;
  
  mode.growth_ratio = grading; // Growth ratio
  //std::cout << "gr : " << mode.growth_ratio << std::endl;

  if (metric_interp_type != 0)
    mode.metric_interpol_type = mode.GEOMETRIC;

  
  E_Int err = 0;
  if (dT3.metrics.rows() < 3) // Isotropic
  {
    DELAUNAY::T3Mesher<E_Float> mesher(mode);
    err = mesher.run(dT3);
  }
  else // Anisotropic
  {
    DELAUNAY::T3Mesher<DELAUNAY::Aniso2D> mesher(mode);
    err = mesher.run(dT3);
  }

  PyObject* tpl = NULL;
  if (!err)
  {
    dT3.pos->resize(3, dT3.pos->cols(), &zero);
    tpl = K_ARRAY::buildArray(*dT3.pos, varString, dT3.connectM, 
                                      -1, "TRI", false);
  }
  delete f; delete cn;
  return tpl;
}

//========================  Generator/T3mesher2D.cpp ========================
