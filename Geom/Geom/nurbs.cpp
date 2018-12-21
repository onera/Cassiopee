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
# include "geom.h"
# include <vector>
# include <stdio.h>

using namespace K_FLD;
using namespace std;

//===========================================================================
/* Nurbs */
//===========================================================================
PyObject* K_GEOM::nurbs(PyObject* self, PyObject* args)
{ 
  PyObject* Array;
  E_Int N, M;
  E_Int ordern, orderm;
  PyObject* ArrayW;
  double density;
  if (!PYPARSETUPLE(args,
                    "OOlllld", "OOiiiid",
                    "OOllllf", "OOiiiif",
                    &Array, &ArrayW, &ordern, &N, &orderm, &M, &density))
  {
      return NULL;
  }

  if (ordern < 1 )
  {
    PyErr_SetString(PyExc_ValueError,
                    "nurbs: order N must be >= 1.");
    return NULL;
  }
  if (orderm < 1 )
  {
    PyErr_SetString(PyExc_ValueError,
                    "nurbs: order M must be >= 1.");
    return NULL;
  }

  E_Int im, jm, km;

  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  E_Int imw, jmw, kmw;

  FldArrayF* fw; FldArrayI* cnw;
  char* varStringw; char* eltTypew;

  E_Int res = K_ARRAY::getFromArray(Array, varString, f, im, jm, km, 
                                    cn, eltType);
  if (res != 1)
  {
    if (res == 2 ) { delete f; delete cn; }
    PyErr_SetString(PyExc_TypeError,
                    "nurbs: (control points) input array not valid.");
    return NULL;
  }

  E_Int resw = K_ARRAY::getFromArray(ArrayW, varStringw, fw, imw, jmw, kmw, 
                                     cnw, eltTypew);
  if (resw != 1)
  {
    delete f;
    if (resw == 2) { delete fw; delete cnw; }
    PyErr_SetString(PyExc_TypeError,
                    "nurbs: (weights) input array not valid.");
    return NULL;
  }
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    delete f; delete fw;
    PyErr_SetString(PyExc_TypeError,
                    "nurbs : coordinates not found in control points array.");
    return NULL;
  }
  posx++; posy++; posz++;

  if (im < 2) 
  {
    delete f; delete fw;
    PyErr_SetString(PyExc_ValueError,
                    "nurbs: minimum 2 control points required.");
    return NULL;
  }

  if (im != imw || jm != jmw || km != kmw)
  {
    delete f; delete fw;
    PyErr_SetString(PyExc_ValueError,
                    "nurbs: number of control points differ from number of weights points.");
    return NULL;
  }

  if (im < ordern)
  {
    delete f; delete fw;
    PyErr_SetString(PyExc_ValueError,
                    "nurbs: number of control points must be greater than order.");
    return NULL;
  }
  
  if (jm > 1 && jm < orderm) // ij-array 
  {
    delete f; delete fw;
    PyErr_SetString(PyExc_ValueError,
                    "nurbs: number of control points must be greater than order.");
    return NULL;
  }
  
  if (km > 1)
  {
    delete f; delete fw;
    PyErr_SetString(PyExc_ValueError,
                    "nurbs: array must have one or two dimensions.");
    return NULL;
  }

  if (N < im || M < jm)
  {
    delete f; delete fw;
    PyErr_SetString(PyExc_ValueError,
                    "nurbs: N and M must be greater than the number of control points.");
    return NULL;
  }
  E_Int taille = fw->getSize();
  E_Int test = 0;

  for (E_Int i = 0 ; i < taille ; i++)
  {
    if (fw[0][i] < 0)
    {
      delete f; delete fw;
      PyErr_SetString(PyExc_ValueError,
                      "nurbs: all the weight must be positive.");
      return NULL;
    }
    
    if (K_FUNC::fEqualZero(fw[0][i]) == true) test += 1;
  }
  
  if (test == taille)
  {
    delete f; delete fw;
    PyErr_SetString(PyExc_ValueError,
                    "nurbs: at least one weight must be greater than 0.");
    return NULL;
  }
  /*Fin des Tests*/
  
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz); 
  E_Float* W = fw->begin();

  if (im != 1)
  {
    if (jm == 1) 
    { 
      K_FLD::FldArrayF PF;
      K_COMPGEOM::regularNurbs(im, ordern, N, density, xt, yt, zt, W, PF);
      delete f; delete fw;
      PyObject* tpl = K_ARRAY::buildArray(PF, "x,y,z", PF.getSize(), 1, 1);
      return tpl;
    }
    else 
    {
      K_FLD::FldArrayF PF;
      E_Int niout, njout;
      K_COMPGEOM::regularNurbs2D(im, jm, ordern, N, orderm, M, 
                                 density, xt, yt, zt, W, PF, niout, njout);
      delete f; delete fw;
      PyObject* tpl = K_ARRAY::buildArray(PF, "x,y,z", niout, njout, 1);
      return tpl;
    }
  }
  else
  {
    delete f; delete fw;
    PyErr_SetString(PyExc_TypeError,
                    "nurbs: control points array must be 1D, 2D.");
    return NULL;
  }
}
