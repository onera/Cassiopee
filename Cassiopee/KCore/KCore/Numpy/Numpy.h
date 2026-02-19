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
// Interface to numpy arrays

#ifndef _KCORE_NUMPY_H_
#define _KCORE_NUMPY_H_
#include "kPython.h"
#include "Def/DefTypes.h"
#include "Def/DefCplusPlusConst.h"
#include "Numpy/importNumpy.h"
#include <numpy/arrayobject.h>
#include "Fld/FldArray.h"
#include <vector>

#define FldArrayF K_FLD::FldArrayF
#define FldArrayI K_FLD::FldArrayI

// Release shared numpy array
#define RELEASESHAREDN(array, f) {delete f; Py_DECREF(array);}

namespace K_NUMPY
{
  /* Extrait les donnees d'un numpy array
     IN: o: numpy array
     OUT: f: FldArray alloue et rempli.
     Le tableau pointe sur le numpy et la reference sur
     o est incrementee.
     For FldArray, use RELEASESHAREDN to deallocate.
     For E_Int/E_Float, use Py_DECREF to deallocate.
     Retourne 0 (FAIL), 1 (SUCCESS) */
  E_Int getFromNumpyArray(PyObject* o, FldArrayI*& f);
  E_Int getFromNumpyArray(PyObject* o, FldArrayF*& f);
  E_Int getFromNumpyArray(PyObject* o, E_Int*& f, E_Int& size, E_Int& nfld);
  E_Int getFromNumpyArray(PyObject* o, E_Int*& f, E_Int& size);
  E_Int getFromNumpyArray(PyObject* o, E_Float*& f, E_Int& size, E_Int& nfld); 
  E_Int getFromNumpyArray(PyObject* o, E_Float*& f, E_Int& size);

  E_Int getFromNumpyArray(PyObject* o, FldArrayI*& f, E_Bool shared); //#OBSOLETE
  E_Int getFromNumpyArray(PyObject* o, FldArrayF*& f, E_Bool shared); //#OBSOLETE
  E_Int getFromNumpyArray(PyObject* o, E_Int*& f, E_Int& size, E_Int& nfld, E_Bool shared); //#OBSOLETE
  //E_Int getFromNumpyArray(PyObject* o, E_Int*& f, E_Int& size, E_Bool shared); //#OBSOLETE
  E_Int getFromNumpyArray(PyObject* o, E_Float*& f, E_Int& size, E_Int& nfld, E_Bool shared); //#OBSOLETE 
  //E_Int getFromNumpyArray(PyObject* o, E_Float*& f, E_Int& size, E_Bool shared); //#OBSOLETE

  // identical to getFromNumpy but return a flat array in case of (1,nb) arrays
  E_Int getFromPointList(PyObject* o, FldArrayI*& f);
  E_Int getFromPointList(PyObject* o, E_Int*& f, E_Int& size, E_Int& nfld);

  /* Construit un numpy array a partir d'un FldArray (copie) 
     IN: field: fld array 
     IN: fortran: si 1, le numpy a le flag fortran
     OUT: numpy array retourne */
  PyObject* buildNumpyArray(FldArrayF& field, E_Int fortran=0);
  PyObject* buildNumpyArray(FldArrayI& field, E_Int fortran=0);
  PyObject* buildNumpyArray(E_Float* field, E_Int size, E_Int nfld, E_Int fortran=0);
  PyObject* buildNumpyArray(E_Int* field, E_Int size, E_Int nfld, E_Int fortran=0);
  
  /* Construit un numpy array vide.
     Si type=0 (E_Float), si type=1 (E_Int) */
  PyObject* buildNumpyArray(E_Int size, E_Int nfld, E_Int type, E_Int fortran=0);

  /* Recupere le ptr sur le numpy */
  E_Float* getNumpyPtrF(PyObject* o);
  E_Int* getNumpyPtrI(PyObject* o);
}

#undef FldArrayF
#undef FldArrayI
#endif
