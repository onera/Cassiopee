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

# include "connector.h"
# include <cassert>
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Init data in parallel openmp to improve data placement on numa machine */
//=============================================================================
PyObject* K_CONNECTOR::initNuma(PyObject* self, PyObject* args)
{
  PyObject* sourceArray;  PyObject* targetArray; E_Int ideb;  E_Int size_bc; E_Int vartype; E_Float val;

  if (!PYPARSETUPLE_(args, OO_ III_ R_, &sourceArray, &targetArray, &ideb, &size_bc, &vartype, &val ))
  {
    return NULL;
  }

  E_Int iDeb = E_Int(ideb);
  E_Int varType = E_Int(vartype);
  E_Int size_BC = E_Int(size_bc);
 
  E_Int option = 0;
  if (sourceArray == Py_None) option =1;
 
  if (varType == 0)
  { // gestion des tableau Float
    FldArrayF* source; FldArrayF* target;
    E_Float* iptsource = NULL;
    if (option == 0) 
    {
      K_NUMPY::getFromNumpyArray(sourceArray, source); iptsource = source->begin();
    }
    K_NUMPY::getFromNumpyArray(targetArray, target); E_Float* ipttarget = target->begin();

    # pragma omp parallel default(shared)
    {
     #ifdef _OPENMP
     E_Int  ithread           = omp_get_thread_num()+1;
     E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
     #else
     E_Int ithread = 1;
     E_Int Nbre_thread_actif = 1;
     #endif
     E_Int pt_deb, pt_fin;

     // Calcul du nombre de champs a traiter par chaque thread
     E_Int chunk   =  size_BC/Nbre_thread_actif;
     E_Int r       =  size_BC - chunk*Nbre_thread_actif;
     // pts traitees par thread
     if (ithread <= r)
          { pt_deb = iDeb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); }
     else { pt_deb = iDeb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk; } 

     if (option == 0)
     {   for ( E_Int i = pt_deb; i < pt_fin; i++)
         { ipttarget[ i ] = iptsource[i-iDeb]; }
     }
     else
     {   for ( E_Int i = pt_deb; i < pt_fin; i++)
         { ipttarget[ i ] = val; }
     }
    }// omp
    if (option == 0) RELEASESHAREDN( sourceArray  , source  );
    RELEASESHAREDN( targetArray  , target  );
  }
  else
  {// gestion des tableau Int
    FldArrayI* source = NULL; FldArrayI* target;
    E_Int* iptsource = NULL;
    if (option == 0) 
    {
      K_NUMPY::getFromNumpyArray(sourceArray, source); 
      assert(source != NULL);
      iptsource = source->begin();
    }
    K_NUMPY::getFromNumpyArray(targetArray, target); E_Int* ipttarget = target->begin();

    # pragma omp parallel default(shared)
    {
     #ifdef _OPENMP
     E_Int  ithread           = omp_get_thread_num()+1;
     E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
     #else
     E_Int ithread = 1;
     E_Int Nbre_thread_actif = 1;
     #endif
     E_Int pt_deb, pt_fin;

     // Calcul du nombre de champs a traiter par chaque thread
     E_Int chunk   =  size_BC/Nbre_thread_actif;
     E_Int r       =  size_BC - chunk*Nbre_thread_actif;
     // pts traitees par thread
     if (ithread <= r)
          { pt_deb = iDeb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); }
     else { pt_deb = iDeb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk; } 

     if (option == 0)
     {   for ( E_Int i = pt_deb; i < pt_fin; i++)
         { ipttarget[ i ] = iptsource[i-iDeb]; }
     }
     else
     {   E_Int Val = E_Int(val);
         for ( E_Int i = pt_deb; i < pt_fin; i++)
         { ipttarget[ i ] = Val; }
     }
    }// omp
    if (option == 0) RELEASESHAREDN( sourceArray  , source  );
    RELEASESHAREDN( targetArray  , target  );
  }

  Py_INCREF(Py_None);
  return Py_None;
}
