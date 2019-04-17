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

# include "kcore.h"
# include "cplot.h"
# include "Data.h"

#if defined(_WIN32) || defined(_WIN64)
#  include <winsock.h>
#endif

using namespace K_FLD;

//=============================================================================
/* 
   delete: suppress zones from plotter list.
   IN: l: liste du no des zones a supprimer ou liste du nom des zones.
   IN: si l = 'all', supprime toutes les zones
   Ne lance pas de render.
 */
//=============================================================================
PyObject* K_CPLOT::deletez(PyObject* self, PyObject* args)
{
  PyObject* l;
  if (!PyArg_ParseTuple(args, "O", &l)) return NULL;
  
  Data* d = Data::getInstance();
  FldArrayI deleted;
  
  E_Int size;
  if (PyList_Check(l) == true) // list
  {
    size = PyList_Size(l);
    deleted.malloc(size);
    PyObject* tpl;
    for (E_Int i = 0; i < size; i++)
    {
      tpl = PyList_GetItem(l, i);
      if (PyLong_Check(tpl) || PyInt_Check(tpl))
        deleted[i] = PyLong_AsLong(tpl);
      else if (PyString_Check(tpl))
      {
        char* name = PyString_AsString(tpl);
        for (E_Int j = 0; j < d->_numberOfZones; j++)
        {
          if (strcmp(d->_zones[j]->zoneName, name) == 0) 
          { deleted[i] = j; break; }
        }
      }
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(tpl))
      {
        char* name = PyBytes_AsString(PyUnicode_AsUTF8String(tpl));
        for (E_Int j = 0; j < d->_numberOfZones; j++)
        {
          if (strcmp(d->_zones[j]->zoneName, name) == 0) 
          { deleted[i] = j; break; }
        } 
      }
#endif
      else
      {
        PyErr_SetString(PyExc_TypeError, 
                        "delete: arg must be a list of ints or strings.");
        return NULL;
      }
    }
  }
  else if (PyString_Check(l)) // string -> delete all
  {
    size = d->_numberOfZones;
    deleted.malloc(size);
    for (E_Int i = 0; i < size ; i++) deleted[i] = i;
  }
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(l))
  {
    size = d->_numberOfZones;
    deleted.malloc(size);
    for (E_Int i = 0; i < size ; i++) deleted[i] = i;   
  }
#endif
  else
  {
    PyErr_SetString(PyExc_TypeError, 
               "delete: arg must be a list.");
    return NULL;
  }

  // Suppression DL
  d->ptrState->syncGPURes();
  for (int i = 0; i < size; i++)
  {
    Zone* z = d->_zones[deleted[i]];
    z->freeGPURessources( false, false );
  }

  // Suppression de la zone
  E_Int numberOfStructZones = d->_numberOfStructZones;
  E_Int numberOfUnstructZones = d->_numberOfUnstructZones;
  E_Int numberOfZones = d->_numberOfZones-size;
  if (numberOfZones < 0) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "delete: invalid zones in list.");
    return NULL;
  }
  for (E_Int i = 0; i < size; i++)
  {
    if (deleted[i] < d->_numberOfStructZones) numberOfStructZones--;
    else numberOfUnstructZones--;
  }
  if (numberOfStructZones < 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "delete: invalid zones in list.");
    return NULL;
  }
  if (numberOfUnstructZones < 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "delete: invalid zones in list.");
    return NULL;
  }

  // Sauvegarde pointeurs
  Zone** zonesp = d->_zones;
  StructZone** szonesp = d->_szones;
  UnstructZone** uzonesp = d->_uzones;

  // malloc pointeurs
  Zone** zones = (Zone**)malloc(numberOfZones*sizeof(Zone*));
  StructZone** szones = (StructZone**)malloc(numberOfStructZones*
                                             sizeof(StructZone*));
  UnstructZone** uzones = (UnstructZone**)malloc(numberOfUnstructZones*
                                                 sizeof(UnstructZone*));
  int mustDel;
  int c = 0;
  for (E_Int i = 0; i < d->_numberOfStructZones; i++)
  {
    mustDel = 0;
    for (E_Int j = 0; j < size; j++)
      if (i == deleted[j]) {mustDel = 1; break;}
    if (mustDel == 0) { szones[c] = d->_szones[i]; c++; }
  }
  c = 0;
  for (E_Int i = d->_numberOfStructZones;
       i < d->_numberOfStructZones+d->_numberOfUnstructZones; i++)
  {
    mustDel = 0;
    for (E_Int j = 0; j < size; j++)
      if (i == deleted[j]) {mustDel = 1; break;}
    if (mustDel == 0) { uzones[c] = d->_uzones[i-d->_numberOfStructZones]; c++; }
  }
  for (E_Int i = 0; i < numberOfStructZones; i++)
    zones[i] = szones[i];
  for (E_Int i = 0; i < numberOfUnstructZones; i++)
    zones[i+numberOfStructZones] = uzones[i];

  // Switch - Dangerous zone protegee par _state.lock
  if (d->ptrState->selectedZone >= numberOfZones)
    d->ptrState->selectedZone = 0; // RAZ selected zone
  d->ptrState->syncDisplay();
  d->_zones = zones;
  d->_szones = szones;
  d->_uzones = uzones;
  d->_numberOfStructZones = numberOfStructZones;
  d->_numberOfUnstructZones = numberOfUnstructZones;
  d->_numberOfZones = numberOfZones;

  // Mise a jour des min-max globaux
  globMinMax(d->_zones, d->_numberOfZones, 
             d->xmin, d->xmax, d->ymin, d->ymax, d->zmin, d->zmax, 
             d->epsup, d->epsstrafe, d->dmoy);
  globFMinMax(d->_zones, d->_numberOfZones, d->minf, d->maxf);

  // Free
  for (E_Int i = 0; i < size; i++) delete zonesp[deleted[i]];
  free(szonesp); free(uzonesp); free(zonesp);
  return Py_BuildValue("i", KSUCCESS);
}
