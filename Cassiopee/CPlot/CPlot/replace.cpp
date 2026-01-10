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

# include "kcore.h"
# include "cplot.h"
# include "Data.h"
#if defined(_WIN32) || defined(_WIN64)
#  include <winsock.h>
#endif

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   replace: replace a zone in plotter list.
   IN: array.
   IN: (nzs, nzu, oldType): no de la zone a remplacer
   IN: zoneName: nom de la zone (optionnel)
   IN: renderTag: renderTag (optionnel)
   Ne lance pas de render.
 */
//=============================================================================
PyObject* K_CPLOT::replace(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* l;
  PyObject* zoneNameObject;
  PyObject* renderTagObject;
  if (!PyArg_ParseTuple(args, "OOOO", &array, &l, 
                        &zoneNameObject, &renderTagObject)) return NULL;

  // Recuperation du container de donnees
  Data* d = Data::getInstance();

  // Recuperaton des no (nzs, nzu)
  if (PyTuple_Check(l) == false)
  {
    PyErr_SetString(PyExc_TypeError,
               "replace: arg must be a tuple.");
    return NULL;
  }
  E_Int sizel = PyTuple_Size(l);
  if (sizel != 3)
  {
    PyErr_SetString(PyExc_TypeError, 
               "replace: arg must be a tuple (nzs, nzu, oldType).");
    return NULL;
  }
  E_Int nzs = PyLong_AsLong(PyTuple_GetItem(l, 0));
  E_Int nzu = PyLong_AsLong(PyTuple_GetItem(l, 1));
  E_Int oldType = PyLong_AsLong(PyTuple_GetItem(l, 2));

  // Recuperation de l'array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, 
                                     ni, nj, nk, cn, eltType);

  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "replace: invalid array.");
    return NULL;
  }

  E_Int numberOfStructZones = d->_numberOfStructZones;
  E_Int numberOfUnstructZones = d->_numberOfUnstructZones;
  E_Int numberOfZones = d->_numberOfZones;

  // Element a remplacer
  E_Int nz = 0;
  if (oldType == 1) nz = nzs;
  else nz = numberOfStructZones+nzu;
  E_Int* replaced = new E_Int[1]; replaced[0] = nz;

  // Suppression 
  d->ptrState->syncGPURes();
  Zone* z = d->_zones[replaced[0]];
  z->freeGPURessources(false, false);

  // Recuperation des nouveaux noms de zones (eventuellement)
  char* zoneNameI=NULL;
  getStringFromPyObj(zoneNameObject, zoneNameI);

  // Recuperation des tags de render (eventuellement)
  char* zoneTagI=NULL;
  getStringFromPyObj(renderTagObject, zoneTagI);

  // Here we go
  E_Int posx, posy, posz;
  char zoneName[80];

  // Sauvegarde pointeurs
  Zone** zonesp = d->_zones;
  StructZone** szonesp = d->_szones;
  UnstructZone** uzonesp = d->_uzones;

  Zone* referenceZone = NULL;
  E_Int referenceNfield = -1;
  char** referenceVarNames = NULL;
  if (numberOfZones > 0) 
  {
    referenceZone = zonesp[0];
    referenceNfield = referenceZone->nfield;
    referenceVarNames = new char* [referenceNfield];
    for (E_Int i = 0; i < referenceNfield; i++) 
    {
      referenceVarNames[i] = new char [MAXSTRINGLENGTH];
      strcpy(referenceVarNames[i], referenceZone->varnames[i]);
    } 
  }

  // malloc nouveaux pointeurs (copie)
  Zone** zones = (Zone**)malloc(numberOfZones*sizeof(Zone*));
  for (E_Int i = 0; i < numberOfZones; i++) zones[i] = d->_zones[i];
  
  // Creations d'une nouvelle zone structuree
  if (res == 1)
  {
    if (zoneNameI != NULL) 
    { strcpy(zoneName, zoneNameI); delete [] zoneNameI; }
    else sprintf(zoneName, "S-Zone " SF_D_, nzs);
    
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);

    Zone* zz = 
      d->createStructZone(f, varString,
                          posx+1, posy+1, posz+1,
                          ni, nj, nk,
                          zoneName, zoneTagI, 
                          referenceNfield, referenceVarNames, 1);
    StructZone& z = (StructZone&)*zz;

    // Previous
    Zone* pz;
    if (oldType == 1) // structure
    {
      pz = d->_zones[nz];
      zones[nz] = zz;
      StructZone& zp = *(StructZone*)(pz);
      z.activePlane = zp.activePlane;
      z.iPlane = zp.iPlane;
      z.jPlane = zp.jPlane;
      z.kPlane = zp.kPlane;
      z.iLine = zp.iLine;
      z.jLine = zp.jLine;
      z.kLine = zp.kLine;
      z.blank = zp.blank;
      z.active = zp.active;
      z.selected = zp.selected;
    }
    else // non structure
    {
      pz = d->_zones[nz];
      insertAfterNz(zonesp, numberOfZones, zones, nzs, zz);
      deleteNz(zonesp, numberOfZones, zones, nz+1);
      UnstructZone& zp = *(UnstructZone*)(pz);
      z.blank = zp.blank;
      z.selected = zp.selected;
      z.active = zp.active;
    }
  }
  else // res=2
  {
    if (zoneNameI != NULL) { strcpy(zoneName, zoneNameI); delete [] zoneNameI; }
    else sprintf(zoneName, "U-Zone " SF_D_, nzu);
    
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);

    Zone* zz = 
      d->createUnstrZone(f, varString,
                         posx+1, posy+1, posz+1,
                         cn, eltType,
                         zoneName, zoneTagI, 
                         referenceNfield, referenceVarNames, 1);
    UnstructZone& z = (UnstructZone&)*zz;

    // Previous
    Zone* pz;
    if (oldType == 1) // structure
    {
      pz = d->_zones[nz];
      insertAfterNz(zonesp, numberOfZones, zones, numberOfStructZones+nzu, zz);
      deleteNz(zonesp, numberOfZones, zones, nzs);      
      StructZone& zp = *(StructZone*)(pz);
      z.blank = zp.blank;
      z.selected = zp.selected;
      z.active = zp.active;
    }
    else // non structure
    {
      pz = d->_zones[nz];
      zones[nz] = zz;
      UnstructZone& zp = *(UnstructZone*)(pz);
      z.blank = zp.blank;
      z.selected = zp.selected;
      z.active = zp.active;
    }
  }

  if (zoneTagI != NULL) delete [] zoneTagI;

  // Remise a jour des ptrs
  numberOfStructZones = 0;
  numberOfUnstructZones = 0;
  for (E_Int i = 0; i < numberOfZones; i++)
  {
    Zone* z = zones[i];
    if (dynamic_cast<StructZone*>(z) != NULL) numberOfStructZones++;
    else numberOfUnstructZones++;
  }
  StructZone** szones = (StructZone**)malloc(numberOfStructZones*
                                             sizeof(StructZone*));
  UnstructZone** uzones = (UnstructZone**)malloc(numberOfUnstructZones*
                                                 sizeof(UnstructZone*));
  E_Int p=0, q=0;
  for (E_Int i = 0; i < numberOfZones; i++)
  {
    Zone* z = zones[i];
    if (dynamic_cast<StructZone*>(z) != NULL) 
    { szones[p] = (StructZone*)z; p++; }
    else { uzones[q] = (UnstructZone*)z; q++; }
  }
  for (E_Int i = 0; i < numberOfStructZones; i++)
    zones[i] = szones[i];
  for (E_Int i = 0; i < numberOfUnstructZones; i++)
    zones[i+numberOfStructZones] = uzones[i];

  // update de la liste des deactivatedZones
  if (oldType != res && d->ptrState->deactivatedZones != NULL) // la numerotation change que si la zone change de type
  {
    int* old2new = new int [d->_numberOfZones];
    if (oldType == 1) // structure qui est maintenant a la fin en non-structure
    {
      for (E_Int i = 0; i < nzs; i++) old2new[i] = i;
      for (E_Int i = nzs+1; i < d->_numberOfZones; i++) old2new[i] = i-1;
      old2new[nzs] = d->_numberOfZones-1;
    }
    else if (oldType == 2) // non structure qui est maintenant a la fin des structures
    {
      for (E_Int i = 0; i < d->_numberOfStructZones; i++) old2new[i] = i;
      for (E_Int i = d->_numberOfStructZones; i < d->_numberOfZones; i++) old2new[i] = i+1;
      old2new[nzs+nzu] = d->_numberOfStructZones;
    }
    chain_int* c = d->ptrState->deactivatedZones;
    int oldn, newn;
    while (c != NULL)
    {
      oldn = c->value-1;
      newn = old2new[oldn];
      c->value = newn+1; c = c->next;
    }
    delete [] old2new;
    //d->ptrState->printDeactivatedZones();
    //d->ptrState->clearDeactivatedZones();
  }
  
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

  // Free previous zones
  delete zonesp[nz];
  free(szonesp); free(uzonesp); free(zonesp);
  
  // Free the input array
  RELEASESHAREDB(res, array, f, cn);

  for (E_Int i = 0; i < referenceNfield; i++) delete [] referenceVarNames[i];
  delete [] referenceVarNames;

  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/* Insert a zone after nz
   IN: zonesp: previous zone ptr (non modifie)
   IN/OUT: lzonen: nbre de zones dans zonesn
   IN/OUT: zonesn: new updated zone ptr
   IN: nz: nz d'insertion. 
*/
//=============================================================================
void insertAfterNz(Zone** zonesp, E_Int& lzonesn, Zone**& zonesn, E_Int nz, Zone* z)
{
  /*
  // locate nz in zonesp
  Zone* zz = zonesp[nz];

  // locate z in zonen
  int i = 0;
  for (i = 0; i < lzonesn; i++)
  {
    if (zonesn[i] == zz) break;
  }
  */

  lzonesn++;
  Zone** ntzones = (Zone**)malloc(lzonesn * sizeof(Zone*) );
  E_Int i = nz;
  for (E_Int j = 0; j < i; j++) ntzones[j] = zonesn[j];

  ntzones[i] = z;
  for (E_Int j = i+1; j < lzonesn; j++) ntzones[j] = zonesn[j-1];
  free(zonesn);
  zonesn = ntzones;
}

//=============================================================================
/* Delete a zone no nz
   IN: zonesp: previous zone ptr (non modifie)
   IN/OUT: lzonen: nbre de zones dans zonesn
   IN/OUT: zonesn: new updated zone ptr
   IN: nz: nz delete. 
*/
//=============================================================================
void deleteNz(Zone** zonesp, E_Int& lzonesn, Zone**& zonesn, E_Int nz)
{
  lzonesn--;
  Zone** ntzones = (Zone**)malloc(lzonesn * sizeof(Zone*) );
  for (E_Int j = 0; j < nz; j++) ntzones[j] = zonesn[j];
  for (E_Int j = nz; j < lzonesn; j++) ntzones[j] = zonesn[j+1];
  free(zonesn);
  zonesn = ntzones;
}
