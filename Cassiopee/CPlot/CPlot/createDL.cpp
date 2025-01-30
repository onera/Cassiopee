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

# include "DataDL.h"
# include "CPlotStateDL.h"
# include "ZoneImplDL.h"

//=============================================================================
/*
  Create Display Lists (solid et mesh si USEDLMESH == 1)
  Cette fonction cree une seule display list pour la premiere zone
  telle que:
  - sa display list n'a pas encore ete calculee
  - _DLuse = 1
  Pour creer les DL de toutes les zones, il faut donc l'appeler
  plusieurs fois.
*/
//=============================================================================
void DataDL::createGPURes()
{
  E_Int zone, zonet;
  ptrState->lockDisplay();
  // - Solid -
  zone = 0;
  while (zone < _numberOfStructZones)
  {
    StructZone* z    = _szones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    if (zImpl._GPUResUse == 1 && zImpl._DLsolid == 0)
    { createGPUSSolidZone(z, zone); goto end; }
    zone++;
  }

  zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    UnstructZone* z = _uzones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    if (zImpl._GPUResUse == 1 && zImpl._DLsolid == 0)
    {
      zonet = zone + _numberOfStructZones;
      if ( not z->_is_high_order)
      { createGPUUSolidZone(z, zone, zonet); goto end; }
      else
      {
        createGPUUSolidHOZone(z, zone, zonet); goto end;
      }
    }
    zone++;
  }

  // - Mesh -
#if (USEDLMESH == 1)
  zone = 0;
  while (zone < _numberOfStructZones)
  {
    StructZone* z = _szones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    if (zImpl._GPUResUse == 1 && zImpl._DLmesh == 0)
    { createGPUSMeshZone(z, zone); goto end; }
    zone++;
  }

  zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    UnstructZone* z = _uzones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    if (zImpl._GPUResUse == 1 && zImpl._DLmesh == 0)
    {
      zonet = zone + _numberOfStructZones;
      { createGPUUMeshZone(z, zone, zonet); goto end; }
    }
    zone++;
  }
#else
  zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    UnstructZone* z = _uzones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    if (zImpl._GPUResUse == 1 && zImpl._DLmesh == 0 && z->_is_high_order == true )
    {
      zonet = zone + _numberOfStructZones;
      { createGPUUMeshZone(z, zone, zonet); goto end; }
    }
    zone++;
  }
#endif
 end: ;
 ptrState->unlockDisplay();
}

//=============================================================================
// Cree les DL pour le champ scalaire nofield
//=============================================================================
void DataDL::createIsoGPURes(E_Int nofield)
{
  E_Int zone = 0;
  E_Int zonet;
  while (zone < _numberOfStructZones)
  {
    StructZone* z = _szones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    
    if (zImpl._GPUResUse == 1)
    {
      if (zImpl._DLisoField != nofield)
      {
        if (zImpl._DLiso != 0) 
        { glDeleteLists(zImpl._DLiso, 1); zImpl._DLiso = 0; }
        createGPUSIsoSolidZone(z, zone, nofield);
        zImpl._DLisoField = nofield;
        zImpl._DLisoField2 = -1;
        zImpl._DLisoField3 = -1;
      }
    }
    zone++;
  }
  zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    UnstructZone* z = _uzones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    if (zImpl._GPUResUse == 1)
    {
      zonet = zone + _numberOfStructZones;
      if (zImpl._DLisoField != nofield) 
      {
        if (zImpl._DLiso != 0) 
        { glDeleteLists(zImpl._DLiso, 1); zImpl._DLiso = 0; }
        createGPUUIsoSolidZone(z, zone, zonet, nofield);
        zImpl._DLisoField = nofield;
        zImpl._DLisoField2 = -1;
        zImpl._DLisoField3 = -1;
      }
    }
    zone++;
  }
}

//=============================================================================
// Cree les DL pour le mode render et pour les zones le necessitant
//=============================================================================
void DataDL::createIsoGPUResForRender()
{
  E_Int zone = 0;
  E_Int zonet;
  while (zone < _numberOfStructZones)
  {
    StructZone* z = _szones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    int color = (int)(z->colorR);
    if (zImpl._GPUResUse == 1 && color < -1.5)
    {
      int nofield = (int)(-color-2);
      if (zImpl._DLisoField != nofield)
      {
        if (zImpl._DLiso != 0) 
        { glDeleteLists(zImpl._DLiso, 1); zImpl._DLiso = 0; }
        createGPUSIsoSolidZone(z, zone, nofield);
        zImpl._DLisoField = nofield;
        zImpl._DLisoField2 = -1;
        zImpl._DLisoField3 = -1;
      }
    }
    zone++;
  }
  zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    UnstructZone* z = _uzones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    int color = (int)(z->colorR);
    if (zImpl._GPUResUse == 1 && color < -1.5)
    {
      int nofield = (int)(-color-2);
      zonet = zone + _numberOfStructZones;
      if (zImpl._DLisoField != nofield) 
      {
        if (zImpl._DLiso != 0) 
        { glDeleteLists(zImpl._DLiso, 1); zImpl._DLiso = 0; }
        createGPUUIsoSolidZone(z, zone, zonet, nofield);
        zImpl._DLisoField = nofield;
        zImpl._DLisoField2 = -1;
        zImpl._DLisoField3 = -1;
      }
    }
    zone++;
  }
}

//=============================================================================
// Cree les DL pour le champ vector nofield1, nofield2, nofield3
//=============================================================================
void DataDL::createIsoGPURes(E_Int nofield1, E_Int nofield2, E_Int nofield3)
{
  E_Int zone = 0;
  E_Int zonet;
  while (zone < _numberOfStructZones)
  {
    StructZone* z = _szones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    if (zImpl._GPUResUse == 1)
    {
      if (zImpl._DLisoField != nofield1 || zImpl._DLisoField2 != nofield2 
          || zImpl._DLisoField3 != nofield3)
      {
        if (zImpl._DLiso != 0) glDeleteLists(zImpl._DLiso, 1);
        createGPUSIsoSolidZone(z, zone, nofield1, nofield2, nofield3);
        zImpl._DLisoField = nofield1;
        zImpl._DLisoField2 = nofield2;
        zImpl._DLisoField3 = nofield3;
      }
    }
    zone++;
  }
  zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    UnstructZone* z = _uzones[zone];
    ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
    if (zImpl._GPUResUse == 1)
    {
      zonet = zone + _numberOfStructZones;
      if (zImpl._DLisoField != nofield1 || zImpl._DLisoField2 != nofield2
          || zImpl._DLisoField3 != nofield3) 
      {
        if (zImpl._DLiso != 0) glDeleteLists(zImpl._DLiso, 1);
        createGPUUIsoSolidZone(z, zone, zonet, nofield1, nofield2, nofield3);
        zImpl._DLisoField = nofield1;
        zImpl._DLisoField2 = nofield2;
        zImpl._DLisoField3 = nofield3;
      }
    }
    zone++;
  }
}
//=============================================================================
// Free the DL
// mode=0 -> DLmesh
// mode=1 -> DLsolid
// mode=-1 -> DLmesh et DLsolid
// start: debut des zones nettoyees
// end: fin des zones nettoyees (numerotation dans zone)
// permanent=1: met useDL de la zone a 0, la DL ne sera plus calculee
//=============================================================================
void DataDL::freeGPURes(int mode, int start, int end, int permanent)
{
  E_Int i;
  if (mode == 0 || mode == -1)
  {
    for (i = start; i <= end; i++)
    { 
      Zone* z = _zones[i];
      ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
      if (zImpl._DLmesh != 0) {glDeleteLists(zImpl._DLmesh, 1); zImpl._DLmesh = 0;}
      if (permanent == 1) zImpl._GPUResUse = 0;
    }
  }
  if (mode == 1 || mode == -1)
  {
    for (i = start; i <= end; i++)
    { 
      Zone* z = _zones[i];
      ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
      if (zImpl._DLsolid != 0) {glDeleteLists(zImpl._DLsolid, 1); zImpl._DLsolid = 0;}
      if (permanent == 1) zImpl._GPUResUse = 0;
    }
  }
}

//=============================================================================
// Delete les DL d'une liste de zone definis par ptr. 
// Cette fonction delete ptr
//=============================================================================
void DataDL::freeGPURes(int mode, int size, int* ptr, int permanent)
{
  E_Int i;
  if (mode == 0 || mode == -1)
  {
    for (i = 0; i < size; i++)
    {
      if (ptr[i] >= 0 && ptr[i] < _numberOfZones)
      {
        Zone* z = _zones[ptr[i]];
        ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
        if (zImpl._DLmesh != 0) {glDeleteLists(zImpl._DLmesh, 1); zImpl._DLmesh = 0;}
        if (permanent == 1) zImpl._GPUResUse = 0;
      }
    }
  }
  if (mode == 1 || mode == -1)
  {
    for (i = 0; i < size; i++)
    {
      if (ptr[i] >= 0 && ptr[i] < _numberOfZones)
      {
        Zone* z = _zones[ptr[i]];
        ZoneImplDL& zImpl = *static_cast<ZoneImplDL*>(z->ptr_impl);
        if (zImpl._DLsolid != 0) {glDeleteLists(zImpl._DLsolid, 1); zImpl._DLsolid = 0;}
        if (permanent == 1) zImpl._GPUResUse = 0;
      }
    }
  }
  delete [] ptr;
}
