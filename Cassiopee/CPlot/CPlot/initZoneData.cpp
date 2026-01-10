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
# include "Data.h"
# include "cplot.h"
#include <time.h>
#if defined(_WIN32) || defined(_WIN64)
#  include <winsock.h>
#endif

using namespace K_FLD;
using namespace std;

//=============================================================================
// Initialise les donnees des zones a partir des arrays.
// si retourne 0: echec
// si retourne 1: OK
//=============================================================================
E_Int Data::initZoneData(
  vector<FldArrayF*>& structF,
  vector<char*>& structVarString,
  vector<E_Int>& nit,
  vector<E_Int>& njt,
  vector<E_Int>& nkt,
  vector<FldArrayF*>& unstrF,
  vector<char*>& unstrVarString,
  vector<FldArrayI*>& cnt,
  vector<char*>& eltType,
  vector<char*>& zoneNames,
  vector<char*>& zoneTags,
  E_Int referenceNfield,
  char** referenceVarNames)
{
  // Calcul de la position des coords dans chaque arrays
  E_Int posx, posy, posz;
  vector<E_Int> sposx; vector<E_Int> sposy; vector<E_Int> sposz;
  E_Int sSize = structF.size();
  E_Int uSize = unstrF.size();
  for (E_Int i = 0; i < sSize; i++)
  {
    posx = K_ARRAY::isCoordinateXPresent(structVarString[i]);
    posy = K_ARRAY::isCoordinateYPresent(structVarString[i]);
    posz = K_ARRAY::isCoordinateZPresent(structVarString[i]);
    if (posx == -1 || posy == -1 || posz == -1) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "display: coordinates must be present in array.");
      for (E_Int i = 0; i < sSize; i++)
      {
        delete structF[i];
      }
      for (E_Int i = 0; i < uSize; i++)
      {
        delete unstrF[i]; delete cnt[i];
        delete [] eltType[i];
      }
      return 0;
    }
    sposx.push_back(posx+1); sposy.push_back(posy+1); sposz.push_back(posz+1);
  }
  vector<E_Int> uposx; vector<E_Int> uposy; vector<E_Int> uposz;
  for (E_Int i = 0; i < uSize; i++)
  {
    posx = K_ARRAY::isCoordinateXPresent(unstrVarString[i]);
    posy = K_ARRAY::isCoordinateYPresent(unstrVarString[i]);
    posz = K_ARRAY::isCoordinateZPresent(unstrVarString[i]);
    if (posx == -1 || posy == -1 || posz == -1) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "display: coordinates must be present in array.");
      for (E_Int i = 0; i < sSize; i++)
      {
        delete structF[i];
      }
      for (E_Int i = 0; i < uSize; i++)
      {
        delete unstrF[i]; delete cnt[i];
      }
      return 0;
    }
    uposx.push_back(posx+1); uposy.push_back(posy+1); uposz.push_back(posz+1);
  }

  // Recuperation des precedents pointeurs
  E_Int numberOfStructZonesp = _numberOfStructZones;
  E_Int numberOfUnstructZonesp = _numberOfUnstructZones;
  Zone** zonesp = _zones;
  StructZone** szonesp = _szones;
  UnstructZone** uzonesp = _uzones;

  // Remplissage du container de donnees
  E_Int numberOfStructZones = structF.size();
  E_Int numberOfUnstructZones = unstrF.size();
  E_Int numberOfZones = numberOfStructZones + numberOfUnstructZones;
  Zone** zones = (Zone**)malloc(numberOfZones*sizeof(Zone*));
  StructZone** szones = (StructZone**)malloc(numberOfStructZones*
                                             sizeof(StructZone*));
  UnstructZone** uzones = (UnstructZone**)malloc(numberOfUnstructZones*
                                                 sizeof(UnstructZone*));
  E_Int zoneNamesSize = zoneNames.size();
  E_Int zoneTagsSize = zoneTags.size();
  char zoneName[MAXSTRINGLENGTH];
  char* zTags;

  // mise a jour des pointeurs + BB pour les grilles struct
  for (E_Int i = 0; i < sSize; i++)
  {
    if (i < zoneNamesSize) strcpy(zoneName, zoneNames[i]);
    else sprintf(zoneName, "S-Zone " SF_D_, i);
    if (i < zoneTagsSize) zTags = zoneTags[i];
    else zTags = NULL;
    szones[i] = createStructZone(structF[i], structVarString[i],
                                 sposx[i], sposy[i], sposz[i],
                                 nit[i], njt[i], nkt[i],
                                 zoneName, zTags,
                                 referenceNfield, referenceVarNames);
    StructZone& z = *(szones[i]);
 
    // Essai de retrouver les reglages dans la previous zonesp
    // Pour les zones inactives en previous, elles ne doivent pas
    // avoir ete modifiees
    for (E_Int j = 0; j < numberOfStructZonesp; j++)
    {
      StructZone& zp = *(szonesp[j]);
      if (zp.active == 0)
      {
        E_Int ind = z.ni * z.nj * z.nk - 1;
        if (z.ni == zp.ni &&
            z.nj == zp.nj &&
            z.nk == zp.nk &&
            z.x[0] == zp.x[0] &&
            z.y[0] == zp.y[0] &&
            z.z[0] == zp.z[0] &&
            z.x[ind] == zp.x[ind] &&
            z.y[ind] == zp.y[ind] &&
            z.z[ind] == zp.z[ind] &&
            z.x[ind/2] == zp.x[ind/2] &&
            z.y[ind/2] == zp.y[ind/2] &&
            z.z[ind/2] == zp.z[ind/2])
          z.active = 0;
      }
    }

    // Pour le reste, il faut que la grille soit toujours positionnees
    // identiquement
    if (numberOfStructZonesp == numberOfStructZones)
    {
      StructZone& zp = *(szonesp[i]);
      if (z.ni == zp.ni &&
          z.nj == zp.nj &&
          z.nk == zp.nk)
      {
        z.activePlane = zp.activePlane;
        z.iPlane = zp.iPlane;
        z.jPlane = zp.jPlane;
        z.kPlane = zp.kPlane;
        z.iLine = zp.iLine;
        z.jLine = zp.jLine;
        z.kLine = zp.kLine;
        z.blank = zp.blank;
        //z.active = zp.active;
        z.selected = zp.selected;
      }
    }
  }

  // mise a jour des pointeurs + BB pour les grilles unstruct
  for (E_Int i = 0; i < uSize; i++)
  {
    if (i+sSize < zoneNamesSize) strcpy(zoneName, zoneNames[i+sSize]);
    else sprintf(zoneName, "U-Zone " SF_D_, i);
    if (i+sSize < zoneTagsSize) zTags = zoneTags[i+sSize];
    else zTags = NULL;
    uzones[i] = createUnstrZone(
      unstrF[i], unstrVarString[i],
      uposx[i], uposy[i], uposz[i],
      cnt[i], eltType[i],
      zoneName, zTags,
      referenceNfield, referenceVarNames);
    UnstructZone& z = *(uzones[i]);

    // Essai de retrouver les reglages dans la previous zonesp
    // Pour les zones inactives en previous, elles ne doivent pas
    // avoir ete modifiees
    for (E_Int j = 0; j < numberOfUnstructZonesp; j++)
    {
      UnstructZone& zp = *(uzonesp[j]);
      if (zp.active == 0)
      {
        E_Int ind = z.npts-1;
        if (z.ne == zp.ne &&
            z.npts == zp.npts &&
            z.x[0] == zp.x[0] &&
            z.y[0] == zp.y[0] &&
            z.z[0] == zp.z[0] &&
            z.x[ind] == zp.x[ind] &&
            z.y[ind] == zp.y[ind] &&
            z.z[ind] == zp.z[ind] &&
            z.x[ind/2] == zp.x[ind/2] &&
            z.y[ind/2] == zp.y[ind/2] &&
            z.z[ind/2] == zp.z[ind/2])
          z.active = 0;
      }
    }

    // Autres reglages
    if (numberOfUnstructZonesp == numberOfUnstructZones)
    {
      UnstructZone& zp = *(uzonesp[i]);
      if (z.material == 6) // volumetric
      {
        if (z.xmin == zp.xmin &&
            z.xmax == zp.xmax &&
            z.ymin == zp.ymin &&
            z.ymax == zp.ymax &&
            z.zmin == zp.zmin &&
            z.zmax == zp.zmax)
        {
          z.blank = uzonesp[i]->blank;
          z.selected = zp.selected;
        }
      }
      
      if (z.ne == zp.ne && z.npts == zp.npts)
      {
        z.blank = zp.blank;
        //z.active = zp.active;
        z.selected = zp.selected;
      }
    }
  }

  CPlotState& state = *ptrState;
  // Mise a jour du pointeur zone global
  E_Int structFSize = structF.size();
  for (E_Int i = 0; i < structFSize; i++)
  {
    zones[i] = szones[i];
  }
  E_Int unstrFSize = unstrF.size();
  for (E_Int i = 0; i < unstrFSize; i++)
  {
    zones[i+structF.size()] = uzones[i];
  }

  // On assure que le mode est valide
  if (state.mode == SCALARFIELD)
  {
    if (numberOfZones > 0 && zones[0]->nfield <= state.scalarField)
    { state.scalarField = 0; state.mode = MESH; }
  }
  else if (state.mode == VECTORFIELD)
  {
    if (numberOfZones > 0)
    {
      if (zones[0]->nfield <= state.vectorField1) 
      { state.vectorField1 = 0; state.mode = MESH; }
      if (zones[0]->nfield <= state.vectorField2) 
      { state.vectorField2 = 0; state.mode = MESH; }
      if (zones[0]->nfield <= state.vectorField3) 
      { state.vectorField3 = 0; state.mode = MESH; }
    }
  }

  // Switch - Dangerous zone protegee par state.lock
  if (state.selectedZone >= numberOfZones) state.selectedZone = 0; // RAZ selected zone
  state.kcursor = 0; // RAZ clavier
  state.syncDisplay();
  _zones = zones;
  _szones = szones;
  _uzones = uzones;
  _numberOfStructZones = numberOfStructZones;
  _numberOfUnstructZones = numberOfUnstructZones;
  _numberOfZones = numberOfZones;
  
  // Mise a jour des min-max globaux
  globMinMax(_zones, _numberOfZones, 
             xmin, xmax, ymin, ymax, zmin, zmax, epsup, epsstrafe, dmoy);
  globFMinMax(_zones, _numberOfZones, minf, maxf);
  
  // Free the previous
  E_Int i;
  for (i = 0; i < numberOfStructZonesp; i++) delete szonesp[i];
  for (i = 0; i < numberOfUnstructZonesp; i++) delete uzonesp[i];
  if (szonesp != NULL) free(szonesp);
  if (uzonesp != NULL) free(uzonesp);
  if (zonesp != NULL) free(zonesp);

  // Modifie les zones pour le rendu volumetrique
  replaceVolumetricZones();
  
  return 1;
}

//=============================================================================
// Remplace une zone qui est en rendu volumetrique par un parallepipede 
// correspondant a la BBOX
//=============================================================================
void Data::replaceVolumetricZones()
{
  for (E_Int i = 0; i < _numberOfZones; i++)
  {
    Zone& z = *_zones[i];
    if (z.material == 6) // smoke
    { 
      if (i < _numberOfStructZones)
      { // Cree une nouvelle zone structuree et recopie le max de trucs   
        StructZone* znp = new StructZone(ptrState, createZoneImpl());
        StructZone& zn = *znp;
        strcpy(zn.zoneName, z.zoneName);
        zn.npts = 8;
        zn.ni = 2; zn.nj = 2; zn.nk = 2;
        zn.nfield = z.nfield;
        zn.dim = z.dim;
        zn.x = new E_Float[zn.npts];
        double factor = 0.3;
        double deltax = factor*(z.xmax - z.xmin);
        double deltay = factor*(z.ymax - z.ymin);
        double deltaz = factor*(z.zmax - z.zmin);
        zn.x[0] = z.xmin - deltax; zn.x[1] = z.xmax + deltax;
        zn.x[2] = z.xmin - deltax; zn.x[3] = z.xmax + deltax;
        zn.x[4] = z.xmin - deltax; zn.x[5] = z.xmax + deltax;
        zn.x[6] = z.xmin - deltax; zn.x[7] = z.xmax + deltax;
        zn.y = new E_Float[zn.npts];
        zn.y[0] = z.ymin - deltay; zn.y[1] = z.ymin - deltay;
        zn.y[2] = z.ymax + deltay; zn.y[3] = z.ymax + deltay;
        zn.y[4] = z.ymin - deltay; zn.y[5] = z.ymin - deltay;
        zn.y[6] = z.ymax + deltay; zn.y[7] = z.ymax + deltay;
        zn.z = new E_Float[zn.npts];
        zn.z[0] = z.zmin - deltaz; zn.z[1] = z.zmin - deltaz;
        zn.z[2] = z.zmin - deltaz; zn.z[3] = z.zmin - deltaz;
        zn.z[4] = z.zmax + deltaz; zn.z[5] = z.zmax + deltaz;
        zn.z[6] = z.zmax + deltaz; zn.z[7] = z.zmax + deltaz;
        
        for (E_Int n = 0; n < zn.nfield; n++)
        {
          strcpy(zn.varnames[n], z.varnames[n]);
          zn.f[n] = new E_Float [zn.npts];
          memset(zn.f[n], 0, zn.npts*sizeof(E_Float));
        }
        strcpy(zn.renderTag, z.renderTag);
        zn.colorR = z.colorR;
        zn.colorG = z.colorG;
        zn.colorB = z.colorB;
        zn.material = z.material;
        zn.blending = z.blending;
        zn.meshOverlay = z.meshOverlay;
        zn.meshColorR = z.meshColorR;
        zn.meshColorG = z.meshColorG;
        zn.meshColorB = z.meshColorB;
        zn.meshWidth = z.meshWidth;
        zn.shaderParam1 = z.shaderParam1;
        zn.shaderParam2 = z.shaderParam2;

        zn.compNorm();
        zn.blank = z.blank;
        zn.active = z.active;
        zn.selected = z.selected;
        findMinMax(&zn); findFMinMax(&zn);
  
        // Voxelize la zone
        voxelize(zn, (StructZone&)z);

        delete _zones[i]; _szones[i] = znp;
        _zones[i] = znp;
      }
      else
      {
        // Cree une nouvelle zone non structuree et recopie le max de trucs
        UnstructZone* znp = new UnstructZone(ptrState, createZoneImpl());
        UnstructZone& zn = *znp;
        strcpy(zn.zoneName, z.zoneName);
        zn.npts = 8;
        zn.np = zn.npts;
        zn.nfield = z.nfield;
        zn.dim = z.dim;
        double factor = 0.3;
        double deltax = factor*(z.xmax - z.xmin);
        double deltay = factor*(z.ymax - z.ymin);
        double deltaz = factor*(z.zmax - z.zmin);
        zn.x = new E_Float[zn.npts];
        zn.x[0] = z.xmin - deltax; zn.x[1] = z.xmax + deltax;
        zn.x[2] = z.xmax + deltax; zn.x[3] = z.xmin - deltax;
        zn.x[4] = z.xmin - deltax; zn.x[5] = z.xmax + deltax;
        zn.x[6] = z.xmax + deltax; zn.x[7] = z.xmin - deltax;
        zn.y = new E_Float[zn.npts];
        zn.y[0] = z.ymin - deltay; zn.y[1] = z.ymin - deltay;
        zn.y[2] = z.ymax + deltay; zn.y[3] = z.ymax + deltay;
        zn.y[4] = z.ymin - deltay; zn.y[5] = z.ymin - deltay;
        zn.y[6] = z.ymax + deltay; zn.y[7] = z.ymax + deltay;
        zn.z = new E_Float[zn.npts];
        zn.z[0] = z.zmin - deltaz; zn.z[1] = z.zmin - deltaz;
        zn.z[2] = z.zmin - deltaz; zn.z[3] = z.zmin - deltaz;
        zn.z[4] = z.zmax + deltaz; zn.z[5] = z.zmax + deltaz;
        zn.z[6] = z.zmax + deltaz; zn.z[7] = z.zmax + deltaz;

        zn.ne = 1;
        zn.nec.push_back(1);
        
        for (E_Int n = 0; n < zn.nfield; n++)
        {
          strcpy(zn.varnames[n], z.varnames[n]);
          zn.f[n] = new E_Float [zn.npts];
          memset(zn.f[n], 0, zn.npts*sizeof(E_Float));
        }
        strcpy(zn.renderTag, z.renderTag);
        zn.colorR = z.colorR;
        zn.colorG = z.colorG;
        zn.colorB = z.colorB;
        zn.material = z.material;
        zn.blending = z.blending;
        zn.meshOverlay = z.meshOverlay;
        zn.meshColorR = z.meshColorR;
        zn.meshColorG = z.meshColorG;
        zn.meshColorB = z.meshColorB;
        zn.meshWidth = z.meshWidth;
        zn.shaderParam1 = z.shaderParam1;
        zn.shaderParam2 = z.shaderParam2;

        zn.eltType.push_back(7);
        zn.eltSize.push_back(8);
        
        zn.connect.push_back(new E_Int[zn.nec[0]*zn.eltSize[0]]);
        E_Int* cn = zn.connect[0];

        //cn[0] = 1; cn[6] = 2; cn[12] = 3; cn[18] = 4;
        //cn[1] = 5; cn[7] = 6; cn[13] = 7; cn[19] = 8;
        //cn[2] = 1; cn[8] = 2; cn[14] = 6; cn[20] = 5;
        //cn[3] = 4; cn[9] = 3; cn[15] = 7; cn[21] = 8;
        //cn[4] = 1; cn[10] = 4; cn[16] = 8; cn[22] = 5;
        //cn[5] = 2; cn[11] = 3; cn[17] = 7; cn[23] = 6;
        
        cn[0] = 1; cn[1] = 2; cn[2] = 3; cn[3] = 4; cn[4] = 5;
        cn[5] = 6; cn[6] = 7; cn[7] = 8;

        //zn.surf = NULL; 
        zn.compNorm();
        zn.blank = z.blank;
        zn.active = z.active;
        zn.selected = z.selected;
        findMinMax(&zn); findFMinMax(&zn);
  
        // Voxelize la zone
        voxelize(zn, (UnstructZone&)z);

        delete _zones[i]; _uzones[i-_numberOfStructZones] = znp;
        _zones[i] = znp;
      }
    }
  }
}
