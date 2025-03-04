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
#include "../Data.h"

//=============================================================================
/* 
   loadPlugins.
   This function loads plugins functions in the plugins structure.
   Called once at beginning.
*/
//=============================================================================
void Data::loadPlugins()
{
  struct chain_function_void* pv;
  struct chain_function_void2* pv2;
  struct chain_function_void3* pv3;
  struct chain_function_int* pi;
  struct chain_function_double* pd;

  // -- Init of plugins structure -- 
  _plugins.colorMap = NULL;
  _plugins.lookFor = NULL;
  _plugins.addAVariable = NULL;
  _plugins.blanking = NULL;
  _plugins.selectNextZone = NULL;
  _plugins.selectPreviousZone = NULL;
  _plugins.mouseClick = NULL;
  _plugins.mouseRightClick = NULL;
  _plugins.screenDump = NULL;

  // -- Default blanking functions --
  // Blanking by cellN
  _plugins.blanking = 
    (struct chain_function_int*)malloc(sizeof(struct chain_function_int));
  pi = _plugins.blanking;
  strcpy(pi->functionName, "Blanking by cellN");
  strcpy(pi->varName, "cellN");
  pi->f = blankCellN;
  pi->next = NULL;
  // Blanking by cellNF
  pi->next = 
    (struct chain_function_int*)malloc(sizeof(struct chain_function_int));
  pi = pi->next;
  strcpy(pi->functionName, "Blanking by cellNF");
  strcpy(pi->varName, "cellNF");
  pi->f = blankCellNF;
  pi->next = NULL;
  // Blanking by status
  pi->next = 
    (struct chain_function_int*)malloc(sizeof(struct chain_function_int));
  pi = pi->next;
  strcpy(pi->functionName, "Blanking by status");
  strcpy(pi->varName, "status");
  pi->f = blankStatus;
  pi->next = NULL;

//   // -- Default add-a-variable functions --
//   // Add  pressure variable
//   _plugins.addAVariable = 
//     (struct chain_function_void*)malloc(sizeof(struct chain_function_void));
//   pv = _plugins.addAVariable;
//   strcpy(pv->functionName, "Pressure variable");
//   pv->f = addPressureVariable;
//   pv->next = NULL;
// //   pv->next = NULL;
  
  // -- Default colormap functions --
  
  // Blue to red
  _plugins.colorMap = 
    (struct chain_function_double*)malloc(sizeof(struct chain_function_double));
  pd = _plugins.colorMap;
  strcpy(pd->functionName, "Blue to red colormap");
  strcpy(pd->varName, "0");
  pd->f = colBlueToRed;
  pd->next = NULL;

  // col2RGB
  pd->next =
    (struct chain_function_double*)malloc(sizeof(struct chain_function_double));
  pd = pd->next;
  strcpy(pd->functionName, "Bi-color RGB colormap");
  strcpy(pd->varName, "1");
  pd->f = col2RGB;
  pd->next = NULL;

  // col2HSV
  pd->next =
    (struct chain_function_double*)malloc(sizeof(struct chain_function_double));
  pd = pd->next;
  strcpy(pd->functionName, "Bi-color HSV colormap");
  strcpy(pd->varName, "2");
  pd->f = col2HSV;
  pd->next = NULL;

  // col3RGB
  pd->next =
    (struct chain_function_double*)malloc(sizeof(struct chain_function_double));
  pd = pd->next;
  strcpy(pd->functionName, "Tri-color RGB colormap");
  strcpy(pd->varName, "3");
  pd->f = col3RGB;
  pd->next = NULL;

  // col3HSV
  pd->next =
    (struct chain_function_double*)malloc(sizeof(struct chain_function_double));
  pd = pd->next;
  strcpy(pd->functionName, "Tri-color HSV colormap");
  strcpy(pd->varName, "4");
  pd->f = col3HSV;
  pd->next = NULL;

  // colMRGB
  pd->next =
    (struct chain_function_double*)malloc(sizeof(struct chain_function_double));
  pd = pd->next;
  strcpy(pd->functionName, "Multi-color RGB colormap");
  strcpy(pd->varName, "5");
  pd->f = colMRGB;
  pd->next = NULL;

  // colMHSV
  pd->next =
    (struct chain_function_double*)malloc(sizeof(struct chain_function_double));
  pd = pd->next;
  strcpy(pd->functionName, "Multi-color HSV colormap");
  strcpy(pd->varName, "6");
  pd->f = colMHSV;
  pd->next = NULL;

  // Diverging
  pd->next =
    (struct chain_function_double*)malloc(sizeof(struct chain_function_double));
  pd = pd->next;
  strcpy(pd->functionName, "Diverging colormap");
  strcpy(pd->varName, "7");
  pd->f = diverging;
  pd->next = NULL;

  // Zone color map plugins
  _plugins.zoneColorMap = 
    (struct chain_function_double*)malloc(sizeof(struct chain_function_double));
  pd = _plugins.zoneColorMap;
  strcpy(pd->functionName, "Green to red colormap");
  strcpy(pd->varName, "0");
  pd->f = colGreenToRed;
  pd->next = NULL;

  // -- Default screenDump functions -- 
  // RGB ppm screen dump
  _plugins.screenDump = 
    (struct chain_function_void2*)malloc(sizeof(struct chain_function_void2));
  pv2 = _plugins.screenDump;
  strcpy(pv2->functionName, "PPM screen dump");
  strcpy(pv2->extension, "ppm");
  pv2->f = writePPMFile;
  pv2->next = NULL;
  // RGB png screen dump
  pv2->next = 
    (struct chain_function_void2*)malloc(sizeof(struct chain_function_void2));
  pv2 = pv2->next;
  strcpy(pv2->functionName, "PNG screen dump");
  strcpy(pv2->extension, "png");
  pv2->f = writePNGFile;
  pv2->next = NULL;
  // RGB png screen dump + depth dans alpha channel
  pv2->next = 
    (struct chain_function_void2*)malloc(sizeof(struct chain_function_void2));
  pv2 = pv2->next;
  strcpy(pv2->functionName, "PNG + depth screen dump");
  strcpy(pv2->extension, "dpng");
  pv2->f = writeDPNGFile; 
  pv2->next = NULL;
  // mpeg screen dump
  pv2->next = 
    (struct chain_function_void2*)malloc(sizeof(struct chain_function_void2));
  pv2 = pv2->next;
  strcpy(pv2->functionName, "MPEG screen dump");
  strcpy(pv2->extension, "mpeg");
  pv2->f = writeMPEGFrame; 
  pv2->next = NULL;
  // postscript screen dump
  pv2->next =
     (struct chain_function_void2*)malloc(sizeof(struct chain_function_void2));
  pv2 = pv2->next;
  strcpy(pv2->functionName, "Postscript screen dump");
  strcpy(pv2->extension, "ps");
  pv2->f = NULL; // this is handled another way
  pv2->next = NULL;

  // -- Default select functions --
  // Incremental selection
  _plugins.selectNextZone =
     (struct chain_function_void*)malloc(sizeof(struct chain_function_void));
  pv = _plugins.selectNextZone;
  strcpy(pv->functionName, "Incremental selection");
  pv->f = selectNextZoneIncr;
  pv->next = NULL;
  _plugins.selectPreviousZone =
     (struct chain_function_void*)malloc(sizeof(struct chain_function_void));
  pv = _plugins.selectPreviousZone;
  strcpy(pv->functionName, "Incremental selection");
  pv->f = selectPreviousZoneIncr;
  pv->next = NULL;

  // -- Default mouse functions --
  // Click selects
  _plugins.mouseClick =
     (struct chain_function_void3*)malloc(sizeof(struct chain_function_void3));
  pv3 = _plugins.mouseClick;
  strcpy(pv3->functionName, "Shift+Left Click selects");
  pv3->f = mouseClickSelect;
  pv3->next = NULL;

  // Click multiple selects
  _plugins.mouseMultipleClick =
    (struct chain_function_void3*)malloc(sizeof(struct chain_function_void3));
  pv3 = _plugins.mouseMultipleClick;
  strcpy(pv3->functionName, "Shift+Ctrl+Left Click multiple selects");
  pv3->f = mouseClickMultipleSelect;
  pv3->next = NULL;

  // Click accurate selects
  _plugins.mouseAccurateClick =
    (struct chain_function_void3*)malloc(sizeof(struct chain_function_void3));
  pv3 = _plugins.mouseAccurateClick;
  strcpy(pv3->functionName, "Ctrl+Left Click accurate selects");
  pv3->f = mouseClickAccurateSelect;
  pv3->next = NULL;

  // Right click hide 
  _plugins.mouseRightClick = 
    (struct chain_function_void3*)malloc(sizeof(struct chain_function_void3));
  pv3 = _plugins.mouseRightClick;
  strcpy(pv3->functionName, "Shift+Right Click unselect all");
  pv3->f = mouseRightClickSelect;
  pv3->next = NULL;

  // -- Default look-for functions --
  // Find active zone
  _plugins.lookFor =
     (struct chain_function_void*)malloc(sizeof(struct chain_function_void));
  pv = _plugins.lookFor;
  strcpy(pv->functionName, "Find selected zone");
  pv->f = lookForActiveZone;
  pv->next = NULL;

  // Find point of max value
  pv->next = 
    (struct chain_function_void*)malloc(sizeof(struct chain_function_void));
  pv = pv->next;
  strcpy(pv->functionName, "Find point of max value");
  pv->f = lookForMaxValue;
  pv->next = NULL;

  // Find point of min value
  pv->next = 
    (struct chain_function_void*)malloc(sizeof(struct chain_function_void));
  pv = pv->next;
  strcpy(pv->functionName, "Find point of min value");
  pv->f = lookForMinValue;
  pv->next = NULL;

  // -- Add plugins functions by loading external .so -- 
  
}

//=============================================================================
/* 
   loadPluginsPath.
   Plugin directories are found in $HOME/.cplot file
*/
//=============================================================================
void Data::loadPluginsPath()
{
  
}

//=============================================================================
/* 
   autoPlugins.
   Set automatic plugins.
   Those plugins are automatically triggered/computed at beginning or
   at each reload.
*/
//=============================================================================
void Data::autoPlugins()
{
  // if no zone exists, return
  if (_numberOfZones == 0) return;

  if (ptrState->autoblank == 0)
  {
    _pref.blanking = NULL; findBlankedZones(); return;
  };

  // Blanking, if cellN or cellNF or status exists
  if (checkVariable(0, "cellN") != -1)
  {
    _pref.blanking = _plugins.blanking;
    if (_pref.blanking == NULL)
      printf("CPlot: no blanking available.\n");
    else
      findBlankedZones();
  }
  
  else if (checkVariable(0, "cellNF") != -1)
  {
    _pref.blanking = _plugins.blanking;
    if (_pref.blanking != NULL)
    {
      _pref.blanking = _pref.blanking->next;
      if (_pref.blanking == NULL)
        printf("CPlot: no blanking available.\n");
      else
        findBlankedZones();
    }
  }

  else if (checkVariable(0, "status") != -1)
  {
    _pref.blanking = _plugins.blanking;
    if (_pref.blanking != NULL)
    {
      _pref.blanking = _pref.blanking->next;
      _pref.blanking = _pref.blanking->next;
      if (_pref.blanking == NULL)
        printf("CPlot: no blanking available.\n");
      else
        findBlankedZones();
    }
  }

  else
  {
    for (int i = 0; i < _numberOfZones; i++)
    {
      _zones[i]->blank = 0;
    }
  }
}
