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

#include "Zone.h"
#include "ZoneImpl.h"
#include <stdio.h>
#if defined (_WIN32) ||defined(_WIN64)
#  include <windows.h>
#endif
#if defined(__APPLE__)
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

//=============================================================================
Zone::Zone( CPlotState* states, ZoneImpl* impl ) : ptr_impl(impl)
{
  dim = 3;
  nfield = 0;
  npts = 0;
  x = NULL; y = NULL; z = NULL;
  f = NULL;
  varnames = NULL;
  surf = NULL;
  minf = NULL; maxf = NULL;
  xmin = 1.e6; ymin = 1.e6; zmin = 1.e6;
  xmax = -1.e6; ymax = -1.e6; zmax = -1.e6;
  active = 0;
  selected = 0;
  previouslySelected = 0;
  blank = 0;
  shaderParam1 = 1.; shaderParam2 = 1.;
  _voxelArray = NULL;
  ptrState = states;
}

//=============================================================================
Zone::~Zone()
{
  if (x != NULL) delete [] x;
  if (y != NULL) delete [] y;
  if (z != NULL) delete [] z;
  if (surf != NULL) delete [] surf;
  if (f != NULL)
  {
    for (int i = 0; i < nfield; i++) delete [] f[i];
    delete [] f;
  }
  if (varnames != NULL)
  {
    for (int i = 0; i < nfield; i++) delete [] varnames[i];
    delete [] varnames;
  }
  if (minf != NULL) delete [] minf;
  if (maxf != NULL) delete [] maxf;
  
  if (regtexu != NULL) delete [] regtexu;
  if (regtexv != NULL) delete [] regtexv;
  
  delete [] _voxelArray;
  ptr_impl->freeGPURes(ptrState);
  delete ptr_impl;
}
// ------------------------------------------------------------------------
void
Zone::freeGPURessources( bool useGPURessources, bool freeIso )
{
  ptr_impl->freeGPURes( ptrState, freeIso );
  if ( useGPURessources ) setUseGPURessources();
  else unsetUseGPURessources();
}
